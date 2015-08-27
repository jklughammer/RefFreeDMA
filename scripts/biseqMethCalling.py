#! /usr/bin/env python
'''
biseqMethCalling.py
(C) Christoph Bock, Natalie Jaeger, Fabian Mueller
version: 2011-07-25

----------------------------------------------------------------------------------------------
Methylation calling for bisulfite sequencing methods: 
- Whole genome bisulfite sequencing (WGBS)
- Reduced-representation bisulfite sequencing (RRBS)
- Hybrid-selection bisulfite sequencing
- Other variants of region-specific bisulfite sequencing (e.g. MethylCap-biseq and ChIP-biseq)

By default CpG are analyzed, so in the descriptions, we often refer to CpGs as the general Idea of Methylated 'patterns'. 
However, keep in mind that the sequence patterns are configurated and may well be CpA, CpH, CpHpG, etc.

INPUT:  - (multiple) BAM files containing containing read alignments
        - parameters (see --help for more info)
OUTPUT: - statistics (see the description in StatisticsCollection.writeStatistics2files() for details)
        - CpGs with genomic coordinates and methylation ratios as a bed file
        - Fragment coverage: coordinates, readcounts and CpG methylation ratios of each generated fragment (RRBS: MspI digestion, WGBSS: Paired end reads)
        - Reads: coords (read and fragment the read is on), CpG (possibly other patterns) of __processed__ reads
        - ReadDetails: more detailed read info for processed reads including base quality details, alignment/mismatch details, conversion rates, ...
        - runtime profile files for each processed region, if desired
        - if the verbose option is set: information on the processed windows/regions
        - log files
        - a vast amount of temporary files for each processed region

SCRIPT OUTLINE:
The entire genome is processed in the order

(A) for chrom in chromosomes:
(B)   for region in sliding_windows_across_chromosome:
(C)        for read in all_reads_in_all_input_lanes_in_current_region:
(D)            analyze_read_for_methylation
(E) summarize_regions

Steps A and B are distributed by manageProcesses() which iterates over all chromosomes and chops them up into subtasks
Step C is taken care of by sequenceDataProcessing() which retrieves all the reads for the current region, and processes them,
generates region specific outputs for Reads,ReadDetails,CpG,Fragments and Statistics
Step D is perfomed by methylationCalling_bitarray() which identifies patterns/CpGs in each read and updates the region specific
summary of all CpGs in the region. The main class which is operated upon is the bsAlignedRead class which stores all necassary read
associated information.
Step E is taken care of by performAnalysis_wrapup() which combines the output files generated for each region, converts them
(to the desired genome version on compressed format, if required), and unifies the statistics into a global statistics (objects of the class
StatisticsCollection). Also deletes temporary files and performs other cleanup options

For parallelization purposes, the script can be called as one of three tasks:
'main'        : basically levels A and B form the structure above
                prepares the parameters and invokes manageProcesses() to devide up the genome into subtasks (ready to be processes in a 'region' tasked
                session of the script)
                invokes 'region' tasks and one 'wrapup' task which depends upon completion of the region tasks
'region'      : basically Step C and D from the structure above
'wrapup'      : Step E
'diagnostics' : Print diagnostics about processed regions, submitted jobs and errors, warnings, notes in the log files

If parallelization option is set the 'region' and 'wrapup' tasks are submitted to the compute cluster (as specified in the parameters). Otherwise they are invoked as
system commands on the same machine as the 'main' task. Communication between the tasks is done via pickle files.
'''
#--------Changes by Johanna Klughammer 25.07.2015----------#
#change x.tostring() to str(x)




#TODOS:
# - whole genome scale!
# - nonCpG conversion rates currently are not precisely defined
# - sort errors in readDetails files (e.g. see nonCG methylation project), would probably also resolve memory problems in wrapup
# - pairedEnd for SpikeInControls
# - aligned Reads matching to spike ins?
# - nonCpG methylation (nonsymmetrical patterns: adjust pattern coordinates): i.e. remove motifIsSymmetrical. implemented, but left the comment till the final version of the nonCpG methylation project goes through;

#IDEAS:
# - SNP invariance by taking into account the complementary strand
#Known Issues:
#- in paired end processing, the mate is not checked for quality criteria (mismatches, base qualities, restriction site, ...) for determining the bases which overlap with the current read. Therefore,
#  some overlapping bases between the read pairs will be excluded from the analysis. Implementing this feature is postponed, as it takes a significant amount of time and computing resources
#- in the statistics of medians, although not mathematically sound medians of medians are computed as an approximation. Precise computation would require a significant amount of computing resources
# core libraries
import cPickle
import math
import os
import shutil
import string
import re
import subprocess
import sys
import time
import random
# additional libraries
#import Bio
import Bio.Seq
import Bio.SeqIO

import pysam
from bitarray import bitarray
#profiling
import guppy
import cProfile


# custom libraries
#import downloadGenome

# operating parameters
genomePath = "" # directory containing the genomic DNA sequence ordered by chromosomes -> global variable initialized in performAnalysis()

# biological parameters
chromSizeHash = {} # hashtable with chromosome names and lengths -> global variable initialized in performAnalysis()
seqMotifList = [] # base(s) adjacent to the cytosine that will be analyzed for DNA methylation, e.g. "G" for CpG methylation (global variable initialized in performAnalysis)
numberOfLanes = 0 # total number of sequencing lanes to be combined into the current sample (global variable initialized in performAnalysis)

# biological constants
ONE_MB = 1000000
#MIN_READ_FRAGMENT_LENGTH = 15 # minimum readlength after adapter clipping or fragment length for a read to be processed, see minimumReadLength()
BAD_QUALITY_SCORE = 3 # everything below will be regarded as a really bad read quality score
#WINDOW_AROUND_NONSTRAND_MATE = 5 # which windowsize (in multiples of readlength) should be looked up in the reference sequence if the mate is not on the same strand as the original read
SEQ_MOTIF_IS_SYMMETRICAL = {'G':True,'A':False,"C":False,"T":False,'HpG':True,'HpH':False,'H':False}
COMPLEMENTARY_BASES = {'G':'C','A':'T','C':'G','T':'A','H':'D'}

#lookup tables
manualGenomeCorrection = {"1":"chr1", "2":"chr2", "3":"chr3", "4":"chr4", "5":"chr5", "6":"chr6", "7":"chr7", "8":"chr8", "9":"chr9", "10":"chr10", "11":"chr11", "12":"chr12", "13":"chr13", "14":"chr14", "15":"chr15", "16":"chr16", "17":"chr17", "18":"chr18", "19":"chr19", "20":"chr20", "21":"chr21", "22":"chr22", "X":"chrX", "Y":"chrY", "chrMT":"chrM", "gi|149288852|ref|NC_000067.5|NC_000067 Mus musculus chromosome 1, reference assembly (C57BL/6J)":"chr1","gi|149338249|ref|NC_000068.6|NC_000068 Mus musculus chromosome 2, reference assembly (C57BL/6J)":"chr2","gi|149352351|ref|NC_000069.5|NC_000069 Mus musculus chromosome 3, reference assembly (C57BL/6J)":"chr3","gi|149354223|ref|NC_000070.5|NC_000070 Mus musculus chromosome 4, reference assembly (C57BL/6J)":"chr4","gi|149354224|ref|NC_000071.5|NC_000071 Mus musculus chromosome 5, reference assembly (C57BL/6J)":"chr5","gi|149361431|ref|NC_000072.5|NC_000072 Mus musculus chromosome 6, reference assembly (C57BL/6J)":"chr6","gi|149361432|ref|NC_000073.5|NC_000073 Mus musculus chromosome 7, reference assembly (C57BL/6J)":"chr7","gi|149361523|ref|NC_000074.5|NC_000074 Mus musculus chromosome 8, reference assembly (C57BL/6J)":"chr8","gi|149361524|ref|NC_000075.5|NC_000075 Mus musculus chromosome 9, reference assembly (C57BL/6J)":"chr9","gi|149288869|ref|NC_000076.5|NC_000076 Mus musculus chromosome 10, reference assembly (C57BL/6J)":"chr10","gi|149288871|ref|NC_000077.5|NC_000077 Mus musculus chromosome 11, reference assembly (C57BL/6J)":"chr11","gi|149292731|ref|NC_000078.5|NC_000078 Mus musculus chromosome 12, reference assembly (C57BL/6J)":"chr12","gi|149292733|ref|NC_000079.5|NC_000079 Mus musculus chromosome 13, reference assembly (C57BL/6J)":"chr13","gi|149292735|ref|NC_000080.5|NC_000080 Mus musculus chromosome 14, reference assembly (C57BL/6J)":"chr14","gi|149301884|ref|NC_000081.5|NC_000081 Mus musculus chromosome 15, reference assembly (C57BL/6J)":"chr15","gi|149304713|ref|NC_000082.5|NC_000082 Mus musculus chromosome 16, reference assembly (C57BL/6J)":"chr16","gi|149313536|ref|NC_000083.5|NC_000083 Mus musculus chromosome 17, reference assembly (C57BL/6J)":"chr17","gi|149321426|ref|NC_000084.5|NC_000084 Mus musculus chromosome 18, reference assembly (C57BL/6J)":"chr18","gi|149323268|ref|NC_000085.5|NC_000085 Mus musculus chromosome 19, reference assembly (C57BL/6J)":"chr19","gi|149361525|ref|NC_000086.6|NC_000086 Mus musculus chromosome X, reference assembly (C57BL/6J)":"chrX","gi|149361526|ref|NC_000087.6|NC_000087 Mus musculus chromosome Y, reference assembly (C57BL/6J)":"chrY" }
# addition for canFam3 (has 38 chromosomes, X)
manualGenomeCorrection.update({"23":"chr23", "24":"chr24", "25":"chr25", "26":"chr26", "27":"chr27", "28":"chr28", "29":"chr29", "30":"chr30", "31":"chr31", "32":"chr32", "33":"chr33", "34":"chr34", "35":"chr35", "36":"chr36", "37":"chr37", "38":"chr38"})

ASCII2PHRED = {} # conversion of ASCII-encoded into numeric quality scores -> global variable initialized in performAnalysis()
#readcount thresholds for statistics
SEQMOTIF_THRES = [0,1,2,5,10]
READ_THRES_POSITIONS = [1,5,10,50,100]
READ_THRES_WINDOWS = [0,1,5,10,50,100]

exceptionOccurred = False


def setGlobals(options):
    '''
    initialize (some) global variables
    others are set in performAnalysis()
    '''
    global debug
    debug = options.verbose
    global ASCII2PHRED
    ASCII2PHRED = {}
    for i in range(33,123): #end coordinate chosen arbitrarily. adapt if problems occurr
        ASCII2PHRED[chr(i)] = i - 33
        
def setGlobalsMaxMismatches(options):
    if options.mismatchTruncate > 0:
        MAX_ALLOWED_READLEN_FOR_MAX_MISMATCHES = 1000
        global DYNAMIC_MAX_MISMATCHES
        if options.maxMismatches > 0 and options.maxMismatches < 1:
            DYNAMIC_MAX_MISMATCHES = map(lambda x: int(round(options.maxMismatches*x)),range(1,MAX_ALLOWED_READLEN_FOR_MAX_MISMATCHES+1))
        else:
            DYNAMIC_MAX_MISMATCHES = MAX_ALLOWED_READLEN_FOR_MAX_MISMATCHES * [int(options.maxMismatches)]

def getDebugString(s):
    '''
    converts and formats the string s to a status string for output to console
    '''
    u1 = ((len(s) + 4) - 21) / 2
    u2 = ((len(s) + 4) - 21) / 2 + (len(s)+1)%2
    ss =  u1*"#" + " " + str(time.strftime("%Y-%m-%d %H:%M:%S")) + " " + u2*"#" + "\n"
    ss += "# " + s + " #\n"
    ss += (len(s)+4)*"#" + "\n"
    return ss

def isValidFile(fn):
    '''
    check if the file specified under fn causes an exception while trying to read it. Eg because the file does not exist
    '''
    try:
        f = open(fn,"r")
        f.close()
        return True
    except:
        return False
    
def listAsString(l):
    '''
    converts [a,b,c,d] into "a,b,c,d"
    '''
    return str(l)[1:-1].replace(" ","")

def getSepStringFromList(l,sep="\t",addquotes=False):
    '''
    return a string containing the elements of the list as tab separated
    '''
    s = ''
    qtstr = ''
    if addquotes:
        qtstr = '"'
    for elem in l:
        s += qtstr + str(elem) + qtstr + sep
    return s.rstrip(sep)

def strRevComp(s):
    '''
    transforms a DNA sequence given in a string to its reverse complement using Biopython functions
    '''
    return(str(Bio.Seq.reverse_complement(Bio.Seq.Seq(s))))

def strComp(s):
    '''
    transforms a DNA sequence given in a string to its complement using Biopython functions
    '''
    return(strRevComp(s)[::-1])

def mean(list):
    '''
    calculated the mean of a list
    '''
    length = len(list)
    if length == 0:
        return "NA"   
    sum = 0
    for l in list:
        sum += l
    return float(sum) / length
   
def stddev(list):
    '''
    calculated the standard deviation of a list
    '''
    if len(list) == 0 or len(list) == 1:
        return "NA" 
    meanV = mean(list)
    length = len(list)
    sum = 0
    for l in list:
        sum += pow(l-meanV,2)
    return math.sqrt(float(sum) / (length-1)) 

def median(theValues):
    '''
    calculated the median of a list, special NA handling
    '''
    # takes a sorted list as input
    if len(theValues) == 0:
        return "NA"    
    if len(theValues) % 2 == 1:
        return theValues[(len(theValues)+1)/2-1]
    else:
        lower = theValues[len(theValues)/2-1]
        upper = theValues[len(theValues)/2]
        if lower == "NA" and upper == "NA":
            return "NA" 
        elif lower == "NA":
            return upper 
        elif upper == "NA":
            return lower 
        else:
            return (float(lower + upper)) / 2

class genomeLocation:
    '''
    parent class for bsAlignedRead, UniqueCpgsPerLane, UniqueFragmentsPerLane
    important is the cmp method
    '''
    def __init__(self):
        self.chrom = ""
        self.chromstart = 0
        self.chromend = 0
        self.strand = ""
        
    def __cmp__(self, other):
        #alphanumerical chromosome order        
        cmpVal = cmp(manualGenomeCorrection.get(self.chrom,self.chrom),manualGenomeCorrection.get(other.chrom,other.chrom))
        if cmpVal == 0:
            cmpVal = cmp(self.chromstart,other.chromstart)
            if cmpVal == 0:
                cmpVal = cmp(self.chromend,other.chromend)
        return cmpVal


def getGappedSeqRead(rr,clipping=True,removeIns=True):
    '''
    given the read return the clipped read sequence as well as start and end position of the read in genomic coordinates. Optionally you can disable clipping
    readsequence and cigar string are assumed to be given on the forward strand. removeIns determines whether insertions with respect to the 
    reference should be removed from the read sequence
    '''
    #CIGAR format: [(X_1,l_1),...,(X_n,l_n)]
    #where X_i in {0:M,1:I,2:D,3:N,4:S,5:H,6:P}
    retseq = ""
    qualstr = ""
    i = 0
    start = rr.pos
    end = rr.pos
    cigarstr = rr.cigar
    t0 = cigarstr[0] #first tuple
    for t in cigarstr:
        if t[0] == 0: #M
            iEnd = i + t[1]
            retseq += rr.seq[i:iEnd]
            qualstr += rr.qual[i:iEnd]
            i = iEnd
            end += t[1]
        elif t[0] == 1: #I
            iEnd = i + t[1]
            if not removeIns:
                retseq += rr.seq[i:iEnd]
                qualstr += rr.qual[i:iEnd]
            i = iEnd
        elif t[0] in [2,3]: #D,N
            retseq += t[1] * "-"
            qualstr += t[1] * "!" #! corresponds to 0 quality
            end += t[1]
        elif t[0] == 4: #S
            iEnd = i + t[1]
            if not clipping:
                retseq += rr.seq[i:iEnd]
                qualstr += rr.qual[i:iEnd]
            i = iEnd #just move cursor and do not add the bases to the string if clipping
            end += t[1]
        #do nothing for H,P or others
    te = cigarstr[-1] #last tuple
    if rr.is_reverse:
        if t0[0] == 4 and clipping:
            end -= t0[1]
    else:
        if te[0] == 4 and clipping: #S
            end -= te[1]
    return (retseq,qualstr,start,end)


def printread(rr):
    '''
    helper function to print the key attributes of a pysam read object
    '''
    print "------------------------------"
    print rr.qname
    print "------------------------------"
    print "is_unmapped: " + str(rr.is_unmapped)
    print "pos: " + str(rr.pos)
    print "is_reverse: " + str(rr.is_reverse)
    print "rlen: " + str(rr.rlen)
    print "rname: " + str(rr.rname)
    print "seq: " + rr.seq
    print "cigar: " + str(rr.cigar)
    print "is_paired: " + str(rr.is_paired)
    print "is_proper_pair: " + str(rr.is_proper_pair)
    readpairnumstr = 0
    if rr.is_read1:
        readpairnumstr = 1
    elif rr.is_read2:
        readpairnumstr = 2
    print "read_in_pair: " + str(readpairnumstr)
    print "-------------"
    print "MATE:"
    print "mate_is_unmapped: " + str(rr.mate_is_unmapped)
    print "mpos: " + str(rr.mpos)
    print "mrnm: " + str(rr.mrnm)
    print "mate_is_reverse: " + str(rr.mate_is_reverse)
    print "isize: " + str(rr.isize)

                
def getSeqBitArrays(s,padding=0,padVal=False,addBisulfite=1,addH=0,alphabet=["C","G","A","T","-"],otherChar="N"):
        '''
        given a sequence s return a dicionary with keys correcponding to the letters of a given alphabet and bitarrays
        as indicator verctors for the occurrences of that letter in s
        padding determines the final bitarray length = sequence length + padding. padded positions will only contain 0s
        addBisulfite=1: the letter Y=[CT] will be added, addBisulfite=2: the letter R=[AG] will be added, addBisulfite=0: nothing will be added
        addH=1: the letter H=[ACT] will be added, addH=2: the letter D=[AGT] will be added, addH=0: nothing will be added
        '''
        #padding currently not used. leave it in. it might come in handy
        
        bitArrays = {}
        l = len(s)
        ll = l
        if padding > 0:
            ll = l + padding
        for b in alphabet + [otherChar]:
            bitArrays[b] = bitarray(ll)
            bitArrays[b].setall(False)
            if padVal:
                bitArrays[b][l:] = True
        for i in xrange(l):
            c = s[i]
            if c in alphabet:
                bitArrays[c][i] = True
            else:
                bitArrays[otherChar][i] = True
        if addBisulfite==1:
            bitArrays["Y"] = bitArrays["C"] | bitArrays["T"]
        elif addBisulfite==2:
            bitArrays["R"] = bitArrays["A"] | bitArrays["G"]
        if addH==1:
            bitArrays["H"] = bitArrays["A"] | bitArrays["C"] | bitArrays["T"]
        elif addH==2:
            bitArrays["D"] = bitArrays["A"] | bitArrays["G"] | bitArrays["T"]
        return bitArrays

#debug
def printSeqBitArrays(sba):
    for b in sba.keys():
        print b + str(sba[b])

def getBisulfiteMatchArray(refBitA,querBitA,rl=None):
    if not rl:
        rl = min(len(refBitA["C"]),len(querBitA["C"]))
    matchBitArray = bitarray(rl)
    matchBitArray.setall(False)
    #determine matches (bisulfite converted)
    matchBitArray[0:rl] = (refBitA["C"][0:rl] & querBitA["Y"][0:rl]) | \
                          (refBitA["T"][0:rl] & querBitA["T"][0:rl]) | \
                          (refBitA["A"][0:rl] & querBitA["A"][0:rl]) | \
                          (refBitA["G"][0:rl] & querBitA["G"][0:rl])
    return(matchBitArray)

def transformQualStr(s):
    '''
    converts a string of BAM-based quality scores into list of numeric PHRED scores
    '''
    try:
        tmp = map(lambda x: ASCII2PHRED[x],s)
    except KeyError:
        print "WARNING: could not convert ASCII to Phred score: '" + s + "' --> set to 0"
        tmp = len(s) * [0]
    return tmp

def calcItemRgb(read, seqMotif):
    '''
    assign RGB colors to methylation ratios
    '''
    itemRgb = [0,0,0]
    # go to http://genome-test.cse.ucsc.edu/admin/rgbItemExamples.html for RGB values
    if read.meanMeth > 900:
        itemRgb = [57,39,140]   # blue for methylated CpGs 
    elif 800 < read.meanMeth <= 900:
        itemRgb = [49,48,206]       
    elif 700 < read.meanMeth <= 800:
        itemRgb = [57,48,198] 
    elif 600 < read.meanMeth <= 700: 
        itemRgb = [80,80,80]      
    elif 500 < read.meanMeth <= 600:
        itemRgb = [112,112,112]
    elif 400 < read.meanMeth <= 500:
        itemRgb = [144,144,144]  # partially methylated CpGs in grey, with darker hues signifying higher methylation levels
    elif 300 < read.meanMeth <= 400:
        itemRgb = [160,160,160] 
    elif 200 < read.meanMeth <= 300:
        itemRgb = [192,192,192]
    elif 100 <= read.meanMeth <=200:
        itemRgb = [214,214,216]
    else: 
        itemRgb = [238,58,58]                      # red for unmethylated CpGs  
        if seqMotif == 'A':  itemRgb = [255,216,22]   # yellow for unmethylated CpAs
        if seqMotif == 'C':  itemRgb = [153,0,255]    # purple for unmethylated CpCs 
        if seqMotif == 'T':  itemRgb = [255,102,204]  # pink for unmethylated CpTs
        if seqMotif == 'H':  itemRgb = [102,102,0]    # dark green for unmethylated CpHs
        if seqMotif == 'HpH':  itemRgb = [45,198,214] # turquoise for unmethylated CpHpHs 
        if seqMotif == 'HpG':  itemRgb = [96,221,73]  # green for unmethylated CpHpGs 
    return itemRgb

def bedSpecificStartEnd(blockStartList, readstart):
    '''
    be sure that the bed file starts with a pattern block and ends with one. otherwise the file would not be genome browser readable
    concrete application: transform the start and end coordinates in bed fitting format (st it starts and ends with a CpG/pattern)
    '''
    # the CpG positions for the BED file of the first output file have to be represented differently, since BED blocks must span chromStart to chromEnd. (chromStart + chromStarts[last] + blockSizes[last]) must equal chromEnd.
    blockStartsGenomeBrowser = []
    blockStartsGenomeBrowserNew = []
    bedSpecificChromstart = int
    bedSpecificChromend = int
    for cpgStart in blockStartList:
        blockStartsGenomeBrowser.append(cpgStart - readstart)
    if blockStartsGenomeBrowser[0] == 0:
        bedSpecificChromstart = readstart 
        bedSpecificChromend = readstart + blockStartsGenomeBrowser[len(blockStartsGenomeBrowser)-1] + 2
    else:
        bedSpecificChromstart = readstart + blockStartsGenomeBrowser[0]
        bedSpecificChromend = readstart + blockStartsGenomeBrowser[len(blockStartsGenomeBrowser)-1] + 2
        for cpgStart in blockStartsGenomeBrowser:
            blockStartsGenomeBrowserNew.append(cpgStart - blockStartsGenomeBrowser[0])
        blockStartsGenomeBrowser = blockStartsGenomeBrowserNew 
    return [bedSpecificChromstart, bedSpecificChromend, blockStartsGenomeBrowser]

class bsAlignedReadMinimal(genomeLocation):
    '''
    a memory and runtime sparse representation of an aligned read
    '''
    def __init__(self,alignedRead,chrom,seqMotifLength=1):
        self.chrom = chrom
        (self.readSequence,qualstr,self.readstart,self.readend) = getGappedSeqRead(alignedRead,clipping=True,removeIns=True)
        self.readlength = self.readend - self.readstart
        
        #fragment length
        if options.alignmentTool == 'rrbsmap':
            if options.rrbsMode:
                self.fragmentLength = int(alignedRead.opt("ZL"))   #ZL: frangment length in BAM files, ZL is not set for wholegenome aligment
            else:
                self.fragmentLength = self.readlength + seqMotifLength#len(fragmentSequence)
                #fragments for paired end sequencing
                if options.pairedEnd and alignedRead.is_proper_pair:
                    self.fragmentLength = abs(alignedRead.isize) + 1
        
    #reads are sorted according to read coordinates, not chrom coordinates
    def __cmp__(self, other):
        #alphanumerical chromosome order        
        cmpVal = cmp(manualGenomeCorrection.get(self.chrom,self.chrom),manualGenomeCorrection.get(other.chrom,other.chrom))
        if cmpVal == 0:
            cmpVal = cmp(self.readstart,other.readstart)
            if cmpVal == 0:
                cmpVal = cmp(self.readend,other.readend)
        return cmpVal


class bsAlignedRead(bsAlignedReadMinimal):
    '''
    main class responsible for manipulating aligned reads
    '''
    def __init__(self,alignedRead,mappedStrand,lane,chrom,refSeqHandler,options,seqMotif,addH=False,ampStrandI=None,statsObj=None,startsInRegion2Process=True):
        '''
        parse all the read information from a pysam read object
        ''' 
        self.chrom = chrom
        #chromstart = start of FRAGMENT in the genome
        # if whole genome alignment the class tries to infer the fragment length from the pairedEnd information
        # if that fails or the information is not available the fragment is basically the read coordinates with the end coordinated shifted by the length of the pattern - 1 (ie 1 for CpG or 2 for CpHpG,...)
        self.chromstart = None
        #chromend == fragment end
        self.chromend = None
        #READ coordinates
        self.readstart = None
        self.readend = None
        self.readlength = None
        self.ampStrand = None #the strand of the amplificate
        self.mappedStrand = None #the strand which the read is mapped to by the alignment tool. != ampStrand for read#2 of paired end reads
        self.mismatches = None
        self.readSequence = None
        self.readSequence_raw = None #unmodified sequence from the input file, required for strand prediction
        
        self.readSequenceBitArrays = None #bitvectors for each base in the alphabet
        self.fragmentSequence = None
        self.fragmentSequenceBitArrays = None #bitvectors for each base in the alphabet (always regarded gapless)
        self.seqMotifLength = 0 # length of the dimer-1, e.g. 1 for CpG, CpA ,...; 2 for CpHpG,... IMPORTANT since the fragment sequence will be retrieved with this offset
        self.matchBitArray = None #bitvector for matching bases between read and reference for each position
        
        self.fragmentLength = None # store the size of the fragment the read is mapped to. 
        self.qualityScoresPerBaseInRead = None
        self.lane = None
        self.correctAlignedPercent = None
        self.isDiscarded = False
        self.strandAgreesWithPrediction = None #if the option is enabled allowing only reads that map to the strand that is confirmed by a prediction model
        self.isDuplicateFlag = False
        self.isDuplicatePosCount = False
        
        
        #Paired End support
        self.qname = alignedRead.qname #the identifier for the pair, assumed to occur exactly twice in a file: once for each mate
        self.isPaired = alignedRead.is_paired == 1
        self.isProperPair = alignedRead.is_proper_pair == 1 #can we safely assume that the two mates capture the start and the end of a fragment? Yes!
        
        #In the Picard Maq pipeline, a two reads are considered a proper pair if:
        #* neither read has the unmapped flag set
        #* both reads map to the same reference sequence, and that sequence is not "*"
        #* the reads do not map to the same strand
        #* for standard (non-jumping) libraries, the start position of the positive end is less than the start position of the negative end plus the read length
        self.mateFound = None
        self.pairNum = 0 #which read is it in the pair (1 or 2). 0 for unpaired or unknown
        #which read in the pair are you?
        #read 1 is assumed to be the forward read (the read next to the start adapter) and read 2 the reverse read (the read sequenced from the end adapter)
        if self.isPaired:
            if alignedRead.is_read1 == 1:
                self.pairNum = 1
            elif alignedRead.is_read2 == 1:
                self.pairNum = 2
        
        self.ZSfield = None
        self.account4beingRead2 = options.pairedEnd and self.isPaired #does one have to account for the reverse read in a pair
        if options.processZSField and self.isPaired:
            try:
                self.ZSfield = alignedRead.opt("ZS")
            except KeyError:
                print "WARNING: read " + str(alignedRead) +" does not contain ZS field"
            self.account4beingRead2 = self.account4beingRead2 and self.ZSfield[1] == "-"
        else:
            self.account4beingRead2 = self.account4beingRead2  and self.pairNum == 2
        
        #information about the mate. Can be set by the setMateInfo() method
        self.overlapWithMate = 0 #how many basepairs does the read overlap with its mate according to the alignment
        self.overlapWithMateStartI = None #which start and stop indices of the overlap with the mate (0-based)
        self.overlapWithMateEndI = None #0-based: first index not in overlap
        self.processInPair = True #will the read be processed? This will be set to false if isProperPair is false and the read is the one in the pair that is of lower matching quality       
        self.undigestedFragment = False #is the fragment end inbetween the paired end read
        
        self.excludedBases = None #a binary list storing for each position in the sequence whether this base should be included in the processing
        self.basesToAnalyze = 0 #the number of bases that are analyzed from the start of the read
        
        # the values of the following variables are computed during methylation calling
        self.methRatioInt = []
        self.internal_cpgs = 0
        

        self.blockStartList = []
        self.isNoCpgConForReadConRate = 0 #for Non-CpG (or other seqMotif) conversion: all converted Cs
        self.isNoCpgForReadConRate = 0 #for Non-CpG (or other seqMotif) conversion: all Cs
        self.conversionRate = 'NA'
        self.meanMeth = 0
        self.bedScore = 0
        self.correctAlignmentRateQ = [] # list that stores the correct alignment rate of the 4 quartiles
        self.medianQualityScoreCpg = []
        self.medianQualityScoreDiscardedCpg = [] 
        self.medianQualityScoreAllBases = 0    
        self.medianQualityScoreNonCpgQ1 = 0
        self.medianQualityScoreNonCpgQ2 = 0
        self.medianQualityScoreNonCpgQ3 = 0
        self.medianQualityScoreNonCpgQ4 = 0
        
        #the strand of the amplificate
        #read#1 is the forward read and should correspond to the mapping
        #the mapped strand of read#2 has to be inverted
        if ampStrandI:
            ampStrand = ampStrandI
        else:
            ampStrand = mappedStrand
        
                    
        #convert complex CIGAR-defined read into gapped and clipped read string
        (readSequence,qualstr,readstart,readend) = getGappedSeqRead(alignedRead,clipping=True,removeIns=True)
        qualityScoresPerBaseInRead = transformQualStr(qualstr)
        readSequence = readSequence.upper()
        self.readSequence_raw = alignedRead.seq.upper()
        

        if mappedStrand == "-": 
            # transform to reverse complement if mapped to the minus strand
            qualityScoresPerBaseInRead = qualityScoresPerBaseInRead[::-1] # reverse the list of the Phred base quality to match reverse complement
            if options.alignmentTool == 'rrbsmap':
                readSequence = strRevComp(readSequence)
#            #TODO: Check here if still correct for rrbsmap
#            elif options.alignmentTool in ['BSMAP','RRBSMAP']: # only for BSMAP and RRBSMAP, where complement (if on minus strand) is already given, but not yet reversed
#                readSequence = readSequence[::-1]
            else: raise Exception("Unsupported alignment tool: "+options.alignmentTool)

        self.seqMotifLength = len(seqMotif.replace('p',''))
        readLength = len(readSequence)            
                
        if options.rrbsMode:
            if options.alignmentTool == 'rrbsmap':
                chromstart = int(alignedRead.opt("ZP")) - 1    # since 1-based leftmost position in additional TAG fields #ZP: position of fragment in reference in BAM files
                fragmentLength = int(alignedRead.opt("ZL"))   #ZL: frangment length in BAM files, ZL is not set for wholegenome aligment
                chromend = chromstart + fragmentLength   
            else: raise Exception("Unsupported alignment tool: "+options.alignmentTool)
        else:
            if options.alignmentTool == 'rrbsmap':
                fragmentLength = readLength + self.seqMotifLength#len(fragmentSequence)
                #FM:?[20100316] for wholegenome alignment, fragment start == readstart ?
                chromstart = alignedRead.pos
                #fragments for paired end sequencing
                if options.pairedEnd and alignedRead.is_proper_pair:
                    fragmentLength = abs(alignedRead.isize) + 1 # +1 as the ISIZE is specified ads the DIFFERENCE between the 5' coord of the mate and the 5' of the current read. that would leave 1 base out
                    if mappedStrand == '-':
                        chromstart = alignedRead.mpos
                    else:
                        chromstart = alignedRead.pos
                chromend = chromstart + fragmentLength
            else: raise Exception("Unsupported alignment tool: "+options.alignmentTool)     
        
        if fragmentLength < readLength:
            if mappedStrand == '+':
                readend = readstart + fragmentLength
            elif mappedStrand == '-':
                readstart = readend - fragmentLength
            else:
                raise Exception("Unknown Strand: "+mappedStrand)
            readLength = fragmentLength
            readSequence = readSequence[0:readLength]
            qualityScoresPerBaseInRead = qualityScoresPerBaseInRead[0:readLength]
        
        # retrieve genomic sequence from reference file handler
        try:
            if mappedStrand == '+':
                fragmentSequence = refSeqHandler.getRefSeq(readstart-self.seqMotifLength, readend+self.seqMotifLength,reverse=False).upper()
            else: # transform to reverse complement if on minus strand
                fragmentSequence = refSeqHandler.getRefSeq(readstart-self.seqMotifLength, readend+self.seqMotifLength,reverse=True).upper()
        except Exception, ex:
            sys.stdout.flush()
            print 'WARNING: Could not retrieve sequence due to the following exception: '+str(type(ex))+": "+str(ex)    
            self.isDiscarded = True
            exceptionOccurred = True
            return
        
        self.readlength = readend - readstart
        self.readSequence = readSequence
        self.fragmentSequence = fragmentSequence
        #bitarray conversion
        bisulfiteMode = 1 #if paired end sequncing and the read is the reverse read, other bisulfite matching criteria are used (since we have G->A conversion instead of C->T)
        if self.account4beingRead2:
            bisulfiteMode = 2
        self.readSequenceBitArrays = getSeqBitArrays(readSequence,padding=0,addBisulfite=bisulfiteMode,addH=0,alphabet=["C","G","A","T","-"],otherChar="N")
        #remark: you can save the "Y" in the read and just work with the match array, but the code reads better this way.
        self.fragmentSequenceBitArrays = getSeqBitArrays(fragmentSequence,padding=0,addBisulfite=False,addH=bisulfiteMode,alphabet=["C","G","A","T"],otherChar="N")
        
        #discard very short reads
        if self.readlength < options.minReadLength:
            if not self.isDiscarded:
                self.isDiscarded = True
                if startsInRegion2Process: #only count the read towards the statistics if it starts in the processing window that is currently processed
                    statsObj.updateStatsVal(seqMotif,lane,options.isStrandSpecific,ampStrand,'discard_minReadLength')
            if debug:
                print "NOTICE: discarded short read"
                print alignedRead.qname + " on " + self.chrom + "[" + str(readstart) +","+ str(readend) + ")" + mappedStrand
                
        #just an assertion that everything has the correct format
        rl = len(self.readSequence)
        
        if rl > len(self.readSequenceBitArrays["C"]) or rl > (len(self.fragmentSequenceBitArrays["C"])-2*self.seqMotifLength):
            print "ERROR: readlen does not match fragment length"
            print alignedRead.qname + " on " + self.chrom + "[" + str(readstart) +","+ str(readend) + ")" + mappedStrand
            print self.readSequence_raw + " , pos: "+str(alignedRead.pos)+"-- raw"
            print self.readSequence + " " + str(rl) + "bp [" +  str(len(self.readSequenceBitArrays["C"])) + " bitarray] -- processed"
            print self.fragmentSequence + " " + str(len(self.fragmentSequence)) + "bp [" +  str(len(self.fragmentSequenceBitArrays["C"])) + " bitarray]"
            print "Reference handler object:"
            print "chrom :" + refSeqHandler.chromName + "[" + str(refSeqHandler.chromStart) + "," + str(refSeqHandler.chromEnd) +"]"
            print "seqlength: " + str(len(refSeqHandler.seq))
            self.isDiscarded = True
            exceptionOccurred = True
            return

        mismatches = self.__countBsMismatches__()
        maxMismatches = int(options.maxMismatches)
        self.basesToAnalyze = self.readlength
        if options.mismatchTruncate > 0:
            #which bases will be analyzed
            mmPos = self.matchBitArray[self.seqMotifLength:(self.seqMotifLength+rl)].search('0')
            if len(mmPos) <= DYNAMIC_MAX_MISMATCHES[rl-1]: #full length read fulfills the criterion (-1 as the array is 0 based but the readlength is 1 based)
                self.basesToAnalyze = rl
            else:
                #determine the maximal length of bases to analyze where a maximum of --maxMismatches (in percent or absolute value) occurred
                for i in range(len(mmPos)):
                    curPos = mmPos[i] #up to (= not including) read position mmPos[i] (0-based), exactly i mismatches occured
                    if DYNAMIC_MAX_MISMATCHES[curPos] <= i: 
                        self.basesToAnalyze = curPos #-1 not necassary: index is zero based and basesToAnalyze describes the length
                if self.basesToAnalyze < options.mismatchTruncate:
                    if not self.isDiscarded:
                        if startsInRegion2Process: #only count the read towards the statistics if it starts in the processing window that is currently processed
                            statsObj.updateStatsVal(seqMotif,lane,options.isStrandSpecific,ampStrand,'discard_mismatches')
                        self.isDiscarded = True
        else:
            #is the maximum mismatch criterion fulfilled. otherwise discard the read
            if options.maxMismatches > 0 and options.maxMismatches < 1:
                maxMismatches = int(round(options.maxMismatches*self.readlength))
            if mismatches > maxMismatches:
                if not self.isDiscarded:
                    if startsInRegion2Process: #only count the read towards the statistics if it starts in the processing window that is currently processed
                        statsObj.updateStatsVal(seqMotif,lane,options.isStrandSpecific,ampStrand,'discard_mismatches')
                    self.isDiscarded = True
        minMismatches = int(options.minMismatches)
        if options.minMismatches > 0 and options.minMismatches < 1:
            minMismatches = int(round(options.minMismatches*self.readlength))
        if mismatches < minMismatches:
            if not self.isDiscarded:
                if startsInRegion2Process: #only count the read towards the statistics if it starts in the processing window that is currently processed
                    statsObj.updateStatsVal(seqMotif,lane,options.isStrandSpecific,ampStrand,'discard_mismatches')
                self.isDiscarded = True
            
        #actually set the attributes
        self.chromstart = chromstart
        self.chromend = chromend
        self.readstart = readstart
        self.readend = readend
        self.mappedStrand = mappedStrand
        self.ampStrand = ampStrand
        self.mismatches = mismatches
        self.qualityScoresPerBaseInRead = qualityScoresPerBaseInRead
        self.lane = lane
        self.correctAlignedPercent = 100 - ((float(self.mismatches) / len(readSequence))*100)
        self.fragmentLength = fragmentLength
        self.excludedBases = bitarray(self.readlength)
        self.excludedBases.setall(False)
        
    
    def __countBsMismatches__(self):
        '''
        count the alignement mismatches and set a bitarray indicating the read sequence positions where the bases are matched to the reference
        '''
        rl = len(self.readSequence)
        
        #initialize the matcharray with the length of the padded read bitarrays (C chosen arbitrarily)
        self.matchBitArray = bitarray(len(self.fragmentSequenceBitArrays["C"]))
        self.matchBitArray.setall(False)
        #determine matches (bisulfite converted)
        if self.account4beingRead2:
            self.matchBitArray[self.seqMotifLength:(self.seqMotifLength+rl)] = (self.fragmentSequenceBitArrays["C"][self.seqMotifLength:(rl+self.seqMotifLength)] & self.readSequenceBitArrays["C"]) | \
                                       (self.fragmentSequenceBitArrays["T"][self.seqMotifLength:(rl+self.seqMotifLength)] & self.readSequenceBitArrays["T"]) | \
                                       (self.fragmentSequenceBitArrays["A"][self.seqMotifLength:(rl+self.seqMotifLength)] & self.readSequenceBitArrays["A"]) | \
                                       (self.fragmentSequenceBitArrays["G"][self.seqMotifLength:(rl+self.seqMotifLength)] & self.readSequenceBitArrays["R"])
        else:
            self.matchBitArray[self.seqMotifLength:(self.seqMotifLength+rl)] = (self.fragmentSequenceBitArrays["C"][self.seqMotifLength:(rl+self.seqMotifLength)] & self.readSequenceBitArrays["Y"]) | \
                                       (self.fragmentSequenceBitArrays["T"][self.seqMotifLength:(rl+self.seqMotifLength)] & self.readSequenceBitArrays["T"]) | \
                                       (self.fragmentSequenceBitArrays["A"][self.seqMotifLength:(rl+self.seqMotifLength)] & self.readSequenceBitArrays["A"]) | \
                                       (self.fragmentSequenceBitArrays["G"][self.seqMotifLength:(rl+self.seqMotifLength)] & self.readSequenceBitArrays["G"])

        #count the number of bits set to 0 in the matching bitarray
        return self.matchBitArray[self.seqMotifLength:(self.seqMotifLength+rl)].count(0)
    
    def __correctAlignmentRateInQuartiles__(self, quartileNo):
        '''
        calculate the correct alignment rate in 4 segments of the read
        '''
        lenReadSeq = self.readlength/4
        start = (lenReadSeq * quartileNo) - (lenReadSeq)
        end = start + lenReadSeq
        matchesInQNr = self.matchBitArray[(self.seqMotifLength+start):(self.seqMotifLength+end)].count()
        return ((float(matchesInQNr) / lenReadSeq)*100)

    def CpPs(self,p="G",p_compl="C"):
        '''
        returns a bitvector with 1s for the positions the pattern p occurs: A pattern is defined by C (potentially bisulfite converted) followed by the bases given in p. p_rc denotes the reverse complement of the pattern
        both the read sequence and the reference sequence have to fulfill the pattern. the bitarrays of the readsequence is assumed to padded by 0s, to match the length of the reference 
        bitarrays. Corrolary: last positions in a read do NOT fulfill the pattern even if the refernce sequence fulfills it
        '''
        if self.account4beingRead2:
            isC = self.readSequenceBitArrays["R"] & self.fragmentSequenceBitArrays["G"][self.seqMotifLength:(self.readlength+self.seqMotifLength)]
            ret = isC
            for i in range(len(p_compl)):
                #reference fits the pattern
                ret &= self.fragmentSequenceBitArrays[p_compl[i]][(self.seqMotifLength-i-1):(self.seqMotifLength-i-1+self.readlength)]
                #and the read fits the reference (after bisulfite conversion)
                ret &= self.matchBitArray[(self.seqMotifLength-i-1):(self.seqMotifLength-i-1+self.readlength)]
        else:
            isC = self.readSequenceBitArrays["Y"] & self.fragmentSequenceBitArrays["C"][self.seqMotifLength:(self.readlength+self.seqMotifLength)]
            ret = isC
            for i in range(len(p)):
                #reference fits the pattern
                ret &= self.fragmentSequenceBitArrays[p[i]][(self.seqMotifLength+i+1):(self.seqMotifLength+i+1+self.readlength)]
                #and the read fits the reference (after bisulfite conversion)
                ret &= self.matchBitArray[(self.seqMotifLength+i+1):(self.seqMotifLength+i+1+self.readlength)]
        return ret
    
    def nonCpPs(self,p="G",p_compl="C"):
        '''
        analogous to CpPs() but returns positions where a C occurs and the pattern is not matched but read and reference agree on the sequence
        returns a bitvector with 1s for the positions the pattern p does not occur but the there is a C and the read and reference sequences agree for the next len(pattern) bases: A pattern is defined by C (potentially bisulfite converted) followed by the bases given in p
        '''
        if self.account4beingRead2:
            isC = self.readSequenceBitArrays["R"] & self.fragmentSequenceBitArrays["G"][self.seqMotifLength:(self.readlength+self.seqMotifLength)]
            ret = bitarray(self.readlength)
            ret.setall(False)
            for i in range(len(p_compl)):
                #reference does not fit the pattern
                ret |= ~self.fragmentSequenceBitArrays[p_compl[i]][(self.seqMotifLength-i-1):(self.seqMotifLength-i-1+self.readlength)]
                #and the read fits the reference (after bisulfite conversion)
                ret &= self.matchBitArray[(self.seqMotifLength-i-1):(self.seqMotifLength-i-1+self.readlength)]
        else:
            isC = self.readSequenceBitArrays["Y"] & self.fragmentSequenceBitArrays["C"][self.seqMotifLength:(self.readlength+self.seqMotifLength)]
            ret = bitarray(self.readlength)
            ret.setall(False)
            for i in range(len(p)):
                #reference does not fit the pattern
                ret |= ~self.fragmentSequenceBitArrays[p[i]][(self.seqMotifLength+i+1):(self.seqMotifLength+i+1+self.readlength)]
                #and the read fits the reference (after bisulfite conversion)
                ret &= self.matchBitArray[(self.seqMotifLength+i+1):(self.seqMotifLength+i+1+self.readlength)]
        return ret & isC


    def getMateNum(self):
        '''
        retrun the number of the mate in the pair (1 or 2 or 0 for unpaired)
        '''
        if self.pairNum == 1:
            return 2
        elif self.pairNum == 2:
            return 1
        else:
            return 0
    
    #setters
    def setIsDuplicate(self,isDuplicateFlag,isDuplicatePosCount):
        self.isDuplicateFlag = isDuplicateFlag
        self.isDuplicatePosCount = isDuplicatePosCount
    def setMateFound(self,b):
        self.mateFound = b
    def setMethRatioInt(self,x):
        self.methRatioInt = x
    def addSize2BlockSizeList(self,x):
        self.blockSizeList.append(x)
    def setInternal_cpgs(self,x):
        self.internal_cpgs = x                     
    def setBlockStartList(self,bsl):
        self.blockStartList = bsl
    def setChromstart(self,x):
        self.chromstart = x
    def setChromend(self,x):
        self.chromend = x
    
    def setBEDspecificStartEnd(self):
        '''
        transform the start and end coordinates in bed fitting format (st it starts and ends with a CpG/pattern)
        '''
        newChrmStartEnd = bedSpecificStartEnd(self.blockStartList, self.readstart)
        self.chromstart = newChrmStartEnd[0]
        self.chromend   = newChrmStartEnd[1]
        self.blockStartList =newChrmStartEnd[2]
    
    def setMateInfo(self,mate):
        '''
        sets the information about the mate which is given as bsAlignedRead. takes care of assigning positions where read and mate overlap
        and setting the respective bases to not be analyzed if applicable (i.e. if readnumber in pair <> 1)
        also identifies whether for non WGBSS a restriction site is overlooked
        '''        
        #there are the following scenarios how the paired reads in a proper pair can lie:
        #(A) +++++>           (B)    ++>+++>  (C) ++++++>     (D)    +++>     (E) ++++++>++>
        #            <------      <--<--             <-----       <--<------         <---
        #
        #     +++++>
        #     <-----  is handled in (C)
        #REMARK: this only depends on the mappedstrand and is not influenced by being an amplicon of the genomic + or - strand or the readnumber (forward or reverse read).
        #   <-------
        # +++++>
        if self.mappedStrand == "+":
            if mate.readstart >= self.readend: #scenario (A)
                self.overlapWithMate = 0
            elif mate.readstart >= self.readstart:
                if mate.readend >= self.readend: #(C)
                    self.overlapWithMate = self.readend - mate.readstart
                    self.overlapWithMateStartI = self.readlength - self.overlapWithMate
                    self.overlapWithMateEndI = self.readlength
                    #mateOverlapStartI = mate.readend - self.readend #overlap indeices for the mate sequence
                    #mateOverlapEndI = mate.readlength
                else: #(E)
                    self.overlapWithMate = mate.readlength
                    self.overlapWithMateStartI = mate.readstart - self.readstart
                    self.overlapWithMateEndI = mate.readend - self.readstart
                    #mateOverlapStartI = 0
                    #mateOverlapEndI = mate.readlength
                    self.excludedBases[self.overlapWithMateEndI:] = True # truncate the processed positions as case (E) here can only mean that the sequencing ran into the adaper
            else:
                if mate.readend <= self.readend: #(B)
                    self.overlapWithMate = mate.readend - self.readstart
                    self.overlapWithMateStartI = 0
                    self.overlapWithMateEndI = self.overlapWithMate
                    #mateOverlapStartI = 0
                    #mateOverlapEndI =  self.overlapWithMate
                    self.excludedBases[self.overlapWithMateEndI:] = True# truncate the processed positions as case (B) here can only mean that the sequencing ran into the adaper
                else: #(D)
                    self.overlapWithMate = self.readlength
                    self.overlapWithMateStartI = 0
                    self.overlapWithMateEndI = self.readlength
                    #mateOverlapStartI = mate.readend -self.readend
                    #mateOverlapEndI = mate.readend - self.readstart
        else:
            if mate.readend <= self.readstart: #(A)
                self.overlapWithMate = 0
            elif mate.readstart <= self.readstart:
                if self.readend >= mate.readend: #(C)
                    self.overlapWithMate = mate.readend -self.readstart
                    self.overlapWithMateStartI = self.readlength - self.overlapWithMate
                    self.overlapWithMateEndI = self.readlength
                    #mateOverlapStartI = self.readstart - mate.readstart
                    #mateOverlapEndI = mate.readlength
                else: #(E)
                    self.overlapWithMate = self.readlength
                    self.overlapWithMateStartI = 0
                    self.overlapWithMateEndI = self.readlength
                    #mateOverlapStartI = self.readstart - mate.readstart
                    #mateOverlapEndI = self.readend - mate.readstart
            else:
                if self.readend <= mate.readend: #(B)
                    self.overlapWithMate = self.readend - mate.readstart
                    self.overlapWithMateStartI = 0
                    self.overlapWithMateEndI = self.overlapWithMate 
                    #mateOverlapStartI = 0
                    #mateOverlapEndI = self.overlapWithMate
                    self.excludedBases[self.overlapWithMateEndI:] = True # truncate the processed positions as case (B) here can only mean that the sequencing ran into the adaper
                else: # (D)
                    self.overlapWithMate = mate.readlength
                    self.overlapWithMateStartI = self.readend - mate.readend
                    self.overlapWithMateEndI = self.readend - mate.readstart
                    #mateOverlapStartI = 0
                    #mateOverlapEndI = mate.readlength
                    self.excludedBases[self.overlapWithMateEndI:] = True # truncate the processed positions as case (D) here can only mean that the sequencing ran into the adaper
        
        #if the current read number is 1 process the entire read
        #else do not process the overlap
        if self.overlapWithMate > 0 and self.pairNum <> 1:
            self.excludedBases[self.overlapWithMateStartI:self.overlapWithMateEndI] = True
            
        #check if there has been a digestion error
        if options.rrbsMode:
            if self.mappedStrand == "+":
                if self.readstart + self.fragmentLength < mate.readend:
                    self.undigestedFragment = True
            else:
                if mate.readstart + mate.fragmentLength < self.readend:
                    self.undigestedFragment = True
                    
    def setBluntEndRepairNeglectedBases(self,restrSite):
        '''
        set the bases not to be analyzed that correspont to blunt end fill up
        those sites do not contain any information, as they are filled in with methylated Cs
        (only makes sense in RRBS mode)
        '''
        if self.account4beingRead2:
            self.excludedBases[:restrSite.fillInsEnd] = True
            self.excludedBases[(self.fragmentLength-restrSite.fillInsStart):] = True #if the index is too large nothing happens
        else:
            self.excludedBases[:restrSite.fillInsStart] = True
            self.excludedBases[(self.fragmentLength-restrSite.fillInsEnd):] = True #if the index is too large nothing happens
            
    def setAgreesWithStrandPrediction(self,model):
        '''
        use a statistical model to classify the strand of the read and see whether it agrees with the strand assignment
        '''
        self.strandAgreesWithPrediction = model.modelAgreesWithStrand(self)
    
    def coordsFitFragmentPlacement(self):
        '''
        used to discard reads in which the sequencing start does not fit their respective fragment placement, i.e.
        check whether the aligned read start position fits the start of a fragment (genomic + strand) or the end coordinate fits the fragment's end (- strand)
        '''
        #TODO: check whether the fragment coordinates actually account for the fill-in-bases
        if self.mappedStrand == "+":
            return self.readstart == self.chromstart
        else:
            return self.readend == self.chromend
    
    #setter function used in methylationCalling()
    def addMedianQualityScoreCpgFromQualityScoresPerBaseInRead(self,index):
        '''
        add the base quality of a position to the list of Cpg qualities for later median calculation
        '''
        self.medianQualityScoreCpg.append(self.qualityScoresPerBaseInRead[index])
    def addMedianQualityScoreDiscardedCpgFromQualityScoresPerBaseInRead(self,index):
        '''
        like addMedianQualityScoreCpgFromQualityScoresPerBaseInRead but just for discarded CpG positions
        '''
        self.medianQualityScoreDiscardedCpg.append(int(self.qualityScoresPerBaseInRead[index])) 
        
    def updateNoCpgInfo(self,isNoCpg,isNoCpgConv):
        '''
        update the counts for the C that do not fall into the pattern but still agree with the reference
        '''
        if isNoCpg == 1:
            self.isNoCpgForReadConRate += 1        
        if isNoCpgConv == 1: 
            self.isNoCpgConForReadConRate += 1
    def updateNoCpgInfoX(self,x_isNoCpg,x_isNoCpgConv):
        '''
        update the counts for the C that do not fall into the pattern but still agree with the reference.
        important for conversion rate calculation
        '''
        self.isNoCpgForReadConRate += x_isNoCpg   
        self.isNoCpgConForReadConRate += x_isNoCpgConv

    def setConversionRate4nonCpG(self):
        '''
        set the reads conversion rate: calculated from the Cs that are out of CpG context
        '''
        if self.isNoCpgForReadConRate > 0:
            self.conversionRate = int(round((float(self.isNoCpgConForReadConRate)/self.isNoCpgForReadConRate)*100)) #FM:? Why int()?
    
    def discardMe(self):
        '''
        set the discard flag
        '''
        self.isDiscarded = True    
    
    def calcStats(self):
        '''
        calculate single read statistics: quality scores (in segments of the read), correct alignemtn rate,...
        '''
        if len(self.methRatioInt) > 0: self.meanMeth = int(round(float(sum(self.methRatioInt)) / len(self.methRatioInt) * 1000.0))
        if options.readDetailsFile:
            for i in range(4):
                self.correctAlignmentRateQ.append(self.__correctAlignmentRateInQuartiles__(i+1))    
        # from here on SORTED qualityScoresPerBaseInRead --> do NOT USE them anymore, only for the following statistics
        quarter = len(self.qualityScoresPerBaseInRead)/4
        q1 = self.qualityScoresPerBaseInRead[0:quarter]
        q1.sort()
        self.medianQualityScoreNonCpgQ1 = median(q1)
        q2 = self.qualityScoresPerBaseInRead[quarter:(quarter*2)]
        q2.sort()
        self.medianQualityScoreNonCpgQ2 = median(q2)
        q3 = self.qualityScoresPerBaseInRead[(quarter*2):(quarter*3)]
        q3.sort()
        self.medianQualityScoreNonCpgQ3 = median(q3)
        q4 = self.qualityScoresPerBaseInRead[(quarter*3):len(self.qualityScoresPerBaseInRead)]
        q4.sort()
        self.medianQualityScoreNonCpgQ4 = median(q4)          
        self.qualityScoresPerBaseInRead.sort() #FM maybe do a deepcopy here or is this too slow and we really not use this one anymore? NJ: would leave this as is, should not hurt
        self.medianQualityScoreAllBases = median(self.qualityScoresPerBaseInRead)      
        self.medianQualityScoreCpg.sort() #FM? maybe do a deepcopy here or is this too slow and we really not use this one anymore?
        self.medianQualityScoreCpg = median(self.medianQualityScoreCpg) # e.g. CpG quality score >20
        self.medianQualityScoreDiscardedCpg.sort()
        self.medianQualityScoreDiscardedCpg = median(self.medianQualityScoreDiscardedCpg) # e.g. CpG quality score <20
    
    def getDebugString(self):
        s = self.qname + "\n"
        s += "Coords (read): " + self.chrom + "[" + str(self.readstart) + ":" + str(self.readend) + ") ampStrand:" + self.ampStrand+ " mappedStrand" + self.mappedStrand + " length:"+str(self.readlength)+"\n"
        s += "Coords (frag): " + self.chrom + "[" + str(self.chromstart) + ":" + str(self.chromend) + ") length:" + str(self.fragmentLength)+ "\n"
        s += "__Seq__\n"
        s += "read:      " + self.readSequence + "\n"
        s += "read(raw): " + self.readSequence_raw + "\n"
        s += "qual     : " + str(self.qualityScoresPerBaseInRead) + "\n"
        s += "matches  : " + str(self.matchBitArray) + " (" + str(self.mismatches)+" mismatches)\n"
        s += "frag     :"  + self.fragmentSequence + '\n'
        s += "isDiscarded: " + str(self.isDiscarded) + "\n"
        s += "isDuplicate(Flag,PositionCount): " + str(self.isDuplicateFlag) + "," + str(self.isDuplicatePosCount) + "\n"
        s += "basesToAnalyze: " + str(self.basesToAnalyze) +"\n"
        s += "excludedBases: " + str(self.excludedBases) +"\n"
        s += "__paired end__\n"
        s += "isPaired: " + str(self.isPaired) + "\n"
        s += "isProperPair: " + str(self.isProperPair) + "\n"
        s += "pairNum: " + str(self.pairNum) + "\n"
        s += "account4beingRead2: " + str(self.account4beingRead2) + "\n"
        s += "__Patterns__\n"
        s += "internal_cpgs: " + str(self.internal_cpgs) + "\n"
        return s
        
    def __str__(self):
        return "%s\t%i\t%i\t%i\t%i\t%s\t%i\t%s\t%s\t%s\t%i\n" % (manualGenomeCorrection.get(self.chrom), self.chromstart, self.chromend, self.readstart, self.readend, self.mappedStrand, self.mismatches, self.readSequence, self.fragmentSequence[self.seqMotifLength:(self.readlength+self.seqMotifLength)], str(self.qualityScoresPerBaseInRead).strip('[]').translate(string.maketrans("", ""), string.whitespace), self.lane)


class bsAlignedReadCollection:
    '''
    container class for bsAlignedRead
    '''
    def __init__(self):
        self.reads = []
    def addRead(self,read):
        self.reads.append(read)
    
    def getNumReads(self):
        return len(self.reads)
    
    def sort(self):
        self.reads.sort()
        
    def updateChromStartEnds(self):
        '''
        for every read, transform the start and end coordinates in bed fitting format (st it starts and ends with a CpG/pattern)
        '''
        for read in self.reads:
            newChrmStartEnd = bedSpecificStartEnd(read.blockStartList, read.readstart)
            read.setChromstart(newChrmStartEnd[0])
            read.setChromend(newChrmStartEnd[1])
            read.setBlockStartList(newChrmStartEnd[2])
#        self.sort() # it's possible that .bed file is not sorted correctly anymore after the .bed format specific chromstart and chromend have been determined
    
    def writeReadOutFile(self,seqMotif,readOutfileObj):
        '''
        write all the reads to file
        '''
        dimLen = len(seqMotif.replace('p',''))+1
        dimLenStr = str(dimLen)+','
        
        for read in self.reads:
            blockSizeListString = (read.internal_cpgs * dimLenStr)[:-1] #:-1: dont include the last ','; before: str(read.blockSizeList).strip('[]').translate(string.maketrans("", ""), string.whitespace)
            itemRgb = calcItemRgb(read, seqMotif)
            readOutfileObj.write(manualGenomeCorrection.get(read.chrom,read.chrom) + "\t" +str(read.chromstart) + "\t"  + str(read.chromend) + "\t" +  str(read.meanMeth / 10.0)+"%" + "\t" + str(read.meanMeth) + "\t" + read.mappedStrand+ "\t" + \
                                 str(read.chromstart) + "\t" + str(read.chromend)+ "\t" + str(itemRgb).strip('[]').translate(string.maketrans("", ""), string.whitespace) + "\t" + str(read.internal_cpgs) + "\t" + blockSizeListString + "\t" + str(read.blockStartList).strip('[]').translate(string.maketrans("", ""), string.whitespace) + "\t" + str(read.pairNum) + "\n")

    def writeReadDetailsFile(self,readQualOutfileObj):
        '''
        write all the reads to file in extended format
        '''
        for read in self.reads:
            readQualOutfileObj.write(manualGenomeCorrection.get(read.chrom,read.chrom) + "\t" +str(read.chromstart) + "\t"  + str(read.chromend) + "\t" + listAsString(read.methRatioInt) + "\t" + str(read.meanMeth) + "\t" + read.mappedStrand+ "\t" + str(read.readstart) + "\t" + str(read.readend)+ "\t" + str(read.internal_cpgs) + "\t" + listAsString(read.blockStartList) + "\t" + str(read.lane) + "\t" + str(read.mismatches) + "\t" + str(read.correctAlignedPercent) +"\t" + \
                                     listAsString(read.correctAlignmentRateQ)  + "\t" +  str(read.medianQualityScoreAllBases)+ "\t" + str(read.medianQualityScoreNonCpgQ1) + \
                                     "," +str(read.medianQualityScoreNonCpgQ2) + "," + str(read.medianQualityScoreNonCpgQ3) + "," + str(read.medianQualityScoreNonCpgQ4) + "\t" + str(read.medianQualityScoreCpg)+ "\t" +  str(read.medianQualityScoreDiscardedCpg) + "\t" + str(read.conversionRate) + "\t" + str(read.isNoCpgConForReadConRate) + '/' + str(read.isNoCpgForReadConRate) + "\t" + str(read.basesToAnalyze) + "\t" + str(read.pairNum) + "\n") # str(read.maxSeedWithout_mismatches) + "\t" + str(read.mismatches_completelyMethylated) + "\t" + str(read.maxSeedWithout_mismatches_completelyMethylated) + "\t" + str(read.mismatches_completelyUnmethylated)  + "\t" + str(read.maxSeedWithout_mismatches_completelyUnmethylated)
            
            
class ReadOutputFileHandler:
    '''
    handles the output of reads to an outputfile
    after adding reads up to the buffer size, it writes them to the output file
    '''
    def __init__(self,fileName, bufferSize, seqMotif, detailsFile=None):
        self.fileName = fileName
        self.detailsFileName = detailsFile
        self.seqMotif = seqMotif
        self.f = open(fileName,"w")
        self.df = None
        if detailsFile:
            self.df = open(detailsFile,"w")
        self.numReads = 0
        self.curReads = 0
        self.bufferSize = bufferSize
        self.readCollection = bsAlignedReadCollection()
        
    def __del__(self):
        self.f.close()
        if self.detailsFileName:
            self.df.close()
        
    def sort(self):
        '''
        sort the BED file according to chrom and start coordinate
        sorting closes the file!
        the return value should be 0 for successful sorting
        '''
        self.f.close()
        unsortedFile = self.fileName + ".unsorted"
        os.system("mv " + self.fileName + " " + unsortedFile)
        sortCmd = "sort -k1,1 -k2,2n " + unsortedFile + " > " + self.fileName
        os.system(sortCmd)
        os.system("rm " + unsortedFile)
        
        if self.detailsFileName:
            self.df.close()
            unsortedFile = self.detailsFileName + ".unsorted"
            os.system("mv " + self.detailsFileName + " " + unsortedFile)
            sortCmd = "sort -k1,1 -k2,2n " + unsortedFile + " > " + self.detailsFileName
            os.system(sortCmd)
            os.system("rm " + unsortedFile)
            return(isValidFile(self.fileName) and isValidFile(self.detailsFileName))
        
        return(isValidFile(self.fileName))
        
    def flush(self):
        if self.detailsFileName:
            self.readCollection.writeReadDetailsFile(self.df)

        self.readCollection.updateChromStartEnds()
        self.readCollection.writeReadOutFile(self.seqMotif, self.f)
        self.readCollection = bsAlignedReadCollection()
        self.curReads = 0
            
    def addRead(self,read):  
        self.readCollection.addRead(read)
        self.numReads += 1
        self.curReads += 1
        #flush if the read buffer is full
        if self.curReads >= self.bufferSize:
            self.flush()


class UniqueCpgsPerLane(genomeLocation):
    '''
    handles a single CpG position and its methylation ratio
    used to generate the list of all CpGs for a region in sequenceDataProcessing()
    '''  
    def __init__(self, fragmentID,strand, cpgMethylated, cpgTotal, seqMotif):        
        fragmentID_details = fragmentID.split("/")
        self.chrom = fragmentID_details[0]            
        self.chromstart = int(fragmentID_details[1])
        self.chromend = self.chromstart + len(seqMotif.split("p")) + 1
        self.cpgMethylated = cpgMethylated
        self.cpgTotal = cpgTotal
        self.meanMethVal = float(cpgMethylated) / float(cpgTotal)
        self.meanMethBEDscore = int(round(self.meanMethVal*1000))#bed score is in the range 0-1000 and resemble the color intensity for the genome browser. here used to represent mthylation ratio
        if seqMotif == 'G': self.bedScore = self.meanMethBEDscore
        else: 
            if self.meanMethBEDscore == 0: self.bedScore = 0
            elif (self.meanMethBEDscore > 0) and (self.meanMethBEDscore <= 100): self.bedScore = 500 #FM: what does this do? what is a bedScore? #NJ bed format specific score in range 0-1000 #FM:? but what does it indicate?
            else: 
                if int(500 + (self.meanMethBEDscore /2.0)) <= 1000: self.bedScore = int(500 + (self.meanMethBEDscore /2.0))
                else: self.bedScore = 1000
        self.strand = strand

    def __str__(self):
        return "%s\t%i\t%i\t'%i/%i'\t%i\t%s\n" % (manualGenomeCorrection.get(self.chrom,self.chrom), self.chromstart, self.chromend, self.cpgMethylated, self.cpgTotal, self.bedScore, self.strand) 


class UniqueFragmentsPerLane(genomeLocation):
    '''
    handles a genomic fragment (as determined by digestion or paired end sequencing of fragmented genome)
    used to generate the list of all CpGs for a region in sequenceDataProcessing()
    '''       
    def __init__(self, fragmentID, methList, strand):    
        maxCoverage = 100
        self.fragmentID = fragmentID
        fragmentID_details = fragmentID.split(":")
        self.chrom = fragmentID_details[0]            
        self.chromstart = int(fragmentID_details[1])
        self.chromend = int(fragmentID_details[2])
        for methRate in range(len(methList)):
            methList[methRate] = int(round(methList[methRate] / 10))
        self.methList = methList
        coverageScore = ((math.log(len(methList))/math.log(2)) / (math.log(maxCoverage)/math.log(2))) * 1000
        if coverageScore > 1000: coverageScore = 1000
        self.fragmentCoverage = coverageScore
        self.strand = strand
        
    def __add__(self,other):
        '''
        for adding 2 fragments resembling the same coordinates.
        if the strand is + on either instance, the new strand will be + 
        '''
        ss = self.strand
        newMethList = self.methList+other.methList
        newMethList10 = map(lambda x: x *10,newMethList) #since the constructor of UniqueFragmentsPerLane expects the bed specific ratios: 100% corresponds to 1000
        #set the strand to + if either one of the strands is +
        if ss == "-" and other.strand == "+":
            ss = "+"
        return UniqueFragmentsPerLane(self.fragmentID,newMethList10,ss)
    
    def __str__(self):
        return "%s\t%i\t%i\t'%i:%s'\t%i\t%s\n" % (manualGenomeCorrection.get(self.chrom,self.chrom), self.chromstart, self.chromend, len(self.methList), str(self.methList).strip('[]').translate(string.maketrans("", ""), string.whitespace), self.fragmentCoverage, self.strand)
        
class IncludedRegions:
    '''
    handles the genomic regions that will be processed in the script
    in essence, sliding window that overlap with the user specified included regions are used for analysis
    '''
    def __init__(self,s):
        '''
        parse from a string with the format #[CHROM]:[START]-[STOP]#[CHROM]:[START]-[STOP]#...
        '''
        self.regions = {}
        if s == "":
            return
        try:
            ss = s.strip().split("#")
            for elem in ss:
                chromSplit = elem.strip().split(":")
                chrom = chromSplit[0]
                posSplit = chromSplit[1].strip().split("-")
                start = int(posSplit[0])
                stop = int(posSplit[1])
                try:
                    self.regions[chrom].append((start,stop))
                except KeyError:
                    self.regions[chrom] = [(start,stop)]
        except IndexError:
            self.regions = {}
            print "WARNING: could not parse includedRegion string: '" + s + "' --> resetted to '' (all regions)"
            
    def overlapsWithRegion(self,chrom,start,end):
        '''
        does the given genomic interval overlap with any of the stored regions
        '''
        if len(self.regions) == 0:
            return True
        try:
            chromList = self.regions[chrom]
        except KeyError:
            return False
        for region in chromList:
            if (start >= region[0] and start < region[1]) or (end > region[0] and end <= region[1]) or (start <= region[0] and end > region[1]):
                return True
        return False
    
class FragmentsInfo:
    '''
    stores read counts for the given fragment sizes
    important for statistics
    '''
    def __init__(self,fragmentSizeRangeStr):
        self.coverage4fragmentSizes = {}
        self.fragmentSizeRanges = fragmentSizeRangeStr.split(",")
        self.ranges = [] #lower and upper bounds for ranges each element is assumed to be a 2 tuple t where t[0] contains the lower and t[1] the upper bound
        self.ranges4header = ""
        
        #prepare lower and upper bound lists
        for s in self.fragmentSizeRanges:
            try:
                ss = s.split("-")                 
                self.ranges.append([int(ss[0]),int(ss[1])])
            except:
                raise Exception("Invalid Fragment Range: " + s )
        self.ranges.sort()
            
        # prepare header for global statistics output file ???
        for size in self.fragmentSizeRanges:
            self.ranges4header +='fragmentSizeRange'+ size.replace("-","to") +'bases_readCount\t'
            self.ranges4header +='fragmentSizeRange'+ size.replace("-","to") +'bases_fragmentCount\t'

class Statistics:
    '''
    main container for storing the statistics
    (a) lane, strand and region wise
    (b) a global statistics object
    see also StatisticsCollection which is the corresponding container class
    
    see the constructor/str() methods and the extensive description in StatisticsCollection.writeStatistics2files() for details on the stored information
    '''            
    def __init__(self, SampleName, inputFileName, strand, seqMotif, fragmentInfo):       
        self.SampleName = SampleName
        self.inputFileName = inputFileName
        self.strand = strand     
        self.seqMotifName = seqMotif
        self.cytosineStatsVal = {'total_c':{'base':0, 'firstBase':0}, 'total_c_conv':{'base':0, 'firstBase':0}, 'total_cpg':{'base':0, 'firstBase':0}, 'total_cpg_conv':{'base':0, 'firstBase':0}, 'total_c_nocpg':{'base':0, 'firstBase':0}, 'total_c_nocpg_conv':{'base':0, 'firstBase':0} }
        self.statsVal = {'mappingError':0, 'totalMappedReads':0, 'totalPairedReads':0, 'totalProperPairedReads':0, 'readsPassedQc':0, 'total_unique_cpg':0, 'processedDuplicatesFlag':0,'processedDuplicatesPosCount':0,'pairedReads':0,'properPairedReads':0,'mateNotFoundInRegion':0,'overlappingPairedReads':0,'pairedReadsUndigestedFragments':0,'discard_disagreeWithStrandPrediction':0,'discard_PFstatus':0,'discard_minReadLength':0,'discard_mismatches':0,'discard_restrSite':0,'discard_improperPairs':0,'discard_noSeqMotif':0,'discard_duplicatesFlag':0,'discard_duplicatesPosCount':0,'discard_noSeqMotif':0,'discard_fragmentsOutsideSizeRange':0,'mappingOrderDoesNotMatchStrand':0}
        #add statistics values that depend on MULTIPLE specific thresholds
        for t in SEQMOTIF_THRES:
            self.statsVal["readsWithAtLeast"+str(t)+"seqMotifs"] = 0
        for t in READ_THRES_POSITIONS:
            self.statsVal["coordsWithAtLeast"+str(t)+"reads"] = 0
        #non-additive statistics where the value for both strands cant be calculated from the sum of both strands
        self.statsValStrandSpecific = {}
        for t in READ_THRES_WINDOWS:
            self.statsValStrandSpecific["windowsWithAtLeast"+str(t)+"reads"] = 0

        self.allReadsMethSum = 0
        self.allReadsMethSumOfSquares = 0
        self.allReadsMethCount = 0
        
        self.uniqueCpgMethSum = 0
        self.uniqueCpgMethSumOfSquares = 0
        self.uniqueCpgMethCount = 0
        self.medianLists = {'correctAlignmentRate':[], 'conversionRate':[], 'medianQualityScoreAllBases':[], 'medianQualityScoreCpg':[], 'medianQualityScoreDiscardedCpg':[]}
        self.Q1 = {}; self.Median = {}; self.Q3 = {}
        
        #set the fragment information
        self.fragmentRanges = fragmentInfo.ranges
        #prepare dictionary for read counts
        self.coverage4fragmentSizes = {}
        for size in fragmentInfo.fragmentSizeRanges: 
            self.coverage4fragmentSizes[size] = 0  #readCount in fragment size range
        
        # attributes set by calculateStats
        self.acceptanceRate = None
        self.meanMethTotal = None
        self.meanMethAtStart = None
        self.conversionRate = None
        self.adjustedMeanMeth = None
        self.meanMethAllUniqueCpgs = None
        self.meanMethAllReads = None
        
        # attributes set by calcMeanSD()
        self.allReadsMethMean = None
        self.allReadsMethStdev = None
        self.uniqueCpgMethMean = None
        self.uniqueCpgMethStdev = None
    
    def printDebugAll(self):
        print "SampleName: " + str(self.SampleName)
        print "inputFileName: " + str(self.inputFileName)
        print "strand: " + str(self.strand)
        print "seqMotifName: " + str(self.seqMotifName)
        print "cytosineStatsVal: " + str(self.cytosineStatsVal)
        print "statsVal: " + str(self.statsVal)
        print "allReadsMethSum: " + str(self.allReadsMethSum)
        print "allReadsMethSumOfSquares: " + str(self.allReadsMethSumOfSquares)
        print "allReadsMethCount: " + str(self.allReadsMethCount)
        print "uniqueCpgMethSum: " + str(self.uniqueCpgMethSum)
        print "uniqueCpgMethSumOfSquares: " + str(self.uniqueCpgMethSumOfSquares)
        print "uniqueCpgMethCount: " + str(self.uniqueCpgMethCount)
        print "medianLists: " + str(self.medianLists)
        print "Q1: " + str(self.Q1)
        print "Median: " + str(self.Median)
        print "Q3: " + str(self.Q3)
        print "fragmentRanges: " + str(self.fragmentRanges) 
        print "coverage4fragmentSizes: " + str(self.coverage4fragmentSizes)
        print "acceptanceRate: " + str(self.acceptanceRate)
        print "meanMethTotal: " + str(self.meanMethTotal)
        print "meanMethAtStart: " + str(self.meanMethAtStart)
        print "conversionRate: " + str(self.conversionRate)
        print "adjustedMeanMeth: " + str(self.adjustedMeanMeth)
        print "meanMethAllUniqueCpgs: " + str(self.meanMethAllUniqueCpgs)
        print "meanMethAllReads: " + str(self.meanMethAllReads)
        print "allReadsMethMean: " + str(self.allReadsMethMean)
        print "allReadsMethStdev: " + str(self.allReadsMethStdev)
        print "uniqueCpgMethMean: " + str(self.uniqueCpgMethMean)
        print "uniqueCpgMethStdev: " + str(self.uniqueCpgMethStdev)
        
    # update the main (process) Statistic object with the stats values from the genome fraction's Statistic object             
    def updateWith(self, genomeFracStats):
        '''
        take the values of another Statistics objects and update your own values accordingly
        corresponds to addition of 2 Statistics objects
        '''
        self.allReadsMethSum += genomeFracStats.allReadsMethSum
        self.allReadsMethSumOfSquares += genomeFracStats.allReadsMethSumOfSquares
        self.allReadsMethCount += genomeFracStats.allReadsMethCount

        self.uniqueCpgMethSum += genomeFracStats.uniqueCpgMethSum
        self.uniqueCpgMethSumOfSquares += genomeFracStats.uniqueCpgMethSumOfSquares
        self.uniqueCpgMethCount += genomeFracStats.uniqueCpgMethCount
        
        for cStatus in genomeFracStats.cytosineStatsVal.keys():
            for pos in genomeFracStats.cytosineStatsVal[cStatus].keys():
                self.cytosineStatsVal[cStatus][pos] += genomeFracStats.cytosineStatsVal[cStatus][pos]
                        
        for val in genomeFracStats.statsVal.keys():
            self.statsVal[val] += genomeFracStats.statsVal[val]
        for val in genomeFracStats.statsValStrandSpecific.keys():
            self.statsValStrandSpecific[val] += genomeFracStats.statsValStrandSpecific[val]
        
        for medianList in genomeFracStats.medianLists.keys():
            for quartile in [genomeFracStats.Q1[medianList], genomeFracStats.Median[medianList], genomeFracStats.Q3[medianList]]:
                self.medianLists[medianList].append(quartile)
        for kk in genomeFracStats.coverage4fragmentSizes.keys():
            self.coverage4fragmentSizes[kk] += genomeFracStats.coverage4fragmentSizes[kk]
    
    #setters
    def setConversionRate(self,cr):
        self.conversionRate = cr
    def addOneStatsVal(self,statsValKey):
        self.statsVal[statsValKey] += 1
    def setStatsVal(self,statsValKey,x,isStrandSpecificStat=False):
        if isStrandSpecificStat:
            self.statsValStrandSpecific[statsValKey] = x
        else:
            self.statsVal[statsValKey] = x
    def addOneCytosineStatsVal(self,seqMotifStatus,posVal='base'):
        self.cytosineStatsVal[seqMotifStatus][posVal] += 1
    def addXCytosineStatsVal(self,x,seqMotifStatus,posVal='base'):
        self.cytosineStatsVal[seqMotifStatus][posVal] += x 
    def setUniqueCpgMethSum(self,x):
        self.uniqueCpgMethSum = x
    def setUniqueCpgMethSumOfSquares(self,x):
        self.uniqueCpgMethSumOfSquares = x
    def setUniqueCpgMethCount(self,x):
        self.uniqueCpgMethCount = x
    
    def updatePlusStrandSpecific(self, statsBoth, statsMinus):
        '''
        given statistics for the minus strand a combined statistics object representing both strands,
        calculate the plus strand specific properties.
        only applicable for a plus strand statistics object
        '''
        # most statistic values for the plus strand are obtained by simply subtracting the values for minus strand from the values for both strands 
        for statsValue in self.statsVal.keys():
            self.statsVal[statsValue] += statsBoth.statsVal[statsValue] - statsMinus.statsVal[statsValue]
        for methStatsValue in self.cytosineStatsVal.keys():
            for basePos in self.cytosineStatsVal[methStatsValue].keys():
                self.cytosineStatsVal[methStatsValue][basePos] += statsBoth.cytosineStatsVal[methStatsValue][basePos] - statsMinus.cytosineStatsVal[methStatsValue][basePos]
    
    def updateUniqueCpgMeth(self,methVal):
        '''
        update count, methylation status and sum of squeres for methylation (for later standard deviation calculation)
        for a CpG postion (unique position in the genome)
        '''
        self.uniqueCpgMethSum += methVal
        self.uniqueCpgMethSumOfSquares += pow(methVal, 2)
        self.uniqueCpgMethCount += 1
    
    def __getFragmentWindowKeys4Length__(self,l):
        '''
        for a length l returns a list of corresponding fragment windows
        '''
        retlist =[]
        for i in range(len(self.fragmentRanges)):
            if l >= self.fragmentRanges[i][0] and l < self.fragmentRanges[i][1]:
                retlist.append(str(self.fragmentRanges[i][0]) + "-" + str(self.fragmentRanges[i][1]))
        return retlist
        
    def __addReadFragmentInfo__(self,read):
        '''
        determine the window of fragment sizes corresponding to the fragment the read maps to
        and increment the read counts
        '''
        fragmentWindowKeys = self.__getFragmentWindowKeys4Length__(read.fragmentLength)
        for kk in fragmentWindowKeys:
            self.coverage4fragmentSizes[kk] += 1 
                
    def calcReadStats(self,read,options):
        '''
        update for a given processed read
        needed for calcReadStats method of StatisticsCollection class
        '''
        # for statistics output file: for each of these lists the Median, Q1 and Q3 are computed in the Statistics object
        if read.correctAlignedPercent != 'NA': self.medianLists['correctAlignmentRate'].append(read.correctAlignedPercent)
        if read.conversionRate != 'NA': self.medianLists['conversionRate'].append(read.conversionRate)
        if read.medianQualityScoreAllBases != 'NA': self.medianLists['medianQualityScoreAllBases'].append(read.medianQualityScoreAllBases)
        if read.medianQualityScoreCpg != 'NA': self.medianLists['medianQualityScoreCpg'].append(read.medianQualityScoreCpg)
        if read.medianQualityScoreDiscardedCpg != 'NA': self.medianLists['medianQualityScoreDiscardedCpg'].append(read.medianQualityScoreDiscardedCpg)  
        self.allReadsMethSum += read.meanMeth/10.0 #/10 since bedscores run to 1000
        self.allReadsMethSumOfSquares += pow(read.meanMeth/10.0,2)
        self.allReadsMethCount += 1

        self.__addReadFragmentInfo__(read)
    
    def calcQ1Q3Median(self,debuginfo=""):
        '''
        calculate quartile and median statistics for statistics indicated by the keys in self.medianLists
        Note: empties the lists after computation
        '''
        for listM in self.medianLists.keys():
            lenList = len(self.medianLists[listM])
            if lenList != 0:
                self.medianLists[listM].sort()
                
                self.Q1[listM] = self.medianLists[listM][int(lenList*0.25)]
                self.Median[listM] = self.medianLists[listM][int(lenList*0.5)]
                self.Q3[listM] = self.medianLists[listM][int(lenList*0.75)]
            else:
                if debug:
                    print "Could not calculate Q1, median and Q3 for statistics [" + listM + "]" + " " + debuginfo
                self.Q1[listM] = "NA"; self.Median[listM] = "NA"; self.Q3[listM] = "NA"
                
            self.medianLists[listM] = [] # delete lists to save memory (and not pickle it), after extracting median and quartiles
                
    def calculateStats(self,debuginfo=""):
        '''
        calculated multiple statistics values from its own attributes
        final lane statistics calculation in wrapup task
        '''
        self.acceptanceRate = 0
        if self.statsVal['totalMappedReads'] > 0: self.acceptanceRate = (float(self.statsVal['readsPassedQc']) / self.statsVal['totalMappedReads'])
        self.meanMethTotal = 0
        self.meanMethAtStart = 0
        self.conversionRate = 'NA'
        if self.cytosineStatsVal['total_cpg']['base'] > 0: #there are 2 cases: stats for the First Base Of Read (usually the first CpG) names ['firstBase'], and just ['base'] means all bases, including first base
            self.meanMethTotal = (float(self.cytosineStatsVal['total_cpg']['base'] - self.cytosineStatsVal['total_cpg_conv']['base'])/self.cytosineStatsVal['total_cpg']['base'] *100)        
        if self.cytosineStatsVal['total_cpg']['firstBase'] > 0: self.meanMethAtStart = (float(self.cytosineStatsVal['total_cpg']['firstBase'] - self.cytosineStatsVal['total_cpg_conv']['firstBase'])/self.cytosineStatsVal['total_cpg']['firstBase'] *100)        
        if self.cytosineStatsVal['total_c_nocpg']['base'] > 0:
            print str(float(self.cytosineStatsVal['total_c_nocpg_conv']['base'])) + " / " + str(self.cytosineStatsVal['total_c_nocpg']['base'])
            self.conversionRate = (float(self.cytosineStatsVal['total_c_nocpg_conv']['base'])/self.cytosineStatsVal['total_c_nocpg']['base']) *100            
        if self.conversionRate != 'NA': self.adjustedMeanMeth = self.meanMethTotal * (self.conversionRate/100.0)
        else: self.adjustedMeanMeth = 'NA'
        if self.uniqueCpgMethCount > 0: self.meanMethAllUniqueCpgs = self.uniqueCpgMethSum / self.uniqueCpgMethCount
        else: self.meanMethAllUniqueCpgs = 0            
        if self.allReadsMethCount > 0: self.meanMethAllReads = self.allReadsMethSum / self.allReadsMethCount 
        else: self.meanMethAllReads = 0

        self.calcQ1Q3Median(debuginfo=debuginfo)

    def calcMeanSD(self):
        '''
        compute methylation means and standard deviation from the accumulated values for each CpG
        final lane statistics calculation in wrapup task
        '''
        if self.allReadsMethCount > 0:
            self.allReadsMethMean = self.allReadsMethSum / self.allReadsMethCount
            self.allReadsMethStdev = pow((self.allReadsMethSumOfSquares - pow(self.allReadsMethMean,2)*self.allReadsMethCount) / (self.allReadsMethCount - 1), 0.5) # cf. http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
        else: self.allReadsMethStdev = "NA"
        if self.uniqueCpgMethCount > 0:
            self.uniqueCpgMethMean = self.uniqueCpgMethSum / self.uniqueCpgMethCount
            self.uniqueCpgMethStdev = pow((self.uniqueCpgMethSumOfSquares - self.uniqueCpgMethCount*pow(self.uniqueCpgMethMean,2)) / (self.uniqueCpgMethCount - 1), 0.5)
        else:self.uniqueCpgMethStdev = "NA"
        
    def __str__(self):
        outputStringFragCov = ''
        sortListOfSizeRanges = []
        for sizeRange in self.coverage4fragmentSizes.keys():
            sortListOfSizeRanges.append((int(sizeRange.split('-')[0]), self.coverage4fragmentSizes[sizeRange]))
        sortListOfSizeRanges.sort()
        for item in sortListOfSizeRanges:
            outputStringFragCov += str(item[1]) + "\t"  + "0" +"\t"
        outputStringFragCov = outputStringFragCov[:-1]
            
        #statistics that depend on multiple thresholds
        multipleThresStatsStr = ''
        for t in SEQMOTIF_THRES:
            multipleThresStatsStr += str(self.statsVal["readsWithAtLeast"+str(t)+"seqMotifs"]) + "\t"
        for t in READ_THRES_POSITIONS:
            multipleThresStatsStr += str(self.statsVal["coordsWithAtLeast"+str(t)+"reads"]) + "\t"
        for t in READ_THRES_WINDOWS:
            multipleThresStatsStr += str(self.statsValStrandSpecific["windowsWithAtLeast"+str(t)+"reads"]) + "\t"
        multipleThresStatsStr = multipleThresStatsStr[:-1] #cut off trailing '\t'
            
        # paired end total number and accepted number of read = number of mate pairs investigated (not number of single end reads)      
        totalReadNum = self.statsVal['totalMappedReads']
        acceptedReadNum = self.statsVal['readsPassedQc']
        acceptanceRate = self.acceptanceRate 
        if options.pairedEnd:
            totalSingleReads = float(self.statsVal['totalMappedReads'] - self.statsVal['totalPairedReads'])
            totalReadNum = float(self.statsVal['totalPairedReads']) / 2.0 + totalSingleReads
            acceptedSingleReads = float(self.statsVal['readsPassedQc'] - self.statsVal['pairedReads'])
            acceptedReadNum = float(self.statsVal['pairedReads']) / 2.0 +  acceptedSingleReads
            acceptanceRate = acceptedReadNum/totalReadNum
        
        return str(self.SampleName + "\t" + self.inputFileName + "\t" + 'Cp' +self.seqMotifName +"\t" + str(self.strand)+ "\t" + str(totalReadNum) +"\t" + str(acceptedReadNum)+ "\t" + str(acceptanceRate) + "\t" \
                   + str(self.statsVal['total_unique_cpg'])+ "\t" + str(self.conversionRate) + "\t" +str(self.meanMethTotal) + "\t" + str(self.adjustedMeanMeth) + "\t" + str(self.meanMethAtStart) + "\t" + str(self.meanMethAllUniqueCpgs) + "\t" + \
                   str(self.uniqueCpgMethStdev) + "\t" + str(self.meanMethAllReads)  + "\t" + str(self.allReadsMethStdev)+ "\t" + str(self.cytosineStatsVal['total_c_conv']['base']) + "\t" + str(self.cytosineStatsVal['total_c']['base']) + "\t" + str(self.cytosineStatsVal['total_cpg_conv']['base']) + "\t" + \
                   str(self.cytosineStatsVal['total_cpg']['base']) + "\t" + str(self.cytosineStatsVal['total_c_nocpg_conv']['base']) + "\t" + str(self.cytosineStatsVal['total_c_nocpg']['base']) + "\t" + str(self.Median['correctAlignmentRate'])  + "\t" + str(self.Q1['correctAlignmentRate']) + "\t" + str(self.Q3['correctAlignmentRate']) + "\t" + \
                    str(self.Median['conversionRate']) + "\t" + str(self.Q1['conversionRate']) + "\t" + str(self.Q3['conversionRate']) + "\t" + str(self.Q1['medianQualityScoreAllBases']) + "\t" + str(self.Median['medianQualityScoreAllBases']) + "\t" + str(self.Q3['medianQualityScoreAllBases']) + "\t" + \
                    str(self.Q1['medianQualityScoreCpg']) + "\t" + str(self.Median['medianQualityScoreCpg']) + "\t" + str(self.Q3['medianQualityScoreCpg']) + "\t" + str(self.Q1['medianQualityScoreDiscardedCpg']) + "\t" + str(self.Median['medianQualityScoreDiscardedCpg']) + "\t" + \
                    str(self.Q3['medianQualityScoreDiscardedCpg'])  + "\t" + str(self.statsVal['totalMappedReads']) + "\t" + str(self.statsVal['readsPassedQc']) + "\t" + str(self.statsVal['totalPairedReads'])  + "\t" + str(self.statsVal['pairedReads']) + "\t" + str(self.statsVal['totalProperPairedReads']) + "\t" + str(self.statsVal['properPairedReads'])  + "\t" + str(self.statsVal['mateNotFoundInRegion']) + "\t" + str(self.statsVal['overlappingPairedReads']) + "\t" + str(self.statsVal['pairedReadsUndigestedFragments']) + "\t" + \
                    str(self.statsVal['processedDuplicatesFlag'])+'\t'+str(self.statsVal['processedDuplicatesPosCount'])+'\t'+\
                    str(self.statsVal['discard_PFstatus'])+'\t'+str(self.statsVal['discard_minReadLength'])+'\t'+str(self.statsVal['discard_mismatches'])+'\t'+str(self.statsVal['discard_restrSite'])+'\t'+str(self.statsVal['discard_fragmentsOutsideSizeRange']) + '\t' + str(self.statsVal['discard_disagreeWithStrandPrediction']) +'\t'+str(self.statsVal['discard_improperPairs']) + "\t" + \
                    str(self.statsVal['discard_duplicatesFlag']) +'\t'+str(self.statsVal['discard_duplicatesPosCount']) +'\t' + str(self.statsVal['discard_noSeqMotif'])+'\t'+ str(self.statsVal['mappingOrderDoesNotMatchStrand']) +'\t'+ \
                    str(self.statsVal['mappingError']) + "\t" + outputStringFragCov + "\t" + multipleThresStatsStr + "\n")


class StatisticsCollection:
    '''
    stores Statistics objects for lanes and a global statistic
    indexed by seqMotif and strand
    provides the framework for other functions to update lane statistics and compute global statistics
    Class instances are generated for regions (sequenceDataProcessing)
    or on a global scale in performAnalysis()
    '''
    def __init__(self,filelist,seqMotiflist,strands,fragmentInfo):
        self.laneStats = {}
        #helpful attributes, not really necassary, but nice to have
        self.numberOfLanes = 0
        self.seqMotifs = seqMotiflist
        self.strands = strands
        #combines laneStats into a global statistic
        self.globalStats = {}
        self.fragmentInfo = fragmentInfo # an object of type FragmentInfo to compute fragment size coverages, etc.
        self.liftOverMappingErrorDict = {} #dictionary containing the number of mapping errors (liftOver) for each seqMotif and strand
        self.regionStats = {} #statistics applying to a genomic region accross all lanes, dictionary for seqMotifs and strands
        self.spikeInStats = {}
        self.spikeInStatsHeaders = []
        for seqMotif in self.seqMotifs:
            self.regionStats[seqMotif] = {}
            for s in self.strands:
                self.regionStats[seqMotif][s] = {"uniqueseqMotifCount": 0,"uniqueseqMotifMethSum":0.0, "uniqueseqMotifMethSumSquares":0.0}
                #multiple threshold based statistics
                for t in READ_THRES_POSITIONS:
                    self.regionStats[seqMotif][s]["coordsWithAtLeast"+str(t)+"reads"] = 0

        # prepare lane-specific statistics
        for seqMotif in self.seqMotifs:
            statsLaneSpecificTable = [] # list that stores all lane-specific statistics objects
            for inputFile in filelist: 
                statsLaneSpecific = {}  # dictionary that stores the lane-specific statistics objects. statsLaneSpecific[0]= statsLaneSpecific for both strands, statsLaneSpecific[1]= statsLaneSpecific for - strand, statsLaneSpecific[2]= statsLaneSpecific for + strand
                for strand in self.strands:
                    statsLaneSpecific[strand] = Statistics(options.sampleName, inputFile, strand, seqMotif, self.fragmentInfo)       
                statsLaneSpecificTable.append(statsLaneSpecific)           
            self.laneStats[seqMotif] = statsLaneSpecificTable
            self.numberOfLanes = len(filelist)
        
        for seqMotif in self.seqMotifs:
            self.liftOverMappingErrorDict[seqMotif] = {}
            for strand in self.strands:
                self.liftOverMappingErrorDict[seqMotif][strand] = 0
        for s in self.strands:
            self.spikeInStats[s] = {}
                
    #interface for the region wide statistics over all lanes
    #used for all the regionwise StatisticsCollection objects (in the region pickle files)
    def addOneRegionStats(self,seqMotif,strand,k):
        self.regionStats[seqMotif][strand][k] += 1
    def setRegionStats(self,seqMotif,strand,statKey,x):
        self.regionStats[seqMotif][strand][statKey] = x
    def updateRegionStats4process(self,seqMotif,strand,d): #d should be a dictionary containing a subset of keys as self.regionStats
        for k in d.keys():
            self.regionStats[seqMotif][strand][k] += d[k]
            
    def updateStatsVal(self,seqMotif,lanekey,isStrandSpecific,strand,statsValKey):
        '''
        basically an update wrapper:
        for a given lane and strand, add 1 to a specific statsVal in that StatisticsObjects
        also set the minus strand if the strand specific option is enabled
        called from sequenceDataProcessing() : statsValKey: 'totalMappedReads','duplicates'
        called from calcReadStatistics(): 'readsPassedQc'
        '''
        self.laneStats[seqMotif][lanekey]["both"].addOneStatsVal(statsValKey)
        # the + strand is not treated beacuse the values can be derived by subtracting the - values from the both strand values (at the end when calculating stats)
        if isStrandSpecific and strand == '-': self.laneStats[seqMotif][lanekey]["-"].addOneStatsVal(statsValKey)
    
    def setLaneStatsVal(self,seqMotif,lane,strand,statsValKey,x,isStrandSpecificStat=False):
        '''
        basically a setter wrapper:
        for a given lane and strand, set a specific statsVal in that StatisticsObjects
        '''
        self.laneStats[seqMotif][lane][strand].setStatsVal(statsValKey,x,isStrandSpecificStat=isStrandSpecificStat)
    
    def updateCytosineStatsVal_posVal(self,seqMotif,lane,seqMotifStatusList,isStrandSpecific,strand,posVal='base'):
        '''
        basically a setter wrapper:
        for a given lane and strand, update a specific cytosineStatsVal in that StatisticsObjects by 1
        called from methylationCalling()
        '''
        for seqMotifStatus in seqMotifStatusList:
            self.laneStats[seqMotif][lane]["both"].addOneCytosineStatsVal(seqMotifStatus,posVal)
            if isStrandSpecific and strand == '-': self.laneStats[seqMotif][lane]["-"].addOneCytosineStatsVal(seqMotifStatus,posVal)
    
    def updateCytosineStatsVal(self,seqMotif,lane,seqMotifStatusList,valuesToUpdate,isStrandSpecific,strand):
        '''
        wrapper for updateCytosineStatsVal_posVal to update multiple values
        called from methylationCalling()
        '''
        for posVal in valuesToUpdate:
            self.updateCytosineStatsVal_posVal(seqMotif,lane,seqMotifStatusList,isStrandSpecific,strand,posVal)
            
    def updateXCytosineStatsVal_seqMotifStatus_posVal(self,x,seqMotif,lane,seqMotifStatus,isStrandSpecific,strand,posVal='base'):
        '''
        basically a setter wrapper:
        for a given lane and strand, update a specific cytosineStatsVal in that StatisticsObjects by x
        called from methylationCalling()
        '''
        self.laneStats[seqMotif][lane]["both"].addXCytosineStatsVal(x,seqMotifStatus,posVal)
        if isStrandSpecific and strand == '-': self.laneStats[seqMotif][lane]["-"].addXCytosineStatsVal(x,seqMotifStatus,posVal)
            
    def calcReadStats(self,read,seqMotif,options):
        '''
        feed the read to the corresponding lanes object to be able to calculate its statistics
        basically median lists and methylation means
        needed for calcReadStatistics() function
        ''' 
        self.laneStats[seqMotif][read.lane]["both"].calcReadStats(read,options)
        if options.isStrandSpecific:
            self.laneStats[seqMotif][read.lane][read.ampStrand].calcReadStats(read,options)
    
    def updateUniqueCpgMeth(self,seqMotif,lane,strand,methVal):
        '''
        for a discovered CpG update the corresponding lane statistics
        is called from performAnalysis()
        '''    
        self.laneStats[seqMotif][lane][strand].updateUniqueCpgMeth(methVal)                                           
        #self.laneStats[seqMotif][lane][strand].meanMethAllUniqueCpgsList.append(methVal)# for meanMethAllUniqueCpgs in statistics file
        self.laneStats[seqMotif][lane][strand].addOneStatsVal('total_unique_cpg')
    
    #TODO: calculate Q1, Q3 and median correctly and not as medians of medians
    def setQ1Q3Median(self,seqMotif,debuginfo=""):
        '''
        calculate the quariles and median from the values in the median lists
        '''
        for lane in range(self.numberOfLanes):
            for strand in self.strands:
                self.laneStats[seqMotif][lane][strand].calcQ1Q3Median(debuginfo=debuginfo)
    
    #interface for the global object summarizing the region statisticsCollection objects
    def updateProcessStats(self,seqMotif,procStatsCol):
        '''
        take a StatisticsCollection object from a processed region and update the corresponding Statistics object in the global collection
        is called from performAnalysis_wrapup()
        '''
        if self.numberOfLanes <> procStatsCol.numberOfLanes or seqMotif not in self.seqMotifs or seqMotif not in procStatsCol.seqMotifs or self.strands <> procStatsCol.strands:
            raise Exception("Incompatible StatisticsCollections")
        for lane in range(self.numberOfLanes):
            for strand in self.strands:    
                self.laneStats[seqMotif][lane][strand].updateWith(procStatsCol.laneStats[seqMotif][lane][strand])
        for s in self.strands:
            for k in self.regionStats[seqMotif][s].keys():
                self.regionStats[seqMotif][s][k] += procStatsCol.regionStats[seqMotif][s][k]
                
    def updateStrandSpecific(self,seqMotif):
        '''
        for each lane compute the + strand specific statistics from the combined object and the minus strand object
        is called from performAnalysis()
        '''
        for lane in range(self.numberOfLanes):
            self.laneStats[seqMotif][lane]["+"].updatePlusStrandSpecific(self.laneStats[seqMotif][lane]["both"],self.laneStats[seqMotif][lane]["-"])            
                
    def setConversionRate(self,seqMotif,lane,strand):
        try:  # take BS conversionRate only for CpGs                  
            self.laneStats[seqMotif][lane][strand].setConversionRate(self.laneStats['G'][lane][strand].conversionRate) 
        except: # CpG was not calculated in this analysis
            self.laneStats[seqMotif][lane][strand].setConversionRate('NA')
                
    def calcStats_setConversionRate(self,seqMotif):
        '''
        a multifuntion that does all the things in its name. used in performAnalysis
        final lane statistics calculation in wrapup task
        '''
        for lane in range(self.numberOfLanes):
            for s in self.strands:
                self.laneStats[seqMotif][lane][s].calculateStats()
                self.setConversionRate(seqMotif,lane,s)
                
    def computeGlobalStats(self):
        '''
        compute the global statistics from lane statistics and store the result in self.globalStats
        Note some values are just simple additions, some are more complex (eg. the unique CpGs)
        also set non region specific statistics (as liftover mapping errors)
        '''
        
        # GLOBAL STATISTICS over all lanes of sample, if more than one lane was given as input file
        for seqMotif in self.seqMotifs:
            #Initialize
            #for global statistics object over all lanes (= input files)
            statsGlobal = {}  #dictionary that stores the global statistics objects ## list that stores the global statistics objects. statsGlobal[0]= statsGlobal for both strands, statsGlobal[1]= statsGlobal for - strand, statsGlobal[2]= statsGlobal for + strand
            for strand in self.strands:
                statsGlobal[strand] = Statistics(options.sampleName, 'all lanes (input files)', strand, seqMotif, self.fragmentInfo)        
            self.globalStats[seqMotif] = statsGlobal
            
            if self.numberOfLanes > 1:
                for lane in range(self.numberOfLanes):
                    for s in self.strands:
                        #update all statistics
                        #self.globalStats[seqMotif][s].updateWithStatsVal(self.laneStats[seqMotif][lane][s])
                        self.globalStats[seqMotif][s].updateWith(self.laneStats[seqMotif][lane][s])
            else:
                for s in self.strands:
                    self.globalStats[seqMotif][s] = self.laneStats[seqMotif][0][s]
        #calculate statistics
        for seqMotif in self.seqMotifs: 
            for s in self.strands:
                if self.numberOfLanes>1:        
                    self.globalStats[seqMotif][s].calculateStats()
        
        #set the LiftOver mapping errors
        for seqMotif in self.seqMotifs: 
            for s in self.strands:
                self.globalStats[seqMotif][s].setStatsVal('mappingError',self.liftOverMappingErrorDict[seqMotif][s])
            
        #set the uniqueseqMotif, mean and sd values from the region stats
        for seqMotif in self.seqMotifs: 
            for s in self.strands:
                self.globalStats[seqMotif][s].setUniqueCpgMethSum(self.regionStats[seqMotif][s]["uniqueseqMotifMethSum"])
                self.globalStats[seqMotif][s].setUniqueCpgMethSumOfSquares(self.regionStats[seqMotif][s]["uniqueseqMotifMethSumSquares"])
                self.globalStats[seqMotif][s].setUniqueCpgMethCount(self.regionStats[seqMotif][s]["uniqueseqMotifCount"])
                self.globalStats[seqMotif][s].setStatsVal('total_unique_cpg',self.regionStats[seqMotif][s]["uniqueseqMotifCount"])

 
    def calcMeanSDs(self):
        '''
        calculate the mean and standard deviation methylation for every lane, strand and pattern
        '''
        for seqMotif in self.seqMotifs: 
            for s in self.strands:
                self.globalStats[seqMotif][s].calcMeanSD()
                for lane in range(self.numberOfLanes):
                    self.laneStats[seqMotif][lane][s].calcMeanSD()
                   
    def updateLiftOverMappingErrors(self,unmappable,seqMotif,outfile,options):
        '''
        updating unmappable after liftover counts 
        '''
        try:
            if options.isStrandSpecific:
                if outfile.find('OutfileM') >= 0:
                    self.liftOverMappingErrorDict[seqMotif]["-"] += unmappable # for - strand
                elif outfile.find('OutfileP') >= 0:
                    self.liftOverMappingErrorDict[seqMotif]["+"] += unmappable # for + strand
            else:
                self.liftOverMappingErrorDict[seqMotif]["both"] += unmappable
        except KeyError:
            print "WARNING: unable to update LiftOver mapping error due to KeyError (strand)"
            exceptionOccurred = True
    
    def readSpikeInCtrlAnalysis(self,sic):
        '''
        add the information contained in a SpikeInControls object
        '''
        processPstrand = "+" in self.strands
        processMstrand = "-" in self.strands
        processBothStrands = "both" in self.strands
        statsHeader_minQual = []
        for motif in sic.motifs.keys():
            sumMotifUnmethConvP = 0
            sumMotifUnmethConvM = 0
            sumMotifUnmethAllP = 0
            sumMotifUnmethAllM = 0
            sumMotifMethConvP = 0
            sumMotifMethConvM = 0
            sumMotifMethAllP = 0
            sumMotifMethAllM = 0
            sumMotifUnmethConvPminQual = 0
            sumMotifUnmethConvMminQual = 0
            sumMotifUnmethAllPminQual = 0
            sumMotifUnmethAllMminQual = 0
            sumMotifMethConvPminQual = 0
            sumMotifMethConvMminQual = 0
            sumMotifMethAllPminQual = 0
            sumMotifMethAllMminQual = 0
            for (cid,cPos) in  sic.motifs[motif]:
                if sic.controls[cid]["methStatus"] == "methylated":
                    if sic.controls[cid]["strand"] == "-":
                        sumMotifMethConvM += sic.controls[cid]["Ccounts"][cPos]["conv"]
                        sumMotifMethConvMminQual += sic.controls[cid]["CcountsQual"][cPos]["conv"]
                        sumMotifMethAllM += sic.controls[cid]["Ccounts"][cPos]["total"]
                        sumMotifMethAllMminQual += sic.controls[cid]["CcountsQual"][cPos]["total"]
                    else:
                        sumMotifMethConvP += sic.controls[cid]["Ccounts"][cPos]["conv"]
                        sumMotifMethConvPminQual += sic.controls[cid]["CcountsQual"][cPos]["conv"]
                        sumMotifMethAllP += sic.controls[cid]["Ccounts"][cPos]["total"]
                        sumMotifMethAllPminQual += sic.controls[cid]["CcountsQual"][cPos]["total"]
                elif sic.controls[cid]["methStatus"] == "unmethylated":
                    if sic.controls[cid]["strand"] == "-":
                        if motif=="CTT" and sic.controls[cid]["Ccounts"][cPos]["total"] > 0:
                            print cid + " pos" + str(cPos) + " counted " + str(sic.controls[cid]["Ccounts"][cPos]["total"])
                        sumMotifUnmethConvM += sic.controls[cid]["Ccounts"][cPos]["conv"]
                        sumMotifUnmethConvMminQual += sic.controls[cid]["CcountsQual"][cPos]["conv"]
                        sumMotifUnmethAllM += sic.controls[cid]["Ccounts"][cPos]["total"]
                        sumMotifUnmethAllMminQual += sic.controls[cid]["CcountsQual"][cPos]["total"]
                    else:
                        sumMotifUnmethConvP += sic.controls[cid]["Ccounts"][cPos]["conv"]
                        sumMotifUnmethConvPminQual += sic.controls[cid]["CcountsQual"][cPos]["conv"]
                        sumMotifUnmethAllP += sic.controls[cid]["Ccounts"][cPos]["total"]
                        sumMotifUnmethAllPminQual += sic.controls[cid]["CcountsQual"][cPos]["total"]                       
            
            if processPstrand:
                self.spikeInStats["+"]["unmethCtrl_"+motif+"_converted"] = sumMotifUnmethConvP
                self.spikeInStats["+"]["unmethCtrl_"+motif+"_all"] = sumMotifUnmethAllP
                if sumMotifUnmethAllP > 0:
                    self.spikeInStats["+"]["unmethCtrl_"+motif+"_convRate"] = float(sumMotifUnmethConvP)/float(sumMotifUnmethAllP)
                else:
                    self.spikeInStats["+"]["unmethCtrl_"+motif+"_convRate"] = None
                self.spikeInStats["+"]["methCtrl_"+motif+"_converted"] = sumMotifMethConvP
                self.spikeInStats["+"]["methCtrl_"+motif+"_all"] = sumMotifMethAllP
                if sumMotifMethAllP > 0:
                    self.spikeInStats["+"]["methCtrl_"+motif+"_convRate"] = float(sumMotifMethConvP)/float(sumMotifMethAllP)
                else:
                    self.spikeInStats["+"]["methCtrl_"+motif+"_convRate"] = None
                self.spikeInStats["+"]["unmethCtrl_"+motif+"_converted_minQual"] = sumMotifUnmethConvPminQual
                self.spikeInStats["+"]["unmethCtrl_"+motif+"_all_minQual"] = sumMotifUnmethAllPminQual
                if sumMotifUnmethAllPminQual > 0:
                    self.spikeInStats["+"]["unmethCtrl_"+motif+"_convRate_minQual"] = float(sumMotifUnmethConvPminQual)/float(sumMotifUnmethAllPminQual)
                else:
                    self.spikeInStats["+"]["unmethCtrl_"+motif+"_convRate_minQual"] = None
                self.spikeInStats["+"]["methCtrl_"+motif+"_converted_minQual"] = sumMotifMethConvPminQual
                self.spikeInStats["+"]["methCtrl_"+motif+"_all_minQual"] = sumMotifMethAllPminQual
                if sumMotifMethAllPminQual > 0:
                    self.spikeInStats["+"]["methCtrl_"+motif+"_convRate_minQual"] = float(sumMotifMethConvPminQual)/float(sumMotifMethAllPminQual)
                else:
                    self.spikeInStats["+"]["methCtrl_"+motif+"_convRate_minQual"] = None
            if processMstrand:
                self.spikeInStats["-"]["unmethCtrl_"+motif+"_converted"] = sumMotifUnmethConvM
                self.spikeInStats["-"]["unmethCtrl_"+motif+"_all"] = sumMotifUnmethAllM
                if sumMotifUnmethAllM > 0:
                    self.spikeInStats["-"]["unmethCtrl_"+motif+"_convRate"] = float(sumMotifUnmethConvM)/float(sumMotifUnmethAllM)
                else:
                    self.spikeInStats["-"]["unmethCtrl_"+motif+"_convRate"] = None
                self.spikeInStats["-"]["methCtrl_"+motif+"_converted"] = sumMotifMethConvM
                self.spikeInStats["-"]["methCtrl_"+motif+"_all"] = sumMotifMethAllM
                if sumMotifMethAllM > 0:
                    self.spikeInStats["-"]["methCtrl_"+motif+"_convRate"] = float(sumMotifMethConvM)/float(sumMotifMethAllM)
                else: 
                    self.spikeInStats["-"]["methCtrl_"+motif+"_convRate"] = None
                self.spikeInStats["-"]["unmethCtrl_"+motif+"_converted_minQual"] = sumMotifUnmethConvMminQual
                self.spikeInStats["-"]["unmethCtrl_"+motif+"_all_minQual"] = sumMotifUnmethAllMminQual
                if sumMotifUnmethAllMminQual > 0:
                    self.spikeInStats["-"]["unmethCtrl_"+motif+"_convRate_minQual"] = float(sumMotifUnmethConvMminQual)/float(sumMotifUnmethAllMminQual)
                else:
                    self.spikeInStats["-"]["unmethCtrl_"+motif+"_convRate_minQual"] = None
                self.spikeInStats["-"]["methCtrl_"+motif+"_converted_minQual"] = sumMotifMethConvMminQual
                self.spikeInStats["-"]["methCtrl_"+motif+"_all_minQual"] = sumMotifMethAllMminQual
                if sumMotifMethAllMminQual > 0:
                    self.spikeInStats["-"]["methCtrl_"+motif+"_convRate_minQual"] = float(sumMotifMethConvMminQual)/float(sumMotifMethAllMminQual)
                else:
                    self.spikeInStats["-"]["methCtrl_"+motif+"_convRate_minQual"] = None
            if processBothStrands:
                self.spikeInStats["both"]["unmethCtrl_"+motif+"_converted"] = sumMotifUnmethConvP + sumMotifUnmethConvM
                self.spikeInStats["both"]["unmethCtrl_"+motif+"_all"] = sumMotifUnmethAllP + sumMotifUnmethAllM
                if sumMotifUnmethAllP + sumMotifUnmethAllM > 0:
                    self.spikeInStats["both"]["unmethCtrl_"+motif+"_convRate"] = float(sumMotifUnmethConvP + sumMotifUnmethConvM)/float(sumMotifUnmethAllP + sumMotifUnmethAllM)
                else:
                    self.spikeInStats["both"]["unmethCtrl_"+motif+"_convRate"] = None
                self.spikeInStats["both"]["methCtrl_"+motif+"_converted"] = sumMotifMethConvP + sumMotifMethConvM
                self.spikeInStats["both"]["methCtrl_"+motif+"_all"] = sumMotifMethAllP + sumMotifMethAllM
                if sumMotifMethAllP + sumMotifMethAllM > 0:
                    self.spikeInStats["both"]["methCtrl_"+motif+"_convRate"] = float(sumMotifMethConvP + sumMotifMethConvM)/float(sumMotifMethAllP + sumMotifMethAllM)
                else:
                    self.spikeInStats["both"]["methCtrl_"+motif+"_convRate"] = None
                self.spikeInStats["both"]["unmethCtrl_"+motif+"_converted_minQual"] = sumMotifUnmethConvPminQual + sumMotifUnmethConvMminQual
                self.spikeInStats["both"]["unmethCtrl_"+motif+"_all_minQual"] = sumMotifUnmethAllPminQual + sumMotifUnmethAllMminQual
                if sumMotifUnmethAllPminQual + sumMotifUnmethAllMminQual > 0:
                    self.spikeInStats["both"]["unmethCtrl_"+motif+"_convRate_minQual"] = float(sumMotifUnmethConvPminQual + sumMotifUnmethConvMminQual)/float(sumMotifUnmethAllPminQual + sumMotifUnmethAllMminQual)
                else:
                    self.spikeInStats["both"]["unmethCtrl_"+motif+"_convRate_minQual"] = None
                self.spikeInStats["both"]["methCtrl_"+motif+"_converted_minQual"] = sumMotifMethConvPminQual + sumMotifMethConvMminQual
                self.spikeInStats["both"]["methCtrl_"+motif+"_all_minQual"] = sumMotifMethAllPminQual + sumMotifMethAllMminQual
                if sumMotifMethAllPminQual + sumMotifMethAllMminQual > 0:
                    self.spikeInStats["both"]["methCtrl_"+motif+"_convRate_minQual"] = float(sumMotifMethConvPminQual + sumMotifMethConvMminQual)/float(sumMotifMethAllPminQual + sumMotifMethAllMminQual)
                else:
                    self.spikeInStats["both"]["methCtrl_"+motif+"_convRate_minQual"] = None
            self.spikeInStatsHeaders += ["unmethCtrl_"+motif+"_converted","unmethCtrl_"+motif+"_all","unmethCtrl_"+motif+"_convRate","methCtrl_"+motif+"_converted","methCtrl_"+motif+"_all","methCtrl_"+motif+"_convRate"]
            statsHeader_minQual += ["unmethCtrl_"+motif+"_converted_minQual","unmethCtrl_"+motif+"_all_minQual","unmethCtrl_"+motif+"_convRate_minQual","methCtrl_"+motif+"_converted_minQual","methCtrl_"+motif+"_all_minQual","methCtrl_"+motif+"_convRate_minQual"]
        self.spikeInStatsHeaders.sort()
        statsHeader_minQual.sort()
        self.spikeInStatsHeaders = ["acceptedReadsMethCtrl","acceptedReadsUnmethCtrl"] + self.spikeInStatsHeaders + statsHeader_minQual
        for s in self.spikeInStats.keys():
            self.spikeInStats[s]["acceptedReadsMethCtrl"] = 0
            self.spikeInStats[s]["acceptedReadsUnmethCtrl"] = 0
        
        for ctrl in sic.controls.keys():
            if sic.controls[ctrl]["methStatus"] == "methylated":
                if sic.controls[ctrl]["strand"] == "-":
                    if processMstrand:
                        self.spikeInStats["-"]["acceptedReadsMethCtrl"] += sic.controls[ctrl]["acceptedReadCount"]
                else:
                    if processPstrand:
                        self.spikeInStats["+"]["acceptedReadsMethCtrl"] += sic.controls[ctrl]["acceptedReadCount"]
                if processBothStrands:
                    self.spikeInStats["both"]["acceptedReadsMethCtrl"] += sic.controls[ctrl]["acceptedReadCount"]
            elif sic.controls[ctrl]["methStatus"] == "unmethylated":
                if sic.controls[ctrl]["strand"] == "-":
                    if processMstrand:
                        self.spikeInStats["-"]["acceptedReadsUnmethCtrl"] += sic.controls[ctrl]["acceptedReadCount"]
                else:
                    if processPstrand:
                        self.spikeInStats["+"]["acceptedReadsUnmethCtrl"] += sic.controls[ctrl]["acceptedReadCount"]
                if processBothStrands:
                    self.spikeInStats["both"]["acceptedReadsUnmethCtrl"] += sic.controls[ctrl]["acceptedReadCount"]            

                
    def writeStatistics2files(self,path,options):
        '''
        writes the global and if options required also the lanewise statistics
        returns the filename of the global statistics file
        optionally also appends the statistics to a given user specified summary statistics files for multiple samples
        '''
        statisticsHeader = "sampleName\tdataFiles\tseqMotif\tstrand\ttotalMappedReads\tinformativeReads\tinformativeRate\tuniqueSeqMotifCount\tbisulfiteConversionRate\t" + \
                           "globalMethylationMean\tglobalMethylationMean_conversionAdjusted\tfirstSeqMotifMethylationMean\tuniqueSeqMotifMethylationMean\tuniqueSeqMotifMethylationStdev\t" + \
                           "readMethylationMean\treadMethylationStdev\tconvertedCytosineCount\ttotalCytosineCount\tconvertedSeqMotifCount\ttotalSeqMotifCount\tconvertedCytosineCount_nonSeqMotif\t" + \
                           "totalCytosineCount_nonSeqMotif\tcorrectAlignmentRate_median\tcorrectAlignmentRate_q1\tcorrectAlignmentRate_q3\tbisulfiteConversionRate_median\tbisulfiteConversionRate_q1\tbisulfiteConversionRate_q3\t" + \
                           "medianSequenceQualityScoreAllBases_readFirst25p\tmedianSequenceQualityScoreAllBases_readCenter50p\tmedianSequenceQualityScoreAllBases_readLast25p\tmedianSequenceQualityScoreAcceptedSeqMotif_readFirst25p\t" + \
                           "medianSequenceQualityScoreAcceptedSeqMotif_readCenter50p\tmedianSequenceQualityScoreAcceptedSeqMotif_readLast25p\tmedianSequenceQualityScoreDiscardedSeqMotif_readFirst25p\tmedianSequenceQualityScoreDiscardedSeqMotif_readCenter50p\t" + \
                           "medianSequenceQualityScoreDiscardedSeqMotif_readLast25p\ttotalMappedReads\tinformativeReads\ttotalPairedReadCount\tacceptedPairedReadCount\ttotalProperPairedReadCount\tacceptedProperPairReadCount\tmateNotFoundInRegionCount\toverlappingPairedReads\tpairedReadsUndigestedFragments\tprocessedDuplicatesFlag\tprocessedDuplicatesPosCount\t" + \
                           "discard_PFstatus\tdiscard_minReadLength\tdiscard_mismatches\tdiscard_restrSite\tdiscard_fragmentsOutsideSizeRange\tdiscard_disagreeWithStrandPrediction\tdiscard_improperPairs\tdiscard_duplicatesFlag\tdiscard_duplicatesPosCount\tdiscard_noSeqMotif\tmappingOrderDoesNotMatchStrand\tliftOverMappingErrorCount\t"+ self.fragmentInfo.ranges4header
        #add headings that depend on multiple thresholds
        for t in SEQMOTIF_THRES:
            statisticsHeader += "readsWithAtLeast"+str(t)+"seqMotifs\t"
        for t in READ_THRES_POSITIONS:
            statisticsHeader += "coordsWithAtLeast"+str(t)+"reads\t"
        for t in READ_THRES_WINDOWS:
            statisticsHeader += "windowsWithAtLeast"+str(t)+"reads\t"
        statisticsHeader = statisticsHeader[:-1] + "\n" # cut off trailing '\t' and add linebreak
        '''
        HEADER LEGEND:
        sampleName:
        tdataFiles:
        seqMotif:    sequence pattern to be investigated for methylation
        strand:    
        totalMappedReads:    number of reads processed
                        for pairedEnd Sequencing: #paired_reads/2 + #single_reads according to pysam !Careful! not every paired read is also in a proper pair and the partner is not guaranteed to be found)
        informativeReads:    number of reads passing the quality checks and containing at least 1 CpG/seqMotif
                            for pairedEnd Seqeuncing: #accepted_paired_reads/2 + #accepted_single_reads according to pysam !Careful! not every paired read is also in a proper pair and the partner is not guaranteed to be found)
        informativeRate:    informativeReads/totalMappedReads
        uniqueSeqMotifCount:    number of unique CpG positions (a CpG is counted only once for forward and reverse strand)
        bisulfiteConversionRate: convertedCytosineCount_nonSeqMotif/totalCytosineCount_nonSeqMotif; for individual lanes in not CpG context the Conversion rate wrt CpG is reported
        globalMethylationMean:    (totalSeqMotifCount - convertedSeqMotifCount) /totalSeqMotifCount * 100
        globalMethylationMean_conversionAdjusted:    globalMethylationMean * bisulfiteConversionRate
        firstSeqMotifMethylationMean:    fraction of converted SeqMotifs where the first base of a read is the first base of the seqMotif
        uniqueSeqMotifMethylationMean, uniqueSeqMotifMethylationStdev: Mean methylation ratio of all unique seqMotifs and its Standard Deviation 
        readMethylationMean, readMethylationStdev:    Mean and Standard Deviation of #meth_seqMotifs_read/#seqMotifs_read across all accepted reads (in CpG sites that pass the quality criteria)
        convertedCytosineCount:    number of Cs that have undergone methylation (sum across ALL cytosines, not only the ones with high base quality)
        totalCytosineCount:    number of Cs total (sum across ALL cytosines, not only the ones with high base quality)
        convertedSeqMotifCount:    number of seqMotifs in which the C has been converted
        totalSeqMotifCount:    number of times the seqMotif pattern has been detected
        convertedCytosineCount_nonSeqMotif:    number of converted Cs out of seqMotif context
        totalCytosineCount_nonSeqMotif:    total number of Cs out of seqMotif context
        correctAlignmentRate_median(*): median accross all reads of the following fraction: 1 - mismatches_alignedRead / length_alignedRead
        correctAlignmentRate_q1, correctAlignmentRate_q3(*): 25 and 75 percent quartiles of the fraction above
        bisulfiteConversionRate_median(*): median conversion rate across all reads: #converted_Cs_nonseqMotif/#Cs_nonSeqMotif
        bisulfiteConversionRate_q1, bisulfiteConversionRate_q3(*): 25 and 75 percent quartiles of the fraction above
        medianSequenceQualityScoreAllBases_readFirst25p(*): 25 percent quartile of the median quality among all bases from each read (median is calculated for each read)
        medianSequenceQualityScoreAllBases_readCenter50p(*): median of the median quality among all bases from each read (median is calculated for each read)
        medianSequenceQualityScoreAllBases_readLast25p(*): 75 percent quartile of the median quality among all bases from each read (median is calculated for each read)
        medianSequenceQualityScoreAcceptedSeqMotif_readFirst25p(*): 25 percent quartile of the median quality among all Cs where the C and the subsequent base in the seqMotif fulfill the quality criterion in seqMotif context from each read (median is calculated for each read)
        medianSequenceQualityScoreAcceptedSeqMotif_readCenter50p(*): median of the median quality among all Cs where the C and the subsequent base in the seqMotif fulfill the quality criterion in seqMotif context from each read (median is calculated for each read)
        medianSequenceQualityScoreAcceptedSeqMotif_readLast25p(*): 75 percent quartile of the median quality among all Cs where the C and the subsequent base in the seqMotif fulfill the quality criterion in seqMotif context from each read (median is calculated for each read)
        medianSequenceQualityScoreDiscardedSeqMotif_readFirst25p(*): 25 percent quartile of the median quality among all Cs where the C fulfills the quality criterion but the subsequent base does not in seqMotif context from each read (median is calculated for each read)
        medianSequenceQualityScoreDiscardedSeqMotif_readCenter50p(*): median of the median quality among all Cs where the C fulfills the quality criterion but the subsequent base does not in seqMotif context from each read (median is calculated for each read)
        medianSequenceQualityScoreDiscardedSeqMotif_readLast25p(*): 75 percent quartile of the median quality among all Cs where the C fulfills the quality criterion but the subsequent base does not in seqMotif context from each read (median is calculated for each read)
        totalMappedReads:  number of reads processed (this time this is the absolute number, i.e. not grouped by pairs)
        informativeReads: number of reads passing the quality checks and containing at least 1 CpG/seqMotif (this time this is the absolute number, i.e. not grouped by pairs)
        totalPairedReadCount: if the pairedEnd option is enabled: number of paired reads (according to the SAM/BAM flag)
        acceptedPairedReadCount: if the pairedEnd option is enabled: number of ACCEPTED paired reads (according to the SAM/BAM flag)
        totalProperPairedReadCount: if the pairedEnd option is enabled: number of properly paired reads (according to the SAM/BAM flag)
                                #In the Picard Maq pipeline, a two reads are considered a proper pair if:
                                #* neither read has the unmapped flag set
                                #* both reads map to the same reference sequence, and that sequence is not "*"
                                #* the reads do not map to the same strand
                                #* for standard (non-jumping) libraries, the start position of the positive end is less than the start position of the negative end plus the read length
                                #(for jumping libraries, we assert that the negative end's start position is less than the positive end's start position plus the read length)
        properPairReadCount:    if the pairedEnd option is enabled: number of ACCEPTEDproperly paired reads (according to the SAM/BAM flag)

        mateNotFoundInRegionCount: if the pairedEnd option is enabled: number of Properly Paired reads for which the mate was not found in the same processed region
        overlappingPairedReads: if the pairedEnd option is enabled: number of proper paired Sequences (where the partner sequence is found and that pass the alignment mismatch criterion (and the CGG criterion if required)) that run into each other. Note: its the number of paired sequences, not the number of pairs (this number can be approximated by deviding by 2)
        pairedReadsUndigestedFragments: if the pairedEnd option is enabled and the rrbsMode option is enabled: the number of reads indicating undigested fragments:
                                        case: + strand: fragment_end_read < mate_end, case: - strand: fragment_end_mate < read_end
        processedDuplicatesFlag: number of reads processed but marked with the duplicate flag in the input file
        processedDuplicatesPosCount: number of reads processed but which map to genomic position (start,end,AMPLIFIED strand,lane=input file) that is already occupied
        
        
        The following filter statistics are successive and exclusive. I.e. a read filtered out by a previous step can not be subject to filtering of a successive step
        Filter statistics are given in the order the filters are applied
        discard_PFstatus:      number of reads discarded due to quality check flag in input file (options.pfStatus)
        discard_minReadLength: nuber of reads discarded due too short reads (shorter than the minReadLength parameter)
        discard_mismatched: number of reads discarded due to mismatch criteria (see the following parameters for details: minMismatches,maxMismatches,mismatchTruncate)
        discard_restrSite:        number of reads discarded due to restriction site start criterion
        discard_fragmentsOutsideSizeRange: number of reads discarded due them being matched to a fragment that does not fit the [min|max]FragmentLength criteria
        discard_disagreeWithStrandPrediction: if the strandPredict option is enabled: number of reads where the logistic regression model strand prediction disagrees with the reads aligned strand
        discard_improperPairs:   number of reads discarded due to proper pair Paired End read handling (if the parameter discardImproperPairs is set, all reads not flagged with proper pairs will be counted as discarded, otherwise read#2s of improper pairs will be counted)
        discard_duplicatesFlag: number of reads discarded due to the duplicate flag in the input file
        discard_duplicatesPosCount: number of reads discarded due to mapping to a genomic position that already has a mapped read (!!!reads discarded and counted in discard_duplicatesFlag can also appear here!!!)
        discard_noSeqMotif:      number of reads discarded because they do not contain a CpG (or different sequence motif looked for)
                
        liftOverMappingErrorCount: number of mapping errors in the liftover call
        mappingOrderDoesNotMatchStrand: assuming readnumbers correspond to the order of mapping coordinates (on genomic + strand), this reports the number of reads where the readnumber does not correspond to the mapped strand. I.e. it is assumed that read#1s map to the genomic + strand while read#2s map to the genomic - strand.
        self.fragmentInfo.ranges4header: Fragment info/counts
        readsWithAtLeastXseqMotifs: number of reads containing at least X CpGs/SequenceMotifs. Only processed reads and reads containing no CpG/motif are accounted for
        coordsWithAtLeastXreads: number of unique genomic coordinates (start,end,strand) occupied by at least X reads. Only processed reads, duplicate reads and reads containing no CpG/motif are accounted for
        windowsWithAtLeastXreads: number of tiling windows (size: generally: 50bp (differs for processed region boundaries)) covered by at least X reads. Only processed reads and reads containing no CpG/motif are accounted for
        
        (*) in order to enhance computational performance, Medians and Quartiles are calculated for each subprocess window. then medians and quartiles are calculated from those medians for the overall statistics  
        
        unless otherwise indicated, for paired end sequencing the read numbers correspont to actual number of sequences, NOT sequenced fragments
        '''
        
        
        
        # prepare global statistics output file(s)
        statsFileName = path + options.methodPrefix +"_statistics_"+ options.sampleName + ".txt"

        outfileGlobalStats = open(statsFileName, 'w')
        outfileGlobalStats.write(statisticsHeader)
        if len(options.appendStatisticsOutput)>1:
            outfileAppendStats = open(options.appendStatisticsOutput, 'a')
            if os.path.getsize(options.appendStatisticsOutput) == 0: outfileAppendStats.write(statisticsHeader)     
        
        for seqMotif in self.seqMotifs: 
            for s in self.strands:
                if self.numberOfLanes>1:
                    for lane in range(self.numberOfLanes): 
                        #FM: why does this have to be set again? (it is set twice in the original script)?
                        #NJ: don't really know what you mean, since code looks so different now, but maybe it was set for global and lane separately. Or seqMotif specific, since only the conversion rate for CpG is of interest
                        #FM:? I mean, it is already set in the call of calcStats_setConversionRate() and here again?
                        self.setConversionRate(seqMotif,lane,s)
                outfileGlobalStats.write(str(self.globalStats[seqMotif][s])) 
                if len(options.appendStatisticsOutput)>1: outfileAppendStats.write(str(self.globalStats[seqMotif][s]))
            for lane in range(self.numberOfLanes):
                for s in self.strands:
                    if options.laneSpecificStatistics and (self.numberOfLanes>1):
                        outfileGlobalStats.write(str(self.laneStats[seqMotif][lane][s]))
                        if len(options.appendStatisticsOutput)>1: outfileAppendStats.write(str(self.laneStats[seqMotif][lane][s]))
        outfileGlobalStats.close()
        if len(options.appendStatisticsOutput)>1: outfileAppendStats.close()
        return statsFileName
    
    def writeSpikeInStatistics(self,path,options):
        '''
        process spike in control statistics if required
        the output will have the following format:
        sampleName
        strand
        acceptedReadsUnmethCtrl:       number of reads processed for the unmethylated Control (after filtering)
        acceptedReadsMethCtrl:         number of reads processed for the methylated Control (after filtering)
        for every sequence motif in which a C occurs in the spike in control sequences the following 12 statistisc are written:
        unmethCtrl[motif]converted:    number of bisulfite converted Cs found in all reads matched to all unmethylation controls (the reads have to pass the --spikeInMinMatchLen and the --spikeInMaxMismatches criteria)
        unmethCtrl[motif]all:          number of all Cs found in all reads matched to all unmethylation controls 
        unmethCtrl[motif]convRate:     unmethCtrl[motif]converted/unmethCtrl[motif]all
        methCtrl[motif]converted:      number of bisulfite converted Cs found in all reads matched to all methylation controls (the reads have to pass the --spikeInMinMatchLen and the --spikeInMaxMismatches criteria)
        methCtrl[motif]all:            number of all Cs found in all reads matched to all methylation controls 
        methCtrl[motif]convRate:       methCtrl[motif]converted/methCtrl[motif]all
        the analogous definitions hold for unmethCtrl[motif]converted_minQual,unmethCtrl[motif]all_minQual,unmethCtrl[motif]convRate_minQual, methCtrl[motif]converted_minQual,methCtrl[motif]all_minQual,methCtrl[motif]convRate_minQual
        except for the fact that only read positions where the C has at least base quality --baseQualityScoreC
        '''
        statsFileName = path + options.methodPrefix +"_statistics_spikeInCtrls_"+ options.sampleName + ".txt"
        
        header = "sampleName\tstrand"
        for hh in self.spikeInStatsHeaders:
            header += "\t" + hh
        
        ofile = open(statsFileName, 'w')
        ofile.write(header + "\n")
        for strand in self.strands:
            ofile.write(options.sampleName+"\t"+strand)
            for hh in self.spikeInStatsHeaders:
                ofile.write("\t"+str(self.spikeInStats[strand][hh]))
            ofile.write("\n")
        ofile.close()
        return statsFileName
    
    def printMe(self):
        '''
        method just used for debugging
        print the datastructure as a debug string
        '''
        for seqMotif in self.seqMotifs:
            for s in self.strands:
                print getDebugString("Global Statistics (" + str(seqMotif) + "," + str(s) + "):") 
                self.globalStats[seqMotif][s].printDebugAll()
                for lane in range(self.numberOfLanes):
                    print getDebugString("Lane Statistics (" + str(seqMotif) + "," + str(s) + "," + str(lane) +  "):")
                    self.laneStats[seqMotif][lane][s].printDebugAll()


def removeLinefeed(s):
    '''
    remove trailing \n and \t from a string
    '''
    try:
        while s[-1] in ["\n","\r"]: s = s[0:-1]
        return s
    except IndexError: return ""
                    
class OutFileHandler:
    '''
    handles the output filenames for performAnalysis and manageProcesses, i.e. the general combined outputfiles, and the combined files for seqMotifs
    Pattern and strand specific! dictionaries
    files handled: Reads, ReadDetails, CpGs, Fragments,
    output files generated for each region are stored in a list of files that require later concatenation
    also contains methods to get a filename dictionary for the genomic regions (i.e. the subprocesses) and combining those files after all the processes finished using cat 
    '''
    def __init__(self,path):
        self.outputFilenames = [] #final output files to be (zipped and) moved later. read details file and either bigbed or plain bed files
        self.bedFileList = {} # collect all bed file names for easy zipping and conversion to bigBed
        self.bedMethylationFileList = {} # collect all bed methylation file names for easy conversion to bigBed and prior format changes of the cpg methylation file
        self.seqMotifOutfileNames = {} # store all output files and make them accessible to the rest of the script for each seqMotif
        self.seqMotifOutfilePaths = {} #dictionary that stores just the base paths. should contain the same keys as self.seqMotifOutfileNames
        self.seqMotifFiles2cat = {} # for each seqMotif store strings of concatenated filenames for each genomic region process for later use of cat in manageProcesses()
        self.numFiles2cat = 0 # counter for the number of files to concatenate. important for controlling that a cat command does not get too many arguments
        self.numConcats4cat = 0 #counts the number of times the method __collapseCatFileStrings__ has been applied
        self.basePath = path
        if not os.path.exists(path): os.mkdir(path)
    
    def addFile(self,fn):
        '''
        add a final outputfile to be (zipped and) moved later
        '''
        self.outputFilenames.append(fn)
        
    def addFiles4seqMotif(self,seqMotif,options):
        '''
        add the filenames for a seqMotif and return a dictionary containing all relevant output files for a seqMotif
        if strand strand specificity is required, this is also handled
        '''
        self.seqMotifOutfileNames[seqMotif] = {}
        self.seqMotifOutfilePaths[seqMotif] = {}
        
        #renamed cph to noncpg in outputfilenames and headers
        seqMotifLower = "cp" + seqMotif.lower()
        if seqMotifLower == "cph":
            seqMotifLower = "noncpg"
        
        seqMotifPath = self.basePath + seqMotifLower + os.sep
        if not os.path.exists(seqMotifPath): os.mkdir(seqMotifPath)
        
        #Prepare filenames
        if not options.statisticsOnly:  
            if options.readDetailsFile:
                readDetailsPath = seqMotifPath + "readDetails" + os.sep
                if not os.path.exists(readDetailsPath): os.mkdir(readDetailsPath)
                detailsFilename = "%s%s_%sReadDetails_%s.txt" % (readDetailsPath, options.methodPrefix, seqMotifLower, options.sampleName)       
                self.outputFilenames.append(detailsFilename)  
                self.seqMotifOutfileNames[seqMotif]['readQualOutfile'] = detailsFilename
                self.seqMotifOutfilePaths[seqMotif]['readQualOutfile'] = readDetailsPath
                headerReadDetailsFile = "chrom\tchromstart\tchromend\tmethRatios\tmeanMeth\tstrand\treadstart\treadend\t" + seqMotifLower +"Count\t" + seqMotifLower +"Positions_inGenome" + "\tlaneNumber\tmismatches\tcorrectAlignmentRate\tcorrectAlignmentRateBySegment\tmedianQualityScoreAllBases\tmedianQualityScoreBySegment\tmedianQualityScoreForAcceptedCp" + seqMotif +"s\tmedianQualityScoreForDiscardedCp" + seqMotif +"s\tconversionRate\tconversionRateValues\tbasesAnalyzed\tPEnum\n" # maxSeedWithout_mismatches\tmismatches_completelyMethylated\tmaxSeedWithout_mismatches_completelyMethylated\tmismatches_completelyUnmethylated\tmaxSeedWithout_mismatches_completelyUnmethylated\t
                tempF = open(self.basePath + 'headerForMergeReadDetailsFile.txt','w') 
                tempF.write(headerReadDetailsFile)
                tempF.close()
                
            #create subdirectories if not present already:
            readsPath = seqMotifPath + "reads" + os.sep
            if not os.path.exists(readsPath): os.mkdir(readsPath)
            fragmentPath = seqMotifPath + "fragments" + os.sep
            if not os.path.exists(fragmentPath): os.mkdir(fragmentPath)
            methPath = seqMotifPath + "meth" + os.sep
            if not os.path.exists(methPath): os.mkdir(methPath)
            
            #the others are bed or bigbed
            genericOutfileName = "%s" + options.methodPrefix + "_" + seqMotifLower +"%s_"+ options.sampleName + ".bed"
            outfileName = genericOutfileName % (readsPath,"Reads")
            self.bedFileList[outfileName] = outfileName
            self.seqMotifOutfileNames[seqMotif]['readOutfile'] = outfileName
            self.seqMotifOutfilePaths[seqMotif]['readOutfile'] = readsPath
            if options.isStrandSpecific:
                outfileName = genericOutfileName % (fragmentPath,"FragmentCoverageMinus")
                self.bedFileList[outfileName] = outfileName
                self.seqMotifOutfileNames[seqMotif]['fragmentOutfileM'] = outfileName
                outfileName = genericOutfileName % (fragmentPath,"FragmentCoveragePlus")
                self.bedFileList[outfileName] = outfileName
                self.seqMotifOutfileNames[seqMotif]['fragmentOutfileP'] = outfileName
                self.seqMotifOutfilePaths[seqMotif]['fragmentOutfileP'] = fragmentPath
                self.seqMotifOutfilePaths[seqMotif]['fragmentOutfileM'] = fragmentPath
            else:   
                outfileName = genericOutfileName % (fragmentPath,"FragmentCoverage")  
                self.bedFileList[outfileName] = outfileName
                self.seqMotifOutfileNames[seqMotif]['fragmentOutfile'] = outfileName
                self.seqMotifOutfilePaths[seqMotif]['fragmentOutfile'] = fragmentPath
            if options.isStrandSpecific:
                outfileName = genericOutfileName % (methPath,"MethylationMinus")
                self.bedFileList[outfileName] = outfileName
                self.seqMotifOutfileNames[seqMotif]['cpgOutfileM'] = outfileName
                self.bedMethylationFileList[outfileName] = outfileName
                outfileName = genericOutfileName % (methPath,"MethylationPlus")
                self.bedFileList[outfileName] = outfileName
                self.seqMotifOutfileNames[seqMotif]['cpgOutfileP'] = outfileName
                self.bedMethylationFileList[outfileName] = outfileName
                self.seqMotifOutfilePaths[seqMotif]['cpgOutfileP'] = methPath
                self.seqMotifOutfilePaths[seqMotif]['cpgOutfileM'] = methPath
            else:   
                outfileName = genericOutfileName % (methPath,"Methylation")
                self.bedFileList[outfileName] = outfileName
                self.bedMethylationFileList[outfileName] = outfileName
                self.seqMotifOutfileNames[seqMotif]['cpgOutfile'] = outfileName
                self.seqMotifOutfilePaths[seqMotif]['cpgOutfile'] = methPath
    
    def initializeseqMotifFiles2cat(self,seqMotif):
        '''
        set the output file list for a pattern to the empty list
        '''
        self.seqMotifFiles2cat[seqMotif] = {}
        for fileType in self.seqMotifOutfileNames[seqMotif].keys():
            self.seqMotifFiles2cat[seqMotif][fileType] = []
    
    def __preConcatCatFiles__(self,seqMotif): #TODO: delete source files after concat?
        '''
        pre concatenate files in seqMotifFiles2cat in order to prevent the argumentlist from getting too long
        Potential speedup: if linear function is too slow, do tree recursion
        '''
        MAX_CAT_FILE_NUM = 100   
        for fileType in self.seqMotifOutfileNames[seqMotif].keys():
            running_index = 0
            while len(self.seqMotifFiles2cat[seqMotif][fileType]) > MAX_CAT_FILE_NUM:
                running_index += 1
                collapsedFileName = self.seqMotifOutfileNames[seqMotif][fileType] + ".concat" + str(running_index)
                #concatenate MAX_CAT_FILE_NUM into a preconcatenated file
                catCmd = 'cat ' + getSepStringFromList(self.seqMotifFiles2cat[seqMotif][fileType][:MAX_CAT_FILE_NUM]," ",addquotes=True) + ' > ' + collapsedFileName
                if debug:
                    print "pre-concatenate: " + catCmd
                os.system(catCmd)
                self.seqMotifFiles2cat[seqMotif][fileType] = [collapsedFileName] + self.seqMotifFiles2cat[seqMotif][fileType][MAX_CAT_FILE_NUM:]
                #delete previous concatenated file
                if running_index > 1:
                    collapsedFileName_prev = self.seqMotifOutfileNames[seqMotif][fileType] + ".concat" + str(running_index -1)
                    if debug:
                        print "  removing previous file: " + collapsedFileName_prev
                    os.system('rm ' + collapsedFileName_prev)
                
    def getRegionFileNameDictAndUpdateseqMotifFiles2cat(self,seqMotif,chrom,processUnit):
        '''
        generate a dictionary with the processes output filenames. This dictionary tells seqDataProcessing where to put its output files
        add the region specific output files to the lists of files to be concatenated later
        '''
        fileExt =  '_'+ chrom +'_'+ str(processUnit)
        processFileNames = {}
        self.numFiles2cat += 1
        for fileType in self.seqMotifOutfileNames[seqMotif].keys():
            chromFiletypePath = self.seqMotifOutfilePaths[seqMotif][fileType] + chrom + os.sep
            if not os.path.exists(chromFiletypePath): os.mkdir(chromFiletypePath)
            curFn = chromFiletypePath + os.path.split(self.seqMotifOutfileNames[seqMotif][fileType])[1] + fileExt
            processFileNames[fileType] = curFn
            self.seqMotifFiles2cat[seqMotif][fileType].append(curFn)

        return processFileNames
    
    def catRegionFiles_old(self,seqMotif):
        '''
        concatenate the output files generated for each region
        '''
        self.__preConcatCatFiles__(seqMotif)
        for finalOutfile in self.seqMotifOutfileNames[seqMotif].keys(): 
            if self.seqMotifFiles2cat[seqMotif][finalOutfile] == []:
                continue
            catCmd = 'cat ' + getSepStringFromList(self.seqMotifFiles2cat[seqMotif][finalOutfile]," ",addquotes=True) + ' > ' + self.seqMotifOutfileNames[seqMotif][finalOutfile]
            sys.stdout.flush()
            if debug:
                print 'Concatenate temp files: ' + catCmd[:160] + '[...] > ' + self.seqMotifOutfileNames[seqMotif][finalOutfile]
            os.system(catCmd)
            
    def catRegionFiles(self,seqMotif):
        '''
        concatenate the output files generated for each region
        '''
        
        for finalOutfile in self.seqMotifOutfileNames[seqMotif].keys(): 
            if self.seqMotifFiles2cat[seqMotif][finalOutfile] == []:
                continue
            #output the filenames to concat to a file
            fileNames2catFileName = self.seqMotifOutfileNames[seqMotif][finalOutfile] + ".files2concat"
            fileNames2catFile = open(fileNames2catFileName,"w")
            for fn in self.seqMotifFiles2cat[seqMotif][finalOutfile]:
                fileNames2catFile.write(fn + "\n")
            fileNames2catFile.close()
                
            catCmd = 'cat ' + fileNames2catFileName + ' | xargs cat > ' + self.seqMotifOutfileNames[seqMotif][finalOutfile]
            
            if debug:
                print 'Concatenate temp files: ' + catCmd
                sys.stdout.flush()
            os.system(catCmd)
    
    def addHeaderReadQualOutfile(self,seqMotif,options):
        '''
        add a header to the ReadDetails files
        '''
        if options.readDetailsFile: 
            sys.stdout.flush() 
            os.system('cat ' + self.basePath + 'headerForMergeReadDetailsFile.txt '+ self.seqMotifOutfileNames[seqMotif]['readQualOutfile'] + ' > ' + self.seqMotifOutfilePaths[seqMotif]['readQualOutfile'] +'mergedWheaderTemp.txt')
            os.system('mv ' + self.seqMotifOutfilePaths[seqMotif]['readQualOutfile'] +'mergedWheaderTemp.txt ' + self.seqMotifOutfileNames[seqMotif]['readQualOutfile'])
            
    def processBedFiles(self,options,chromSizesFile):
        '''
        add the final output files for later zipping and moving,
        perform bigBed conversion
        '''
        for bedFile in self.bedFileList.keys():
            if options.bigBedFormat:
                pathSplit = os.path.split(bedFile)
                bedPath = pathSplit[0] + os.sep
                bedFileS = pathSplit[1]   # returns filename only (without path)
                bedFileToConvert = bedFile
                isMethBedFile = False
                if bedFile in self.bedMethylationFileList.keys(): # format changes are necessary for the CpX methylation files prior to bigBed conversion
                    isMethBedFile = True
                    cpgMethBedFile = open(bedFile,"r")
                    bbFormatFilename = bedPath + "tempForBigBedFormat_" + bedFileS 
                    cpgMethBigBedFile = open(bbFormatFilename,"w")
                    for line in cpgMethBedFile:
                        line = removeLinefeed(line)
                        splitLine = line.split("\t")
                        cpgMethBigBedFile.write(splitLine[0] +"\t"+ splitLine[1] +"\t"+ splitLine[2] +"\t'"+ str(int(splitLine[4])/10)+"%["+ splitLine[3].split("/")[1].split("'")[0] +"]'\t"+splitLine[4] +"\t"+ splitLine[5]+ "\n")
                    cpgMethBedFile.close()
                    cpgMethBigBedFile.close()    
                    bedFileToConvert = bbFormatFilename
                    
                if options.keepBedFiles == "all" or (options.keepBedFiles == "CpGonly" and isMethBedFile):
                    # keep file and mark for final gzip
                    self.addFile(bedFile)
                    
                bedFileBigBed = os.path.splitext(bedFileS)[0]
                
                bigBedFormatCmd = options.toolsPath + 'bedToBigBed '+ bedFileToConvert + ' ' + chromSizesFile + ' ' + options.webOutputDir + os.sep + bedFileBigBed +'.bb'

                print "Converting to bigBed format: " + bigBedFormatCmd
                sys.stdout.flush()
                bigBedReturnCode = os.system(bigBedFormatCmd)
                try:
                    if bigBedReturnCode <> 0:
                        raise Exception("BigBed conversion failed")
                except:
                    print "WARNING: BigBed conversion failed. Command used: " + bigBedFormatCmd
                    exceptionOccurred = True
            else: self.addFile(bedFile)
            
    def zipAndMoveOutfiles(self):
        '''
        (zip and) move final output files
        '''
        for outputFilename in self.outputFilenames:
            if options.gzip:
                print "Zipping file: " + outputFilename
                sys.stdout.flush()
                os.system("gzip -f "+outputFilename)
                outputFilename = outputFilename + ".gz"
            print "Moving file: " + outputFilename
            sys.stdout.flush()
            if (not options.bigBedFormat) and (outputFilename.find('.bed.') > 0):
                os.system("mv "+ outputFilename+" "+ options.webOutputDir)   
            else: 
                os.system("mv "+ outputFilename+" "+ options.outputDir)

def pysamReadCmp(r,s):
    '''
    compare function for comparing pysam reads
    first criterion: start pos, second: readlength (Note that this corresponds to end pos for ungapped alignment like in rrbsmap)
    '''
    cmpVal = cmp(r.pos,s.pos)
    if cmpVal == 0:
        cmpVal = cmp(r.rlen,s.rlen)
        if cmpVal == 0:
            # - strand < + strand
            cmpVal = cmp(not r.is_reverse, not s.is_reverse)
    return cmpVal

def pysamReadCmp_withFileIndex(r,s):
    '''
    compare function for comparing tuples of pysam reads and file indices
    first criterion: read start pos, second: readlength (Note that this corresponds to end pos for ungapped alignment like in rrbsmap)
    '''
    return pysamReadCmp(r[0],s[0])

def getLengthOnReferenceFromCIGAR(read):
    '''
    determine in the reference sequence from an aligned read's CIGAR string
    '''
    #CIGAR format: [(X_1,l_1),...,(X_n,l_n)]
    #where X_i in {0:M,1:I,2:D,3:N,4:S,5:H,6:P}
    if read.cigar:
        return reduce(lambda x, y: x+y, map(lambda x:x[1], filter(lambda x: x[0] in [0,2,3,4] ,read.cigar))) # filter discards insertions and unkown flags
    else:
        return read.rlen

class ReadfileCollection:
    '''
    manages pysam bamfiles from multiple lanes
    '''
    def __init__(self,fileNameList,format="BAM"):
        self.readfiles = []
        self.readfileNames = fileNameList
        #open all the read/alignment files
        for rf in self.readfileNames:
            if debug: print "Loading file: " + rf
            if format == "BAM":
                self.readfiles.append(pysam.Samfile( rf, "rb" ))
            else:
                raise Exception("Unrecognized file format")

    def getEstimatedReadCount4Region(self,chrom,startPos,endPos):
        '''
        estimate the number of reads in a given region based on the first lane in the collection
        '''
        n = self.readfiles[0].count(chrom, startPos, endPos)
#FM[20110603]: installed pysam 0.4.1 which has a count function
#        n = 0
#        for read in self.readfiles[0].fetch(chrom, startPos, endPos):
#            n += 1
        
        return n * len(self.readfiles)
    
    def getExactReadCount4Region(self,chrom,startPos,endPos):
        '''
        count the exact number of reads in a given region based on all lanes in the collection
        '''
        n = 0
        for rf in self.readfiles:
            n += rf.count(chrom, startPos, endPos)
#FM[20110603]: installed pysam 0.4.1 which has a count function
#            for read in self.readfiles[0].fetch(chrom, startPos, endPos):
#                n += 1
        return n
            
    def getReadsInRegion(self,chrom,startPos,endPos):
        '''
        for given pysam samfile objects returns a SORTED list of all reads in that region using the pysam fetch() method
        '''
        #readList is a dynamic list initialized with length readListSize. readCount-1 determines the highest index in which a meaningful entry can be found
        #when the number of meaningful elements reaches the size limit the size of the list is doubled
        readListSize = 100000
        readList = readListSize * [None]
        readCount = 0
        for rf in self.readfiles:
            #FM: pysam.fetch() returns a sorted iterator when fetching (see pysam API: http://wwwfgu.anat.ox.ac.uk/~andreas/documentation/samtools/api.html#pysam.Samfile.fetch)
            #syntax: samfile.fetch(CHROMOSOME, FROMPOS, TOPOS) FROMPOS is the first position included, TOPOS is the first position not included (0-based indexing)
            for read in rf.fetch(chrom, startPos, endPos):
                #double the size of the dynamic list as it fills up
                if readListSize <= readCount:
                    readList = readList + readListSize * [None]
                    readListSize *= 2
                readList[readCount] = read
                readCount += 1
        
        #sort according to start and end pos
        readList = readList[0:readCount]
        readList.sort(pysamReadCmp)    
        return readList
    
    def getReadsInRegion_withFileIndex(self,chrom,startPos,endPos):
        '''
        does the same as getReadsInRegion(), but also returns the index of the file a read came from, i.e. a list of tuples: (pysam.AlignedRead,fileIndex)
        '''
        #readList is a dynamic list initialized with length readListSize. readCount-1 determines the highest index in which a meaningful entry can be found
        #when the number of meaningful elements reaches the size limit the size of the list is doubled
        readListSize = 100000
        readList = readListSize * [None]
        readCount = 0
        for i in range(len(self.readfileNames)): 
            rf = self.readfiles[i]
            #FM: pysam.fetch() returns a sorted iterator when fetching (see pysam API: http://wwwfgu.anat.ox.ac.uk/~andreas/documentation/samtools/api.html#pysam.Samfile.fetch)
            #syntax: samfile.fetch(CHROMOSOME, FROMPOS, TOPOS) FROMPOS is the first position included, TOPOS is the first position not included (0-based indexing)
            for read in rf.fetch(chrom, startPos, endPos):
                #double the size of the dynamic list as it fills up
                if readListSize <= readCount:
                    readList = readList + readListSize * [None]
                    readListSize *= 2
                readList[readCount] = (read,i)
                readCount += 1
        #sort according to start and end pos
        readList = readList[0:readCount]
        readList.sort(pysamReadCmp_withFileIndex)  
        return readList
    
    def getReadIteratorList(self,chrom,startPos,endPos):
        '''
        returns a list of iterators over all reads in all files of the read file collection
        '''
        ril = [None] * len(self.readfileNames)
        for i in range(len(self.readfileNames)): 
            rf = self.readfiles[i]
            ril[i] = rf.fetch(chrom, startPos, endPos)
        return ril
     
    def __getLastReadEnd4region__(self,readIterList):
        '''
        returns the last end position on the reference
        readIterList is a list containing read iterators
        '''
        maxEndPos = -1
        numReadsConsidered = 0
        for ri in readIterList:
            for r in ri:
                numReadsConsidered += 1
                lenOnRef = getLengthOnReferenceFromCIGAR(r)
                if lenOnRef == 0:
                    continue
                curEndPos = r.pos + lenOnRef
                if curEndPos > maxEndPos:
                    maxEndPos = curEndPos
        return((maxEndPos,numReadsConsidered))

    def findNextGap(self,pos,chrom,chromSize):
        '''
        Sweeping algorithm to
        deterministically find the next gap position downstream of a genomic position (on a given chrom/scaffold)
        '''
        if debug:
            print "trying unit boundary: " + str(pos)
         
        runPos = pos
        minGapSize = 1 # minimum size of the gap
        while runPos < chromSize:
            readIter = self.getReadIteratorList(chrom,runPos,runPos+minGapSize)

            prevPos = runPos
            (runPos,numReads) = self.__getLastReadEnd4region__(readIter)
            if numReads < 1:
                runPos = prevPos  #avoid -1 value coming from self.__getLastReadEnd4region__()
                break #found a region with no reads
            if prevPos == runPos:
                raise Exception("Encountered endless loop in finding the next window gap")
            if debug:
                print "  --> " + str(numReads) + " reads found in region --> trying unit boundary: " + str(runPos)
            if runPos == -1:
                raise Exception("could not compute last read end position for chrom " + chrom + " start position " + str(pos) )
            del readIter
                
        #ran across end of the chromosome/genomic region    
        if runPos > chromSize:
            runPos = chromSize
        
        return runPos
    
    def getMinStartMaxEnd4region(self,chrom,startPos,endPos):
        '''
        find the minimum start coordinate and the maximum end coordinate for a genomic region
        i.e. the minimum position of a read overlapping the region start site and the maximum read end position overlapping with the end site
        if no reads overlap with start and/or end, the original coordinates will be returned
        '''
        minPos = startPos
        readIterMin = self.getReadIteratorList(chrom,startPos,startPos+1)
        for ri in readIterMin:
            for r in ri:
                curStartPos = r.pos
                if curStartPos < minPos:
                    minPos = curStartPos
        
        maxPos = endPos
        readIterMax = self.getReadIteratorList(chrom,endPos,endPos+1)
        for ri in readIterMax:
            for r in ri:
                lenOnRef = getLengthOnReferenceFromCIGAR(r)
                if lenOnRef == 0:
                    lenOnRef = r.qlen
                curEndPos = r.pos + lenOnRef
                if curEndPos > maxPos:
                    maxPos = curEndPos
                    
        return((minPos,maxPos))
                    
    def getReferenceName(self,refIndex,fileIndex):
        return self.readfiles[fileIndex].references[refIndex]
    
    def checkForZSfield(self):
        '''
        check whether the costum ZS field is present in all read files as determined by the presence in the first aligned, paired read of each file
        '''
        allZS = True
        for i in range(len(self.readfiles)):
            rf = self.readfiles[i]
            rfname = self.readfileNames[i]
            reads = rf.fetch()
            read2check = reads.next()
            #get the first properly paired (and therefor also aligened in a pair) read
            try:
                while read2check and read2check.is_proper_pair == 0:
                    read2check = reads.next()
            except StopIteration: #ran through all reads in the file and found no suitable one
                read2check = None
            if not read2check:
                print "WARNING: no properly paired read found for checking for ZS field in file " + rfname
                sys.stdout.flush()
                allZS = False
                continue
            try:
                v = read2check.opt("ZS")
            except KeyError:
                print "NOTE: detected no ZS flag in file " + rfname
                print "  read: " + str(read2check)
                allZS = False
                
        return allZS
    
    def __del__(self):
        for rfo in self.readfiles:
            rfo.close()

class PysamMateHash:
    '''
    stores mate pairs (as pysam read objects) for a genomic region (as identified by the same qname) for each fileindex (eHarmony for reads in a region ;-))
    '''
    def __init__(self,readList,readfileNum):
        self.mateHash = readfileNum * [None]
        for i in range(readfileNum):
            self.mateHash[i] = {}
        for (rr,ii) in readList:
            if rr.is_paired:
                try:
                    self.mateHash[ii][rr.qname].append(rr)
                except KeyError:
                    self.mateHash[ii][rr.qname] = [rr]
                    
    def getReadfromPair(self,qname,readNum,fileIndex):
        '''
        retrieve a pysam read object as specified by qname and readnumber (for a given lane)
        '''
        try:
#            if len(self.mateHash[fileIndex][qname]) > 2:
#                if debug:
#                    print "WARNING: read " + qname + " has more than 2 occurrences"
#                return None
#            else:
            for rr in self.mateHash[fileIndex][qname]:
                if (rr.is_read1 and readNum == 1) or (rr.is_read2 and readNum == 2):
                    return rr
            return None                
        except KeyError:
            return None
    
    def addReads(self,readList):
        '''
        add all reads from a readlist
        !might lead to duplicates in the dictionary for an entry (i.e. identifiers with more than 2 reads)
        '''
        for (rr,ii) in readList:
            if rr.is_paired:
                try:
                    self.mateHash[ii][rr.qname].append(rr)
                except KeyError:
                    self.mateHash[ii][rr.qname] = [rr]
            
class chromFileCollection:
    '''
    handles chromosome files in an adaptive datastructure:
    each chromosomeis assigned a file OBJECT
    if there are too many open chromosome files, files will be closed at random until used again
    useful for handling genomes with huge chromosome number or simply scaffolds
    '''
    def __init__(self):
        self.chromFileObjects = {}
        self.chromFileNames = {}
        self.openFiles = 0
        self.MAXOPENFILES = 500
    
    def addFile(self,chromKey,fileName):
        self.chromFileNames[chromKey] = fileName
        
    def getFileObject(self,chromKey):
        '''
        get the fileobject corresponding to chromosome chromKey
        '''
        try:
            return self.chromFileObjects[chromKey]
        except KeyError:
            while self.openFiles >= self.MAXOPENFILES: #kick out one file object at random to make space
                k = self.chromFileObjects.keys()[random.randint(0, len(self.chromFileObjects.keys())-1)]
                del self.chromFileObjects[k]
                self.openFiles -= 1
            self.chromFileObjects[chromKey] = open(self.chromFileNames[chromKey],"r")
            self.openFiles += 1
            return self.chromFileObjects[chromKey]
    def __del__(self):
        for k in self.chromFileObjects.keys():
            self.chromFileObjects[k].close()
        

class strandClassifierLogReg:
    '''
    a logistic regression classifier for strand determination: class 0: "+" strand, class 1: "-" strand
    it was trained in R using a random subset of 0.1 percent of the samples of BiSeq_alignment_HUES3_WGBS_28_61EV0_7p using only proper pairs, on the "+" strand only sequences which did not contain >= 3 consecutive C and at least one G (vice versa for the reverse strand).
        Coefficients:
                        Estimate Std. Error z value Pr(>|z|)    
        (Intercept)      2.31527    5.45137   0.425    0.671    
        frac.c           0.70421    5.46871   0.129    0.898    
        frac.t          -2.50403    5.45496  -0.459    0.646    
        frac.a          -2.07207    5.45607  -0.380    0.704    
        frac.g          -5.47688    5.46804  -1.002    0.317    
        log10.ratio.c.t  0.01110    0.03825   0.290    0.772    
        log10.ratio.c.g  0.45502    0.02737  16.627   <2e-16 ***
        ---
        Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

    '''
    def __init__(self):
        #                           fraction #X/length               log10(X/Y)
        #              inter    c        t         a         g     c/t       c/g
        self.coeffs = [2.31527,0.70421,-2.50403,-2.07207,-5.47688,0.01110,0.45502]
        self.epsilon = 1e-10 #this was also the epsilon parameter for the trained model  
    
    def __determineLogFrac__(self,cx,cy):
        return math.log((cx+self.epsilon)/(cy+self.epsilon),10)
            
    def classifyRead(self,bsRead):
        '''
        determine parameters and and classify the strand of the read
        '''
        ll = len(bsRead.readSequence_raw)
        countC = float(bsRead.readSequence_raw.count("C"))
        countT = float(bsRead.readSequence_raw.count("T"))
        countA = float(bsRead.readSequence_raw.count("A"))
        countG = float(bsRead.readSequence_raw.count("G"))
        intercept = 1.0
        fracC = countC/ll
        fracT = countT/ll
        fracA = countA/ll
        fracG = countG/ll
#        logFracAC = self.__determineLogFrac__(countA, countC)
#        logFracAT = self.__determineLogFrac__(countA, countT)
#        logFracAG = self.__determineLogFrac__(countA, countG)
        logFracCT = self.__determineLogFrac__(countC, countT)
        logFracCG = self.__determineLogFrac__(countC, countG)
#        logFracTG = self.__determineLogFrac__(countT, countG)
        vec = [intercept,fracC,fracT,fracA,fracG,logFracCT,logFracCG]
        scalarProd = reduce(lambda x, y: x+y,map(lambda i: self.coeffs[i]*vec[i],range(len(self.coeffs)))) # compute the scalar product of the coefficients and the vector
        expScalarProd = math.exp(scalarProd)
        probMstrand = expScalarProd/(1.0+expScalarProd) #assigned probability of the model to classify to reverse strand
        return probMstrand
        
    def modelAgreesWithStrand(self,bsRead,confidenceM=0.7,confidenceP=0.3):
        '''
        does the assigned strand of a read agree with the predicted strand
        minus strand model confidence thresholds are used to determine agreement or not
        '''
        probM = self.classifyRead(bsRead)
        if bsRead.mappedStrand == "+":
            return probM <= confidenceP
        elif bsRead.mappedStrand == "-":
            return probM >= confidenceM
        else:
            raise Exception("Unknown Strand: " +bsRead.mappedStrand)


def lookupSequence(chromFile, chromstart, chromend, offset):
    '''
    obtain the DNA sequence for a given genomic region (0-based)
    lookup a sequence in a fasta file. offset should be the number of characters in the identifier line (including the > and \n)
    '''
    # retrieve genomic sequence
    pos = chromstart / 50 * 51 + chromstart % 50
    length = (chromend - chromstart + 50) / 50 * 51
    chromFile.seek(offset+pos,0)
    seq = chromFile.read(length)
    seq = seq.replace("\n","")
    seq = seq[0:(chromend-chromstart)]
    return seq                
                
class ReferenceGenomeRegion:
    '''
    responsible for reference sequence handling for a certain genomic region
    mainly used for retrieving the reference sequencs
    '''
    def __init__(self,chromFile,start,end,regionName,offset=-1):
        self.isWeirdChrom = False
        self.seq = ""
        self.chromFileObj = chromFile
        self.chromName = regionName
        self.offset = offset
        if offset < 0:
            self.offset = self.__qualityCheckChromFileAndGetOffset__()
        self.chromStart = start
        self.chromEnd = end
        
        self.seq = lookupSequence(self.chromFileObj,start,end, self.offset)
        self.len = len(self.seq)
        #just an assertion:
        ls = self.len
        lr = end - start
        if ls <> lr:
            print "WARNING: in reference sequence handling: length of retrieved sequence  (" + str(ls) + ") does not match genomic coordinates ("+regionName +"["+str(start)+","+str(end)+") - length: " +str(lr)+")\nfile used: " + str(chromFile)
        
    def __qualityCheckChromFileAndGetOffset__(self):
        '''
        check the offset (number of characters in the identifier line of fasta) and the format (50 columns)
        chromFileObj needs to be set at position 0
        '''
        self.chromFileObj.seek(0)
        offset = 0
        # determine offset arising from sequence name (first line in FASTA format)
        firstByte = self.chromFileObj.read(1)
        if firstByte == ">": offset = len(self.chromFileObj.readline()) + 1 # +1 for first byte
        # make sure that the row length is 50 bp and the line separator is "\n"
        self.chromFileObj.seek(offset+49,0)
        seq = self.chromFileObj.read(3)
        validLetters = ["A","C","G","T","N"]
        if seq == '':
            print "WARNING(severe): small chromosome detected: " + self.chromName
            self.isWeirdChrom = True
            return offset
        if not seq[0].upper() in validLetters or not seq[1] == "\n" or not seq[2].upper() in validLetters:
            print "Incorrect file format detected for "+ self.chromName +". Sequence around position 50: "+seq
            raise SystemExit
        return offset
    
    def getRefSeq(self,chromStart,chromEnd,reverse=False):
        '''
        retrieve the reference sequence
        0-based, end: first base not in sequence
        '''
        iStart_org = chromStart - self.chromStart
        iEnd = chromEnd - self.chromStart
        iStart = max(0,iStart_org)
        retseq = self.seq[iStart:iEnd]
        #if out of bounds, fill up the reference with 'N's
        if iStart_org < 0:
            retseq = abs(iStart_org)*"N" + retseq
        if iEnd > self.len:
            retseq = retseq + (iEnd-self.len)*"N"
        if reverse:
            return strRevComp(retseq)
        else:
            return retseq


class Windows4region:
    '''
    sliding window statistics for a genomic region (for each lane)
    keeps track of read counts for each window within a certain genomic region
    '''
    def __init__(self,start,end,numLanes,size=50):
        '''
        inputs should be integers
        '''
        self.start = start
        self.end = end
        self.numLanes = numLanes
        self.size = size
        self.numWindows = (end-start)/size + 1 #implicit int concersion assumed and the / operator should result in rounding down
        self.readCounts = []
        for i in range(numLanes):
            self.readCounts.append(self.numWindows * [0])
        
    def addReadStats(self,read,lane):
        '''
        calculate the windows overlapping with the read and increment the readcounts for those windows
        ! also counts reads outside the region!
        '''
        affectedWindowsStart = max(0,(read.readstart-self.start)/self.size)
        affectedWindowsEnd = min((read.readend-1-self.start)/self.size,self.numWindows - 1) #-1 because readend is the first position not spanned by the read
        for i in range(affectedWindowsStart,affectedWindowsEnd+1):
            self.readCounts[lane][i] += 1
    
    def getReadCounts4thres(self,thres,lane):
        '''
        given a SORTED list of threshold values, returns a dictionary containing the number of windows with a readcount >= the thresholds 
        '''
        retD = {}
        for t in thres:
            retD[t] = 0
        for rc in self.readCounts[lane]:
            for t in thres:
                if rc >= t:
                    retD[t] +=1
                else:
                    break
        return retD

class ReadPositions4region:
    '''
    keeps track of the unique read positions (start and end coordinate (remark: strand is handled in seperate class instances)) for a genomic region (for each lane)
    and their respective readcounts
    '''
    def __init__(self,start,end,numLanes):
        '''
        inputs should be integers
        '''
        self.start = start
        self.end = end
        self.numLanes = numLanes
        self.readCounts = numLanes * [None] # double dictionaries: first index: start position, second index: end position; contains tuples: first element: total reads seen at this position, second element: #reads that are not flagged as PCR duplicates seen at this position
        for i in range(numLanes):
            self.readCounts[i] = {}
        
    def addReadStats(self,read,lane):
        '''
        identify if the reads position is already in the datastructure and increment the readcount or insert the coordinates into the datastructure
        '''
        incNonFlagged = 0
        if not read.isDuplicateFlag:
            incNonFlagged = 1
        try:#does the index [read.readstart][read.readend] exist
            self.readCounts[lane][read.readstart][read.readend][0] += 1
            self.readCounts[lane][read.readstart][read.readend][1] += incNonFlagged       
        except KeyError:
            try: #does the index [read.readstart] exist
                self.readCounts[lane][read.readstart][read.readend] = [1,0]
                self.readCounts[lane][read.readstart][read.readend][1] += incNonFlagged
            except KeyError:
                self.readCounts[lane][read.readstart] = {}
                self.readCounts[lane][read.readstart][read.readend] = [1,0]
                self.readCounts[lane][read.readstart][read.readend][1] += incNonFlagged

    def getReadCounts4thres(self,thres,lane,startPos,endPos):
        '''
        given a SORTED list of threshold values, returns a dictionary containing the number of positions occupied by >= threshold reads
        for a given region
        '''
        retD = {}
        for t in thres:
            retD[t] = 0
        #looks more terrible than it is: eg self.readCounts.keys() is maximally ["+","-"]
        #the worst case iteration number is the number of reads in the region (if each read has pairwise distinct coordinates)
        for start in self.readCounts[lane].keys():
            if start >= startPos and start < endPos:
                for end in self.readCounts[lane][start].keys():
                    for t in thres:
                        if self.readCounts[lane][start][end][0] >= t:
                            retD[t] +=1
                        else:
                            break
        return retD
    
    def containsPosition(self,start,end,lane):
        '''
        has the position given in the arguments been recorded in the readCounts. useful for duplicate marking
        '''
        try:
            a = self.readCounts[lane][start][end]
            return True
        except KeyError:
            return False
        
    def containsPosition_nonFlagged(self,start,end,lane):
        '''
        has the position given in the arguments been recorded in the readCounts as a non-flagged read. useful for duplicate marking
        '''
        try:
            a = self.readCounts[lane][start][end]
            if a[1] > 0:
                return True
            else:
                return False
        except KeyError:
            return False
    
                    
class PatternBuckets4region:
    '''
    keeps track of reads with a MINIMUM of X CpGs/patterns for a genomic region (for each lane) the thresholds are called countBuckets here
    '''
    def __init__(self,start,end,numLanes,countBuckets=[0,1,2,5,10]):
        '''
        inputs should be integers. countBuckets specifies the buckets in which the reads are put into according to minimum CpG/pattern counts. Note: smaller buckets are subsets of larger buckets
        countBuckets are assumed to be sorted
        '''
        self.start = start
        self.end = end
        self.numLanes = numLanes
        self.countBuckets = countBuckets
        self.numBuckets = len(countBuckets)
        self.readCounts = numLanes * [None]
        for i in range(numLanes):
            self.readCounts[i] = {}
        for l in range(numLanes):
            for cb in countBuckets:
                self.readCounts[l][cb] = 0
        
    def addReadStats(self,read,lane):
        '''
        check if the reads contains more than thres (for thres in buckets) motifs and increment the count
        '''
        for cb in self.countBuckets: #count buckets assumed to be sorted
            if read.internal_cpgs >= cb:
                self.readCounts[lane][cb] += 1
            else:
                break
            
            
class ProcessShellScriptHandler:
    '''
    for parallelization: handles the shellscript to be executed in the LSF queue (or on a single node if --parallelize is not set)
    shellscript contains 1 or more calls of the 'region' tasks
    seperated by '&'
    if --parallelize is set, the LSF prepends the queue submission comman to the tasks
    also collects the pickle filenames containing the region statistics objects to be unified in the 'wrapup' task
    '''
    def __init__(self,path,seqMotif,maxPnum,submitCmd="bsub"):
        self.path = path + "cp" + seqMotif.lower() + os.sep + "shell" + os.sep#temp path + seq motif subdirectory
        if not os.path.exists(self.path): os.mkdir(self.path)
        self.seqMotif = seqMotif
        self.maxPnum = maxPnum #max number of processes in a shell script
        self.curNumP = 0
        self.submitCmd = submitCmd
        
        self.jobnames = [] # all the jobnames for submitted jobs (later important for dependencies)
        self.logFileNames = [] #all filenames for logfiles (for concatenating later into single logfile)
        self.pickleFileNames = []
        self.regionCoods = []
        self.curFileNameBase = "process4regions_" + self.seqMotif
        self.curFileName = self.curFileNameBase
        self.curProcessStr = ""
    
    def addPKLfile(self,curPKLfn):
        '''
        update the picklefiles to be unified in the 'wrapup' task. used in re-runs, when a region does not need to be processed as the picklefile
        is already existent
        '''
        self.pickleFileNames.append(curPKLfn)
        
    def addProcess(self,pString,regionStr,curPKLfn,regionCoordTupel):
        '''
        adds a 'region' task to the shellscript. if the number of processes per shelscript is reached it writes the shellscript and resets the variables
        for a new shellscript to start
        '''
        self.pickleFileNames.append(curPKLfn)
        self.regionCoods.append(regionCoordTupel)
        self.curFileName += "_" + regionStr
        self.curProcessStr += pString + " &\n"
        self.curNumP += 1
        if self.curNumP >= self.maxPnum:
            return self.writeShellScript()
        return ""
        
    def writeShellScript(self):
        '''
        output the shellscript
        '''
        #Write the content of the shellscript
        curPname = options.sampleName + "_" + self.curFileName + "_rand" + str(random.randint(0, 10000000))
        self.curFileName += ".sh"
        shellScript = open(self.path + self.curFileName,"w")
        shellScript.write("#!/bin/sh\n")
        shellScript.write(self.curProcessStr + "wait") #write a final wait for preventing the job from terminating when the last process has finished but previous processes have not
        shellScript.close()
        
        self.jobnames.append(curPname)
        curLogFile = self.path + curPname + ".log"
        self.logFileNames.append(curLogFile)
        jobString = "sh " + self.path + self.curFileName
        if options.parallelize:
#            jobString = self.submitCmd + " -o " + curLogFile + " -e " + curLogFile + ".err" + " -J " + curPname + " " + jobString
            jobString = self.submitCmd + " -o " + curLogFile + " -J " + curPname + " " + jobString
        #N#print "-------------------------------------------------------------------"
        #N#print "REGION job:"
        #N#print jobString
        
        #reset
        self.curFileName = self.curFileNameBase
        self.curProcessStr = ""
        self.curNumP = 0
        return jobString
        
class jobHandler:
    '''
    handles the jobs submitted to the cluster or executed on a single node
    1. collecting via addJob, then execute
    finally output to file
    collects all the summary shellscript executions  for the 'region' tasks (see class: ProcessShellScriptHandler) and the 'wrapup' call
    and executes them 
    '''
    def __init__(self):
        self.jobs = []
    
    def addJob(self,j):
        '''
        add a job to be executed
        '''
        self.jobs.append(j)
    
    def execute(self,wait=5):
        '''
        execute the jobs in the collection using os.system and wait inbetween
        '''
        for j in self.jobs:
            os.system(j)
            time.sleep(wait)
    
    def write2file(self,fn):
        '''
        write all the submitted jobs to an output file
        '''
        f = open(fn,"w")
        f.write("#!/bin/sh\n")
        for j in self.jobs:
            f.write(j + "\n")
        f.close()

class SpikeInControls:
    '''
    takes care of everything related to Spike In Controls
    '''
    def __init__(self,filename,motifLen):
        '''
        input: a fasta file containing spike in control sequences (forward strand only)
        'methylated' or 'unmethylated' should be last word of every description line
        '''
        self.controls = {}
        self.ids = []
        self.matchedReadFiles = {}
        self.motifs = {} #dictionary with motif as the key and the control ids and positions as tupel lists
        
        f = open(filename,"r")
        for record in Bio.SeqIO.parse(f, "fasta"):
            self.ids.append(record.id)
            curSeqFW=str(record.seq).upper()
            curSeqFWbitarrays = getSeqBitArrays(curSeqFW,padding=motifLen-1,addBisulfite=True,addH=True,alphabet=["C","G","A","T"],otherChar="N")
            curSeqRC=str(Bio.Seq.reverse_complement(record.seq)).upper()
            curSeqRCbitarrays = getSeqBitArrays(curSeqRC,padding=motifLen-1,addBisulfite=True,addH=True,alphabet=["C","G","A","T"],otherChar="N")
            methStatus = record.description.split()[-1]
            if not methStatus in ["methylated","unmethylated"]:
                raise Exception("Unknown Spike In Control methylation status:"+methStatus)
            #construct dictionaries storing the coverage and conversion information for each C
            countsFW = {}
            countsFWqual = {}
            for i in curSeqFWbitarrays["C"].search('1'): #go through all C positions
                countsFW[i] = {"conv":0,"total":0}
                countsFWqual[i] = {"conv":0,"total":0}
            countsRC = {}
            countsRCqual = {}
            for i in curSeqRCbitarrays["C"].search('1'): #go through all C positions
                countsRC[i] = {"conv":0,"total":0}
                countsRCqual[i] = {"conv":0,"total":0}
            
            self.controls[record.id+"+"]   = {"id":record.id,"methStatus":methStatus,"strand":"+","seq":curSeqFW,"seqBitArrays":curSeqFWbitarrays,"Ccounts":countsFW,"CcountsQual":countsFWqual,"acceptedReadCount":0}
            self.controls[record.id+"-"]   = {"id":record.id,"methStatus":methStatus,"strand":"-","seq":curSeqRC,"seqBitArrays":curSeqRCbitarrays,"Ccounts":countsRC,"CcountsQual":countsRCqual,"acceptedReadCount":0}
            
            seqLen = len(curSeqFW)
            for cPos in countsFW.keys():
                i = 1
                while i <= motifLen and (cPos+i) < seqLen:
                    mm = curSeqFW[cPos:(cPos+i)]
                    try:
                        self.motifs[mm].append((record.id+"+",cPos))
                    except KeyError:
                        self.motifs[mm] = [(record.id+"+",cPos)]
                    i += 1
            for cPos in countsRC.keys():
                i = 1
                while i <= motifLen and (cPos+i) < seqLen:
                    mm = curSeqRC[cPos:(cPos+i)]
                    try:
                        self.motifs[mm].append((record.id+"-",cPos))
                    except KeyError:
                        self.motifs[mm] = [(record.id+"-",cPos)]
                    i += 1
        f.close()
    
    def getRegExpr4id(self,id,length2check):
        '''
        returns a regular expression encoding the sequence pattern of length length2check (from the beginning of the sequence) of spike in control 'id' 
        bisulfite conversion is taken into account: both C and T fit Cs in the original sequence
        forward and reverse complement are combined into one regular expression
        '''
        length2check = min(length2check,len(self.controls[id+"+"]["seq"]))
        exprFW = self.controls[id+"+"]["seq"][0:length2check].replace("C","[CT]")
        exprRC = self.controls[id+"-"]["seq"][0:length2check].replace("C","[CT]")
        return "(("+exprFW+")|("+exprRC+"))"

    
    def extractSICreads(self,options):
        '''
        given alignment files in BAM format, extract all unaligned reads whose beginning matches exactly (disregarding C-T mismatches)
        to a Spike In Control Seqence on a length of length2check
        to seperate output SAM files
        '''
        self.matchedReadFiles = {}
        #-f 4: only unaligned reads
        callStringBase = "samtools view -f 4 %s | grep -i -P '\t" #-P for pearl regexpr dialect (to grep for tab). grep for tab in order to filter in the sequence field (heuristic)
        for id in self.ids:
            print "  extracting spike in control reads for " + id + "..."
            curRegEx = self.getRegExpr4id(id, options.spikeInMinMatchLen)
            print "    pattern to find: " + curRegEx
            curMatchedReadFile = options.tempDir + os.sep + options.methodPrefix + "_" + options.sampleName +"_spikeInControl_" + id + "_matchedReads.sam"
            self.matchedReadFiles[id] = curMatchedReadFile
            tmpHeaderFile = options.tempDir + os.sep + options.methodPrefix + "_" + options.sampleName +"_spikeInControl_" + id + "_samheader.txt"
            os.system("samtools view -H -o " + tmpHeaderFile + " " + options.alignmentFiles[0])
            callString = (callStringBase) % (options.alignmentFiles[0]) + curRegEx + "' > " + curMatchedReadFile + "_withoutHeader"
            if debug:
                print "Executing extraction command:"
                print callString
            os.system(callString)
            callString = "cat " + tmpHeaderFile + " " + curMatchedReadFile + "_withoutHeader >" + curMatchedReadFile
            if debug:
                print "concatenating SAM header and matched reads: "
                print callString
            os.system(callString)
            if len(options.alignmentFiles) > 1:
                for i in range(1,len(options.alignmentFiles)):
                    callString = (callStringBase) % (options.alignmentFiles[i])+ curRegEx + "' >> " + curMatchedReadFile
                    if debug:
                        print "Executing extraction command:"
                        print callString
                    os.system(callString)
        
    def processSICreads(self,options):
        '''
        process the reads found associated with SIC
        '''
        for id in self.ids:
            curMethCtrlStatus = self.controls[id+"+"]["methStatus"]
            curSam = pysam.Samfile(self.matchedReadFiles[id],"r")
            for alignedRead in curSam.fetch():
                #quality controls
                if options.pfStatus == 'PF':
                    if alignedRead.is_qcfail:
                        continue
                elif options.pfStatus == 'NonPF':
                    if not alignedRead.is_qcfail:
                        continue
                readBitArrays = getSeqBitArrays(alignedRead.seq,padding=0,padVal=False,addBisulfite=True,addH=False,alphabet=["C","G","A","T"],otherChar="N")
                readLength = len(readBitArrays["C"])
                matchArrayP = getBisulfiteMatchArray(self.controls[id+"+"]["seqBitArrays"],readBitArrays,rl=readLength)
                mismatchesP = matchArrayP.count(0)
                matchArrayM = getBisulfiteMatchArray(self.controls[id+"-"]["seqBitArrays"],readBitArrays,rl=readLength)
                mismatchesM = matchArrayM.count(0)
                #mismatches
                if min(mismatchesP,mismatchesM) > options.spikeInMaxMismatches:
                    continue
                #take the strand with fewest mismatches to the reference
                isReverse = False
                refControl = None
                if mismatchesM < mismatchesP:
                    isReverse = True
                    refControl = self.controls[id+"-"]
                else:
                    refControl = self.controls[id+"+"]
                refControl["acceptedReadCount"] += 1
                #look for Cs and check which ones are converted and which ones fulfill the quality criteria
                isConvertedC = refControl["seqBitArrays"]["C"][0:readLength] & readBitArrays["T"]
                isC = refControl["seqBitArrays"]["C"][0:readLength] & readBitArrays["Y"]
                readQuals = transformQualStr(alignedRead.qual)
                isQual =  bitarray(map(lambda x: x>=options.baseQualityScoreC,readQuals))
                for i in range(readLength):
                    if isC[i]:
                        refControl["Ccounts"][i]["total"] += 1
                        qq = isQual[i]
                        if qq:
                            refControl["CcountsQual"][i]["total"] += 1
                        if isConvertedC[i]:
                            refControl["Ccounts"][i]["conv"] += 1
                            if qq:
                                refControl["CcountsQual"][i]["conv"] += 1 
                    
            curSam.close()

class RestrictionSite:
    '''
    takes care of restriction site parsing and checking if read fulfills the criteria
    e.g. MspI: "C-CG_G"
    '''
    def __init__(self,resStr,bisulfite=True):
        '''
        parses the restriction site from a string.
        use lower case c to specify required UNmethylated Cs
        e.g. useful to specify that Cs out of CG context are unmethylated and thus bisulfite converted. Useful for read filtering by start sequence
        '''
        if resStr.find("_") < 0 or resStr.find("-") < 0 or resStr.count("_") > 1 or resStr.count("-") > 1:
            print "misspecified restriction site: " + resStr
            raise SystemExit
        
        resStr = resStr.strip()
        self.pattern = resStr.replace("_","").replace("-","")

        self.patternRC = strRevComp(self.pattern)
        self.siteFW = resStr.find("-")
        self.siteRConFWstrand = resStr.find("_")
        self.RCsiteBeforeFWsite = self.siteRConFWstrand < self.siteFW
        
        minSiteFW = min(self.siteFW,self.siteRConFWstrand)
        minSiteRC = len(resStr) - max(self.siteFW,self.siteRConFWstrand) - 1
        self.patStart = self.pattern[minSiteFW:]
        self.patEnd = self.patternRC[minSiteRC:]
        self.fillInsStart = 0 # the number of blunt end repair bases, they are methylated and thus bisulfite converted
        self.fillInsEnd = 0
        distRestrSites = abs(self.siteFW - self.siteRConFWstrand) - 1
        if self.RCsiteBeforeFWsite:
            self.fillInsStart = distRestrSites
        else:
            self.fillInsEnd = distRestrSites
        if bisulfite: # blunt end repair bases are methylated and thus bisulfite converted
            self.patStart = self.patStart[0:self.fillInsStart] + self.patStart[self.fillInsStart:].replace("C","Y").replace("c","T")
            self.patStart = self.patStart.upper()
            self.patEnd = self.patEnd[0:self.fillInsEnd] + self.patEnd[self.fillInsEnd:].replace("G","R").replace("g","A") #paired end second read
            self.patEnd = self.patEnd.upper()
        self.patStartLen = len(self.patStart)
        self.patEndLen = len(self.patEnd)
                  
    def readPassesRestrictionPattern(self,read):
        '''
        check if a given reads startsequnce fits the restriction site
        '''
        if read.account4beingRead2:
            for i in range(self.patEndLen):
                if not read.readSequenceBitArrays[self.patEnd[i]][i]:
                    return False
        else:
            for i in range(self.patStartLen):
                if not read.readSequenceBitArrays[self.patStart[i]][i]:
                    return False
        return True


def methylationCalling_bitarray(read, seqMotif, pattern, pattern_compl, statisticsLanes, allCpgsOfSampleGlobal,regionStartPos,regionEndPos):
    '''
    sets the reads C and methylation attributes and updates the statistics (given in statisticsLanes)
    based on the reads bitarrays encoding the read and reference sequence
    identifies the patterns in the read fulfilling the quality criteria
    updates the CpG/pattern information given in allCpgsOfSampleGlobal
    the many variants of seqMotif, pattern and pattern_rc are just in the signuture to improve runtime, as this function is called a lot
    '''
    patLen = len(pattern)
    readIsFW = read.mappedStrand == '+'
    readAmpRev = read.ampStrand == "-"
    # METHYLATION CALLING
    basesToAnalyze = read.basesToAnalyze
    if read.account4beingRead2:
        isC = read.readSequenceBitArrays["R"][:basesToAnalyze] & read.fragmentSequenceBitArrays["G"][patLen:(basesToAnalyze+patLen)] #bitarray for positions in read sequence being a Cytosine
        isCconv = read.readSequenceBitArrays["A"][:basesToAnalyze] & read.fragmentSequenceBitArrays["G"][patLen:(basesToAnalyze+patLen)] #C at current base at index has been converted/daminated to T    
    else:
        isC = read.readSequenceBitArrays["Y"][:basesToAnalyze] & read.fragmentSequenceBitArrays["C"][patLen:(basesToAnalyze+patLen)] #bitarray for positions in read sequence being a Cytosine
        isCconv = read.readSequenceBitArrays["T"][:basesToAnalyze] & read.fragmentSequenceBitArrays["C"][patLen:(basesToAnalyze+patLen)] #C at current base at index has been converted/daminated to T 
    
    isPattern = read.CpPs(pattern,pattern_compl)[:basesToAnalyze] #the seqMotif pattern is fit
    isPatternConv = isPattern & isCconv#the seqMotif pattern is fit and the C has been deaminated
    isNoCpg = read.nonCpPs(pattern,pattern_compl)[:basesToAnalyze]#the seqMotif pattern is not fit, eg CpA if looking for CpG methylation
    isNoCpgConv = isNoCpg & isCconv#the seqMotif pattern is not fit and the C has been deaminated, eg TpA if looking for CpG methylation

    internal_patterns_max = isPattern.count() # number of internal Patterns
    
    #find starting and end coordinates for bases to analyze for reads overlapping the processing window boundaries:
    #         ################
    #        [a              [b
    #         |               |
    #     ----p###>       ####q--->
    #     |   |    |      |   |    |
    #     s  qp    e      s  pq    e
    #     |   |    |      |   |    |
    #     <---q####       <###p----
    #     <---l--->       <---l--->
    #where a is the first base of the region to be processed, b is the first base not in the region anymore
    #s and e are readstart and readend respectively (defined as the first base included in the bases and the first base not included anymore)
    #and p and q are the starting and ending indecis of the read (relative to read orientation) of readbases within the region (###)
    # --> corresponds to a forward, <-- to a reversely aligned read
    #l designates the length of the read
    #if the region boundaries happen to cut through a CpG, i.e. a or b are the G in a CpG context, this CpG is counted towards the region where the C is located (wrt forward strand):
    #i.e. a CpG starting in in a-1 would be counted towards the previous region while a CpG in b-1 would be counted towards the current region. The example above shows a pattern length of one (e.g. CpG).
    # for non symmetrical patterns it is the same: when looking for CpA methylation and finding it on the genomic - strand, the reverse complement on the genomic + strand is TpG. If the T has coordinate a-1, it counts towards the previous region, if the T has coordinates b-1 it counts towards the current region
    # hence we also have to take into account which strand was amplified
    #For longer pattern length (eg 2 for CHG), the reverse read's p and q shift to the right
    ampStrandAdd = 0
    if readAmpRev:
        ampStrandAdd = patLen
    if readIsFW: # -->
        inRegionStart = max(0,regionStartPos - read.readstart + ampStrandAdd) #p = max(0,a-s) #+ampStrandAdd to account for amplified genomic - strand: exclude patterns  starting in a-patLen
        inRegionEnd = min(basesToAnalyze,regionEndPos - read.readstart + ampStrandAdd) #q = min(l,b-s) #+ampStrandAdd to account for amplified genomic - strand: include patterns  starting in b-patLen
    else: # <--
        inRegionStart = max(0,read.readend - regionEndPos - ampStrandAdd) #p = max(0,e-b) #-ampStrandAdd to include patterns starting in b-patLen (on the genomic + strand) if the genomic - strand is amplified
        inRegionEnd = min(basesToAnalyze, read.readend - regionStartPos - ampStrandAdd) #q = min(l,e-a)  #-ampStrandAdd to exclude patterns starting in a-patLen (on the genomic + strand) if the genomic - strand is amplified
    isInRegion = bitarray(basesToAnalyze)
    isInRegion.setall(False)
    isInRegion[inRegionStart:inRegionEnd] = True

    positions = [] #for each position in the read, calculate the actual genomic position
    if readIsFW:
        positions = [read.readstart + i for i in xrange(basesToAnalyze)]  # readStart plus index equals position of C if on + strand 
    else:
#TODO: leave this in till the nonCpG meth project is done. In this project all the regions have been generated using those modified coords
#        if seqMotifIsSymmetrical:
#            positions = [read.readend - i -1 - patLen for i in xrange(basesToAnalyze)] # since the - strand is given as the reverse complement, subtract the base positions (= index) from the readEnd (since it is reverse) and then subtract an additional len(seqMotif) (since it is a complement, too, e.g. if the seqMotif is actually a trimer subtract 2) and an additional -1 for counting in reverse direction
#        else:
#            positions = [read.readend - i -1 for i in xrange(basesToAnalyze)]
        positions = [read.readend - i - 1 for i in xrange(basesToAnalyze)]
    #remove an additional pattern length to end up on the genomic + strand coordinates of the start. i.e. coordinate of the T (+ strand) for CpA (+ strand motif TpG) methylation or C (+ strand) for CpG
    if readAmpRev:
        positions = map(lambda x: x-patLen,positions)
            
    #initialize and setting has much better runtime than appending
    blockStartList = internal_patterns_max * [0]
    methRatioInt   = internal_patterns_max * [0]
    
    isQCandIncluded = bitarray(basesToAnalyze)
    isQCandIncluded.setall(False)
    internal_patterns_qc = 0 #how many of the internal patterns actually pass the qualtiy criteria
    
    for i in xrange(basesToAnalyze):
        # only nucleotides with sufficient quality scores are used for further analysis
        if (read.qualityScoresPerBaseInRead[i] >= options.baseQualityScoreC) and not read.excludedBases[i]:
            if i > 0: readQualityScoreNextToC = options.baseQualityScoreNextToC
            else: readQualityScoreNextToC = max(BAD_QUALITY_SCORE,options.baseQualityScoreNextToC) #to exclude guessed reads with base-specific quality scores like: 20,1,1,1,1,1,1,1,...
            baseNextToCfitsQC = True
            try:
                indNext = i+1
                if read.account4beingRead2:
                    indNext = i-1
                baseNextToCfitsQC = read.qualityScoresPerBaseInRead[indNext] >= readQualityScoreNextToC
            except IndexError:
                #neglect the basequality next to C for the last base of a read
                #if (i + 1) == truncLength:
                baseNextToCfitsQC = False #set the quality criterion unfulfilled if the last base fits the quality criterion as the quality criterion for neighbor base cannot be determined
            
            if baseNextToCfitsQC:
                isQCandIncluded[i] = True
                if isPattern[i]:
                    if isInRegion[i]:
                        # the hashFragmentKey for hash table 'allCpgsOfSample' has the following format: chrom/CpGposition/strand
                        hashFragmentKey = "%s/%i" % (read.chrom, positions[i])                               
                        # increment totalCpG count for current CpG position in global (over all lanes) hashtable of unique CpGs  
                        try:
                            allCpgsOfSampleGlobal[hashFragmentKey][read.lane][read.ampStrand][1] +=1  
                        except KeyError:   #strand or lane or hashFragmentKey not in keys
                            try:
                                allCpgsOfSampleGlobal[hashFragmentKey][read.lane][read.ampStrand] = [0, 1]
                            except KeyError: #lane or hashFragmentKey not in keys
                                try:
                                    allCpgsOfSampleGlobal[hashFragmentKey][read.lane] = {read.ampStrand: [0,1]}
                                except KeyError: #hashFragmentKey not in allCpgsOfSampleGlobal
                                    allCpgsOfSampleGlobal[hashFragmentKey] = {read.lane:{read.ampStrand: [0,1]}}
                        # Bisulfite sequencing: unconverted Cs at investigated sites are assumed to be methylated --> increment methCpG count at current CpG position 
                        if not isCconv[i]:
                            allCpgsOfSampleGlobal[hashFragmentKey][read.lane][read.ampStrand][0] +=1   # hashCpgMethRatio[0] +=1 to count only methylated Cs
                    if not isCconv[i]:
                        methRatioInt[internal_patterns_qc] = 1
                    blockStartList[internal_patterns_qc] = positions[i]
                    internal_patterns_qc += 1
                    read.addMedianQualityScoreCpgFromQualityScoresPerBaseInRead(i)
            elif isPattern[i] and isInRegion[i]: # for quality score statistics of discarded CpGs
                read.addMedianQualityScoreDiscardedCpgFromQualityScoresPerBaseInRead(i)
    
    c_count = (isC & isQCandIncluded & isInRegion).count()
    if c_count > 0:
        statisticsLanes.updateXCytosineStatsVal_seqMotifStatus_posVal(c_count,seqMotif,read.lane,'total_c',options.isStrandSpecific,read.ampStrand,posVal='base')
    c_conv_count = (isCconv & isQCandIncluded & isInRegion).count()
    if c_conv_count > 0:
        statisticsLanes.updateXCytosineStatsVal_seqMotifStatus_posVal(c_conv_count,seqMotif,read.lane,'total_c_conv',options.isStrandSpecific,read.ampStrand,posVal='base')
    cpg_count = (isPattern & isQCandIncluded & isInRegion).count()
    if cpg_count > 0:
        statisticsLanes.updateXCytosineStatsVal_seqMotifStatus_posVal(cpg_count,seqMotif,read.lane,'total_cpg',options.isStrandSpecific,read.ampStrand,posVal='base')
    cpg_conv_count = (isPatternConv & isQCandIncluded & isInRegion).count()
    if cpg_conv_count > 0:
        statisticsLanes.updateXCytosineStatsVal_seqMotifStatus_posVal(cpg_conv_count,seqMotif,read.lane,'total_cpg_conv',options.isStrandSpecific,read.ampStrand,posVal='base')
    c_nocpg = (isNoCpg & isQCandIncluded)
    c_nocpg_count = c_nocpg.count()
    c_nocpg_count_region = (c_nocpg & isInRegion).count()
    if c_nocpg_count_region > 0:
        statisticsLanes.updateXCytosineStatsVal_seqMotifStatus_posVal(c_nocpg_count_region,seqMotif,read.lane,'total_c_nocpg',options.isStrandSpecific,read.ampStrand,posVal='base')
    c_nocpg_conv = (isNoCpgConv & isQCandIncluded)
    c_nocpg_conv_count = c_nocpg_conv.count()
    c_nocpg_conv_count_region = (c_nocpg_conv & isInRegion).count()
    if c_nocpg_conv_count_region > 0:
        statisticsLanes.updateXCytosineStatsVal_seqMotifStatus_posVal(c_nocpg_conv_count_region,seqMotif,read.lane,'total_c_nocpg_conv',options.isStrandSpecific,read.ampStrand,posVal='base')
    
    #first base
    if isQCandIncluded[0] and isInRegion[0]:
        if isC[0]:
            statisticsLanes.updateXCytosineStatsVal_seqMotifStatus_posVal(1,seqMotif,read.lane,'total_c',options.isStrandSpecific,read.ampStrand,posVal='firstBase')
        if isCconv[0]:
            statisticsLanes.updateXCytosineStatsVal_seqMotifStatus_posVal(1,seqMotif,read.lane,'total_c_conv',options.isStrandSpecific,read.ampStrand,posVal='firstBase')
        if isPattern[0]:
            statisticsLanes.updateXCytosineStatsVal_seqMotifStatus_posVal(1,seqMotif,read.lane,'total_cpg',options.isStrandSpecific,read.ampStrand,posVal='firstBase')
        if isPatternConv[0]:
            statisticsLanes.updateXCytosineStatsVal_seqMotifStatus_posVal(1,seqMotif,read.lane,'total_cpg_conv',options.isStrandSpecific,read.ampStrand,posVal='firstBase')
        if isNoCpg[0]:
            statisticsLanes.updateXCytosineStatsVal_seqMotifStatus_posVal(1,seqMotif,read.lane,'total_c_nocpg',options.isStrandSpecific,read.ampStrand,posVal='firstBase')
        if isNoCpgConv[0]:
            statisticsLanes.updateXCytosineStatsVal_seqMotifStatus_posVal(1,seqMotif,read.lane,'total_c_nocpg_conv',options.isStrandSpecific,read.ampStrand,posVal='firstBase')
    
    # for conversion rate per read                  
    read.updateNoCpgInfoX(c_nocpg_count,c_nocpg_conv_count)
    
    #remove the rest of the list that was estimated too much due to quality checks
    blockStartList = blockStartList[:internal_patterns_qc]
    methRatioInt = methRatioInt[:internal_patterns_qc]
    # invert order of methRatios since for the BED file the positions of the CpGs have to be given in ascending order 
    # but on the minus strand these values are obtained in a descending manner (since meth calls start at the readend), so they have to be inverted 
    # (that is why the "blockStartList" in this command block is also sorted afterwards)
    if not readIsFW: 
        if len(methRatioInt) > 0: methRatioInt = methRatioInt[::-1]
        if len(blockStartList) > 0: blockStartList = blockStartList[::-1]
    read.setMethRatioInt(methRatioInt)                       
                      
    read.setInternal_cpgs(internal_patterns_qc)
    read.setBlockStartList(blockStartList)
    read.setConversionRate4nonCpG()


def calcReadStatistics(read, seqMotif, statisticsLanes,options):
    '''
    update the statistics object with the statistics obtained from processing the read
    calculate the per read statistics for the read output files
    '''
    statisticsLanes.updateStatsVal(seqMotif,read.lane,options.isStrandSpecific,read.ampStrand,'readsPassedQc')
    if read.isDuplicateFlag:
        statisticsLanes.updateStatsVal(seqMotif,read.lane,options.isStrandSpecific,read.ampStrand,'processedDuplicatesFlag')
    if read.isDuplicatePosCount:
        statisticsLanes.updateStatsVal(seqMotif,read.lane,options.isStrandSpecific,read.ampStrand,'processedDuplicatesPosCount')
    if options.pairedEnd:        
        if read.isPaired:
            statisticsLanes.updateStatsVal(seqMotif,read.lane,options.isStrandSpecific,read.ampStrand,'pairedReads')
        if read.isProperPair:
            statisticsLanes.updateStatsVal(seqMotif,read.lane,options.isStrandSpecific,read.ampStrand,'properPairedReads')
        if read.mateFound == False:
            statisticsLanes.updateStatsVal(seqMotif,read.lane,options.isStrandSpecific,read.ampStrand,'mateNotFoundInRegion')
        if read.undigestedFragment:
            statisticsLanes.updateStatsVal(seqMotif,read.lane,options.isStrandSpecific,read.ampStrand,'pairedReadsUndigestedFragments')
        if read.overlapWithMate > 0:
            statisticsLanes.updateStatsVal(seqMotif,read.lane,options.isStrandSpecific,read.ampStrand,'overlappingPairedReads')
    
    read.calcStats()
    statisticsLanes.calcReadStats(read,seqMotif,options)


def seqMotif2pattern(seqMotif):
    '''
    transforms the seqMotif (trailing the C) into a list: G --> ['G"]; HpG --> ['H','G']
    '''
    pat = "".join(seqMotif.split("p"))
    return pat


def sequenceDataProcessing(genomicRegion, startPos, endPos, statsPickledFile, outfileNames, strands, seqMotif, fragmentsInfo, chromSizeHash,chromFileHandler,restrSite):
    '''
    retrieves all reads for the current region and iterates over them
    takes care of quality assessment (e.g.alignment mismatches, MspI start site, mate pair partner found, strand prediction, ...), invokes the methylationCalling for each read
    outputs Read,ReadDetails,CpG and FragmentCoverage informations (in the temp directory. Moreover a pickle file containing region specific statistics is written
    '''    
    readfileCollection_subp = ReadfileCollection(options.alignmentFiles)
    chromString = manualGenomeCorrection.get(genomicRegion,genomicRegion) #string for chromosome name e.g. chr21 or scaffold131289

    seqMotifContainsH = seqMotif.find("H") > -1
    seqMotifLength = len(seqMotif.replace('p',''))
    pattern = seqMotif2pattern(seqMotif)
    pattern_compl = strComp(pattern)

    #N#print "Processing genomic region "+genomicRegion +' '+ str(startPos)+' '+ str(endPos)
    #N#sys.stdout.flush()

    chromFile = chromFileHandler.getFileObject(chromString)
    
    (lb,ub) = readfileCollection_subp.getMinStartMaxEnd4region(genomicRegion, startPos, endPos)
    lb = max(0,lb-seqMotifLength)
    ub = min(chromSizeHash[genomicRegion],ub+seqMotifLength)
    refSeqHandler = ReferenceGenomeRegion(chromFile,lb,ub,genomicRegion)
    
    numberOfLanes = len(options.alignmentFiles)
    statisticsLanes = StatisticsCollection(range(numberOfLanes),[seqMotif],strands,fragmentsInfo)
    
    allCpgsOfSampleGlobal = {} # dictionary for the identified patterns/CpGs in all the reads for the region. first index: hashKey (made up of chromosome and position of the C) second index: lane, third index: strand;  set in methylationCalling_bitarray()
    listCpgMethSummary = [] 
    fragmentCoverage = {}  # hashtable for fragment coverage summary, global (process-specific) version (no lane-specific version necessary)
#    tempReadUnit = bsAlignedReadCollection() # container storing the reads for writing to files later
    
    #handler for writing read file output
    if not options.statisticsOnly:
        if options.readDetailsFile:
            readFileH = ReadOutputFileHandler(outfileNames['readOutfile'],options.readOutputBuffer,seqMotif,detailsFile=outfileNames['readQualOutfile'])
        else:
            readFileH = ReadOutputFileHandler(outfileNames['readOutfile'],options.readOutputBuffer,seqMotif,detailsFile=None)
    
    countDiscardedAlignments = 0
    
    regionWindows = Windows4region(startPos,endPos,numberOfLanes,size=50)
    readPositions = ReadPositions4region(startPos,endPos,numberOfLanes)
    readPositionsSumLanes = ReadPositions4region(startPos,endPos,1)#one summary that captures the combined positions of all lanes
    #+/- strand specific: used for statistics and duplicate marking
    readPositionsSumLanesP = ReadPositions4region(startPos,endPos,1)
    readPositionsSumLanesM = ReadPositions4region(startPos,endPos,1)
    patternBuckets = PatternBuckets4region(startPos,endPos,numberOfLanes,countBuckets=SEQMOTIF_THRES)
    if options.isStrandSpecific:
        regionWindowsM = Windows4region(startPos,endPos,numberOfLanes,size=50)
        regionWindowsP = Windows4region(startPos,endPos,numberOfLanes,size=50)
        readPositionsM = ReadPositions4region(startPos,endPos,numberOfLanes)
        patternBucketsM = PatternBuckets4region(startPos,endPos,numberOfLanes,countBuckets=SEQMOTIF_THRES)
    
    if debug:
        print "fetching reads..."
        sys.stdout.flush()
    try:
        pysamReadList = readfileCollection_subp.getReadsInRegion_withFileIndex(genomicRegion, startPos, endPos)
    except ValueError:
        print '\tChromosome ' + genomicRegion + "[" + str(startPos) + "," + str(endPos) + ') does not exist in current alignment file'
        sys.stdout.flush()
    if len(pysamReadList) == 0:
        print 'No reads found for Chromosome ' + genomicRegion + "[" + str(startPos) + "," + str(endPos) + ')'

    count = 0
    processedReads = 0
    unalignedCount = 0
    if debug:
        if count == 0: print "Non-overlapping unit borders: "+str(startPos)+'  '+str(endPos)

    #build a mate pair dictionary if Paired end handling is enabled
    if options.pairedEnd:
        mateReadHash = PysamMateHash(pysamReadList,len(readfileCollection_subp.readfileNames))
        if ub > endPos:#extend the mate dictionary to the end of the last read in the region to cover all overlapping pairs
            readListExtension = readfileCollection_subp.getReadsInRegion_withFileIndex(genomicRegion, endPos, ub) #all reads covering the bases of the last reads extending over the region boundary
            mateReadHash.addReads(readListExtension) #extend the mate lookup table
        
    if debug:
        print "processing reads |"+genomicRegion + "[" + str(startPos) + "," + str(endPos) + ')'+"| ..."
    #N#print "number of reads: " + str(len(pysamReadList))
    sys.stdout.flush()
    
#DEBUG
    i = 0

    for readAndFileIndex in pysamReadList:
        (alignedRead,inputFileNo) = readAndFileIndex
        sys.stdout.flush()
        if count % 100000 == 0:
            #N#print "Processing " + genomicRegion + ":" + str(startPos) + " to " + str(endPos) + " record "+str(count)
            if options.debugEnableMemoryProfiling:
                h = guppy.hpy()
                print "Local time: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
#                print "MEMORY USAGE (heap): " + genomicRegion + ":" + str(startPos) + " to " + str(endPos) + " record "+str(count) + " (" + str(len(tempReadUnit.reads)) + " in read collection)"
                print "MEMORY USAGE (heap): " + genomicRegion + ":" + str(startPos) + " to " + str(endPos) + " record "+str(count) + " (" + str(processedReads) + " reads accepted)"
                print "Heap: "
                print h.heap()
            sys.stdout.flush()
        count = count + 1
        if alignedRead.is_unmapped <> 0:
            countDiscardedAlignments += 1
            continue
        
        mappedStrand = "+"
        if alignedRead.is_reverse == 1: mappedStrand = "-"
        ampStrand = mappedStrand
        
        if options.pairedEnd and alignedRead.is_paired:
            if options.processZSField: #use the custom ZS field to determine the amplification strand ["+","-"]
                ampStrand = alignedRead.opt("ZS")[0]
            elif alignedRead.is_read2: #use the flag field to determine the amplified strand: revert for read#2
                if mappedStrand == "+":
                    ampStrand = "-"
                elif mappedStrand == "-":
                    ampStrand = "+"

        curReadStartsInWindow = alignedRead.pos >= startPos and alignedRead.pos < endPos
        
        if curReadStartsInWindow:
            statisticsLanes.updateStatsVal(seqMotif,inputFileNo,options.isStrandSpecific,ampStrand,'totalMappedReads')
            if alignedRead.is_paired:
                statisticsLanes.updateStatsVal(seqMotif,inputFileNo,options.isStrandSpecific,ampStrand,'totalPairedReads')
            if alignedRead.is_proper_pair:
                statisticsLanes.updateStatsVal(seqMotif,inputFileNo,options.isStrandSpecific,ampStrand,'totalProperPairedReads')
              
        if options.pfStatus == 'PF':
            if alignedRead.is_qcfail:
                if curReadStartsInWindow:
                    statisticsLanes.updateStatsVal(seqMotif,inputFileNo,options.isStrandSpecific,ampStrand,'discard_PFstatus')
                countDiscardedAlignments += 1
                continue
        elif options.pfStatus == 'NonPF':
            if not alignedRead.is_qcfail:
                if curReadStartsInWindow:
                    statisticsLanes.updateStatsVal(seqMotif,inputFileNo,options.isStrandSpecific,ampStrand,'discard_PFstatus')
                countDiscardedAlignments += 1
                continue

        lane = inputFileNo
        read = bsAlignedRead(alignedRead,mappedStrand,lane,genomicRegion,refSeqHandler,options,seqMotif,addH=seqMotifContainsH,ampStrandI=ampStrand,statsObj=statisticsLanes,startsInRegion2Process=curReadStartsInWindow)
        
        #assuming the readnumber corresponds to the order of mapping positions, determine the number of reads where  read#1 does not map to the + strand
        #hints at sequencing into the adapter
        if options.pairedEnd and options.processZSField and read.isPaired and curReadStartsInWindow:
            if (mappedStrand == "+" and read.pairNum <> 1) or (mappedStrand == "-" and read.pairNum <> 2):
                statisticsLanes.updateStatsVal(seqMotif,inputFileNo,options.isStrandSpecific,ampStrand,'mappingOrderDoesNotMatchStrand')
                
        if read.isDiscarded: #because of mismatch criteria or the reference sequence could not be retrieved or there was some error in parsing
            countDiscardedAlignments += 1
##DEBUG
#            print 20*"@" +"discarded due to mismatchLenOrReference"+ 20*"@"
#            print "pysam read:"
#            printread(alignedRead)
#            print "parsed read:"
#            print read.getDebugString()
            continue 
         
        #just an assertion
        if read.chromend <= read.chromstart or read.readend <= read.readstart:
            print "WARNING: improper readlength detected:"
            print "readend : " + str(read.readend)
            print "chromend: " + str(read.chromend)
            print "readstart : " + str(read.readstart)
            print "chromstart: " + str(read.chromstart)
            printread(alignedRead)
            
            countDiscardedAlignments += 1
            continue 
        
        #check restriction site
        if options.checkRestriction:
            if not restrSite.readPassesRestrictionPattern(read) or not read.coordsFitFragmentPlacement():
                read.discardMe()
                countDiscardedAlignments += 1
                if curReadStartsInWindow:
                    statisticsLanes.updateStatsVal(seqMotif,inputFileNo,options.isStrandSpecific,ampStrand,'discard_restrSite')
##DEBUG
#                print 20*"@" +"discarded due to restriction site"+ 20*"@"
#                print "pysam read:"
#                printread(alignedRead)
#                print "parsed read:"
#                print read.getDebugString()
                
                continue
            
        #check fragment length
        if options.rrbsMode and (read.fragmentLength < options.minFragmentLength or read.fragmentLength > options.maxFragmentLength):
            read.discardMe()
            countDiscardedAlignments += 1
            if curReadStartsInWindow:
                statisticsLanes.updateStatsVal(seqMotif,inputFileNo,options.isStrandSpecific,ampStrand,'discard_fragmentsOutsideSizeRange')                
            continue
        
        #strand prediction
        if options.strandPredict:
            read.setAgreesWithStrandPrediction(strandModel)
            if not read.strandAgreesWithPrediction:
                if curReadStartsInWindow:
                    statisticsLanes.updateStatsVal(seqMotif,read.lane,options.isStrandSpecific,read.ampStrand,'discard_disagreeWithStrandPrediction')
                countDiscardedAlignments += 1
                read.discardMe()
                continue

        if options.pairedEnd and read.isPaired:
            if read.isProperPair:
                alignedMate = mateReadHash.getReadfromPair(read.qname, read.getMateNum(), inputFileNo)
                if not alignedMate:
                    read.setMateFound(False)
                else:
                    mate = bsAlignedReadMinimal(alignedMate,seqMotifLength)
                    read.setMateFound(True)
                    read.setMateInfo(mate)
            else:
                if options.discardImproperPairs:
                    read.discardMe()
                    countDiscardedAlignments += 1
                    if curReadStartsInWindow:
                        statisticsLanes.updateStatsVal(seqMotif,inputFileNo,options.isStrandSpecific,ampStrand,'discard_improperPairs')
                    continue
                #only take read #1 into account for improperly paired reads
                elif read.pairNum <> 1:                    
                    read.discardMe()
                    countDiscardedAlignments += 1
                    if curReadStartsInWindow:
                        statisticsLanes.updateStatsVal(seqMotif,inputFileNo,options.isStrandSpecific,ampStrand,'discard_improperPairs')
                    continue
        
        manualDuplicateMark = False # for later duplicate marking. needs to be set before the position lookup table is updated
        if ampStrand == "-":
            manualDuplicateMark = readPositionsSumLanesM.containsPosition_nonFlagged(read.readstart,read.readend,0)
        else:
            manualDuplicateMark = readPositionsSumLanesP.containsPosition_nonFlagged(read.readstart,read.readend,0)
        read.setIsDuplicate(alignedRead.is_duplicate,manualDuplicateMark)
        #update position statistics (also include duplicate reads and reads with no CpGs/patterns)

        readPositions.addReadStats(read,lane)
        readPositionsSumLanes.addReadStats(read,0)
        
        if options.isStrandSpecific and ampStrand == "-":
            readPositionsM.addReadStats(read,lane)
                
        if ampStrand == "-":
            readPositionsSumLanesM.addReadStats(read,0)
        else:
            readPositionsSumLanesP.addReadStats(read,0)
            
            
        #duplicate read handling AFTER quality parsing
        if options.removeDuplicates and (read.isDuplicateFlag or read.isDuplicatePosCount):
            if curReadStartsInWindow:
                if alignedRead.is_duplicate:
                    statisticsLanes.updateStatsVal(seqMotif,inputFileNo,options.isStrandSpecific,ampStrand,'discard_duplicatesFlag')
                if manualDuplicateMark:
                    statisticsLanes.updateStatsVal(seqMotif,inputFileNo,options.isStrandSpecific,ampStrand,'discard_duplicatesPosCount')
            read.discardMe()
            countDiscardedAlignments += 1
            continue
        
        #update sliding window statistics (also include reads with no CpGs/patterns)
        #here, reads that dont start within this window are also counted
        regionWindows.addReadStats(read,lane)
        if options.isStrandSpecific:
            if ampStrand == "-":
                regionWindowsM.addReadStats(read,lane)
            else:
                regionWindowsP.addReadStats(read,lane)
        
        #set bases not taken into account due to blunt end repair/base fill in with methylated Cs        
        if options.rrbsMode:
            read.setBluntEndRepairNeglectedBases(restrSite)
            
        # START METHYLATION CALLING      
        methylationCalling_bitarray(read, seqMotif, pattern, pattern_compl, statisticsLanes, allCpgsOfSampleGlobal,startPos,endPos)
              
        #update read CpG count statistics (also include reads with no CpGs/patterns)
        if curReadStartsInWindow:
            patternBuckets.addReadStats(read,lane)
            if options.isStrandSpecific and ampStrand == "-":
                patternBucketsM.addReadStats(read,lane)
        if read.internal_cpgs < 1: # due to minimum pattern requirements
            read.discardMe()
            countDiscardedAlignments += 1
            if curReadStartsInWindow:
                statisticsLanes.updateStatsVal(seqMotif,inputFileNo,options.isStrandSpecific,ampStrand,'discard_noSeqMotif')
            continue
        else: 
            processedReads += 1
#            tempReadUnit.addRead(read)
 
            if curReadStartsInWindow:
                calcReadStatistics(read, seqMotif, statisticsLanes,options)
                # store the read's mean meth in hash table 'fragmentCoverage' for the (unique) fragment's read-coverage summary
                methRatioListPerFragment = [] # list to append all read.meanMeth values for one unique fragment
                plusMinusStrandPerFragment = {}
                uniqueFragmentKey = "%s:%i:%i" % (read.chrom, read.chromstart, read.chromend)
                try:
                    fragmentCoverage[uniqueFragmentKey][read.ampStrand].append(read.meanMeth)
                except KeyError:
                    try:
                        fragmentCoverage[uniqueFragmentKey][read.ampStrand] = methRatioListPerFragment  
                        fragmentCoverage[uniqueFragmentKey][read.ampStrand].append(read.meanMeth)
                    except KeyError:
                        plusMinusStrandPerFragment[read.ampStrand] = methRatioListPerFragment
                        fragmentCoverage[uniqueFragmentKey] = plusMinusStrandPerFragment       
                        fragmentCoverage[uniqueFragmentKey][read.ampStrand].append(read.meanMeth)
                
                if not options.statisticsOnly:
                    #write read to read output files:
                    readFileH.addRead(read)
        

    if debug:
        print "updating CpG methylation |"+genomicRegion + "[" + str(startPos) + "," + str(endPos) + ')'+"| ..."
        sys.stdout.flush()
    del readfileCollection_subp
    #N#print "[" + str(time.strftime("%Y-%m-%d %H:%M:%S")) + "] Region " + genomicRegion + " [" + str(startPos) + "," + str(endPos) + "): " + str(unalignedCount) + " of " + str(count) + " reads unaligned, " + str(countDiscardedAlignments) + " discarded"
    #N#print "[" + str(time.strftime("%Y-%m-%d %H:%M:%S")) + "] Region " + genomicRegion + " [" + str(startPos) + "," + str(endPos) + "): updating statistics..."
    sys.stdout.flush()
        
    # create statistics for CpG methylation summary, list each unique CpG once (per lane and accumulated (not summed!) over all lanes in the region)    
    for hashFragmentKey in allCpgsOfSampleGlobal.keys():        
        # for region specific statistics
        uniqueMethAcrossLanes = {"both": 0.0, "-": 0.0, "+":0.0}

        #update CpG/seqMotif specific statistics
        if options.isStrandSpecific:
            methTupleP = [0, 0]
            methTupleM = [0, 0]
            for lane in allCpgsOfSampleGlobal[hashFragmentKey]: # should the case be handled when the denominator becomes zero?
                for strand in allCpgsOfSampleGlobal[hashFragmentKey][lane].keys():
                    curVal = (float(allCpgsOfSampleGlobal[hashFragmentKey][lane][strand][0]) / allCpgsOfSampleGlobal[hashFragmentKey][lane][strand][1])*100.0 #allCpgsOfSampleGlobal[hashFragmentKey][read.lane][0(mseqMotif),0(totalseqMotif)] 
                    statisticsLanes.updateUniqueCpgMeth(seqMotif, lane, strand, curVal)
                    statisticsLanes.updateUniqueCpgMeth(seqMotif, lane, "both", curVal)
                    if strand == "+":
                        methTupleP[0] += allCpgsOfSampleGlobal[hashFragmentKey][lane]['+'][0]
                        methTupleP[1] += allCpgsOfSampleGlobal[hashFragmentKey][lane]['+'][1]
                    else:
                        methTupleM[0] += allCpgsOfSampleGlobal[hashFragmentKey][lane]['-'][0]
                        methTupleM[1] += allCpgsOfSampleGlobal[hashFragmentKey][lane]['-'][1]            
            #for region specific statistics
            if methTupleP[1] > 0:
                uniqueMethAcrossLanes = float(methTupleP[0])/methTupleP[1]*100.0
                statisticsLanes.addOneRegionStats(seqMotif,"both","uniqueseqMotifCount")
                statisticsLanes.updateRegionStats4process(seqMotif, "both", {"uniqueseqMotifMethSum": uniqueMethAcrossLanes,"uniqueseqMotifMethSumSquares": pow(uniqueMethAcrossLanes,2)})
                listCpgMethSummary.append(UniqueCpgsPerLane(hashFragmentKey,"+",methTupleP[0],methTupleP[1],seqMotif))
                statisticsLanes.addOneRegionStats(seqMotif,"+","uniqueseqMotifCount")
                statisticsLanes.updateRegionStats4process(seqMotif, "+", {"uniqueseqMotifMethSum": uniqueMethAcrossLanes,"uniqueseqMotifMethSumSquares": pow(uniqueMethAcrossLanes,2)})
            if methTupleM[1] > 0:
                uniqueMethAcrossLanes = float(methTupleM[0])/methTupleM[1]*100.0
                statisticsLanes.addOneRegionStats(seqMotif,"both","uniqueseqMotifCount")
                statisticsLanes.updateRegionStats4process(seqMotif, "both", {"uniqueseqMotifMethSum": uniqueMethAcrossLanes,"uniqueseqMotifMethSumSquares": pow(uniqueMethAcrossLanes,2)})
                listCpgMethSummary.append(UniqueCpgsPerLane(hashFragmentKey,"-",methTupleM[0],methTupleM[1],seqMotif))
                statisticsLanes.addOneRegionStats(seqMotif,"-","uniqueseqMotifCount")
                statisticsLanes.updateRegionStats4process(seqMotif, "-", {"uniqueseqMotifMethSum": uniqueMethAcrossLanes,"uniqueseqMotifMethSumSquares": pow(uniqueMethAcrossLanes,2)})
        else:
            methTuple = [0, 0]
            seenPlus = False
            seenMinus = False
            for lane in allCpgsOfSampleGlobal[hashFragmentKey]: # should the case be handled when the denominator becomes zero?
                #if the position has already been processed on one strand, dont process it again on the second strand
                addedMethRatios = [0, 0]
                #add the methylation values from both strands
                for strand in allCpgsOfSampleGlobal[hashFragmentKey][lane].keys():
                    if strand == "-":
                        seenMinus = True
                    else:
                        seenPlus = True
                    addedMethRatios[0] += allCpgsOfSampleGlobal[hashFragmentKey][lane][strand][0]
                    addedMethRatios[1] += allCpgsOfSampleGlobal[hashFragmentKey][lane][strand][1]
                methTuple[0] += addedMethRatios[0]
                methTuple[1] += addedMethRatios[1]
                curVal = (float(addedMethRatios[0]) / addedMethRatios[1])*100
                statisticsLanes.updateUniqueCpgMeth(seqMotif, lane, "both", curVal)
            uniqueMethAcrossLanes["both"] += float(methTuple[0])/methTuple[1]*100.0
            
            # if same CpG was sequenced on both strands (write to output file only as + strand)
            if seenPlus: #only detected reads on the + exclusively or on both strands --> write a + to the output
                listCpgMethSummary.append(UniqueCpgsPerLane(hashFragmentKey,"+", methTuple[0], methTuple[1], seqMotif))
            elif seenMinus: #only detected reads on the - strand for this position
                listCpgMethSummary.append(UniqueCpgsPerLane(hashFragmentKey,"-", methTuple[0], methTuple[1], seqMotif)) 
            #for region specific statistics
            statisticsLanes.addOneRegionStats(seqMotif,"both","uniqueseqMotifCount")
            statisticsLanes.updateRegionStats4process(seqMotif, "both", {"uniqueseqMotifMethSum": uniqueMethAcrossLanes["both"],"uniqueseqMotifMethSumSquares": pow(uniqueMethAcrossLanes["both"],2)})

    #N#print "[" + str(time.strftime("%Y-%m-%d %H:%M:%S")) + "] Region " + genomicRegion + " [" + str(startPos) + "," + str(endPos) + "): creating subprocess output..."
    sys.stdout.flush()
        
    # create process specific output files (output files of all processes will be merged in main process once all processes have finished)
    if not options.statisticsOnly:
        # open all output files
        for outFileName in outfileNames.keys():
            if not outFileName in ['readOutfile','readQualOutfile']: # skip the read output files as they are processed in the handler object
                outfile = open(outfileNames[outFileName], 'w')
                outfileNames[outFileName] = outfile
             
#        print str(tempReadUnit.getNumReads()) + ' of ' + str(len(pysamReadList)) +' reads processed for Chromosome ' + genomicRegion + "[" + str(startPos) + "," + str(endPos) + ')'
        #N#print str(processedReads) + ' of ' + str(len(pysamReadList)) +' reads processed for Chromosome ' + genomicRegion + "[" + str(startPos) + "," + str(endPos) + ')'
        
#        # write read and possibly read details output file
#        if debug:
#            print "writing output file: Reads |"+genomicRegion + "[" + str(startPos) + "," + str(endPos) + ')'+"| ..."
#            print "  sorting reads..."
#            sys.stdout.flush()
#        tempReadUnit.sort()
#        if options.readDetailsFile:
#            tempReadUnit.writeReadDetailsFile(outfileNames['readQualOutfile'])
#        tempReadUnit.updateChromStartEnds()
#        
#        tempReadUnit.writeReadOutFile(seqMotif,outfileNames['readOutfile'])
        
        #write the rest of the reads in the buffer to the output files
        readFileH.flush()

        #N#print "sorting reads in read output file: Reads |"+genomicRegion + "[" + str(startPos) + "," + str(endPos) + ')'+"| ..."
        sys.stdout.flush()
        readFileH.sort()
        
        # create output file for CpG methylation summary, list each CpG once. (if command line parameter "--isStrandSpecific" is set, the same CpG will appear separately for "+" and "-" strand --> 2 output files)          
        if debug:
            print "writing output file: CpGs |"+genomicRegion + "[" + str(startPos) + "," + str(endPos) + ')'+"| ..."
            sys.stdout.flush()
        listCpgMethSummary.sort()
        for cpgMethSummary in listCpgMethSummary:
            if options.isStrandSpecific:
                if cpgMethSummary.strand == '-': outfileNames['cpgOutfileM'].write(str(cpgMethSummary))
                else: outfileNames['cpgOutfileP'].write(str(cpgMethSummary))
            else: outfileNames['cpgOutfile'].write(str(cpgMethSummary))  
    
        if debug:
            print "assessing fragment coverage |"+genomicRegion + "[" + str(startPos) + "," + str(endPos) + ')'+"| ..."
        # create global fragment coverage summary file - each fragment appears only once, reporting the number of times the fragment is covered by reads (if command line parameter "options.isStrandSpecific" is set, the same fragment will appear separately for "+" and "-" strand --> 2 output files)                
        coverageSortingList = []
        for fragment in fragmentCoverage.keys():
            dnaStrand = fragmentCoverage[fragment]    
            if not options.isStrandSpecific:
                methList = []
                strandCount = 0
                for strandInfo in dnaStrand.keys():  # for + or - strand
                    for ratio in fragmentCoverage[fragment][strandInfo]:
                        methList.append(ratio)
                    if strandCount > 0:
                        strandInfo = '+' # if same fragment was sequenced on both strands (write to output file as only +)
                    strandCount += 1                 
                coverageSortingList.append(UniqueFragmentsPerLane(fragment, methList, strandInfo))
            else: 
                for strandInfo in dnaStrand.keys():  # for + or - strand
                    methList = fragmentCoverage[fragment][strandInfo]  
                    coverageSortingList.append(UniqueFragmentsPerLane(fragment, methList, strandInfo))
        coverageSortingList.sort()  
        if debug:
            print "writing output file: fragments |"+genomicRegion + "[" + str(startPos) + "," + str(endPos) + ')'+"| ..."
            sys.stdout.flush()
        for frag in coverageSortingList:      
            if options.isStrandSpecific:        
                if frag.strand == '-': outfileNames['fragmentOutfileM'].write(str(frag)) 
                else: outfileNames['fragmentOutfileP'].write(str(frag))
            else: outfileNames['fragmentOutfile'].write(str(frag))    
    # close all files of process
    if not options.statisticsOnly:
        for filename in outfileNames.keys():
            if not filename in ['readOutfile','readQualOutfile']: # skip the read output files as they are processed in the handler object
                outfileNames[filename].close()    
    chromFile.close()

    statisticsLanes.setQ1Q3Median(seqMotif,debuginfo='@Chromosome ' + genomicRegion + "[" + str(startPos) + "," + str(endPos) + ')')
    
    #threshold specific statistics windowsWithAtLeastXreads, coordsWithAtLeastXreads, readsWithAtLeastXseqMotifs
    for laneNum in range(numberOfLanes):
        #update combined strands
        regionWindowThresCounts = regionWindows.getReadCounts4thres(READ_THRES_WINDOWS,laneNum)
        readPositionsThresCounts = readPositions.getReadCounts4thres(READ_THRES_POSITIONS,laneNum,startPos=startPos,endPos=endPos)
        for t in READ_THRES_WINDOWS:
            statisticsLanes.setLaneStatsVal(seqMotif,laneNum,"both","windowsWithAtLeast"+str(t)+"reads",regionWindowThresCounts[t],isStrandSpecificStat=True)
        for t in READ_THRES_POSITIONS:
            statisticsLanes.setLaneStatsVal(seqMotif,laneNum,"both","coordsWithAtLeast"+str(t)+"reads",readPositionsThresCounts[t])
        for t in SEQMOTIF_THRES:
            statisticsLanes.setLaneStatsVal(seqMotif,laneNum,"both","readsWithAtLeast"+str(t)+"seqMotifs",patternBuckets.readCounts[laneNum][t])
        if options.isStrandSpecific:
            #update the - strand, the + strand will get set later using the information of 'both' and '-'
            regionWindowThresCountsM = regionWindowsM.getReadCounts4thres(READ_THRES_WINDOWS,laneNum)
            regionWindowThresCountsP = regionWindowsP.getReadCounts4thres(READ_THRES_WINDOWS,laneNum)
            readPositionsThresCountsM = readPositionsM.getReadCounts4thres(READ_THRES_POSITIONS,laneNum,startPos=startPos,endPos=endPos)
            for t in READ_THRES_WINDOWS:
                statisticsLanes.setLaneStatsVal(seqMotif,laneNum,"-","windowsWithAtLeast"+str(t)+"reads",regionWindowThresCountsM[t],isStrandSpecificStat=True)
                statisticsLanes.setLaneStatsVal(seqMotif,laneNum,"+","windowsWithAtLeast"+str(t)+"reads",regionWindowThresCountsP[t],isStrandSpecificStat=True)
            for t in READ_THRES_POSITIONS:
                statisticsLanes.setLaneStatsVal(seqMotif,laneNum,"-","coordsWithAtLeast"+str(t)+"reads",readPositionsThresCountsM[t])
            for t in SEQMOTIF_THRES:
                statisticsLanes.setLaneStatsVal(seqMotif,laneNum,"-","readsWithAtLeast"+str(t)+"seqMotifs",patternBucketsM.readCounts[laneNum][t])
    
    readPositionsSumLanesThresCounts = readPositionsSumLanes.getReadCounts4thres(READ_THRES_POSITIONS,0,startPos=startPos,endPos=endPos)
    if options.isStrandSpecific:
        readPositionsSumLanesThresCountsP = readPositionsSumLanesP.getReadCounts4thres(READ_THRES_POSITIONS,0,startPos=startPos,endPos=endPos)
        readPositionsSumLanesThresCountsM = readPositionsSumLanesM.getReadCounts4thres(READ_THRES_POSITIONS,0,startPos=startPos,endPos=endPos)
    for t in READ_THRES_POSITIONS:
        statisticsLanes.setRegionStats(seqMotif,"both","coordsWithAtLeast"+str(t)+"reads",readPositionsSumLanesThresCounts[t])
        if options.isStrandSpecific:
            statisticsLanes.setRegionStats(seqMotif,"+","coordsWithAtLeast"+str(t)+"reads",readPositionsSumLanesThresCountsP[t])
            statisticsLanes.setRegionStats(seqMotif,"-","coordsWithAtLeast"+str(t)+"reads",readPositionsSumLanesThresCountsM[t])
    
    #  pickle the statistics object for the 'wrapup' task
    #pickle to tmp file and move to final file name later for easier error checking
    if debug:
        print "writing output file: Pickle |"+genomicRegion + "[" + str(startPos) + "," + str(endPos) + ')'+"| ..."
        sys.stdout.flush()
    output = open(statsPickledFile + ".tmp", 'wb')
    cPickle.dump(statisticsLanes, output)
    cPickle.dump(exceptionOccurred, output)
    output.close()
    if debug:
        print "  moving pickle file |"+genomicRegion + "[" + str(startPos) + "," + str(endPos) + ')'+"| ..."
        sys.stdout.flush()
#    os.system("cp " + statsPickledFile + ".tmp " + statsPickledFile) #TODO: comment in for debug for WGBSS: strangely some tmp files were not moved
    os.system("mv " + statsPickledFile + ".tmp " + statsPickledFile)

def subprocProfileWrapper_sequenceDataProcessing(genomicRegion, startPos, endPos, queue, outfileNames, strands, seqMotif, path, fragmentsInfo, debugQueue, chromSizeHash,chromFileHandler):
    '''
    wrapper for cProfiling the sequenceDataProcessing function. needs to maintain the same signature
    '''
    cProfile.runctx("sequenceDataProcessing(genomicRegion, startPos, endPos, queue, outfileNames, strands, seqMotif, fragmentsInfo, debugQueue, chromSizeHash, chromFileHandler)",globals(), locals(),options.outputDir + os.sep + "runtime_profile_" + genomicRegion +'_'+ str(startPos) + "to" + str(endPos) + ".cProfile")
    
#def getLastReadEnd4region(readList):
#    '''
#    given a list of pysam reads, returns the last end position on the reference
#    '''
#    endPos = []
#    for r in readList:
#        lenOnRef = getLengthOnReferenceFromCIGAR(r)
#        if lenOnRef == 0:
#            continue
#        endPos.append(r.pos + lenOnRef)
#    if len(endPos) < 1:
#        return -1
#    return(max(endPos))
#
#def findNextGap(pos,chrom,chromSize,readfileCollection):
#    '''
#    Sweeping algorithm to
#    deterministically find the next gap position downstream of a genomic position (on a given chrom/scaffold)
#    '''
#    
#    if debug:
#        print "trying unit boundary: " + str(pos)
#     
#    runPos = pos
#    minGapSize = 1 # minimum size of the gap
#    while runPos < chromSize:
#        readList = readfileCollection.getReadsInRegion(chrom,runPos,runPos+minGapSize)
#        # the position does not hit any reads. 
#        if len(readList) == 0:
#            break
#        #jump to the last end position of all reads in the window
#        else:
#            prevPos = runPos
#            runPos = getLastReadEnd4region(readList)
#            if prevPos == runPos:
#                raise Exception("Encountered endless loop in finding the next window gap")
#            if debug:
#                print "  --> " + str(len(readList)) + " reads found in region --> trying unit boundary: " + str(runPos)
#            if runPos == -1:
#                raise Exception("could not compute last read end position for chrom " + chrom + " start position " + str(pos) )
#        del readList #try to remove, because the loops seems to consume too much memory. is this really necassary????
#            
#    #ran across end of the chromosome/genomic region    
#    if runPos > chromSize:
#        runPos = chromSize
#    
#    return runPos

    
def checkRestart4region(outfileNames,statsPickledFile):
    '''
    for a given output file structure and a pickle file name, this function checks whether the pickle file is valid. If not, it removes all files associated with the given output filename structure
    '''
    if isValidFile(statsPickledFile):
        return False
    #delete existing files
    for outFileName in outfileNames.keys():
        if isValidFile(outfileNames[outFileName]):
            os.system("rm " + outfileNames[outFileName])
    return True

def cleanupChromName(chromName,badChars=["(",")"]):
    '''
    escape weird characters in a chromosome name in order to properly process them later
    (for example "(CG)n" in repeat pseudogenomes which downstream cause trouble in file naming) 
    '''
    cn = chromName
    for c in badChars:
        cn = cn.replace(c,"")
    return cn

def manageProcesses(chromSortList,chromSizeHash,chromFileHandler,includedRegions,outfileH,shellScriptHs,jobH,inPKLfn,strands,seqMotifs,path,options,fragmentsInfo,debug):
        '''
        for the 'main' task:
        this function is responsible for chopping up the reference genome into regions for the 'region' tasks
        it discards regions not included in the analysis and collects the region tasks in the shellScriptHs (ProcessShellScriptHandler) and jobH (JobHandler) objects
        in a region the jobs are created for each seMotif under investigation
        also writes an ouput file containing information about the windows processed: window coordinates, estimated number of reads, isRegionProcessed
        if enabled, the number of reads per region are estimated and window sizes adjusted accordingly        
        '''
        nextChromNo = 0
        numChromUnits = None #theoretically maximum number of units of size options.genomeFraction * ONE_MB that are present
        processUnit = None #parts in which a chromosome is devided for parallel processing purposes
        curChromPos = 0 # current chromosome position for the current processUnit
        chrom = None
        
        extendWindowBy = int(options.genomeFraction * ONE_MB)
        minWindowSize = int(options.minSWlength * ONE_MB)

        for seqMotif in seqMotifList:
            outfileH.initializeseqMotifFiles2cat(seqMotif)

        if debug:
            windowStatsF = open(path + "window_statistics.txt","w")
            windowStatsF.write("chrom\tstart\tend\tprocessUnit\testimatedReads\tprocessed\n")
            
        #create a directory for the pickle files
        picklePath = path + "pickle" + os.sep
        if not os.path.exists(picklePath): os.mkdir(picklePath)

        numProcesses = 0 # total number of submitted processes
        numRegions = 0
        
        readfileCollection = ReadfileCollection(options.alignmentFiles)
        processZSField = None
        #check for custom ZS field in the bam files to determine the read numbers and orientations
        if options.pairedEnd:
            print "checking for ZS field in BAM files..."
            processZSField = readfileCollection.checkForZSfield()
            if processZSField:
                print "  all files contain ZS field --> using ZS to determine read numbers and orientation"
            else:
                print "  not all files contain ZS field --> using flag field to determine read numbers and orientation"      
            
        runtimeProfileDir = options.outputDir + os.sep + "runtime_profiles"

        while nextChromNo <= len(chromSortList):
            if debug: print "Processing nextChromNo: " + str(nextChromNo)
            # identify next valid processUnit and start process
            if (processUnit == None or chrom == None or curChromPos >= chromSizeHash[chrom]):
                #only get next chromosome info if not working on the last chromosome
                if nextChromNo < len(chromSortList):
                    # if out of chromUnits: proceed to next chromosome
                    chrom = chromSortList[nextChromNo]
                    nextChromNo += 1
                    # determine single process genomic unit borders for parallel analysis
                    numChromUnits = int(math.ceil(float(chromSizeHash[chrom]) / float(options.genomeFraction * ONE_MB)))
                    #N#print "estimated number of chromUnits for " + chrom + ": " + str(numChromUnits)
                    processUnit = 0
                    curChromPos = 0
                else:
                    nextChromNo += 1 # do not provoke infinite loop
            #once the last region of the last chromosome is done no jobs need to be submitted
            if curChromPos < chromSizeHash[chrom]:
                #determine the genomic window for the current process
                windowBoundary2try = curChromPos + extendWindowBy#the search for the next boundary should start here
                windowEstimatedReadCount = -1
                if options.smartWindows > 0:
                    smartWindowSize = windowBoundary2try - curChromPos
                    #set the window boundary to half the distance as long as the estimated readnumber is larger than the user specified threshold
                    #as an extra precaution, the number of iterations is limited by MAX_SMART_WINDOW_ITERATIONS
                    if options.SWexactReadCount:
#                        windowBoundary2try = readfileCollection.findNextGap(windowBoundary2try,chrom,chromSizeHash[chrom])
                        windowEstimatedReadCount = readfileCollection.getExactReadCount4Region(chrom, curChromPos, windowBoundary2try)
                    else:
                        windowEstimatedReadCount = readfileCollection.getEstimatedReadCount4Region(chrom, curChromPos, windowBoundary2try)
                    while smartWindowSize > minWindowSize and windowEstimatedReadCount > options.smartWindows:
                        if debug:
                            print "Number of reads estimated in " +chrom+"["+str(curChromPos)+","+str(windowBoundary2try)+"): "+str(windowEstimatedReadCount)+" > "+str(options.smartWindows)+" --> decreasing window size "
                        windowBoundary2try = int((curChromPos + windowBoundary2try)/2)
                        if options.SWexactReadCount:
#                            windowBoundary2try = readfileCollection.findNextGap(windowBoundary2try,chrom,chromSizeHash[chrom])
                            windowEstimatedReadCount = readfileCollection.getExactReadCount4Region(chrom, curChromPos, windowBoundary2try)
                        else:
                            windowEstimatedReadCount = readfileCollection.getEstimatedReadCount4Region(chrom, curChromPos, windowBoundary2try)
                        smartWindowSize = windowBoundary2try - curChromPos
                        
                    if smartWindowSize <= minWindowSize:
                        print "INFO: minimum smart window length reached for window " +chrom+"["+str(curChromPos)+","+str(windowBoundary2try)+") ! --> Proceeding anyway"
                processUnitStart = curChromPos
#                processUnitEnd = readfileCollection.findNextGap(windowBoundary2try,chrom,chromSizeHash[chrom])
                processUnitEnd = windowBoundary2try
                curChromPos = processUnitEnd
                numRegions += 1
                
                # start process for current processUnit
                outfileNameHashTemp = {}
                # if options.debugRandomSubset is set to a value 0 < v < 1 with a probability of 1-v this variable is set to false with and the process will not be run 
                runprocess = (options.debugRandomSubset == 1 or random.random() < options.debugRandomSubset)
                overlapsIncludedRegions = includedRegions.overlapsWithRegion(manualGenomeCorrection.get(chrom,chrom),processUnitStart,processUnitEnd) 
                windowStatsStr = chrom +"\t"+str(processUnitStart)+"\t"+str(processUnitEnd)+"\t"+str(processUnit)+"\t"+str(windowEstimatedReadCount)+"\t"
                if runprocess and overlapsIncludedRegions:
                    windowStatsStr += "True"
                    for seqMotif in seqMotifs:
                        if not options.statisticsOnly:
                            outfileNameHashTemp = outfileH.getRegionFileNameDictAndUpdateseqMotifFiles2cat(seqMotif,chrom,processUnit)
            
                        # launch process
                        numProcesses += 1
                        #N#print "-------------------------------------------------------------------------------------------------"
                        #N#print "Initializing process at " + time.strftime("%Y-%m-%d %H:%M:%S") + ": "+chrom+" "+str(processUnit)+" window: [" + str(processUnitStart) + "," + str(processUnitEnd) + ")"
                        chrom_clean = cleanupChromName(chrom)
                        procName = "P_"+ seqMotif + "_" + chrom_clean+"_"+str(processUnit)+"_" + str(processUnitStart) + "to" + str(processUnitEnd)
                        regionID = chrom_clean + "-" + str(processUnit)
                        regionCoordTupel = (chrom,processUnitStart,processUnitEnd)
                        picklePathChrom = picklePath + chrom + os.sep
                        if not os.path.exists(picklePathChrom): os.mkdir(picklePathChrom)
                        outPKLfn = picklePathChrom + procName + "_out.cPickle"
                        #check whether region is already processed. if so, then nothing needs to be done
                        if not checkRestart4region(outfileNameHashTemp,outPKLfn):
                            shellScriptHs[seqMotif].addPKLfile(outPKLfn)
                            print "Skipping region ["+chrom +'_' + str(processUnit)+"_" + str(processUnitStart) + "to" + str(processUnitEnd) + "] <-- already processed"
                        else:
                            paramPKLfn = picklePathChrom + procName + "_InParams.cPickle"
                            paramPKLfile = open(paramPKLfn, 'wb')
                            cPickle.dump(outPKLfn, paramPKLfile)
                            cPickle.dump(seqMotif, paramPKLfile)
                            cPickle.dump(chrom,paramPKLfile)
                            cPickle.dump(processUnitStart, paramPKLfile)
                            cPickle.dump(processUnitEnd, paramPKLfile)
                            cPickle.dump(outfileNameHashTemp, paramPKLfile)
                            cPickle.dump(processZSField, paramPKLfile)
                            paramPKLfile.close() 
                            
                            pString = "python %s --task region --pickleFile %s --RP_paramPickle %s" % (sys.argv[0],inPKLfn,paramPKLfn) # sys.argv[0]: currently executed script
                            
                            if options.debugEnableRuntimeProfiling:
                                pString += " --debugEnableRuntimeProfiling --runtimeProfileFile " + runtimeProfileDir + os.sep + "runtime_profile_region_" + seqMotif + "_" + chrom + "_" + str(processUnitStart) + "to" + str(processUnitEnd) + ".cProfile"
                            jobStr = shellScriptHs[seqMotif].addProcess(pString, regionID, outPKLfn,regionCoordTupel)
                            if len(jobStr) > 0:
                                jobH.addJob(jobStr)
       
                else:
                    windowStatsStr += "False"
                    if overlapsIncludedRegions:
                        print "Randomly skipping current interval at " + time.strftime("%Y-%m-%d %H:%M:%S") + ": "+chrom+" "+str(processUnit)+" window: [" + str(processUnitStart) + "," + str(processUnitEnd) + ")"
                    else:
                        print "[" + time.strftime("%Y-%m-%d %H:%M:%S") + "] skipped intervall not included in --includedRegions" + ": "+chrom+" "+str(processUnit)+" window: [" + str(processUnitStart) + "," + str(processUnitEnd) + ")"
                
                if debug:
                    windowStatsF.write(windowStatsStr+"\n")
                processUnit += 1     
        #close the pysam read files
        del readfileCollection
        
        if debug:
            windowStatsF.close()
        
        for seqMotif in seqMotifs:    
            #execute the remaining processes if there are any
            if shellScriptHs[seqMotif].curNumP > 0:
                jobStr = shellScriptHs[seqMotif].writeShellScript()
                jobH.addJob(jobStr)
        
        print "Region processing completed at " + time.strftime("%Y-%m-%d %H:%M:%S") + " using " + str(numProcesses) + " of " + str(numRegions) + " processes/chromUnits total" 
            
            
def isBedPlus8applicable(filename):
    '''
    checks if a tab seperated file has > MAX_LIFTOVER_COLS (ie 20) columns and print a warning message
    this is useful as LiftOver with the option -bedPlus returns a maximum of 20 columns and neglects the others
    '''
    MAX_LIFTOVER_COLS = 20
    f = open(filename,"r")
    l = f.readline()
    numcols = len(l.split("\t"))
    f.close()
    if  numcols > MAX_LIFTOVER_COLS:
        print "WARNING: LiftOver using -bedPlus=8: too many columns (" + str(numcols) + ") could not transfer all of them"
        
        
def performLiftOver(infile, inGenome, outGenome, execute, bedPlus8 = False):
    '''
    performs liftOver
    '''
    success = False
    outfile = os.path.splitext(infile)[0]+"."+outGenome
    if not os.path.exists(infile):
        print "liftOver input file does not exist: "+infile
        raise Exception("LiftOver file does not exist")
    chainFilename = options.toolsPath + inGenome + "To" + outGenome.capitalize()+ ".over.chain.gz"
    callString = options.toolsPath+"liftOver.linux.x86_64 "
    if bedPlus8:
        callString = callString + "-bedPlus=8 "
    callString = callString + infile+" "+chainFilename+" "+outfile+" "+outfile+".unmappable"
    if execute:
        print "["+time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())+"] " + "Calling liftOver: " + callString
        sys.stdout.flush()
        if os.system(callString) == 0:
            print "["+time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())+"] " + "liftOver mapping succeeded"
            success = True
        else: print "liftOver mapping failed: " + callString
    else:
        print "Preparing liftOver call: " + callString
    return (outfile, callString, success)

def countHashLines(fileName):
    '''
    counts the number of lines in a file starting with '#'
    '''
    #open a subprocess to count the lines beginning with '#'
    myprocess = subprocess.Popen(['grep','-c','^\#',fileName],stdout=subprocess.PIPE)
    #and get the output
    (sout,serr) = myprocess.communicate()
    #lines are the first number in the string
    lineNum = int(sout.strip())
    return lineNum

def combineDuplicatesInCpGmethFile(infilename,outfilename,seqMotif): #TODO: function is relatively ugly: be sure to be well documented
    '''
    go through a SORTED bed file for methylation (CpG info) and combine subsequent lines that share the same coordinates
    + and - strand will be combined
    called after LiftOver which might introduce CpG duplicates 
    '''
    #does the reverse complement pattern match the forward pattern?
    mergeOppositeStrands = SEQ_MOTIF_IS_SYMMETRICAL[seqMotif] 
    
    infile = open(infilename,"rb")
    outfile = open(outfilename,"wb")
    
    duplicates = 0
    
    meth = 0 #values for methylated and total sites for current C
    total = 0
    strand = "INITIAL_STRING" # iitial string because "" is the caracter that is to be written if CpG from both strands are combined
    llprev = [] #previous splitted line
    for curline in infile:
        llcur = curline.strip().split("\t")
        
        if len(llprev) >= 6 and len(llcur) >= 6:
            curMethInfo = llcur[3].strip("'").split("/")
            #chrom, startPos, endPos and strand are identical
            if llcur[0] == llprev[0] and llcur[1] == llprev[1] and llcur[2] == llprev[2] and (mergeOppositeStrands or llcur[5] == llprev[5]):
                duplicates += 1
                meth += int(curMethInfo[0])
                total += int(curMethInfo[1]) 
                if llcur[5] == llprev[5] and strand <> "b":
                    strand = llcur[5]
                else:
                        strand = "b"                 
            else:
                meanMeth = int(round((float(meth) / total)*1000))
                outfile.write(getSepStringFromList(llprev[:3],"\t") + "\t'" + str(meth) + "/" + str(total) + "'\t" + str(meanMeth) + "\t" + strand + "\n")
                meth = int(curMethInfo[0])
                total = int(curMethInfo[1])
                strand = llcur[5]
        #if a line is not empty but doesnt fit the pattern just copy that line and reset
        elif len(llcur) >= 6: #only the previous row broken (or first row)
            if len(llprev) > 0:
                outfile.write(getSepStringFromList(llprev,"\t") + "\n")
            curMethInfo = llcur[3].strip("'").split("/")
            meth = int(curMethInfo[0])
            total = int(curMethInfo[1])
            strand = llcur[5]
        #reset if it does not fit
        elif len(llprev) >= 6: #only current row broken
            meanMeth = int(round((float(meth) / total)*1000))
            outfile.write(getSepStringFromList(llprev[:3],"\t") + "\t'" + str(meth) + "/" + str(total) + "'\t" + str(meanMeth) + "\t" + strand + "\n")
            meth = 0
            total = 0
            strand = "INITIAL_STRING"
        else: #both rows broken
            if len(llprev) > 0:
                outfile.write(getSepStringFromList(llprev,"\t") + "\n")
            meth = 0
            total = 0
            strand = "INITIAL_STRING"

        llprev = llcur
        
    #do it once more for the last line
    if len(llprev) >= 6:
        meanMeth = int(round((float(meth) / total)*1000))
        outfile.write(getSepStringFromList(llprev[:3],"\t") + "\t'" + str(meth) + "/" + str(total) + "'\t" + str(meanMeth) + "\t" + strand + "\n")
    else:
        outfile.write(getSepStringFromList(llprev,"\t") + "\n")
                
    infile.close()
    outfile.close()
    return duplicates

def combineDuplicatesInFragmentFile(infilename,outfilename):
    '''
    combines duplicates in the fragment file (fragments with same coordinates)
    '''
    duplicates = 0
    infile = open(infilename,"rb")
    outfile = open(outfilename,"wb")
    llprev = [] #previous splitted line
    prevline = ""
    curFragmentObj = None
    prevFragmentObj = None
    for curline in infile:
        llcur = curline.strip().split("\t")
        str2write = ""
        if len(llcur) >= 6: # current line valid?
            curMethInfo = llcur[3].strip("'").split(":")[1].split(",")
            curMethList = map(int,curMethInfo)
            curMethList10 = map(lambda x: x *10,curMethList) #since the constructor of UniqueFragmentsPerLane expects the bed specific ratios: 100% corresponds to 1000
            #chrom, startPos, endPos
            curFragmentObj = UniqueFragmentsPerLane(llcur[0]+":"+llcur[1]+":"+llcur[2], curMethList10, llcur[5])
            if len(llprev) >= 6: #previous line valid ?
                if llcur[0] == llprev[0] and llcur[1] == llprev[1] and llcur[2] == llprev[2]:
                    duplicates += 1
                    prevFragmentObj += curFragmentObj               
                else:
                    #write the previous accumulated fragment information
                    if prevFragmentObj: #check if its the first line
                        outfile.write(str(prevFragmentObj))
                    prevFragmentObj = curFragmentObj
            elif len(llprev) > 0:
                #write the previous line
                outfile.write(prevline + "\n")
                prevFragmentObj = curFragmentObj
            elif not prevFragmentObj: #check if its the first line
                prevFragmentObj = curFragmentObj   
        elif len(llcur) > 0: #current line invalid and not empty
            outfile.write(str(prevFragmentObj))
        #if the current line is empty, do nothing
        
        prevline = curline
        llprev = llcur
    
    #also write the last line
    if len(llprev) >= 6:
        outfile.write(str(prevFragmentObj))
    elif len(llprev) > 0:
        outfile.write(prevline + "\n")
    infile.close()
    outfile.close()
    return duplicates

def writeHeaderFiles(options):
    '''
    creates files for headers for CpG methylation, Read and fragment information (not for readDetails)
    '''
    for fileName in ["header_for_"+options.methodPrefix+"_cpgMethylation_file.txt","header_for_"+options.methodPrefix+"_cpgReads_file.txt","header_for_"+options.methodPrefix+"_cpgFragmentCoverage_file.txt"]:
        headerFilename = options.outputDir + os.sep + fileName
        if not os.path.exists(headerFilename):
            try:
                headerFile = open(headerFilename,"w")
                if fileName == "header_for_"+options.methodPrefix+"_cpgMethylation_file.txt":
                    headerFile.write("chrom\tchromstart\tchromend\tmethRatio\tmeanMeth\tstrand\n")
                if fileName == "header_for_"+options.methodPrefix+"_cpgReads_file.txt":
                    headerFile.write("chrom\tchromstart\tchromend\tmethRatios\tmeanMeth\tstrand\treadstart\treadend\titemRbg\tseqMotif_Count\tseqMotif_Lengths\tseqMotif_Positions\tPEnum\n")
                if fileName == "header_for_"+options.methodPrefix+"_cpgFragmentCoverage_file.txt":
                    headerFile.write("chrom\tchromStartOfFragment\tchromEndOfFragment\tmethRatiosOfAllReads\treadCount\tstrand\n")
                headerFile.close()
            except:
                pass # an Exception can occur when several jobs access the file simultaneously

def getNoticesInLog(logFN):
    '''
    greps for keywords of errors in a given logfile
    '''
    #all messages to search for
    searchPat = re.compile("(error|except|fail|exit|warn|note|notice)",flags=re.IGNORECASE)
    noticeStrings = []
    f = open(logFN,"r")
    for ll in f:
        curSearch = searchPat.search(ll)
        if curSearch:
            noticeStrings.append(curSearch.string)
    f.close()
    return noticeStrings
       
# MAIN ANALYSIS PROCEDURE
def performAnalysis(options):
    '''
    'main' task, initiator function
    prepared options format and collector datastructures
    prepares reference and read input files as well as the regions that will be processed
    pickles the required information for the 'region' and 'warpup' tasks
    invokes the chopping of the genome (manageProcesses)
    invokes the 'region' and 'wrapup' tasks
    '''
    print "MAIN process started at "
    print "Local time: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    
    if debug:
        print "Alignment files: " + str(options.alignmentFiles)
    
    # make sure that all alignment files are ready for processing
    for alignmentFile in options.alignmentFiles:
        minTimeSinceAlignment = time.time() - (options.timeDelay*3600)  # alignment files have to be (by default) at least 6h old to be analyzed (to make sure that alignment processing has finished)
        if (os.path.getmtime(alignmentFile) > minTimeSinceAlignment):
            sys.stdout.flush()
            print "Alignment generation is still running! \nRestart analysis when alignment processing has finished."
            raise SystemExit    
    
    # prepare temporary output directory and save program parameters for debugging
    if not os.path.exists(options.tempDir): os.mkdir(options.tempDir)
    if not os.path.exists(options.outputDir): os.mkdir(options.outputDir)
 
    path = options.tempDir + os.sep
    optionsFile = path + "options_object.bin"
    output = open(optionsFile, 'wb'); cPickle.dump(options, output); output.close()
    print "Load option settings as follows (for debugging): input = open('" + optionsFile + "', 'rb'); import pickle; options = pickle.load(input);"

    # obtain genome sequence if necessary
    global genomePath
    genomePath = options.genomeDir + os.sep + options.inGenome
    if not os.path.exists(genomePath):
        print "Genome assembly not found for '"+options.inGenome+"', trying to obtain genome from UCSC Genome Browser"        
        downloadGenome.downloadGenome(options.inGenome, options.genomeDir)    
    if not os.path.exists(genomePath):
        print "Could not obtain required genome assembly for '"+options.inGenome+"'"
        raise SystemExit
    
    # identify which chromosomes to include or exclude
    global chromSizeHash
    includeList = []
    if options.includedChromosomes <> "":
        includeList = options.includedChromosomes.split(",")
        
    excludeList = options.excludedChromosomes.split(",")
    chromSizeHash = {}
    chromFileHandler = chromFileCollection()
    for alignmentFile in options.alignmentFiles:
        if options.inputFormat == "BAM":
            # index the BAM alignment file to enable fast random access
            if not os.path.exists(alignmentFile+'.bai'): 
                try:
                    pysam.index(alignmentFile)
                except pysam.SamtoolsError, e:
                    print "WARNING: the following exception occurred during BAM file indexing of " + alignmentFile + " :"
                    print str(e)
                    raise SystemExit
            samfile = pysam.Samfile(alignmentFile, "rb")
        else:
            raise Exception("Unrecognized file format") 
        # retrieve all chromosomes (from all lanes) to be analyzed and their respective sizes using the pysam module
        assert(len(samfile.references) == len(samfile.lengths))
        
        for i in range(len(samfile.references)):
            #genomic region retrieved from samfile can have either chromosone numbers as strings (e.g. '21', 'X') or full identifiers (e.g. 'chr21','chrX')
            #--> use manualGenomeCorrection.get(genomicRegion,genomicRegion) to correct for that
            #N#print samfile.references[i]
            genomicRegion = samfile.references[i]
            chrSize = samfile.lengths[i]
            chrom = manualGenomeCorrection.get(genomicRegion,genomicRegion)
            if chrom in excludeList: continue
            if len(includeList) > 0 and not chrom in includeList: continue
            if options.ignore_chrRandom_chrM_hap and chrom == "chrM" or chrom.find("random") >= 0 or chrom.find("hap") >= 0: continue
            chromFilename = genomePath + os.sep + chrom + ".fa"
            if not os.path.exists(chromFilename):
                print "WARNING: No sequence data found, skipping chromosome: " + chrom
            elif genomicRegion not in chromSizeHash.keys():                
                chromSizeHash[genomicRegion] = chrSize
                chromFileHandler.addFile(chrom, chromFilename)
        samfile.close()
        
        
    chromSortList = sorted(chromSizeHash.keys())
    if debug: print 'chromSizeHash  ' + str(chromSizeHash)
    if debug: print 'chromSizeHash  ' + str(len(chromSizeHash))
    
    if debug: print 'chromSortList  ' + str(chromSortList)
    if debug: print 'chromSortList  ' + str(len(chromSortList))
    
    includedRegions = IncludedRegions(options.includedRegions)
   
    # initialize analysis parameters
    global seqMotifList
    seqMotifList = []
    if not options.ignoreCpG: seqMotifList.append('G')
    if options.processCpA: seqMotifList.append('A')
    if options.processCpC: seqMotifList.append('C')
    if options.processCpT: seqMotifList.append('T')
    if options.processCpHpG: seqMotifList.append('HpG')
    if options.processCpHpH: seqMotifList.append('HpH')
    if options.processNonCpG: seqMotifList.append('H')
    numberOfLanes = len(options.alignmentFiles)
    strands = ['both']
    if options.isStrandSpecific: strands = ['both', '-', '+']
    
    #N#print 'numberOfLanes---------------------------------'
    #N#print numberOfLanes
    
    fragmentsInfo = FragmentsInfo(options.fragmentSizeRanges)
    statisticsLanes = StatisticsCollection(options.alignmentFiles,seqMotifList,strands,fragmentsInfo)
    
    restrSite = RestrictionSite(options.restrictionSite)
  
    #output file handler class
    outfns = OutFileHandler(path)
    
    #strand classification
    global strandModel
    strandModel = None
    if options.strandPredict:
        strandModel = strandClassifierLogReg()
    
    #main pickleFile for general variables
    mainPickleFileName4regions = path +'main_data_4regions.pkl' 
    
    #write general information to the main picklefile for the wrapup task 
    mpf1 = open(mainPickleFileName4regions, 'wb')
    cPickle.dump(strands, mpf1)
    cPickle.dump(path, mpf1)
    cPickle.dump(fragmentsInfo, mpf1)
    cPickle.dump(restrSite, mpf1)
    cPickle.dump(chromSizeHash, mpf1)
    cPickle.dump(chromFileHandler, mpf1)
    cPickle.dump(strandModel, mpf1)
    cPickle.dump(options, mpf1)
    cPickle.dump(debug, mpf1)
    mpf1.close()
    
    shellScriptHandlers4seqMotifs = {}
    jobH = jobHandler()
    
    for seqMotif in seqMotifList:
        print "Processing Cp" + seqMotif
        #prepare output filenames
        outfns.addFiles4seqMotif(seqMotif, options)
        #shellScriptHandler
        shellScriptH = ProcessShellScriptHandler(path,seqMotif,options.maxProcesses,submitCmd="bsub " + options.lsfOpts)
        shellScriptHandlers4seqMotifs[seqMotif] = shellScriptH
        
    manageProcesses(chromSortList,chromSizeHash,chromFileHandler,includedRegions,outfns,shellScriptHandlers4seqMotifs,jobH,mainPickleFileName4regions,strands,seqMotifList,path,options,fragmentsInfo,debug)
    
    spikeIns = None
    if options.spikeInControlSeqs <> "":
        print "Reading Spike In controls from " + options.spikeInControlSeqs
        spikeIns = SpikeInControls(options.spikeInControlSeqs,options.spikeInMotifLen)
        
    mainPickleFileName4wrapup= path +'main_data_4wrapup.pkl'    
    #write general information to the main picklefile for the wrapup task 
    mpf2 = open(mainPickleFileName4wrapup, 'wb')
    cPickle.dump(statisticsLanes, mpf2)
    cPickle.dump(outfns, mpf2)
    cPickle.dump(shellScriptHandlers4seqMotifs, mpf2)
    cPickle.dump(path, mpf2)
    cPickle.dump(spikeIns, mpf2)
    cPickle.dump(seqMotifList, mpf2)
    cPickle.dump(options, mpf2)
    cPickle.dump(debug, mpf2)
    cPickle.dump(exceptionOccurred, mpf2)
    mpf2.close()
    
    #start the wrapup task
    pString = "python %s --task wrapup --pickleFile %s" % (sys.argv[0],mainPickleFileName4wrapup)
    if options.debugEnableRuntimeProfiling:
        runtimeProfileDir = options.outputDir + os.sep + "runtime_profiles"
        pString += " --debugEnableRuntimeProfiling --runtimeProfileFile " + runtimeProfileDir + os.sep + options.sampleName + "runtime_profile_wrapup.cProfile"
    wrapupLogFn = options.outputDir + os.sep + options.sampleName + "_wrapup.log"
    if options.parallelize:
        jobString = "bsub " + options.lsfOptsWrapup + " -J " + options.sampleName + "_wrapup" + " -o " + wrapupLogFn
        #add dependencies for every single job
        addStr = " -w \""
        #add dependencies for all the job prefixes
        for seqMotif in seqMotifList:
            psch = shellScriptHandlers4seqMotifs[seqMotif]
            addStr += "done('%s_*') && " % (options.sampleName + "_" + psch.curFileNameBase)
#            for jn in psch.jobnames:
#                addStr += "done('%s') && " % (jn)
        addStr = addStr[:-4] + "\"" #:-4 cut off trailing ' && '
        jobString += addStr + " " + pString
        if len(addStr) < 6:
            print "Sorry, NO region jobs were submitted. Exiting..."
            print "This can have multiple reasons:"
            print "(a) Temporary files are already present for all the regions. Have you run the analysis already?"
            print "(b) You specified an illegal subset of regions"
            print "(c) --debugRandomSubset option resulted in no regions to be processed by chance"
            #print "(d) witchcraft. In this case burn all women with red hair around you! ;-)"
            print "(d) programming error. In this case we would be happy about your report."
            sys.stdout.flush()
            raise SystemExit
    else:
        jobString = pString
    print "-------------------------------------------------------------------"
    print "WRAPUP job:"
    print "Local time: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print jobString
    sys.stdout.flush()
    jobH.addJob(jobString)
    
    print "-------------------------------------------------------------------"
    print "executing jobs and writing to shellscript ..."
    print "Local time: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    sys.stdout.flush()
    if not options.prepareOnly:
        jobH.execute(2)
    jobH.write2file(options.outputDir + os.sep + options.sampleName + "_submittedJobs.sh")
    
    print "-------------------------------------------------------------------"
    print "MAIN process finished at "
    print "Local time: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())        


def performAnalysis_region(options_l):
    '''
    'region' task, initiator function
    reads the required options
    reads the required information from the main tasks pickle file
    invokes sequenceDataProcessing()
    '''      
    #read the input parameter pickle file
    ppf = open(options_l.RP_paramPickle, 'rb')
    outPKLfn = cPickle.load(ppf)
    seqMotif = cPickle.load(ppf)
    chrom = cPickle.load(ppf)
    processUnitStart = cPickle.load(ppf)
    processUnitEnd = cPickle.load(ppf)
    outfileNameHashTemp = cPickle.load(ppf)
    processZSField = cPickle.load(ppf)
    ppf.close()

    #read main picklefile
    mpf = open(options_l.pickleFile, 'rb')
    strands = cPickle.load(mpf)
    path = cPickle.load(mpf)
    fragmentsInfo = cPickle.load(mpf)
    restrSite = cPickle.load(mpf)
    global chromSizeHash
    chromSizeHash = cPickle.load(mpf)
    chromFileHandler =  cPickle.load(mpf)    
    global strandModel
    strandModel = cPickle.load(mpf)
    # WARNING: THE INPUT OPTIONS ARE OVERWRITTEN BY THE OPTIONS FROM THE MAIN TASK HERE !!!
    global options
    options = cPickle.load(mpf)
    global debug
    debug = cPickle.load(mpf)
    mpf.close()
    setGlobalsMaxMismatches(options)
    
    options.processZSField = processZSField
    sequenceDataProcessing(chrom, processUnitStart, processUnitEnd, outPKLfn, outfileNameHashTemp, strands, seqMotif, fragmentsInfo,chromSizeHash,chromFileHandler,restrSite)
    
    #N#print "-------------------------------------------------------------------"
    #N#print "REGION process finished at "
    #N#print "Local time: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    
def performAnalysis_wrapup(options_l):
    '''
    'wrapup' task, initiator function
    reads the required options
    reads the required information from the main tasks pickle file
    checks whether the region tasks created the output files they were supposed to
    combines the region statistics objects and calculate global statistics object
    performs statistics calculations (eg conversion rate,...)
    concatenate the regions list output files (reads, readDetails, CpGs, fragments)
    perfoms Liftover, conversion to bigBed and zipping (if required)
    '''
    print "WRAPUP process started at "
    print "Local time: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())       
    #read main picklefile
    mpf = open(options_l.pickleFile, 'rb')
    statisticsLanes = cPickle.load(mpf)
    outfns = cPickle.load(mpf)
    shellScriptHandlers4seqMotifs = cPickle.load(mpf)
    path = cPickle.load(mpf)
    spikeIns = cPickle.load(mpf)
    global seqMotifList
    seqMotifList = cPickle.load(mpf)
    # THE INPUT OPTIONS ARE OVERWRITTEN BY THE OPTIONS FROM THE MAIN TASK HERE !!!
    global options
    options = cPickle.load(mpf)
    global debug
    debug = cPickle.load(mpf)
    global exceptionOccurred
    exceptionOccurred = cPickle.load(mpf)
    mpf.close()
    #read statistics picklefiles
    
    #check whether all the pickle files are in place
    missingPKLfiles = []
    for seqMotif in seqMotifList:
        for statsPickle in shellScriptHandlers4seqMotifs[seqMotif].pickleFileNames:
            if not isValidFile(statsPickle):
                missingPKLfiles.append(statsPickle)
    if len(missingPKLfiles) > 0:
        print "ERROR: Region picklefile(s) missing. Program will terminate now"
        print str(missingPKLfiles)
        raise SystemExit    
    
    for seqMotif in seqMotifList:
        #update the statistics for each region process
        for statsPickle in shellScriptHandlers4seqMotifs[seqMotif].pickleFileNames:
            pkl_file = open(statsPickle, 'rb')
            processStatCol = cPickle.load(pkl_file)
            exceptionOccurred |= cPickle.load(pkl_file)
            pkl_file.close()
            statisticsLanes.updateProcessStats(seqMotif,processStatCol)
        
        # outfile handler class
        outfns.catRegionFiles(seqMotif)
        
        # LiftOver
        doLiftover = options.inGenome != options.outGenome
        for finalOutfile in outfns.seqMotifOutfileNames[seqMotif].keys():
            print "["+time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())+"] " + "Postprocessing file " + outfns.seqMotifOutfileNames[seqMotif][finalOutfile]
            sys.stdout.flush()
            if doLiftover:
                if finalOutfile <> 'readQualOutfile': 
                    (finalOutfileAfterLiftOver, callString, success) = performLiftOver(outfns.seqMotifOutfileNames[seqMotif][finalOutfile], options.inGenome, options.outGenome, True)
                else:
                    if options.readDetailsFileLiftOver:
                        isBedPlus8applicable(outfns.seqMotifOutfileNames[seqMotif][finalOutfile])
                        (finalOutfileAfterLiftOver, callString, success) = performLiftOver(outfns.seqMotifOutfileNames[seqMotif][finalOutfile], options.inGenome, options.outGenome, True, bedPlus8=True)
                    else:
                        continue
            
                os.system("mv " + outfns.seqMotifOutfileNames[seqMotif][finalOutfile] + " " + outfns.seqMotifOutfileNames[seqMotif][finalOutfile] + ".beforeLiftOver")
                os.system("mv " + finalOutfileAfterLiftOver + " " + outfns.seqMotifOutfileNames[seqMotif][finalOutfile])
            
            doSort = doLiftover or finalOutfile.find('fragmentOutfile') >= 0 or finalOutfile == 'readQualOutfile' or finalOutfile == 'readOutfile'
            
            if doSort:
                print "  ["+time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())+"] " + "Sorting bed file..."
                sys.stdout.flush()
                unsortedFile = outfns.seqMotifOutfileNames[seqMotif][finalOutfile] + ".unsorted"
                os.system("mv " + outfns.seqMotifOutfileNames[seqMotif][finalOutfile] + " " + unsortedFile)
                
                sortCmd = "sort -k1,1 -k2,2n " + unsortedFile + " > " + outfns.seqMotifOutfileNames[seqMotif][finalOutfile]
                os.system(sortCmd)
            if doLiftover:
                if debug:
                    print "  ["+time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())+"] " +"Counting mapping e.r.r.o.r.s.."
                    sys.stdout.flush()
                # count the number of mapping errors == number of lines starting with '#' in the unmappable file
                try:
                    unmappableFile = finalOutfileAfterLiftOver + ".unmappable"
                    unmappableAfterLiftover = countHashLines(unmappableFile)
                    statisticsLanes.updateLiftOverMappingErrors(unmappableAfterLiftover,seqMotif,finalOutfile,options)
                except:
                    print "WARNING: could not count mapping errors from file " + unmappableFile
                    exceptionOccurred = True
            
                #LiftOver might create duplicate CpGs in the CpG Output file: remove them ...   
                if finalOutfile.find('cpgOutfile') >= 0:
                    fileWithDuplicates = outfns.seqMotifOutfileNames[seqMotif][finalOutfile] + ".duplicates"
                    os.system("mv " + outfns.seqMotifOutfileNames[seqMotif][finalOutfile] + " " + fileWithDuplicates)
                    
                    print "  Search and combine duplicates: seqMotif output..."
                    numDuplicates = combineDuplicatesInCpGmethFile(fileWithDuplicates,outfns.seqMotifOutfileNames[seqMotif][finalOutfile],seqMotif)
                    if debug:
                        print "    " + str(numDuplicates) + " duplicates found"
            #...and the fragment output file. Here the duplicates could also be due to the processed Regions cutting through a fragment
            if finalOutfile.find('fragmentOutfile') >= 0:
                fileWithDuplicates = outfns.seqMotifOutfileNames[seqMotif][finalOutfile] + ".duplicates"
                os.system("mv " + outfns.seqMotifOutfileNames[seqMotif][finalOutfile] + " " + fileWithDuplicates)
                
                print "  Search and combine duplicates: fragments..."
                numDuplicates = combineDuplicatesInFragmentFile(fileWithDuplicates,outfns.seqMotifOutfileNames[seqMotif][finalOutfile])
                if debug:
                    print "    " + str(numDuplicates) + " duplicates found"
                    
        if debug:
            print 60*"-"+"\nPostprocessing done"
            sys.stdout.flush()
        
        outfns.addHeaderReadQualOutfile(seqMotif,options)

        # STATISTICS                           
        # most statistic values for the plus strand are obtained by simply subtracting the values for minus strand from the values for both strands   
        if options.isStrandSpecific: 
            statisticsLanes.updateStrandSpecific(seqMotif)
        
        if debug:
            print "[" + str(time.strftime("%Y-%m-%d %H:%M:%S")) + "]: combining and calculating statistics..."
            sys.stdout.flush()
        statisticsLanes.calcStats_setConversionRate(seqMotif)
        
        #log files for submitted region jobs
        if options.parallelize:
            callStringLogs = "cat"
#            callStringErrs = "cat"
            for logFN in shellScriptHandlers4seqMotifs[seqMotif].logFileNames:
                callStringLogs += " " + logFN
#                callStringErrs += " " + logFN + ".err"
            concatRegionLogFN = options.outputDir + os.sep + options.methodPrefix + "_Cp" + seqMotif + "_regionJobLogs_" + options.sampleName +".log"
            callStringLogs += " > " + concatRegionLogFN
#            callStringErrs += " > " + concatRegionLogFN + ".err"
            print "    Concatenating region job logfiles..."
            os.system(callStringLogs)
#            os.system(callStringErrs)

    if options.spikeInControlSeqs <> "":
        print "Analysing Spike In controls from " + options.spikeInControlSeqs
        spikeIns.extractSICreads(options)
        spikeIns.processSICreads(options)
        statisticsLanes.readSpikeInCtrlAnalysis(spikeIns)
        spikeInStatsFileName = statisticsLanes.writeSpikeInStatistics(path,options)
        os.system("mv "+ spikeInStatsFileName +" "+ options.outputDir)
    
    writeHeaderFiles(options)    
    toolsPath = options.toolsPath
    chromSizesFile = None
    if options.bigBedFormat:
        print "[" + str(time.strftime("%Y-%m-%d %H:%M:%S")) + "]: BigBed handling..."
        chromSizesFile = options.outGenome.replace("/","")+'.chrom.sizes'
        if not os.path.exists(chromSizesFile): # check again with above code for redundancy
            chromSizesFileCmd = toolsPath + 'fetchChromSizes ' +options.outGenome+' > ' + chromSizesFile
            sys.stdout.flush()
            print "Calling fetchChromSizes: " + chromSizesFileCmd
            sys.stdout.flush()
            os.system(chromSizesFileCmd)

    #create webOutputDir if it does not exist already
    if not os.path.exists(options.webOutputDir): os.mkdir(options.webOutputDir)     
    outfns.processBedFiles(options,chromSizesFile)
    print "[" + str(time.strftime("%Y-%m-%d %H:%M:%S")) + "]: (zipping and) moving files..."
    outfns.zipAndMoveOutfiles()

    # GLOBAL STATISTICS over all lanes of sample, if more than one lane was given as input file
    print "[" + str(time.strftime("%Y-%m-%d %H:%M:%S")) + "]: Computing statistics from lanes and writing output..."
    statisticsLanes.computeGlobalStats()
    statisticsLanes.calcMeanSDs()
    statsFileName = statisticsLanes.writeStatistics2files(path,options)
    os.system("mv "+ statsFileName +" "+ options.outputDir)
    if options.deleteTemp and not exceptionOccurred: 
        sys.stdout.flush()
        shutil.rmtree(path)

    print "-------------------------------------------------------------------"
    print "WRAPUP process finished at "
    print "Local time: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    

def performAnalysis_diagnostics(options):
    '''
    the tempDir parameter has to be specified in the options
    retrieves summary information on a biseeqMethCalling analysis
    by investigating job processes and checking for error messages
    the temp directory and the temporary files need to be still available for proper functionality
    '''
    print "DIAGNOSTICS process started at "
    print "Local time: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    
    pfn_wrapup = options.tempDir + os.sep + 'main_data_4wrapup.pkl'
    try:
        mpf = open(pfn_wrapup, 'rb')
    except:
        print "could not open " + pfn_wrapup
        print "Program will terminate now"
        raise SystemExit
        
    statisticsLanes = cPickle.load(mpf)
    analysisOutfns = cPickle.load(mpf)
    analysisShellScriptHandlers4seqMotifs = cPickle.load(mpf)
    analysisPath = cPickle.load(mpf)
    analysisSpikeIns = cPickle.load(mpf)
    analysisSeqMotifList = cPickle.load(mpf)
    # THE INPUT OPTIONS ARE OVERWRITTEN BY THE OPTIONS FROM THE MAIN TASK HERE !!!
    analysisOptions = cPickle.load(mpf)
    analysisDebug = cPickle.load(mpf)
    analysisExceptionOccurred = cPickle.load(mpf)
    mpf.close()
    
    statsFileName = analysisOptions.outputDir + os.sep + analysisOptions.methodPrefix +"_statistics_"+ analysisOptions.sampleName + ".txt"
    analysisFinished = isValidFile(statsFileName)
    
    print "-------------------------------------------------------------------"
    print "GENERAL INFORMATION"
    print "-------------------------------------------------------------------"
    print "processed motifs:          " + "Cp" + listAsString(analysisSeqMotifList)
    print "# alignment files:         " + str(len(analysisOptions.alignmentFiles))
    print "analysis finished:         " + str(analysisFinished)
    
    print "-------------------------------------------------------------------"
    print "PROGRESS INFORMATION"
    print "-------------------------------------------------------------------"
    for seqMotif in analysisSeqMotifList:
        print "Cp" + seqMotif
        print "-------------------------------------------------------------------"
        #check which jobs' pickle files are missing
        missingPKLfiles = []
        for statsPickle in analysisShellScriptHandlers4seqMotifs[seqMotif].pickleFileNames:
            if not isValidFile(statsPickle):
                missingPKLfiles.append(statsPickle)
        
        numRegionsTotal = len(analysisShellScriptHandlers4seqMotifs[seqMotif].pickleFileNames)
        numRegionsProcessed = numRegionsTotal - len(missingPKLfiles)
        
        #check which jobs have already created a log file yet --> finished
        finishedJobLogFiles = []
        for logFN in analysisShellScriptHandlers4seqMotifs[seqMotif].logFileNames:
            if isValidFile(logFN):
                finishedJobLogFiles.append(logFN)
                
        numJobsTotal = len(analysisShellScriptHandlers4seqMotifs[seqMotif].logFileNames) 
        #print summary
        print "  # regions to process (total):     " + str(numRegionsTotal)
        print "  # regions already processed:      " + str(numRegionsProcessed) + " (" + str(float(numRegionsProcessed)/numRegionsTotal * 100.0) + "%)"
        print "  # region jobs to process (total): " + str(numJobsTotal)
        if numJobsTotal > 0:
            print "  # region jobs already finished:   " + str(len(finishedJobLogFiles)) + " (" + str(float(len(finishedJobLogFiles))/numJobsTotal * 100.0) + "%)"
        else:
            print "  # region jobs already finished:   NA (region jobs to process = 0)"
        
    print "-------------------------------------------------------------------"
    print "ERROR INFORMATION"
    print "-------------------------------------------------------------------"
    
    print "-------------------------------------------------------------------"
    print "Main Task"
    print "-------------------------------------------------------------------"
    mainLogFn = options.mainTaskLogFile
    if isValidFile(mainLogFn):
        noticeStrings_main = getNoticesInLog(mainLogFn)
        for ns in noticeStrings_main:
            print ns
    else:
        print "not a valid log file (task not finished yet?)"
    print "-------------------------------------------------------------------"
    print "Region Tasks"
    print "-------------------------------------------------------------------"
    if not analysisFinished:
        #look in the temporary log files
        for logFN in finishedJobLogFiles:
            noticeStrings = getNoticesInLog(logFN)
            if len(noticeStrings) > 0:
                print logFN + ":"
                for ns in noticeStrings:
                    print "    " + ns
    print "-------------------------------------------------------------------"
    print "Wrapup Task"
    print "-------------------------------------------------------------------"
    wrapupLogFn = analysisOptions.outputDir + os.sep + analysisOptions.sampleName + "_wrapup.log"
    if isValidFile(wrapupLogFn):
        noticeStrings_wrapup = getNoticesInLog(wrapupLogFn)
        for ns in noticeStrings_wrapup:
            print ns
    else:
        print "not a valid log file (task not finished yet?)"
    
    print "-------------------------------------------------------------------"
    print "DIAGNOSTICS task finished at "
    print "Local time: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    
    
if __name__ == '__main__':
    #N#print "Starting program..."
    # constructing command line parser
    import optparse
    parser = optparse.OptionParser()
    parser.add_option('--methodPrefix',action='store',type='string',dest='methodPrefix',help='Specify a method-specific prefix for the output files (typically set to RRBS or BiSeq, default=RRBS)',default='RRBS')
    parser.add_option('--alignmentFile',action='append',type='string',dest='alignmentFiles',help='Specify the names of input files containing aligned read positions (typically BAM files)',default=[])
    parser.add_option('--timeDelay',action='store',type='int',dest='timeDelay',help='Specify the minimum age of the alignment files (in hours, default=2). This parameter should be set if the alignment file is generated in-place, in order to minimize the risk of using incomplete alignment data',default=2)
    parser.add_option('--inputFormat',action='store',type='string',dest='inputFormat',help='Specify the format of the input alignment files (default=BAM, with is currently the only supported file format)', default="BAM")
    parser.add_option('--pfStatus',action='store',type='string',dest='pfStatus',help='Specify whether to use all reads, just the PF reads or just the non-PF reads (values={PF,NonPF,All}, default=All)', default="All")
    parser.add_option('--alignmentTool',action='store',type='string',dest='alignmentTool',help='Specify the alignment tool used to generate the alignment files (default=rrbsmap, with is currently the only supported alignment output format)', default="rrbsmap")
    
    parser.add_option('--rrbsMode',action='store_true',dest='rrbsMode',help='Specify whether the reads have been obtained from an RRBS experiment and aligned using RRBS-specific tools (default=False, indicating whole-genome bisulifte alignment)',default=False)
    parser.add_option('--genomeDir',action='store',type='string',dest='genomeDir',help='Specify the name of the directory containing the genome sequence data in fasta format with exactly 50 columns per line (will be downloaded from the USCS Genome Browser if not available locally)',default="genomes")
    
    parser.add_option('--pairedEnd',action='store_true',dest='pairedEnd',help='Specify whether the input file contains paired end reads and these should be treated as such',default=False)
    parser.add_option('--fullPE',action='store_true',dest='fullPE',help='Enables full paired end support (takes significantly more time: for each pair of reads that are improperly paired, only the one with fewer mismatches will be processed. (By default only the first read will be processed)',default=False)
    parser.add_option('--discardImproperPairs',action='store_true',dest='discardImproperPairs',help='ONLY valid if pairedEnd option is enabled: discard all reads that are not classified as proper pairs in the input file, including singleton reads (default=False)',default=False)
    parser.add_option('--strandPredict',action='store_true',dest='strandPredict',help='Specify whether only reads mapped to a strand that is confirmed by a logistic regression model should be included in the analysis (default=False, but highly recommended for whole-genome bisulfite sequencing)',default=False)
    parser.add_option('--noReadNumberField',action='store_false',dest='processZSField',help='Specify if the custum ZS field (defining read orientations in paired end sequencing) in the alignemnt BAM files should not be read, but the flag information should be used instead',default=True)
    
    parser.add_option('--toolsDir',action='store',type='string',dest='toolsPath',help='Specify the directory where command line tools such as liftOver, fetchChromSizes and bedToBigBed can be found (default=[current working directory]/tools/)', default='tools/')
    parser.add_option('--outputDir',action='store',type='string',dest='outputDir',help='Specify the name of the output directory', default='')
    parser.add_option('--webOutputDir',action='store',type='string',dest='webOutputDir',help='Specify the name of the web-accessible output directory for UCSC Genome Browser tracks', default='')
    parser.add_option('--tempDir',action='store',type='string',dest='tempDir',help='Specify the name of the directory for temp file storage', default='BiSeq_temp')
    parser.add_option('--sampleName',action='store',type='string',dest='sampleName',help='Specify a descriptive name for the biological sample (will be used as filename for most output files)', default='')
    
    parser.add_option('--genomeFraction',action='store',type='float',dest='genomeFraction',help='Specify the size of the genomic fraction to be analyzed in parallel (in multiples of 1MB, e.g. 0.5 for a half MB; default=5MB)', default=5.0)
    parser.add_option('--smartWindows',action='store',type='int',dest='smartWindows',help='Adapt window sizes based on estimated readnumbers. when determining windows, decrese the window size until the estimated number of reads is below this number. Default: 0 (off)',default=0)
    parser.add_option('--minSWlength',action='store',type='float',dest='minSWlength',help='Minimum Window size criterion to stop deviding up windows due to --smartWindows (in megabases). The actual size of a region might be as small as half of this Default: 0.02',default=0.02)
    parser.add_option('--SWexactReadCount',action='store_true',dest='SWexactReadCount',help='use the exact readcount rather than an estimate based on just 1 lane. Only meaningful if the --smartWindows option is enabled', default=False)
    parser.add_option('--readOutputBuffer',action='store',type='int',dest='readOutputBuffer',help='Specify how many reads should be stored in the buffer before writing to the read output file. Default: 50000',default=50000)
    
    parser.add_option('--maxProcesses',action='store',type='int',dest='maxProcesses',help='Specify maximum number of processes that are being run in parallel', default=8)
    parser.add_option('--inGenome',action='store',type='string',dest='inGenome',help='Specify the UCSC Genome Browser identifier of source genome assembly', default='')
    parser.add_option('--outGenome',action='store',type='string',dest='outGenome',help='Specify the UCSC Genome Browser identifier of destination genome assembly', default='') 
    parser.add_option('--removeDuplicates',action='store_true',dest='removeDuplicates',help='Specify whether reads that are optical or PCR duplicates should be removed', default=False)
    parser.add_option('--ignoreCpG',action='store_true',dest='ignoreCpG',help='Specify whether CpG methylation should be ignored in analysis', default=False)       
    parser.add_option('--processCpA',action='store_true',dest='processCpA',help='Specify whether CpA methylation should also be analyzed', default=False)
    parser.add_option('--processCpC',action='store_true',dest='processCpC',help='Specify whether CpC methylation should also be analyzed', default=False)
    parser.add_option('--processCpT',action='store_true',dest='processCpT',help='Specify whether CpT methylation should also be analyzed', default=False)
    parser.add_option('--processCpHpG',action='store_true',dest='processCpHpG',help='Specify whether CpHpG (H=A,C,T) methylation should also be analyzed', default=False)
    parser.add_option('--processCpHpH',action='store_true',dest='processCpHpH',help='Specify whether CpHpH (H=A,C,T) methylation should also be analyzed (equivalent to all methC - CpG - CpHpG)', default=False)
    parser.add_option('--processNonCpG',action='store_true',dest='processNonCpG',help='Specify whether nonCpG methylation should also be analyzed (does not treat each nonCpG seqMotif differently)', default=False)
    parser.add_option('--maxMismatches',action='store',type='float',dest='maxMismatches',help='Specify a float value between zero and one (e.g. 0.1 for a maximum of 10% mismatches) or an integer value corresponding to the maximum number of mismatches. (default: 0.4)', default=0.4)
    parser.add_option('--minMismatches',action='store',type='float',dest='minMismatches',help='Specify a float value between zero and one (e.g. 0.1 for a minimum of 10% mismatches) or an integer value corresponding to the minimum number of mismatches', default=0)
    
    parser.add_option('--mismatchTruncate',action='store',type='int',dest='mismatchTruncate',help='If set to >0 (0 is the default), this integer value specifies the number of beginning bases to be considered for the --maxMismatches parameter. This length will be extended for each sequencing read until the the maxMismatches number is reached', default=0)
    parser.add_option('--minReadLength',action='store',type='int',dest='minReadLength',help='Specify an integer determining the minimum read length (default: 5)', default=5)
    parser.add_option('--baseQualityScoreC',action='store',type='int',dest='baseQualityScoreC',help='Threshold for the read quality score (default: 20). Cs and Ts in a CpG context are only methylation-called if this threshold is reached or exceeded', default=20)
    parser.add_option('--baseQualityScoreNextToC',action='store',type='int',dest='baseQualityScoreNextToC',help='Threshold for the read quality score (default: 10).  The 3\' neighbors of Cs (or Ts) are only methylation-called if this threshold is reached or exceeded, except at the start position (no threshold exists)', default=10)
    parser.add_option('--restrictionSite',action='store',type='string',dest='restrictionSite',help='Specify the restriction site of the enzyme for --rrbsMode. "-" for cut on forward and "_" for cut on reverse strand. use lower case c to specify required UNmethylated Cs. Default:MspI site: "c-CG_G"', default='c-CG_G')
    parser.add_option('--checkRestriction',action='store_true',dest='checkRestriction',help='Specify whether reads have to start with the specified restriction site in order to be analyzed, e.g. start with (C/T)GG for unpaired MspI runs or strat with (C/T)GG if read#1 of a paired run or with CGA if read#2', default=False)
    parser.add_option('--ignore_chrRandom_chrM_hap',action='store_false',dest='ignore_chrRandom_chrM_hap',help='Specify whether reads with substring random, hap or chrM should be ignored', default=True)
    parser.add_option('--isStrandSpecific',action='store_true',dest='isStrandSpecific',help='Should + and - strand be analyzed seperately (eg useful in search for asymmetric methylation patterns). WARNING: use with caution: not fully tested yet', default=False)
    parser.add_option('--includedChromosomes',action='store',type='string',dest='includedChromosomes',help='Specify in a comma-separated string which chromosomes should be included (default: chr1 to chr30, chrX, chrY, chrZ)', default="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chr23,chr24,chr25,chr26,chr27,chr28,chr29,chr30,chrX,chrY,chrZ")
    parser.add_option('--excludedChromosomes',action='store',type='string',dest='excludedChromosomes',help='Specify in a comma-separated string which chromosomes should be excluded (e.g., it is useful to exclude chrY when aligning to the female genome)', default="")
    
    parser.add_option('--spikeInControlSeqs',action='store',type='string',dest='spikeInControlSeqs',help="Specify a fasta file containing sequences of spike in controls. In the identifier rows the last word should be either 'methylated' or 'unmethylated' depending on the type of control. If this variable is not specified (default) spike in controls will not be handled", default="")
    parser.add_option('--spikeInMinMatchLen',action='store',type='int',dest='spikeInMinMatchLen',help='[Only valid if --spikeInControlSeqs is specified] Length of pattern to be extracted to match reads to the spike in controls (default: 15)', default=15)
    parser.add_option('--spikeInMaxMismatches',action='store',type='int',dest='spikeInMaxMismatches',help='[Only valid if --spikeInControlSeqs is specified] maximum allowed mismatches to spike in control (default: 1)', default=1)
    #TODO: note that for symmetrical patternds like CGs: patterns on forward and reverse strand of the control are counted as individual CGs. is that ok?
    parser.add_option('--spikeInMotifLen',action='store',type='int',dest='spikeInMotifLen',help='[Only valid if --spikeInControlSeqs is specified] maximum length for analysis of bisulfite conversion. choose number carefully and not too large! (default: 2)', default=2)
    
    parser.add_option('--includedRegions',action='store',type='string',dest='includedRegions',help='Specify in a # seperated string which regions should be processed. Format example: chr3:234452-12412346#chrX:235694-1234789. Default:"" i.e. all', default="")
    parser.add_option('--bigBedFormat',action='store_true',dest='bigBedFormat',help='Specify whether to store the output browser tracks in bigBed format; files will have the same name but with .bb extension instead of .bed (for faster UCSC Genome Browser upload); if option --keepBedFiles is set to none, all .bed files will be deleted after generation of .bigBed files', default=False) 
    parser.add_option('--keepBedFiles',action='store',type='string',dest='keepBedFiles',help='Specify which bed files to keep. options are "none","CpGonly" (default) and "all"', default="CpGonly")
    parser.add_option('--readDetailsFile',action='store_true',dest='readDetailsFile',help='Specify whether to produce a file with detailed information for each read, like read-specific mismatches, read-specific conversion rates and quality data as additional columns to the first output file (file extension is .txt). By default this file contains coordinates according to --inGenome (see --readDetailsFileLiftOver for more options)', default=False)
    parser.add_option('--readDetailsFileLiftOver',action='store_true',dest='readDetailsFileLiftOver',help='only valid if --readDetailsFile is enabled: if --inGenome and --OutGenome are not identical, the readDetails file will be liftovered as well', default=False)
    parser.add_option('--fragmentSizeRanges',action='store',type='string',dest='fragmentSizeRanges',help='Specify fragment size ranges for fragment coverage statistics (enter size ranges in this format: 40-60,40-80,80-120)', default='0-50,50-100,100-150,150-200,200-250,250-500')
    parser.add_option('--appendStatisticsOutput',action='store',type='string',dest='appendStatisticsOutput',help='Specify the name of the statistics output file to be appended with sample statistics', default='')
    parser.add_option('--laneSpecificStatistics',action='store_true',dest='laneSpecificStatistics',help='Specify whether the statistics output file should also contain lane-specific statistics', default=False)
    parser.add_option('--statisticsOnly',action='store_true',dest='statisticsOnly',help='Specify whether to skip the data output files and prepare only the statistics files as output', default=False)
    parser.add_option('--deleteTemp',action='store_true',dest='deleteTemp',help='Specify whether to delete or keep the temp folder (with a random number as name in options.tempDir)', default=False)
    parser.add_option('--gzip',action='store_true',dest='gzip',help='Specify whether to store the output files in compressed format (to save disk space)', default=False)       

    parser.add_option('--prepareOnly',action='store_true',dest='prepareOnly',help='Adds all time-consuming analyes to a shell script instead of starting them',default=False)
    parser.add_option('--debugRandomSubset',action='store',type='float',dest='debugRandomSubset',help='Specify whether to run the script only on a random subset of the genome for faster testing (default: 1, corresponding to including each interval with 100% probability ))', default=1)    
    parser.add_option('--debugEnableRuntimeProfiling',action='store_true',dest='debugEnableRuntimeProfiling',help='Specify whether to perform CPU profiling and provide a summary of the runtime performance', default=False)
    parser.add_option('--debugEnableMemoryProfiling',action='store_true',dest='debugEnableMemoryProfiling',help='additional output on memory requirements (heapsize) during read processing', default=False)   
    parser.add_option('-v','--verbose',action='store_true',dest='verbose',help='Print debugging information [default=False]',default=False)
    
    parser.add_option('--minFragmentLength',action='store',type='int',dest='minFragmentLength',help='minimum fragment length reads are assigned to [RRBS protocol only][default:0]', default=0)
    parser.add_option('--maxFragmentLength',action='store',type='int',dest='maxFragmentLength',help='maximum fragment length reads are assigned to [RRBS protocol only][default:MAX_INT]', default=sys.maxint)
    
#    parser.add_option('--customSwitch',action='store_true',dest='customSwitch',help='DO NOT USE UNLESS YOU HAVE A CLUE WHAT IT ACTUALLY DOES AT THE MOMENT. A variable switch that enables a certain debugging option of interest. depending on script status.',default=False)
    #parallelization on the cluster
    parser.add_option('--parallelize',action='store_true',dest='parallelize',help='Subjobs are submitted to the cluster rather than being compute on a single node',default=False)
    parser.add_option('--task',action='store',type='string',dest='task',help='The task for parallelization on the cluster: "main" (Default), "region", "wrapup" or "diagnostics"',default='main')
    parser.add_option('--lsfOpts',action='store',type='string',dest='lsfOpts',help='Only if parallelize: for the main task: options to submit the subjobs to the cluster. Should include the queue. Default:"-q normal_parallel -R \'rusage[mem=15000]\' -r -mig 5" (-oo, -J and dependencies are set according to internal parameters)',default='-q normal_parallel -R "rusage[mem=15000]" -r -mig 5')
    parser.add_option('--lsfOptsWrapup',action='store',type='string',dest='lsfOptsWrapup',help='Only if parallelize and --task=main: options to submit the wrapup job to the cluster. Should include the queue. Default:the paramter will be set to --lsfOpts',default='')
    parser.add_option('--pickleFile',action='store',type='string',dest='pickleFile',help='for the process and the wrapup tasks: the pickle file where variables can be found',default='')
    parser.add_option('--RP_paramPickle',action='store',type='string',dest='RP_paramPickle',help='ONLY FOR REGION TASK: pickle filename for the input file containing the parameters',default='')
    #the following variables are so far only needed for profiling
    parser.add_option('--runtimeProfileFile',action='store',type='string',dest='runtimeProfileFile',help='filename on where to store to runtime profile (only required if debugEnableRuntimeProfiling is enabled)',default='')
    parser.add_option('--mainTaskLogFile',action='store',type='string',dest='mainTaskLogFile',help='filename of the logfile of the main task (only for "diagnostics" task)',default='')
    
    (options,args) = parser.parse_args()
    if options.toolsPath == "":
        options.toolsPath = toolsPath = os.getcwd() + os.sep + "tools" + os.sep
    if not options.toolsPath[-1] == os.sep:
        options.toolsPath += os.sep
    if options.lsfOptsWrapup == '':
        options.lsfOptsWrapup = options.lsfOpts
    
    setGlobals(options)     
    
    if debug:
        print "using pysam version "+ pysam.version.__version__ + " from " + pysam.__file__
    
    if options.task == "main":
        
        if len(options.alignmentFiles) < 1 or len(options.outputDir) < 1 or len(options.sampleName) < 1 or len(options.outGenome) < 1 or len(options.inGenome) < 1:
            print "Mandatory parameters missing. Program will terminate now."
            print "\nYour parameter settings:"
            print options        
            raise SystemExit
        # prepare output directory and filenames
        if not os.path.exists(options.outputDir): os.mkdir(options.outputDir)
        
        if options.debugEnableRuntimeProfiling:
            runtimeProfileDir = options.outputDir + os.sep + "runtime_profiles"
            if not os.path.exists(runtimeProfileDir): os.mkdir(runtimeProfileDir)
            if len(options.runtimeProfileFile) < 1:
                options.runtimeProfileFile = options.outputDir + os.sep + options.sampleName + "runtime_profile_main.cProfile"
            cProfile.run('performAnalysis(options)',options.runtimeProfileFile)
            #in order to analyze the statistics use python commands. example:
            # import pstats
            # stats = pstats.Stats("runtime_profile_overall.cProfile")
            # stats.strip_dirs().sort_stats('time').print_stats()
        else:
            performAnalysis(options)
        print "Program successfully terminating...."  
    elif options.task == "region":
        if options.debugEnableRuntimeProfiling:
            cProfile.run('performAnalysis_region(options)', options.runtimeProfileFile)
        else:
            performAnalysis_region(options)
    elif options.task == "wrapup":
        if options.debugEnableRuntimeProfiling:
            cProfile.run('performAnalysis_wrapup(options)', options.runtimeProfileFile)
        else:
            performAnalysis_wrapup(options)
    elif options.task == "diagnostics":
        performAnalysis_diagnostics(options)
    else:
        raise Exception("invalid task: " + options.task)
    

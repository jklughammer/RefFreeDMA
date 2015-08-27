# RefFreeDMA
Reference Genome Independent Differential DNA Methylation Analysis

RefFreeDMA is a pipeline to perform genome-wide, high-resolution, differential methylation analysis without the need of a reference genome. RefFreeDMA deduces relevant genomic sequences from standard reduced representation bisulfite sequencing (RRBS) reads and concatenates them into a deduced genome which is then used for differential methylation analysis. RefFreeDMA performs all steps of differential DNA methylation analysis from the raw RRBS data over the creation of a deduced reference genome to differential methylation calling.

Quick start
-----------
__Note:__ _The following steps will run RefFreeDMA in linear mode on a small sample data set consisting of severely downsampled RRBS data for human granulocytes (G), lymphocytes (L), and monocytes (M) in four replicates. An exemplary working directory (RefFreeDMA_test) including the raw data, the sample annotation file and the configuration file is provided [here](http://www.biomedical-sequencing.at/bocklab/jklughammer/RefFreeDMA/RefFreeDMA_test.tar.gz). RefFeeDMA should complete within 10 minutes on a desktop computer and produce plots that show clear clustering of the samples by cell type as well as tables reporting differential methylation between granulocytes and lymphocytes. After completion, all output can be found within the [working directory](#reffreedma-results) which in this example is RefFreeDMA_test. RefFreeDMA skips steps if the respective output is already present. Therefore, in order to rerun, RefFreeDMA_test needs to be reset to its original state by running ```reset_test_dir.sh PATH_TO_TESTDIR/RefFreeDMA_test```._ 

####9 steps to test RefFreeDMA:
1\. Download this repository as ZIP or clone it.  
2\. Download the sample [data set](http://www.biomedical-sequencing.at/bocklab/jklughammer/RefFreeDMA/RefFreeDMA_test.tar.gz) and extract it:```tar -xzf RefFreeDMA_test.tar.gz```
3\. Either download and extract the external software bundle ([tools.tar.gz](http://www.biomedical-sequencing.at/bocklab/jklughammer/RefFreeDMA/tools.tar.gz)): ```tar -xzf tools.tar.gz``` or manually install the required [external software](#external-software). 
4\. Set the path to the downloaded working directory (RefFreeDMA_test) in the configuration file (RefFreeDMA_test/meta/RefFreeDMA_test.cfg).
5\. Adjust the [tool paths](#set-tool-paths) in the configuration file (RefFreeDMA_test/meta/RefFreeDMA_test.cfg) or make sure all tools are in your PATH variable.
6\. Install the required [python libraries](#python-libraries) if needed.
7\. Install the required [R packages](#r-packages) if needed.  
```
8\. Run RefFreeDMA.sh.
```
./RefFreeDMA.sh PATH_TO_TESTDIR/RefFreeDMA_test/meta/RefFreeDMA_test.cfg
```
9\. View the most relevant results under PATH_TO_TESTDIR/RefFreeDMA_test/toSelf_filtered_0.08mm_concat/diffMeth

Dependencies
------------

###External software
**SAMtools:** http://www.htslib.org/download/  
**Picard Tools:** http://broadinstitute.github.io/picard/  
**cutadapt:** https://code.google.com/p/cutadapt/  
**trim_galore:** http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/  
**Bowtie2:** http://bowtie-bio.sourceforge.net/bowtie2/index.shtml  
**BSMAP:** https://code.google.com/p/bsmap/  

###System requirements
**Linux 64bit**  
**bash:** http://www.gnu.org/software/bash/  
**Python 2.7:** https://www.python.org/downloads/  
**R/Rscript 3.1.2:** https://cran.r-project.org/  

###Python libraries  
**biopython:** http://biopython.org/DIST/biopython-1.63.zip  
**bitarray:** https://pypi.python.org/packages/source/b/bitarray/bitarray-0.8.1.tar.gz  
**guppy:** https://pypi.python.org/packages/source/g/guppy/guppy-0.1.10.tar.gz  
**pysam:** https://code.google.com/p/pysam/downloads/detail?name=pysam-0.7.5.tar.gz  

```
#install locally (~/.local/lib/) with: 
python2.7 setup.py install --user
```

###R packages
```R
#CRAN
install.packages("VennDiagram")
install.packages("hexbin")
install.packages("data.table")
install.packages("ggplot2")
install.packages("reshape2")
install.packages("MASS")

#Bioconductor
source("http://bioconductor.org/biocLite.R")
biocLite("GenomicRanges",suppressUpdates=TRUE)
biocLite("limma",suppressUpdates=TRUE)
biocLite("Biostrings",suppressUpdates=TRUE)

#Github
install.packages("devtools")
require(devtools)
install_github("sheffien/simpleCache")
```
###Configuration file
The configuration file is a list of key=value pairs that pass mandatory parameters to RefFreeDMA.sh. There are three categories of parameters: [tool paths](#set-toolpaths) that point RefFreeDMA.sh to the required external software if it is not already part of your PATH variable , [default parameters](#adjust-default-parameters-if-required) that might or might not need to be ajusted and [variable parameters](##set-variable-parameters) that are specific to each analysis. An example configuration file comes with the [test data set](#quick-start).

###Reference genomes for cross-mapping
Reference genomes for cross-mapping should be provided as one fasta file containing the chromosomes as separate entries. These genome fasta files can for example be obtained from UCSC Genome Browser or Ensembl. In the [test run](#quick-start) cross-mapping is disabled to avoid this dependency. 

###Sample annotation sheet

The sample annotation sheet has to contain at least two columns: Sample_Name and a collunm that specifies the **two** groups to be compared in the differential methylation analysis. Samples that should not be included in the differential methylation analysis have to be marked with NA in this column. A third column can be specified to indicate groups for plotting purposes. The only predefines column Name is Sample_Name. All other column names can be chosen and passed as paramenters (compCol and groupsCol). The sample annotation sheet is passed as parameter [sample_annotation](#set-command-line-parameters).

| Sample_Name  |comp_gran_lympho|Cell_Type| XXX  | YYY  |
|---|---|---|---|---|
|Sample1|NA|monocyte|   |   |
|Sample2|NA|monocyte|   |   |
|Sample3|granulocyte|granulocyte|   |   |
|Sample4|granulocyte|granulocyte|   |   |
|Sample5|lymphocyte|lymphocyte|   |   |
|Sample6|lymphocyte|lymphocyte|   |   |


###Input files/ working directory
Input files are unmapped BAM files (one per sample). Input files have to be located or linked (ln -s) in a folder called unmapped_bam inside the specified working directory (working_dir). Input file names consist of either the sample name as specified in the sample annotation sheet or prefix (e.g. flowcell information) followed by the sample name. In the latter case prefix and sample name must be separated by a unique character sequence, which is to be specified as nameSeparator ([default parameters](#adjust-default-parameters-if-required)).

Running RefFreeDMA
------------------

###Set tool paths (configuration file)
|path|explanation|
-----|-----------|
|tool_path|Path to a directory containing all the need external software. e.g. the [tools.tar.gz](www.biomedical-sequencing.at/bocklab/jklughammer/RefFreeDMA/tools.tar.gz)|
|picard_path|Needed for sam to fastq conversion|
|trim_galore_path|Needed for read trimming|
|cutadapt_path|Needed for read trimming|
|bowtie2_path|Needed for all agains all alignment (cluster finding)|
|bsmap_path|Needed for bisulfite conversion aware alignment|
|samtools_path|needed for bam/sam conversions, indexing and flagstats|

```bash
tool_path=YOUR_TOOLS_PATH
picard_path=$tool_path/picard-tools_1.118/
trim_galore_path=$tool_path/trim_galore_0.3.3/
cutadapt_path=$tool_path/cutadapt_1.8.3/
bowtie2_path=$tool_path/bowtie_2.2.4/bin/
bsmap_path=$tool_path/bsmap_2.90/
samtools_path=$tool_path/samtools_1.2/bin/
```

###Adjust default parameters if required (configuration file)
|parameter|proposed value|description|
----------|-------|-----------|
|wait_time|10|Check every wait_time minutes weather process is finished (only relevant for parallel mode).|
|nProcesses|1|Maximum number of allowed parallel processes for mapping and methylation calling.|
|nameSeparator|"#"|Unique character(s) that separate flowdell ID from sample name (as indicated in the sample annotation sheet) in the bam-file name. If the bam-file name consists only of the sample name put ""|
|maxReadLen|51|Actual maximum read letngth or if 3' croping is desired for, max length of read after cropping.|
|maxSamples|`ls $bam_dir/*.bam|wc -l`|Per default use all samples that are provided as .bam files in the unmapped_bam folder. Change this parameter to any desired number.|
|filtLim|4|Only consider reads for the generation of the deduced genome that occur at least in 2 out of filtLim samples.|
|cLimit|0.05|Minimum frequency of cytosines at a certain position in order to call a cytosine in the consensus sequence.
|mapToSelf_filter|0.08|Only accept hits as true matches where the mismatch rate is <mapToSelf_filter|
|consensus_dist|0.05|Max accepted mismatch rate between a read and its assigned consensus|
|crossMap_mismatchRate|0.2|Allowed mismatch rate during crossmapping. The parameter is passed to bsmap|
|nTopDiffMeth|500|Maximum number of top differentially methylated fragments to return as fasta sequences|

```bash
wait_time=10
nProcesses=1
nameSeparator="#"
maxReadLen=51
maxSamples=`ls $bam_dir/*.bam|wc -l`
filtLim=4
cLimit=0.05
mapToSelf_filter=0.08
consensus_dist=0.05
crossMap_mismatchRate=0.2
nTopDiffMeth=500
```

###Set variable parameters (configuration file)
|parameter|example|description|
----------|--------|-------|-----------|
working_dir|/home/RefFreeDMA_test|Directory in which analysis is to be perfoormed.|
species|Hum|Identifier that will be part of the differential methylation output files.|
cross_genome_fa|/home/genomes/mm10.fa|Fasta file for genome that should be used for cross-mapping. Set to "" to disable.|
sample_annotation|/home/RefFreeDMA_test/<br>meta/sampleAnnotation.tsv|Sample annotation file.|
compCol|comp_gran_lympho|Column in the sample annotation file that specifies groups of samples for differential methylation analysis.|
groupsCol|Cell_Type|Column in the sample annotation file that specifies groups of samples for plotting. Can be the same as compCol.|
parallel|TRUE or FALSE|Should RefFreeDMA be run in parallel (TRUE) or linear (FALSE) mode|

```bash
working_dir=PATH_TO_TESTDIR/RefFreeDMA_test
species=Hum
cross_genome_fa=-
sample_annotation=$working_dir/meta/sampleAnnotation.tsv
compCol=comp_gran_lympho
groupsCol=Cell_Type
parallel=FALSE
```

###Adapt cluster submission 
If run in parallel mode, RefFreeDMA submits the different tasks to the compute cluster instead of running them locally. The submit commands are valid for SLURM. If you are working with a different resource menagement system you need adjust the submit commands in RefFreeDMA.sh.
The SLURM submit command is of the following structure:
```
sbatch --export=ALL --get-user-env --job-name= --ntasks= --cpus-per-task= --mem-per-cpu= --partition= --time= -e -o
```

###Execute RefFreeDMA
```
./RefFreeDMA.sh configFile.cfg

```

###RefFreeDMA results
```
working_dir
├── crossMapping
|   └── cross_genome
|	    ├──crossMapping.bam
|	    ├──crossMapping.bai
|	    └──crossMapping.flagstat
├── fastq
|   ├── Sample1
|   └── ...
├── log
├── motifAnalysis
|   ├── topXXXGroup1.fa
|   ├── topXXXGroup2.fa
|   └── allDeducedGenomeFragments.fa
├── RCache
|   └── combinedMeth_[species].RData
├── reduced
|   └── consensus
├── toSelf_filtered_XXX_final_concat
|   ├── Sample1
|   |   └── biseqMethcalling
|   ├── ...
|   └── diffMeth
|   	├── diff_meth.tsv
|	    ├── diff_meth_CpG.tsv
|	    ├── mean_meth.tsv
|	    └── plots.pdf
└── unmapped_bam 
    ├── Sample1.bam
    └── ...
```


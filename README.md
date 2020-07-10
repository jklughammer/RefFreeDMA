# RefFreeDMA
Reference Genome Independent Differential DNA Methylation Analysis

RefFreeDMA is a pipeline to perform genome-wide, high-resolution, differential methylation analysis without the need of a reference genome. RefFreeDMA deduces relevant genomic sequences from standard reduced representation bisulfite sequencing (RRBS) reads and concatenates them into a deduced genome which is then used for differential methylation analysis. RefFreeDMA performs all steps of differential DNA methylation analysis from the raw RRBS data over the creation of a deduced reference genome to differential methylation calling.

See also: http://www.ncbi.nlm.nih.gov/pubmed/26673328

Quick start
-----------
__Note:__ The following steps will run RefFreeDMA in linear mode on a small sample data set consisting of severely downsampled RRBS data for human granulocytes (G), lymphocytes (L), and monocytes (M) in four replicates. The exemplary working directory (RefFreeDMA_test) includes the raw data, the sample annotation file, and the configuration file. RefFreeDMA should complete within 10 minutes on a desktop computer and produce plots that show clear clustering of the samples by cell type as well as tables reporting differential methylation between granulocytes and lymphocytes. After completion, all output can be found within the [working directory](#reffreedma-results) which in this example is RefFreeDMA_test. RefFreeDMA skips steps if the respective output is already present. Therefore, in order to rerun, RefFreeDMA_test needs to be reset to its original state by running `reset_test_dir.sh PATH_TO_TESTDIR/RefFreeDMA_test`. 

#### 7 steps to test RefFreeDMA:
1\. Download this repository as ZIP or clone it.  
2\. Download and extract the test data set ([RefFreeDMA_test.tar.gz](http://www.biomedical-sequencing.at/bocklab//papers/klughammer2015/resources/RefFreeDMA_test.tar.gz)): `tar -xzf RefFreeDMA_test.tar.gz`.  
3\. Either download and extract the external software bundle ([tools.tar.gz](http://www.biomedical-sequencing.at/bocklab//papers/klughammer2015/resources/tools.tar.gz)*): `tar -xzf tools.tar.gz` or manually install the required [external software](#external-software).  
4\. Edit the test configuration file (`RefFreeDMA_test/meta/RefFreeDMA_test.cfg`): set YOUR_TOOLS_PATH and PATH_TO_TESTDIR.  
5\. Install the required [R packages](#r-packages) if needed.  
6\. Run RefFreeDMA.sh.  
```
./RefFreeDMA.sh PATH_TO_TESTDIR/RefFreeDMA_test/meta/RefFreeDMA_test.cfg
```
7\. View the most relevant results under `RefFreeDMA_test/toSelf_filtered_0.08mm_concat/diffMeth_cpg`  
(8.) Get summary statistics (`RefFreeDMA_test/summary.txt`) by running: 
```
./scripts/parse_stats.sh PATH_TO_TESTDIR/RefFreeDMA_test/meta/RefFreeDMA_test.cfg
```
(9.) If you encounter problems please check the `RefFreeDMA_test/log` directory. It contains log files for each step.

*The tools included in the external software bundle have been compiled under SL6.

Dependencies
------------
### System requirements
**Linux 64bit**  
**bash:** http://www.gnu.org/software/bash/  
**java:** http://openjdk.java.net/install/  
**Python 2.7:** https://www.python.org/downloads/  
**R/Rscript 3.1.2:** https://cran.r-project.org/  

### External software
__*Tools*__  
**SAMtools:** http://www.htslib.org/download/  
**Picard Tools:** http://broadinstitute.github.io/picard/  
**cutadapt:** https://code.google.com/p/cutadapt/  
**trim_galore:** http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/  
**Bowtie2:** http://bowtie-bio.sourceforge.net/bowtie2/index.shtml  
**BSMAP:** https://code.google.com/p/bsmap/  
**BWA >0.7*:** https://github.com/lh3/bwa  
**bwa-meth*:** https://github.com/brentp/bwa-meth  
__*Python libraries*__  
**biopython:** http://biopython.org/DIST/biopython-1.63.zip  
**bitarray:** https://pypi.python.org/packages/source/b/bitarray/bitarray-0.8.1.tar.gz  
**guppy:** https://pypi.python.org/packages/source/g/guppy/guppy-0.1.10.tar.gz  
**pysam:** https://code.google.com/p/pysam/downloads/detail?name=pysam-0.7.5.tar.gz  
**toolshed*:** https://pypi.python.org/packages/source/t/toolshed/toolshed-0.4.0.tar.gz  

*Optional: These tools and libraries are needed for the optional [decontamination step](#decontamination).  
### R packages
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

#Check if all required packages are installed properly 
library("VennDiagram")
library("hexbin")
library("data.table")
library("ggplot2")
library("reshape2")
library("MASS")
library("GenomicRanges")
library("limma")
library("Biostrings")
library("simpleCache")
```
### Configuration file
The configuration file is a list of key=value pairs that pass mandatory parameters to RefFreeDMA.sh. There are three categories of parameters: [tool paths](#set-tool-paths) that point RefFreeDMA.sh to the required external software if it is not already part of your PATH variable, [default parameters](#adjust-default-parameters-if-required) that might need to be adjusted and [variable parameters](#set-variable-parameters) that are specific to each analysis. An example configuration file comes with the [test data set](#quick-start).

### Reference genomes for cross-mapping
Reference genomes for cross-mapping should be provided as one fasta file containing the chromosomes as separate entries. These genome fasta files can for example be obtained from UCSC Genome Browser or Ensembl. You can run RefFreeDMA on the same working directory for multiple cross-mapping genomes. There will be a separate output folder for each genome. In the [test run](#quick-start) cross-mapping is disabled by default to avoid the genome dependency, but can easily be enabled by passing a genome (cross_genome_fa) as [variable parameter](#set-variable-parameters). 

### Sample annotation sheet
The sample annotation sheet has to contain at least three columns: Sample_Name, a column that specifies the **two** groups to be compared in the differential methylation analysis and a column that specifies which samples are to be used for the generation of the deduced reference (1 = use; 0 = don't use). Samples that should not be included in the differential methylation analysis have to be marked with NA in this column. A fourth column can be specified to indicate groups for plotting purposes. The only two predefined column names are Sample_Name and Select. All other column names can be chosen and passed as paramenters (compCol and groupsCol). The sample annotation sheet is passed as parameter sample_annotation in the [variable parameters](#set-variable-parameters).

| Sample_Name  |comp_gran_lympho|Cell_Type|Select| XXX  | YYY  |
|---|---|---|---|---|---|
|Sample1|NA|monocyte|1|   |   |
|Sample2|NA|monocyte|1|   |   |
|Sample3|granulocyte|granulocyte|0|   |   |
|Sample4|granulocyte|granulocyte|0|   |   |
|Sample5|lymphocyte|lymphocyte|1|   |   |
|Sample6|lymphocyte|lymphocyte|1|   |   |


### Input files/ working directory
Input files are unmapped BAM files (one per sample). Input files have to be located or linked (`ln -s`) in a folder called unmapped_bam inside the specified working directory (working_dir). Input file names consist of either the sample name as specified in the sample annotation sheet or prefix (e.g. flowcell information) followed by the sample name. In the latter case prefix and sample name must be separated by a unique character sequence, which is to be specified as nameSeparator in the [default parameters](#adjust-default-parameters-if-required).

### Analysing DNA methylation in non-CpG context
RefFreeDMA has so far been used and tested extensively only for the analysis of differential DNA methylation in CpG context. However, setting the [nonCpG parameter](#set-variable-parameters) to TRUE will make RefFreeDMA run the same analysis on Cs in CHH and CHG context as for Cs in CpG context. Two more output folders (diffMeth_cphph and diffMeth_cphpg) will be produced. This feature has been tested and validated on simulated data.

### Decontamination
In order to prevent contaminating sequences (originating from microbial species) from being included in the deduced genome and possibly affecting the analysis, RefFreeDMA offers an optional decontamination step. If activated, all reads are mapped to a decoy genome consisting of all sequences in the ncbi blast representative bacterial/archeal genomes database (ftp://ftp.ncbi.nlm.nih.gov/blast/db/Representative_Genomes.*tar.gz) plus the sequence of phiX174. Only reads that don't map to this decoy genome are included in further analysis. For each sample a summary of contaminating sequences as well as the bam file is reported in the fastq directory. To use this functionality the pre-indexed [decoy genome](http://www.biomedical-sequencing.at/bocklab//papers/klughammer2015/resources/decon_reference.tar.gz) (56GB) needs to be downloaded and extracted. The [necessary tools](#external-software) (bwa and bwa-meth) as well as the toolshed python library need to be installed or used from the external software bundle ([tools.tar.gz](http://www.biomedical-sequencing.at//papers/klughammer2015/resources/tools.tar.gz)). And the respective [paths](#set-tool-paths) and [parameters](#adjust-default-parameters-if-required) need to be set in the configuration file.

### Conversion controls (RRBS protocol specific)
Some RRBS protocols include methylated and unmethylated spike-in sequences to be used as conversion controls. They can be analyzed using this bundle of scripts [conversionCrt.tar.gz](http://www.biomedical-sequencing.at/bocklab//papers/klughammer2015/resources/conversionCrt.tar.gz). These scripts expect the environment modules `samtools/1.3`, `bismark/0.12.2`, and `bowtie/1.1.1`. These might need to be adapted to the respective environment. Other versions might not work.

Running RefFreeDMA
------------------

### Set tool paths
**Note:** Tool paths are set in the configuration file.  

|parameter|explanation|
|-----|-----------|
|tool_path|Path to a directory containing all the need external software, e.g. the [tools.tar.gz](http://www.biomedical-sequencing.at/bocklab/papers/klughammer2015/resources/tools.tar.gz)|
|picard_path|Needed for sam to fastq conversion|
|trim_galore_path|Needed for read trimming|
|cutadapt_path|Needed for read trimming|
|bowtie2_path|Needed for all against all alignment (cluster finding)|
|bsmap_path|Needed for bisulfite conversion aware alignment|
|samtools_path|needed for bam/sam conversions, indexing and flagstats|
|bedtools_path|needed for coverage calculation|
|bwa_path|needed for decontamination (mapping)|
|bwameth_path|needed for decontamination (mapping)|
|decon_reference|needed for decontamination (mapping reference)|
```bash
tool_path=YOUR_TOOLS_PATH
picard_path=$tool_path/picard-tools_1.118/
trim_galore_path=$tool_path/trim_galore_0.3.3/
cutadapt_path=$tool_path/cutadapt_1.8.3/
bowtie2_path=$tool_path/bowtie_2.2.4/bin/
bsmap_path=$tool_path/bsmap_2.90/
samtools_path=$tool_path/samtools_1.2/bin/
bedtools_path=$tool_path/bedtools_2.26.0/bin/
bwa_path=$tool_path/bwa_0.7.8/bin/
bwameth_path=$tool_path/bwa-meth-0.10/
decon_reference=YOUR_DECON_PATH/decon_reference/bacterial_extracted_add
```

### Adjust default parameters if required
**Note:** Default parameters are set in the configuration file.  

|parameter|proposed value|description|
|----------|-------|-----------|
|unconv_tag|"_uc"|Unique tag in unmapped bam file names that allows RefFreeDMA to identify samples as unconverted. These samples will be gven priority during reference generation.|
|decon|FALSE|If set to TRUE all reads will be aligned to a microbial decoy genome and only reads that do not align will be used for building the deduced genome/further analysis.|
|restrictionSites|"CGG"|Restriction sites (common 5' motif in RRBS reads). If multiple, specify them separated by pipe e.g. "CGG\|CGA" for MspI\|Taq1
|wait_time|10|Check every wait_time minutes weather process is finished (only relevant for parallel mode).|
|nProcesses|4|Maximum number of allowed parallel processes for mapping and methylation calling.|
|nameSeparator|"#"|Unique character(s) that separate flowcell ID from sample name (as indicated in the sample annotation sheet) in the bam-file name. If the bam-file name consists only of the sample name put ""|
|maxReadLen|51|Actual maximum read length or if 3' cropping is desired for, max length of read after cropping.|
|maxSamples|10-20|Per default use all samples that are provided as .bam files in the unmapped_bam folder. Change this parameter to the number of samples you selected in the sample annotation file (Select column) if the deduced reference is to be produced with only a subset of samples.|
|filtLim|4|Only consider reads for the generation of the deduced genome that occur at least in 2 out of filtLim samples.|
|cLimit|0.05|Minimum frequency of cytosines at a certain position in order to call a cytosine in the consensus sequence.
|mapToSelf_filter|0.08|Only accept hits as true matches where the mismatch rate is < mapToSelf_filter|
|consensus_dist|0.05|Max accepted mismatch rate between a read and its assigned consensus|
|crossMap_mismatchRate|0.2|Allowed mismatch rate during crossmapping. The parameter is passed to bsmap|
|nTopDiffMeth|500|Maximum number of top differentially methylated fragments to return as fasta sequences|

```bash
decon=FALSE
restrictionSites="CGG"
wait_time=10
nProcesses=1
nameSeparator="#"
maxReadLen=51
maxSamples=`ls $working_dir/unmapped_bam/*.bam|wc -l`
filtLim=4
cLimit=0.05
mapToSelf_filter=0.08
consensus_dist=0.05
crossMap_mismatchRate=0.2
nTopDiffMeth=500
```

### Set variable parameters
**Note:** Variable parameters are set in the configuration file.  

|parameter|example|description|
|----------|--------|----------|
|working_dir|/home/RefFreeDMA_test|Directory in which the analysis is to be performed.|
|species|Hum|Identifier that will be part of the differential methylation output files.|
|cross_genome_fa|/home/genomes/mm10.fa|Fasta file for a genome that should be used for cross-mapping. Set to "-" to disable.|
|sample_annotation|/home/RefFreeDMA_test/<br>meta/sampleAnnotation.tsv|Sample annotation file.|
|compCol|comp_gran_lympho|Column in the sample annotation file that specifies groups of samples for differential methylation analysis.|
|groupsCol|Cell_Type|Column in the sample annotation file that specifies groups of samples for plotting. Can be the same as compCol.|
|parallel|TRUE or FALSE|Should RefFreeDMA be run in parallel (TRUE) or linear (FALSE) mode. Parallel mode requires correct cluster submission commands.|
|nonCpG|TRUE or FALSE|Should RefFreeDMA analyse DNA methylation also in CHH and CHG context.|

```bash
working_dir=PATH_TO_TESTDIR/RefFreeDMA_test
species=Hum
cross_genome_fa=-
sample_annotation=$working_dir/meta/sampleAnnotation.tsv
compCol=comp_gran_lympho
groupsCol=Cell_Type
parallel=FALSE
nonCpG=FALSE
```

### Adapt cluster submission 
If run in parallel mode, RefFreeDMA submits the different tasks to the compute cluster instead of running them locally. The submit commands are valid for SLURM. If you are working with a different resource management system you need to adjust the submit commands in RefFreeDMA.sh.
The SLURM submit command is of the following structure:
```
sbatch --export=ALL --get-user-env --job-name= --ntasks= --cpus-per-task= --mem-per-cpu= --partition= --time= -e -o
```

### Execute RefFreeDMA
```
./RefFreeDMA.sh configFile.cfg

```

### RefFreeDMA results
\*\*directories\*\* contain the user-relevant information. The other directories contain intermediate files, that might be needed for reruns and development.
```
working_dir
├── **crossMapping**
|   └── cross_genome
|	    ├──crossMapping.bam
|	    ├──crossMapping.bai
|	    └──crossMapping.flagstat
├── fastq
|   ├── Sample1
|   └── ...
├── log
├── **motifAnalysis**
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
|   └── **diffMeth**
|	    ├── diff_meth.tsv
|	    ├── diff_meth_CpG.tsv
|	    ├── mean_meth.tsv
|	    └── plots.pdf
└── unmapped_bam 
    ├── Sample1.bam
    └── ...
```
### Differential DNA methylation results
RefFreeDMA automatically outputs the nTopDiffMeth differentially methylated deduced reference fragments per group. The default for nTopDiffMeth is set to 500. Because p-values for single CpGs in a deduced reference fragment are combined irrespective of their direction of differential methylation, fragments can display very small or no differences between groups, while still having significant p-values. These heterogeneously differentially methylated regions are likely biologically different from the homogeneously differentially methylated ones and therefore should be analyzed separately.
If the goal is to focus on homogeneously differentially methylated fragments, the number of nTopDiffMeth differentially methylated fragments should be reduced until heterogeneously differentially methylated regions are no longer selected.

### Output statistics
#### Comprehensive statistics
Produced by running: 
```
./scripts/parse_stats.sh PATH_TO_TESTDIR/RefFreeDMA_test/meta/RefFreeDMA_test.cfg
```
**summary.txt** contains all relevant stats concerning each of the samples.\
**ref_summary.txt** contains relevant stats concerning the deduced reference, which is common across all samples. 

#### Header description for summary.txt
header|description
------|-------------------
sample|given sample name
species|given species name
total_reads|number of reads after trimming
mapped_reads|number of mapped reads
mapping_efficiency|% mapped reads
informative_reads|number of qc pass reads with at least 1 motif
CpG_meth|average across all motif methylation values 
avg_meth|read weighted methylation average across all motifs
CpG_measurements|number of times the motif has been detected
coveredCpGs|number of unique motif positions 
conversionRate|based on non-CpG methylation
k1_unmeth|underconversion rate based on spike-in (RRBS protocol specific)
k3_meth|overconversion rate based on spike-in (RRBS protocol specific)
totalMeasurements_k1|number of k1 measurements
totalMeasurements_k3|number of k3 measurements
total_reads_untrimmed|number of reads before trimming
\[CT\]GG|number or reads starting with the restriction site (analogous for other restriction sites)
CGG|% reads starting with a methylated restriction site (analogous for other restriction sites)
TGG|% reads starting with a unmethylated restriction site (analogous for other restriction sites)
others|% reads not starting with a restriction site (indicating pre-fragmentation)
all_bases|number of sequenced bases in the input bam file
As|number of Adenines in the input bam file
Ts|number of Thymines in the input bam file
Cs|number of Cytosines in the input bam file
Gs|number of Guanines in the input bam file
Ns|number of Ns in the input bam file
perc_As|% Adenines in the input bam file
perc_Ts|% Thymines in the input bam file
perc_Cs|% Cytosines in the input bam file
perc_Gs|% Guanines in the input bam file
perc_Ns|% Ns in the input bam file
fragments_ref|number of deduced reference fragments
fragments_uncovered|number of deduced reference fragments that are NOT covered by this sample
fragments_uncovered_perc|% deduced reference fragments that are NOT covered by this sample
max_cont_sp|ncbi ID of species that contributes most to contamination (only if decontamination is run)
max_cont|number of reads assigned to max_cont_sp (only available if decontamination is run)
cont|total number of reads identified as contaminating (only available if decontamination is run)
cont_rat|fraction of reads identified as contaminating (only available if decontamination is run)
blast_count1|number of reads assigned to the most frequent blast hit (only available if [bisulfiteBlast](https://github.com/jklughammer/bisulfiteBlast) is run)
blast_species1|most frequent blast hit species (only available if [bisulfiteBlast](https://github.com/jklughammer/bisulfiteBlast) is run)
blast_count2|number of reads assigned to the second most frequent blast hit (only available if [bisulfiteBlast](https://github.com/jklughammer/bisulfiteBlast) is run)
blast_species2|second most frequent blast hit species (only available if [bisulfiteBlast](https://github.com/jklughammer/bisulfiteBlast) is run)

#### Methylation calling specific statistics
Produced automatically: toSelf_filtered_XXX_final_concat/RRBS_biseq_statistics.txt contains detailed methylation calling statistics produced by the biseqMethCalling.py script. Header descriptions can be found within the script [here](https://github.com/jklughammer/RefFreeDMA/blob/8cfcf83464e6d25e0f7c29e22bb6d7679c2801f3/scripts/biseqMethCalling.py#L2081)



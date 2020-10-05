# HPRep pipeline user manual
## What is HPRep?
HPRep is a methodological framework for determining reproducibility metrics between PLAC-Seq or HiChIP datasets. The pipeline involves three main stages:
1. Preprocessing
2. Regression and normalization
3. Data matrix comparisons
 
The preprocessing of the first stage is borrowed directly from the MAPS pipeline (https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006982). This preprocessing pipeline, called feather, converts aligned, sorted, and merged pair-end reads to long and short .bed/.bedpe files. The short reads are used for determining ChIP efficiency and the long reads are used for reproducibility analysis. The second stage utilizes positive poisson regression to fit contact count predictive models to obtain expected counts for the observed data. Normalization involves reporting log<sub>2</sub>(1 + observed / expected). The last stage requires pairwise comparisons of all samples. Consequently, while the first two stages can be conducted in parallel, the last stage requires all results from stages 1 and 2.

The requirements and details for running the entire pipeline are provided below:

HPRep pipeline runs on Linux and requires several readily available programs and packages. We list the versions used in testing and development with the caveat that (slightly) older and newer versions should work, but caveat utilitor.

**python 3.6.6**
* pandas 0.24.2
* numpy 1.19.0
* pysam 0.15.2\*
* pybedtools 0.8.1

\* pysam 0.16 does not work!

**R 3.6.0**
* MASS 7.3.51.4
* data.table 1.12.8
* VGAM 1.1.3

**samtools 1.11**

**bedtools 2.29**

## Inputs
1. Bam file derived from bwa aligned, sorted, and merged fastq files
2. Genomic features file
3. 1-D ChIP peaks from MACS2

The bam file can be generated using bwa and samtools starting from the individual paired end fastq files using
```
bwa mem -t threads index_file fastq_file1 > sample_R1.bwa
bwa mem -t threads index_file fastq_file2 > sample_R2.bwa
```
where threads in the number of threads available. Consult bwa documentation for building the appropriate index for the species / genomic build being analyzed. Sort both aligned files using
```
samtools sort -o sample_R1.srtn.bam -n -@ threads sample_R1.bwa
samtools sort -o sample_R2.srtn.bam -n -@ threads sample_R2.bwa
```
and finally merge them using
```
samtools merge -n -f sample.srtn.merged.bam sample_R1.bwa sample_R2.bwa
```
Note that while bwa does allow for alignment of paired ends the feather preprocessing pipeline has been optimized for PLAC-Seq and Hi-ChIP data, hence the above procedure. Additionally, while these steps could be incorporated directly into our pipeline, they are left for the user to allow for parallelization as they can be time consuming for large datasets.

Genomic features files can be downloaded from http://enhancer.sdsc.edu/yunjiang/resources/genomic_features/.

## Running stages 1 and 2
The first two stages of the pipeline can be run by copying the bash script run_pipeline_stages_1_2.sh and renaming it run_pipeline_stages_1_2_sample_name.sh, then executing it after editing the following 12 fields:

* python_path - how execute python (python3, e.g.)
* Rscript_path - how you execute R scripts (Rscript, e.g.)
* input_bam_file - full path to input file
* genomic_feat_file - full path to input genomic features file
* macs2_filepath - full path to input MACS2 peaks file
* dataset_name - name of sample
* outdir - name of output directory (recommend using sample name)
* bin_size - bin resolution (typically 5000 or 10000)
* binning_range - maximum bin distance considered (default is 1000000)
* length_cutoff - cutoff length for determining short/long splitting
* chr_count - chromosome count for organism (19 or 22)
* threads - number of theads available for processing

Upon completion there will be two key outputs: 
```
outdir/HPRep_output/dataset_name_date/dataset_name.resolution.anchors.txt
outdir/HPRep_output/dataset_name_date/dataset_name.resolution.normalized.txt
```
When all samples have been run copy both of these files for each sample into a common directory in preparation for stage 3.

## Running stage 3
This stage can be run by copying the bash script run_pipeline_stage_3.sh and renaming it run_pipeline_stage_3_study_name.sh, then executing it after editing the following 7 fileds:

* Rscript_path - same as previous
* dir_name - the directory containing all .normalized.txt and .anchors.txt files
* tune_sample_1 - the sample name of one of the tuning samples, excluding .normalized.txt (see below)
* tune_sample_2 - the sample name of other tuning sample, excluding .normalized.txt (see below)
* chr_count - same as previous
* bin_size - same as previous
* output_name - name of the study, used for final output naming
* seed - sets the tuning seed for reproducibility

The first step of the process tunes the smoothing parameter. The user specifies which samples to use for tuning. It is recommended to use samples that are NOT biological replicates. Note: the sample name should NOT include .normalized.txt.

The final output will be an (n choose 2) x (p + 2) matrix, where n is the number of samples and p is the number of chromosomes. Each row represents a specific pair of samples, the first two columns will be sample names and the remaining columns will be the corresponding reproducibility metric for each chromosome. The final output will be dir_name/ouput_name.results.txt.



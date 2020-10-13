# HPRep user manual
## What is HPRep?
HPRep is a methodological framework to quantify reproducibility between PLAC-Seq or HiChIP datasets. The pipeline involves three main stages:
1. Preprocessing
2. Regression and normalization
3. Data matrix comparisons
 
The first stage, preprocessing stage, is borrowed directly from the [<em>MAPS</em> pipeline](https://github.com/ijuric/MAPS) (detailed in [MAPS paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006982)). This preprocessing pipeline, called <em>feather</em>, converts aligned, sorted, and merged paired-end reads to long-range and short-range .bed/.bedpe files. The second stage utilizes positive Poisson regression estimates expected counts from the observed data, and derives the normalized contact: log<sub>2</sub>(1 + observed / expected). The last stage calculates pairwise reproducibility or similarity statistic among all samples. In general, we recommend processing samples in parallel in the first two stages and processing all samples together in the last stage.

The requirements and details for running the entire pipeline are provided below:

The HPRep pipeline runs on Linux and requires several readily available programs and packages. We list the versions used in testing and development with the caveat that (slightly) older and newer versions should work, but caveat utilitor.

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
1. Bam file: can be obtained from fastq files by applying <em>bwa</em> align, sort, and merge (detailed below)
2. Genomic features file: can be downloaded from [Genomic features](http://enhancer.sdsc.edu/yunjiang/resources/genomic_features/)
3. 1-D ChIP peaks: can be obtained by running <em>MACS2</em> on corresponding ChIP-seq data, or on short-range reads

The bam file can be generated using <em>bwa</em> and <em>samtools</em> starting from the individual paired end fastq files using
```
bwa mem -t threads index_file fastq_file1 > sample_R1.bwa
bwa mem -t threads index_file fastq_file2 > sample_R2.bwa
```
where threads in the number of threads available. Consult bwa documentation for building or finding the appropriate index for the species / genomic build being analyzed. Sort both aligned files using
```
samtools sort -o sample_R1.srtn.bam -n -@ threads sample_R1.bwa
samtools sort -o sample_R2.srtn.bam -n -@ threads sample_R2.bwa
```
and finally merge them using
```
samtools merge -n -f sample.srtn.merged.bam sample_R1.bwa sample_R2.bwa
```
Prior users of the <em>MAPS</em> procedure will realize that these steps could be incorporated directly into our pipeline, however here they are left for the user to perform separately to allow for parallelization as they can be time consuming for large datasets.

## Running stages 1 and 2
The first two stages of the pipeline can be run by copying the bash script <em>run_pipeline_stages_1_2.sh</em> and renaming it <em>run_pipeline_stages_1_2_sample_name.sh</em>, then executing it after editing the following fields:

* python_path - how execute python (python3, e.g.)
* Rscript_path - how you execute R scripts (Rscript, e.g.)
* feather - 1 to run full pipeline, 0 to intersect MAPS (see below)
* input_bam_file - full path to input file
* genomic_feat_file - full path to input genomic features file
* macs2_filepath - full path to input MACS2 peaks file
* sample_name - name of sample
* outdir - name of output directory (recommend using sample name)
* bin_size - bin resolution (typically 5000 or 10000)
* binning_range - maximum bin distance considered (default is 1000000)
* length_cutoff - cutoff length for determining short/long splitting
* chr_count - chromosome count for organism (19 or 22)
* threads - number of theads available for processing

***For users who have previously run the MAPS pipeline***: You can avoid the <em>feather</em> pre-processing steps by setting feather to 0. You will then need to set the long_bedpe_dir and short_bed_dir fields to the directory or directories where these files are located.
  
Upon completion there will be two key outputs for each dataset/sample: 
```
outdir/HPRep_output/dataset_name_date/dataset_name.resolution.anchors.txt
outdir/HPRep_output/dataset_name_date/dataset_name.resolution.normalized.txt
```
When all samples have been run copy both of these files for each sample into a common directory in preparation for stage 3.

## Running stage 3
This stage can be run by copying the bash script <em>run_pipeline_stage_3.sh</em> and renaming it <em>run_pipeline_stage_3_study_name.sh</em>, then executing it after editing the following 9 fileds:

* Rscript_path - same as previous
* dir_name - the directory containing all .normalized.txt and .anchors.txt files
* tune_sample_1 - the sample name of one of the tuning samples, excluding .normalized.txt (see below)
* tune_sample_2 - the sample name of other tuning sample, excluding .normalized.txt (see below)
* chr_count - same as previous
* bin_size - same as previous
* binning_range - same as previous
* output_name - name of the study, used for final output naming
* seed - tuning seed for reproducibility

The first step of the process tunes the smoothing parameter. The user specifies which samples to use for tuning. It is recommended to use samples that are NOT biological replicates. Note: the sample name should NOT include .normalized.txt.

The final output will be an (<em>n</em> choose 2) x (<em>p</em> + 2) matrix, where <em>n</em> is the number of samples and <em>p</em> is the number of chromosomes. Each row represents a specific pair of samples, the first two columns will be sample names and the remaining columns will be the corresponding reproducibility metric for each chromosome. The final output will be dir_name/ouput_name.results.txt.




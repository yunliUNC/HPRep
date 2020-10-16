#!/bin/bash

#######################################
###INPUT VARIABLES
#######################################
Rscript_path="Rscript"
dir_name="../example/"
tune_sample_1="mESC_sample1.10k.normalized.txt.gz"
tune_sample_2="mouse_brain_sample1.10k.normalized.txt.gz"
chr_count="19"
bin_size="10000"
binning_range="1000000"
output_name="example"
seed="12345"
#######################################
$Rscript_path HPRep/HPRep.R $dir_name $tune_sample_1 $tune_sample_2 $chr_count $bin_size $binning_range $output_name $seed





#!/bin/bash

#######################################
###INPUT VARIABLES
#######################################
Rscript_path="Rscript"
dir_name="/nas/longleaf/home/jdrosen/scr/MAPS/files_from_HPRep/scc/"
tune_sample_1="IJ050_merged_trimmedto100_fastp.5k"
tune_sample_2="IJ051_merged_added_trimmed100_fastp.5k"
chr_count="22"
bin_size="5000"
binning_range="1000000"
output_name="test"
seed="1"
#######################################
$Rscript_path HPRep/HPRep.R $dir_name $tune_sample_1 $tune_sample_2 $chr_count $bin_size $binning_range $output_name $seed





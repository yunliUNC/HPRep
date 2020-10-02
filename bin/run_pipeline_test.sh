#!/bin/bash
python_path=python3 #should have pysam, pybedtools installed. bedtools, samtools should be in the path
Rscript_path=Rscript
###################################################################
feather=1 #start from feather or run only MAPS
HPRep=1
dataset_name="test_set1"
input_bam_file="/pine/scr/j/d/jdrosen/MAPS/bin/test_set1/feather_output/test_set1_20201001_133518/tempfiles/test_set1.merged.srtn.bam"
outdir="test_set1"
macs2_filepath="/nas/longleaf/home/jdrosen/scr/MAPS/examples/test_set1/macs2_peaks_final.replicated.narrowPeak"
#organism="mm10"
bin_size=10000
binning_range=1000000
fdr=2 # this is used just for labeling. do not change
filter_file="None"
mapq=30
chr_count=19
length_cutoff=1000
threads=1
model="pospoisson"
sex_chroms_to_process="NA"
genomic_feat_filepath="/pine/scr/j/d/jdrosen/MAPS/MAPS_data_files/mm10/genomic_features/F_GC_M_MboI_10Kb_el.mm10.txt"
optical_duplicate_distance=0
####################################################################
###SET THE VARIABLES AT THIS PORTION ONLY IF 
### number_of_datasets > 1 (merging exisitng datasets)
### specify as many datasets as required
####################################################################
dataset1=""
dataset2=""
dataset3=""
dataset4=""
#...
##################################################################
###SET THESE VARIABLES ONLY IF FEATHER = 0 AND YOU WANT TO RUN
###USING A SPECIFIC FEATHER OUTPUT RATHER THAN $datasetname_Current
###################################################################
feather_output_symlink=""
##################################################################

DATE=`date '+%Y%m%d_%H%M%S'`
#####Armen:
feather_output=$outdir"/feather_output/"$dataset_name"_"$DATE
if [ "$feather_output_symlink" == "" ]; then
	feather_output_symlink=$outdir"/feather_output/"$dataset_name"_current"
fi
resolution=$(bc <<< "$bin_size/1000")
per_chr='True' # set this to zero if you don't want per chromosome output bed and bedpe files
feather_logfile=$feather_output"/"$dataset_name".feather.log"
cwd="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

long_bedpe_dir=$feather_output_symlink"/"
short_bed_dir=$feather_output_symlink"/"
HPRep_output=$outdir"/HPRep_output/"$dataset_name"_"$DATE"/"
HPRep_output_symlink=$outdir"/HPRep_output/"$dataset_name"_current"

if [ $feather -eq 1 ]; then
	mkdir -p $feather_output
	$python_path $cwd/feather/feather_pipe preprocess -o $feather_output -p $dataset_name -i $input_bam_file -q $mapq -l $length_cutoff -t $threads -c $per_chr -a $macs2_filepath -d $optical_duplicate_distance
	qc_filename=$feather_output/$dataset_name".feather.qc"
	temp_qc_file=$feather_output/tempfiles/$dataset_name".feather.qc.modified"
	sed -r 's/  +/\t/g' $qc_filename > $temp_qc_file
	sed -r -i 's/ /\_/g' $temp_qc_file
	cut -f 1-2 $temp_qc_file > $temp_qc_file".cut"
	sed -i 's/\_$//g' $temp_qc_file".cut"
	paste -d"\t" $cwd/feather/qc_template.txt $temp_qc_file".cut" > $temp_qc_file".cut.tmp"
	awk '{print $1,$3,$2}' $temp_qc_file".cut.tmp" >> $qc_filename".tsv"
	sed -i 's/ /\t/g' $qc_filename".tsv"
	cp "$(readlink -f $0)" $feather_output"/execution_script_copy"
	chmod 777 $feather_output
	ln -sfn $feather_output $feather_output_symlink
fi

	
if [ $HPRep -eq 1 ]; then
	mkdir -p $HPRep_output
	echo "$dataset_name $HPRep_output $macs2_filepath $genomic_feat_filepath $long_bedpe_dir $short_bed_dir $bin_size $chr_count $maps_output $sex_chroms_to_process"
	$python_path $cwd/HPRep/make_and_xor_runfile.py $dataset_name $HPRep_output $macs2_filepath $genomic_feat_filepath $feather_output"/" $feather_output"/" $bin_size $chr_count $HPRep_output $sex_chroms_to_process --BINNING_RANGE $binning_range
	#echo "first"
	$python_path $cwd/HPRep/split_and_xor.py $HPRep_output"HPRep_"$dataset_name".split"
	$Rscript_path $cwd/HPRep/find_anchor_bins.R $macs2_filepath $bin_size $chr_count ${HPRep_output}$dataset_name"."$resolution"k.anchors"
	$Rscript_path $cwd/HPRep/normalize_bins.R $HPRep_output reg_raw $dataset_name"."$resolution"k" $bin_size $chr_count "None" "pospoisson" $HPRep_output 
	#echo "second"
	#$Rscript_path $cwd/MAPS/MAPS_regression_and_peak_caller.r $maps_output $dataset_name"."$resolution"k" $bin_size $chr_count$sex_chroms $filter_file $model 
	#$Rscript_path $cwd/MAPS/MAPS_peak_formatting.r $maps_output $dataset_name"."$resolution"k" $fdr $bin_size
	#echo "third"
	#cp "$(readlink -f $0)" $maps_output"/execution_script_copy"
	#chmod 777 $maps_output
	#ln -sfn $maps_output $maps_output_symlink
fi

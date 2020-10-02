import subprocess
import argparse
import os.path
import pandas
import pysam
from io import StringIO
import itertools
import numpy as np
import random
import re
from tempfile import TemporaryFile
import time
import copy
import logger
import sys
from shutil import copyfile

sam_columns = {"query_name" : 0, "flag" : 1, "chr_name" : 2, "pos" : 3, "mapq" : 4,
		"cigar" : 5, "mate_name" : 6, "mate_pos" : 7, "tlen" : 8, "seq" : 9, "quality" : 10}

flag_table_proper={(0, 16): (99, 147), (2048, 16): (99, 147), (0, 2064):  (99, 147), (2048, 2064): (99, 147),     # + -
			(16, 0): (83, 163), (16, 2048): (83, 163), (2064, 0):  (83, 163), (2064, 2048): (83, 163),     # - +
			(0,  0): (65, 129), (0,  2048): (65, 129), (2048, 0):  (65, 129), (2048, 2048): (65, 129),     # + +
			(16, 16):(113, 177), (16, 2064):(113, 177), (2064, 16):(113, 177), (2064, 2064):(113, 177),    # - -
		}
qseq = "K00168:100:HF52YBBXX:4:1101:27336:1103"

def filter_main(sorted_merged_bam, mapq, outdir, prefix, threads, optical_duplicate_distance, to_file = False):
	sys.stdout = logger.Logger(outdir + "/" + prefix + ".feather.log")
	paired_filename, qc_filename = set_filenames(outdir, prefix)
	#get read count for qc file
	proc = subprocess.Popen("samtools view " + sorted_merged_bam + " | awk ' $1 !~ /@/ {print $1}' " + "| uniq -c|wc -l", stdout = subprocess.PIPE, shell = True)
	read_count = proc.stdout.read().decode("utf-8")

	#pairing and filtering alignments for chimeric reads
	print(time.ctime() + " filtering and pairing reads")
	filter_pair_reads(sorted_merged_bam, mapq, paired_filename, qc_filename)
	print(time.ctime() + " paired bam file generated. Calling fixmate")
	print(paired_filename)
	print(os.path.isfile(paired_filename + ".bam")) 
	pysam.fixmate("-m", "-@", str(threads), (paired_filename + ".bam"), (paired_filename + ".fixmated.bam"))
	print(time.ctime() + " sorting by coordinates.")
	pysam.sort("-o", paired_filename + ".srt.bam", "-@", str(threads), paired_filename + ".fixmated.bam")
	print(time.ctime() + " calling samtools markdup")
	proc = subprocess.Popen(" ".join(["samtools markdup", "-r -m s -s -f", paired_filename + ".fixmated.markdup.stats", "-@", str(threads), "-d", str(optical_duplicate_distance), paired_filename + ".srt.bam", paired_filename + ".rmdup.bam"]), shell = True)
	proc.communicate()
	print(time.ctime() + " calling samtools flagstat on mapped file")
	proc = subprocess.Popen("samtools flagstat " + paired_filename + ".srt.bam > " + paired_filename + ".srt.bam.flagstat", 
				shell = True)
	proc.communicate()
	with open(paired_filename + ".srt.bam.flagstat") as flag_file:
		lines = flag_file.readlines()
		uniquely_mapped_count = lines[7].split()[0]
	print(time.ctime() + " calling samtools flagstat on mapped and duplicate-removed file")
	proc = subprocess.Popen("samtools flagstat " + paired_filename + ".rmdup.bam > " + paired_filename + ".rmdup.flagstat", 
				shell = True)
	proc.communicate()
	with open(paired_filename + ".rmdup.flagstat") as flag_file:
		lines = flag_file.readlines()
		duprmd_count = lines[7].split()[0]
		intra_count = lines[11].split()[0]
		intra_count = str(int(float(intra_count)) / 2)
	print(time.ctime() + " calling samtools sort for sorting by query names")
	pysam.sort("-o", paired_filename + ".srtn.rmdup.bam", 
				"-@", str(threads), "-n", paired_filename + ".rmdup.bam")
	print(time.ctime() + " finishing filtering")
	qc_filename = outdir + "/" + prefix + ".feather.qc"
	with open(qc_filename, 'w') as outfile:
		outfile.write("{0:70} {1}".format("number of sequencing pairs", str(read_count)))
		outfile.write("{0:70} {1} ".format("number of uniquely mapped pairs (MAPQ >= " + str(mapq) + ")", str(uniquely_mapped_count)))
		outfile.write("\t({0:.2f}%)\n".format(100 * (int(float(uniquely_mapped_count)) / int(float(read_count)))))
		outfile.write("{0:70} {1} ".format("number of pairs after duplicate removal", str(duprmd_count)))
		outfile.write("\t({0:.2f}%)\n".format(100 * (int(float(duprmd_count)) / int(float(read_count)))))
	return (paired_filename + ".srtn.rmdup.bam")

def is_sorted_queryname(header):
	if("HD" in header):
		if("SO" in header["HD"]):
			if(header["HD"]["SO"] == "queryname"):
				return True
	return False


def extract_chr_list(header):
	chrs = []
	for chr_line in header["SQ"]:
		if "SN" in chr_line:
			chrs.append(chr_line["SN"])
	return(chrs)

def filter_pair_reads(bwa_filename, mapq, paired_filename, qc_filename):
	merged_sam = pysam.AlignmentFile(bwa_filename, "rb")
	with pysam.AlignmentFile(paired_filename + ".bam", "wb", header = merged_sam.header) as outfile:
		current = ""
		counter = 0
		aligned_reads = []
		paired_reads = []
		unmapped_read_count = 0
		total_read_count = 0
		chimeric_read_count = 0
		chrs = set()
		for read in merged_sam.fetch(until_eof = True):
			total_read_count += 1
			if read.flag == 4 or read.mapq < mapq:
				unmapped_read_count += 1
				continue
			else:
				if read.query_name != current:
					if current != "":
						if counter == 2:
							paired_reads.append(pair_up(aligned_reads))
						elif counter == 3:
							chimeric_read_count += 1
							chr_count = len(chrs)
							if chr_count == 1 or chr_count == 2:
								selected_pair = select_valid_pair(aligned_reads, counter, chr_count)
								paired_reads.append(pair_up([aligned_reads[i] for i in selected_pair]))
						for i, j in paired_reads:
							outfile.write(i)
							outfile.write(j)
						paired_reads = []
					aligned_reads = [read]
					counter = 1
					chrs = set()
					chrs.add(read.reference_name)
					current = read.query_name
				else:
					counter += 1
					chrs.add(read.reference_name)
					aligned_reads.append(read)
		if counter == 2:
			paired_reads.append(pair_up(aligned_reads))
		elif counter == 3:
			chimeric_read_count += 1
			selected_pair = select_valid_pair(aligned_reads, counter, len(chrs))
			paired_reads.append(pair_up([aligned_reads[i] for i in selected_pair]))
		for i, j in paired_reads:
			outfile.write(i)
			outfile.write(j)
	with open(qc_filename, "a") as qc_fout:
		output = " ".join(["Total reads:", str(total_read_count)])
		output += "\n" + " ".join(["Unmapped/low quality reads:", str(unmapped_read_count)])
		output += "\n" + " ".join(["Chimeric reads (after low quality read removal):", str(chimeric_read_count)])
		qc_fout.write(output)

def set_filenames(outdir, prefix):
	tempdir = outdir + "/tempfiles"
	if not os.path.exists(tempdir):
		os.makedirs(tempdir)
	paired_filename = outdir + "/" + prefix + ".paired"
	qc_filename = outdir + "/" + prefix + ".qc"
	return(paired_filename, qc_filename)

def set_tempfile(input_content = None, output_content = None, binary = True):
	tfile = TemporaryFile("wb") if (binary) else TemporaryFile("w")
	if (input_content):
		tfile.write(input_content)
		tfile.seek(0)
	return(tfile)
		

def bwa_mem(fastq, bwa_index, threads, output_filename):
	print(time.ctime() + " calling bwa for " + fastq)
	output_file = open(output_filename, "w")
	proc = subprocess.Popen(["bwa", "mem", "-t", str(threads), bwa_index, fastq], stdout = output_file, stderr = open(output_filename + ".log", 'w'))
	proc.wait()
	output_file.close()

def select_valid_pair(aligned_reads, count, chr_count):
	if (chr_count == 1):
		combinations = [(0, 1), (0, 2), (1, 2)]
		dists = []
		for i, j in combinations:
			dists.append(abs(aligned_reads[i].reference_start - aligned_reads[j].reference_start))
		arg_order = np.argsort(dists)
		second_largest_pair = combinations[arg_order[1]]
		pair_indices = second_largest_pair
	elif (chr_count == 2):
		if aligned_reads[0].reference_name == aligned_reads[1].reference_name:
			base_index = 2
			matched_indices = [0, 1]
		elif aligned_reads[0].reference_name == aligned_reads[2].reference_name:
			base_index = 1
			matched_indices = [0, 2]
		else:
			base_index = 0
			matched_indices = [1, 2]
		pair_indices = sorted([base_index, min(matched_indices)])
	return (pair_indices)
		
def pair_up(read_pair):
	r1_cp = pysam.AlignedSegment()
	r2_cp = pysam.AlignedSegment()
	r1_cp = copy.deepcopy(read_pair[0])
	r2_cp = copy.deepcopy(read_pair[1])
	if r1_cp.query_name != r2_cp.query_name:
		print("Error: read name unmathced.\n")
		sys.exit(1);
	# change flag
	flag_swag = flag_table_proper[(r1_cp.flag, r2_cp.flag)]
	r1_cp.flag = flag_swag[0]
	r2_cp.flag = flag_swag[1]
	# now change RNEXT and PNEXT
	if(r1_cp.reference_name  == r2_cp.reference_name):
		r1_cp.next_reference_start = r2_cp.reference_start
		r2_cp.next_reference_start = r1_cp.reference_start
		r1_cp.template_length  = r1_cp.next_reference_start -  r1_cp.reference_start
		r2_cp.template_length  = -r1_cp.template_length
		r1_cp.next_reference_name = r2_cp.next_reference_name = "="
	else:
		r1_cp.next_reference_name, r2_cp.next_reference_name = r2_cp.reference_name, r1_cp.reference_name
		r1_cp.next_reference_start = r2_cp.reference_start
		r2_cp.next_reference_start = r1_cp.reference_start
		r1_cp.template_length  = r1_cp.next_reference_start -  r1_cp.reference_start
		r2_cp.template_length  = -r1_cp.template_length
	return (r1_cp, r2_cp)

def check_requirements():
	for requirement in requires:
		status_which, result = subprocess.getstatusoutput("which " + requirement)
		status_ls = os.path.exists("utils/" + requirement)
		if (status_which % 256 != 0 and not status_ls):
			exit("ERROR: " + requirement + " cannot be found. Please make sure it is " +
			"installed and is in the system path. Exiting!")

def check_file_existance(files):
	for file_name in files:
		if not os.path.exists(file_name):
			exit("ERROR: Input file " + file_name + " does not exist. Please make sure " +
			"the correct path is provided. Exiting!")


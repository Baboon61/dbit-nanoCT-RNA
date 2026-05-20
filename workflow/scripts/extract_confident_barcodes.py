#The script is calculating a penalty score for each barcpde and extracting reads names with lower penalty score than a given cutoff
#The penalty score of each barcode/read is calculated by summing the Qscore differences between the highest base quality Qscore and each base Qscore

import sys
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import gzip

def open_maybe_gzip(path, mode):
  if path.endswith(".gz"):
    return gzip.open(path, mode)
  return open(path, mode)

#Check integer arguments
def check_positive(value):
  ivalue = int(value)
  if ivalue <= 0:
    raise argparse.ArgumentTypeError("%s is an invalid positive int value" % value)
  return ivalue

#Decode Phred score
def decode(c):
  return ord(c) - 33

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--input", required=True, help="input file")
ap.add_argument("-o", "--output", required=True, help="output file")
ap.add_argument("-p", "--primer_length", required=True, help="PCR primer sequence length", type=check_positive)
ap.add_argument("-l1", "--linker1_length", required=True, help="linker1 sequence length", type=check_positive)
ap.add_argument("-l2", "--linker2_length", required=True, help="linker2 sequence length", type=check_positive)
ap.add_argument("-s1", "--spatial_barcode1_length", required=True, help="spatial barcode1 sequence length", type=check_positive)
ap.add_argument("-s2", "--spatial_barcode2_length", required=True, help="spatial barcode2 sequence length", type=check_positive)
ap.add_argument("-c", "--cutoff", nargs='?', default=10, help="percentage of reads extracted (Default : 10)", type=check_positive)
args = vars(ap.parse_args())

input_file = args["input"]
output_file = args["output"]
primer_length = int(args["primer_length"])
linker1_length = int(args["linker1_length"])
linker2_length = int(args["linker2_length"])
spatial_barcode1_length = int(args["spatial_barcode1_length"])
spatial_barcode2_length = int(args["spatial_barcode2_length"])
cutoff = int(args["cutoff"])

#Initialize Phred score directory
phred={}
for x in range(0,94):
  phred.update({x : (chr(x+33).encode('ascii'))})

#Estimate the highest base Qcore of the file
nb_line=0
top_qscore=0
with open_maybe_gzip(input_file, "rt") as in_handle:
  for title, seq, qual in FastqGeneralIterator(in_handle):
    tmp_top_qscore = max(list(map(decode, qual)))
    if tmp_top_qscore > top_qscore:
      top_qscore = tmp_top_qscore
    nb_line=nb_line+1
    if nb_line == 10000:
      break

print("The top Qscore of the fastq file is : Q%i" % (int(top_qscore)))

#Set up spatial barcodes positions
seq_start=primer_length+spatial_barcode2_length+linker2_length+spatial_barcode1_length+linker1_length

bc2_start=primer_length+1
bc2_end=bc2_start+spatial_barcode2_length-1

bc1_start=bc2_end+linker2_length+1
bc1_end=bc1_start+spatial_barcode1_length-1

#Extract the spatial barcodes sequences and qualities and calculate the penalty scores
count=0
read_penalty_dict = {}
with open_maybe_gzip(input_file, "rt") as in_handle:
  for title, seq, qual in FastqGeneralIterator(in_handle):
    pxl_bcd_seq = seq[bc2_start-1:bc2_end] + seq[bc1_start-1:bc1_end] # !!! BC2 + BC1
    pxl_bcd_qual = qual[bc2_start-1:bc2_end] + qual[bc1_start-1:bc1_end]

    score_list = list(map(decode, pxl_bcd_qual))
    penalty_score = sum([(top_qscore-i) for i in score_list])
    read_penalty_dict[title.split(" ")[0]] = penalty_score

quantile_cutoff_penalty = np.quantile(list(read_penalty_dict.values()), (cutoff/100))
best_read_penalty_dict = {key:val for key, val in read_penalty_dict.items() if val <= quantile_cutoff_penalty}

cutoff_penalty = max(best_read_penalty_dict.values())
percent_read_kept = len(best_read_penalty_dict.keys())*100/len(read_penalty_dict.keys())
print("A barcode penalty score of maximum %i was set to extract %f of the reads" % (cutoff_penalty, percent_read_kept))

print("A violin plot of the barcodes penalty will be output at the output file location")
plt.violinplot(read_penalty_dict.values(), showmedians=True, quantiles=[0.25, 0.75])
plt.ylabel("Barcode penalty score", fontsize=10)
plt.title("Distribution of barcodes penalty score over a fastq file", fontsize=10)
plt.suptitle('Selected Cutoff : ' + str(cutoff_penalty) + ' extracting ' + str(percent_read_kept) + ' percent of the reads', fontsize=8)
plt.savefig(output_file + '.png', dpi=300)

with open(output_file, "wt") as out_handle:
  out_handle.write('\n'.join(best_read_penalty_dict.keys()))
  out_handle.write('\n')

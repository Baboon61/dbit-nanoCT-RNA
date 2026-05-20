from Bio.SeqIO.QualityIO import FastqGeneralIterator
import gzip
import sys

import argparse

def open_maybe_gzip(path, mode):
  if path.endswith(".gz"):
    return gzip.open(path, mode)
  return open(path, mode)

def check_positive(value):
  ivalue = int(value)
  if ivalue <= 0:
    raise argparse.ArgumentTypeError("%s is an invalid positive int value" % value)
  return ivalue

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--input", required=True, help="input file")
ap.add_argument("-o", "--output", required=True, help="output file")
ap.add_argument("-p", "--primer_length", required=True, help="PCR primer sequence length", type=check_positive)
ap.add_argument("-l1", "--linker1_length", required=True, help="linker1 sequence length", type=check_positive)
ap.add_argument("-l2", "--linker2_length", required=True, help="linker2 sequence length", type=check_positive)
ap.add_argument("-s1", "--spatial_barcode1_length", required=True, help="spatial barcode1 sequence length", type=check_positive)
ap.add_argument("-s2", "--spatial_barcode2_length", required=True, help="spatial barcode2 sequence length", type=check_positive)
ap.add_argument("-n", "--umi_barcode_length", required=True, help="UMI barcode sequence length", type=check_positive)
args = vars(ap.parse_args())

input_file = args["input"]
output_file = args["output"]
primer_length = int(args["primer_length"])
linker1_length = int(args["linker1_length"])
linker2_length = int(args["linker2_length"])
spatial_barcode1_length = int(args["spatial_barcode1_length"])
spatial_barcode2_length = int(args["spatial_barcode2_length"])
umi_barcode_length = int(args["umi_barcode_length"])

seq_start=primer_length+spatial_barcode2_length+linker2_length+spatial_barcode1_length+linker1_length+umi_barcode_length

bc2_start=primer_length+1
bc2_end=bc2_start+spatial_barcode2_length-1

bc1_start=bc2_end+linker2_length+1
bc1_end=bc1_start+spatial_barcode1_length-1

bc_umi_start=bc1_end+linker1_length+1
bc_umi_end=bc_umi_start+umi_barcode_length-1

with open_maybe_gzip(input_file, "rt") as in_handle, open_maybe_gzip(output_file, "wt") as out_handle:
    for title, seq, qual in FastqGeneralIterator(in_handle):
        if len(seq) < seq_start or len(qual) < seq_start:
            sys.exit("*** Error: read '{}' is too short for configured barcode structure.\n"
                     "Need at least {} bases but found sequence length {} and quality length {}.\n".format(
                         title, seq_start, len(seq), len(qual)
                     ))
    
        pxl_bcd_seq = seq[bc2_start-1:bc2_end] + seq[bc1_start-1:bc1_end] # !!! BC2 + BC1
        pxl_bcd_qual = qual[bc2_start-1:bc2_end] + qual[bc1_start-1:bc1_end]

        umi_bcd_seq = seq[bc_umi_start-1:bc_umi_end]
        umi_bcd_qual = qual[bc_umi_start-1:bc_umi_end]

        pix_full_seq = pxl_bcd_seq+umi_bcd_seq
        pix_full_qual = pxl_bcd_qual+umi_bcd_qual
             
        out_handle.write("@%s\n%s\n+\n%s\n" % (title, pix_full_seq, pix_full_qual))

from Bio.SeqIO.QualityIO import FastqGeneralIterator
import gzip

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
ap.add_argument("-o1", "--output_R1", required=True, help="output file R1")
ap.add_argument("-o2", "--output_R2", required=True, help="output file R2")
ap.add_argument("-p", "--primer_length", required=True, help="PCR primer sequence length", type=check_positive)
ap.add_argument("-l1", "--linker1_length", required=True, help="linker1 sequence length", type=check_positive)
ap.add_argument("-l2", "--linker2_length", required=True, help="linker2 sequence length", type=check_positive)
ap.add_argument("-a", "--adapt_tn5_length", required=True, help="Tn5 adapter sequence length", type=check_positive)
ap.add_argument("-s1", "--spatial_barcode1_length", required=True, help="spatial barcode1 sequence length", type=check_positive)
ap.add_argument("-s2", "--spatial_barcode2_length", required=True, help="spatial barcode2 sequence length", type=check_positive)
ap.add_argument("-n", "--nanoCT_barcode_length", required=True, help="nanoCutTag barcode sequence length", type=check_positive)
args = vars(ap.parse_args())

input_file = args["input"]
output_file_R1 = args["output_R1"]
output_file_R2 = args["output_R2"]
primer_length = int(args["primer_length"])
linker1_length = int(args["linker1_length"])
linker2_length = int(args["linker2_length"])
adapt_tn5_length = int(args["adapt_tn5_length"])
spatial_barcode1_length = int(args["spatial_barcode1_length"])
spatial_barcode2_length = int(args["spatial_barcode2_length"])
nanoCT_barcode_length = int(args["nanoCT_barcode_length"])

seq_start=primer_length+spatial_barcode2_length+linker2_length+spatial_barcode1_length+linker1_length+nanoCT_barcode_length+adapt_tn5_length

bc2_start=primer_length+1
bc2_end=bc2_start+spatial_barcode2_length-1

bc1_start=bc2_end+linker2_length+1
bc1_end=bc1_start+spatial_barcode1_length-1

bc_nano_start=bc1_end+linker1_length+1
bc_nano_end=bc_nano_start+nanoCT_barcode_length-1

adapt_tn5_start=bc_nano_end+1
adapt_tn5_end=adapt_tn5_start+adapt_tn5_length-1

with open_maybe_gzip(input_file, "rt") as in_handle, open_maybe_gzip(output_file_R1, "wt") as out_handle_R1, open_maybe_gzip(output_file_R2, "wt") as out_handle_R2:
    for title, seq, qual in FastqGeneralIterator(in_handle):
    
        new_seq = seq[seq_start:]
        new_qual = qual[seq_start:]
        
        pxl_bcd_seq = seq[bc2_start-1:bc2_end] + seq[bc1_start-1:bc1_end] # !!! BC2 + BC1
        pxl_bcd_qual = qual[bc2_start-1:bc2_end] + qual[bc1_start-1:bc1_end]

        nano_bcd_seq = seq[bc_nano_start-1:bc_nano_end]
        nano_bcd_qual = qual[bc_nano_start-1:bc_nano_end]

        adapt_tn5_seq = seq[adapt_tn5_start-1:adapt_tn5_end]
        adapt_tn5_qual = qual[adapt_tn5_start-1:adapt_tn5_end]

        nano_full_seq = nano_bcd_seq+adapt_tn5_seq+pxl_bcd_seq
        nano_full_qual = nano_bcd_qual+adapt_tn5_qual+pxl_bcd_qual
             
        out_handle_R1.write("@%s\n%s\n+\n%s\n" % (title, new_seq, new_qual))
        out_handle_R2.write("@%s\n%s\n+\n%s\n" % (title, nano_full_seq, nano_full_qual))

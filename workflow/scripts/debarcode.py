import os
import argparse
from glob import glob
import sys
from re import split
from re import compile
import regex
import gzip
from contextlib import ExitStack
from collections import defaultdict
import time

import yaml
from pysam import FastxFile
import Levenshtein

def log(message):
    sys.stderr.write(time.strftime("%Y-%m-%d %H:%M:%S") + " " + message + "\n")

def open_maybe_gzip(path, mode):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)

class bcdCT:
    def __init__(self,args):
        self.detect_input(args.input)
        self.detect_reads()
        self.single_cell=args.single_cell
        self.out_prefix=args.out_prefix
        if self.single_cell:
            self.out_reads = ['R1','R2','R3']
        else:
            self.out_reads = ['R1','R3']

        if args.name:
            self.name = args.name
        else:
            self.autodetect_name()

        self.autodetect_lane()
        self.autodetect_barcodes(args)
        self.prep_out_filenames()

    def detect_input(self,input):
        Error_message="*** Error: Wrong input files specified. The input must be either folder with _R1_*.fastq[.gz] _R2_*.fastq[.gz] _R3_*.fastq[.gz] files or paths to the files themselves ***" +\
        "The files should be placed in the same folder" +\
        "e.g. /data/path_to_my_files/*L001*.fastq.gz or /data/path_to_my_files/"

        input = [os.path.abspath(x) for x in input]
        if len(input) == 1 and os.path.isdir(input[0]):     # Case input is single directory
            self.input_dir = input[0]
            self.input_files = []
            self.input_files.extend(glob(self.input_dir + "/*.fastq.gz"))
            self.input_files.extend(glob(self.input_dir + "/*.fq.gz"))
            self.input_files.extend(glob(self.input_dir + "/*.fastq"))
            self.input_files.extend(glob(self.input_dir + "/*.fq"))

        elif len(input) > 1:                                  # Case input are multiple files
            self.input_files = input
            self.input_dir = list(set([os.path.dirname(x) for x in self.input_files]))
            if not len(self.input_dir) == 1:
                log(Error_message)
                sys.exit(1)
            if not sum([x.endswith('.fastq.gz') or x.endswith('.fq.gz') or x.endswith('.fastq') or x.endswith('.fq') for x in self.input_files]) == len(self.input_files):
                sys.exit(1)
                log(Error_message)
        else:
            sys.exit(1)
            log(Error_message)

    def detect_reads(self):
        Error_message="*** Error: Please specify exactly one _R1_ _R2_ and _R3_ file or folder with exactly one of each files ***" + \
                      "e.g. /data/path_to_my_files/*L001*.fastq.gz or /data/path_to_my_files/"
        self.path_in = {}
        self.path_in['R1'] = [x for x in self.input_files if "_R1_" in x]
        self.path_in['R2'] = [x for x in self.input_files if "_R2_" in x]
        self.path_in['R3'] = [x for x in self.input_files if "_R3_" in x]

        if len(self.path_in['R1']) != 1 or len(self.path_in['R2']) != 1 or len(self.path_in['R3']) != 1:
            log(Error_message)
            sys.exit(1)

        self.path_in = {key:self.path_in[key][0] for key in self.path_in.keys()}



    def in_handles(self,stack):
        in_stack = {x: stack.enter_context(FastxFile(self.path_in[x],'r')) for x in ['R1','R2','R3']}
        return in_stack

    def prep_out_filenames(self):
        self.path_out = {barcode: {} for barcode in self.picked_barcodes}
        self.path_out  = {barcode: {read: "{0}/barcode_{1}/{2}".format(self.out_prefix,barcode,os.path.basename(self.path_in[read])) for read in self.out_reads} for barcode in self.picked_barcodes}
        # If args.name is specified, replace the sample_id prefix with the one specified in args.name
        # e.g. nanoCT_MB22_001_S1_L001_R1_001.fastq.gz is input --name is test
        # Change to test_S1_L001_R1_001.fastq.gz
        if args.name:
            for barcode in self.path_out:
                for read in self.path_out[barcode]:
                    sample_id = split('_S[0-9]+_', os.path.basename(self.path_out[barcode][read]))[0].strip("_")
                    self.path_out[barcode][read] = self.path_out[barcode][read].replace(sample_id,args.name)

    def autodetect_name(self):
        Error_message = "*** Error: Prefix for R1 R2 R3 files not the same. Please use the same prefix for all the files or specify experiment name ***"

        self.name = [split("_R[0-9]_", str(x)) for x in self.path_in.values()]
        self.name = [x[0] for x in self.name]

        if len(list(set(self.name))) > 1:
            log(Error_message)
            sys.exit(1)

        self.name = self.name[0].split("/")[-1]

    def autodetect_lane(self):
        Error_message = "*** Error: Prefix for R1 R2 R3 files not the same. Please use the same prefix for all the files or specify experiment name ***"

        l = compile('S[0-9]+_(.*)_R[0-9]')
        self.lane = [l.findall(str(x)) for x in self.path_in.values()][0]

        if len(list(set(self.lane))) > 1:
            log(Error_message)
            sys.exit(1)

        self.lane = self.lane[0]

    def create_out_handles(self,stack):
        for bcd in self.picked_barcodes:
            os.makedirs(self.out_prefix + "/barcode_" + bcd, exist_ok=True)
        self.out_stack = {barcode: {read: stack.enter_context(open_maybe_gzip(self.path_out[barcode][read],'wt'))for read in self.out_reads} for barcode in self.picked_barcodes}


    def __iter__(self):
        with FastxFile(self.path_in['R1']) as f1, FastxFile(self.path_in['R2']) as f2, FastxFile(self.path_in['R3']) as f3:
            for r1,r2,r3 in zip(f1,f2,f3):
                yield r1, r2, r3

    def autodetect_barcodes(self,args):
        barcodes = defaultdict(int)
        n=0

        hit_patterns = compile_patterns(args.pattern, 0)
        MeA_patterns = compile_patterns(args.no_barcode_seq, 2)

        for read1,read2,read3 in self:
            
            hit = find_seq_compiled(hit_patterns, read2.sequence)
            MeA_hit = find_seq_compiled(MeA_patterns, read2.sequence)
            
            n += 1
            if n == 100000:
                break
            if MeA_hit is not None:
                barcodes['MeA'] += 1
            elif hit is None or hit == 'Multiple':
                barcodes['no_spacer'] += 1
                continue
            else:
                hit = int(hit)
                read_barcode = get_read_barcode(read2, hit)
                try:
                    barcodes[read_barcode] += 1
                except KeyError:
                    barcodes[read_barcode] = 1

        top_barcodes = sorted(barcodes, key=barcodes.get, reverse=True)[:args.Nbarcodes]
        picked_barcodes = {key: barcodes[key] for key in top_barcodes}
        log("Detected following most abundant barcodes out of first {} barcodes:\n{}".format(n, picked_barcodes))
        if args.barcode != "None":
            self.picked_barcodes = args.barcode
            log("Barcode specified for demultiplexing [{barcode}] in top found barcodes: {bool} ".format(bool = [(x,x in picked_barcodes.keys()) for x in args.barcode], barcode = args.barcode))
        else:
            self.picked_barcodes = [i for i in picked_barcodes.keys()]

        if args.report_MeA:
            self.picked_barcodes.append('MeA')
            log("MeA sequence will be reported in the output files due to --report_MeA flag")
        if args.report_no_hit:
            self.picked_barcodes.append('no_spacer')
            log("No spacer sequence will be reported in the output files due to --report_no_hit flag")
        
        print('final barcodes used for demultiplexing:')
        print(self.picked_barcodes)
        

def get_read_barcode(string,index):
    read_barcode = string.sequence[index - len(args.barcode[0]):index].upper()  # Get the barcode sequence
    return read_barcode

def extract_cell_barcode(read,index):
    read.sequence = read.sequence[index:index + 16]  # Get the cell barcode
    read.quality  = read.quality[index:index + 16]   # Get corresponding Quality score
    return read

RC_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")

def revcompl(seq):
    return seq.translate(RC_TABLE)[::-1].upper()

def rev(seq):
    return seq[::-1]

def compile_patterns(pattern, max_mismatch):
    return [
        regex.compile(f"({pattern}){{e<={n}}}")
        for n in range(max_mismatch + 1)
    ]


def find_seq_compiled(compiled_patterns, DNA_string):
    for r in compiled_patterns:
        hits = [x.start() for x in r.finditer(DNA_string)]
        if len(hits) == 0:
            continue
        if len(hits) > 1:
            return None
        return hits[0]
    return None

def flush_buffers(buffers, out_stack):
    for bcd in buffers:
        for read in buffers[bcd]:
            if buffers[bcd][read]:
                out_stack[bcd][read].write("".join(buffers[bcd][read]))
                buffers[bcd][read].clear()

def count_selected_barcode_reads(statistics, picked_barcodes):
    return sum(statistics[barcode] for barcode in picked_barcodes)

def main(args):
    exp = bcdCT(args)
    statistics = defaultdict(int)
    log("Creating file output handles ")
    
    spacer_patterns = compile_patterns(args.pattern, 2)
    MeA_patterns = compile_patterns(args.no_barcode_seq, 2)
    
    with ExitStack() as stack:
        exp.create_out_handles(stack)

        buffers = {
            barcode: {read: [] for read in exp.out_reads}
            for barcode in exp.picked_barcodes
        }

        # Buffered writes are much faster than writing each read immediately.
        # The first flush also acts as an early sanity check for wrong barcodes.
        flush_every = args.flush_every

        n = 0
        sys.stderr.write("Starting demultiplexing \n")

        for read1,read2,read3 in exp:
            n+=1
            if n % 5000000 == 0:
                log("{} reads processed".format(n))
            assert (read1.name == read2.name == read3.name)                                                 # Make sure the fastq files are ok
           
            spacer_hit = find_seq_compiled(spacer_patterns, read2.sequence)
            MeA_hit    = find_seq_compiled(MeA_patterns, read2.sequence)
            
            if spacer_hit is None and MeA_hit is not None:
                read_barcode = get_read_barcode(read2, MeA_hit)                                               # Returns only barcode e.g. ACTGACTG
                hit_barcode  = 'MeA'
                if exp.single_cell:
                    read2 = extract_cell_barcode(read2, MeA_hit-16)     # The cell barcode is 16bp long and is positioned before the MeA spacer

            elif spacer_hit is not None:
                read_barcode = get_read_barcode(read2, spacer_hit)                                               # Returns only barcode e.g. ACTGACTG
                matches = []
                for barcode in exp.picked_barcodes:
                    d = Levenshtein.distance(read_barcode, barcode)
                    if d <= args.mismatch:
                        matches.append(barcode)
                        if len(matches) > 1:
                            break
                if len(matches) == 0:
                    # Spacer hit but no barcode match
                    statistics["no_barcode_match"] += 1
                    continue
                if len(matches) > 1:
                    # Spacer hit but multiple barcode matches
                    statistics["multiple_barcode_matches"] += 1
                    continue

                hit_barcode = matches[0]

                if exp.single_cell:
                    read2 = extract_cell_barcode(read2, spacer_hit + len(args.pattern))     # The cell barcode is 16bp long and is positioned after the spacer
            
            elif spacer_hit is None and MeA_hit is None:
                statistics["no_spacer_found"] += 1
                # No hit, no spacer not nothing found
                if args.report_no_hit:
                    hit_barcode = 'no_spacer'
                else:
                    continue        
            
            # Now continue in the loop
            if len(read2.sequence) < 16:
                    statistics["too_short_read"] += 1
                    continue
            
            statistics[hit_barcode] += 1
            if hit_barcode in exp.picked_barcodes:
                # Write the outputs
                buffers[hit_barcode]["R1"].append(str(read1) + "\n")
                buffers[hit_barcode]["R3"].append(str(read3) + "\n")
                if args.single_cell:
                    buffers[hit_barcode]["R2"].append(str(read2) + "\n")

            if n % flush_every == 0:
                flush_buffers(buffers, exp.out_stack)
                if n == flush_every:
                    selected_reads = count_selected_barcode_reads(statistics, exp.picked_barcodes)
                    selected_ratio = selected_reads / n
                    # Fail early if the selected antibody barcode is essentially absent.
                    if selected_ratio < args.min_first_flush_ratio:
                        sys.exit("*** Error: only {} / {} reads ({:.6f}) matched the selected barcode(s) in the first flush: {}.\n"
                                 "Minimum required by --min_first_flush_ratio is {:.6f}.\n".format(
                            selected_reads,
                            n,
                            selected_ratio,
                            ", ".join(exp.picked_barcodes),
                            args.min_first_flush_ratio
                        ))

        flush_buffers(buffers, exp.out_stack)

    log("Finished demultiplexing {} reads. Statistics:\n{}".format(n, dict(statistics)))

    # Write the statistics file
    with open("{0}/{1}_{2}_statistics.yaml".format(exp.out_prefix,exp.name,exp.lane), 'w') as f:
        yaml.dump(statistics, f)

    if count_selected_barcode_reads(statistics, exp.picked_barcodes) == 0:
        sys.exit("*** Error: no reads matched the selected barcode(s): {}\n".format(
            ", ".join(exp.picked_barcodes)
        ))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="DESCRIPTION: \n\nThis script demultiplexes Nano-CT sequencing data by extracting and matching modality barcodes from the R2 read based on a specified spacer sequence.\nThe script supports both bulk and single-cell data and writes sorted reads into separate output files for each detected barcode""" + 
                                     """The script works with a standard 36-8-48-36 read structure but may also be compatible with other setups where the spacer sequence is similarly arranged.\n""" + 
                                     """AATGATACGGCGACCACCGAGATCTACAC-NNNNNNNNNNNNNNNN-TCGTCGGCAGCGTCTCCACGC-NNNNNNNN-GCGATCGAGGACGGCAGATGTGTATAAGAGACAG"""+
                                     """            P5                |  sc-barcode   |  Linker sequence   | Modality |        Mosaic end               \n """ + 
                                     """""" +
                                     """ Note: If demultiplexing multiple lanes, run for each lane separately and then merge the output files before or after alignment""",
                                     usage=""
                                           "python debarcode.py -i /path/to/input_R1.fastq.gz /path/to/input_R2.fastq.gz /path/to/input_R3.fastq.gz -o /path/to/output_folder --single_cell --barcode ATAGAGGC                      # One specific barcode from single-cell data " 
                                           "python debarcode.py -i /path/to/input_R1.fastq.gz /path/to/input_R2.fastq.gz /path/to/input_R3.fastq.gz -o /path/to/output_folder --single_cell --barcode ATAGAGGC TATAGCCT             # Two specific barcodes from single-cell data "
                                           "python debarcode.py -i /path/to/input_R1.fastq.gz /path/to/input_R2.fastq.gz /path/to/input_R3.fastq.gz -o /path/to/output_folder --single_cell --Nbarcodes 3                           # Top 3 barcodes from single-cell data without specifying the barcodes - use carefully and double check"
                                           "python debarcode.py -i /path/to/input_R1.fastq.gz /path/to/input_R2.fastq.gz /path/to/input_R3.fastq.gz -o /path/to/output_folder --Nbarcodes 3                                         # Top 3 barcodes from bulk data ", 
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--input',
                        required=True,
                        type=str,
                        nargs='+',
                        help='path to input R1,R2,R3 .fastq.gz files [3 files required]')

    parser.add_argument('-o', '--out_prefix',
                        type=str,
                        required=True,
                        help='Prefix to where to put the output files; Directory will be created')

    parser.add_argument('-p', '--pattern',
                        type=str,
                        default="GCGTGGAGACGCTGCCGACGA",
                        help='Pattern that follows the antibody barcode \n \
                                  (Default: %(default)s)')

    parser.add_argument('--single_cell',
                        default=False,
                        action='store_true',
                        help='Data is single cell CUT&Tag (Default: %(default)s)')

    parser.add_argument('--name',
                        type=str,
                        default=None,
                        help='Custom name for the experiment (Default: Autodetect from filename)')

    parser.add_argument('--mismatch',
                        type=int,
                        default=1,
                        help='Maximum mismatches for sample barcode (Default: %(default)s)')

    parser.add_argument('--Nbarcodes',
                        type=int,
                        default=10,
                        help='Number of barcodes in experiment (Default: %(default)s)')

    parser.add_argument('--barcode',
                        type=str,
                        nargs="+",
                        default='None',
                        help='Specific barcode to be extracted [e.g. ATAGAGGC] (Default: All barcodes [see --Nbarcodes])')

    parser.add_argument('--flush_every',
                        type=int,
                        default=100000,
                        help='Number of reads to buffer before writing outputs and checking the first-flush barcode ratio (Default: %(default)s)')

    parser.add_argument('--min_first_flush_ratio',
                        type=float,
                        default=0.00001,
                        help='Minimum selected-barcode read ratio required after the first flush (Default: %(default)s)')

    parser.add_argument('--no_barcode_seq', type=str, 
                        default='GTGTAGATCTCGGTGGTCGCCGTATCATT', 
                        help='Sequence indicating unbarcoded reads')
    
    parser.add_argument('--report_MeA', 
                        action='store_true', 
                        help='Include reads with no barcode and standard MeA sequence in the output')
    
    parser.add_argument('--report_no_hit', 
                        action='store_true', 
                        help='Include reads with no barcode or spacer whatsoever hit in the output')




    args = parser.parse_args()
    if args.flush_every <= 0:
        parser.error("--flush_every must be a positive integer")
    if args.min_first_flush_ratio < 0 or args.min_first_flush_ratio > 1:
        parser.error("--min_first_flush_ratio must be between 0 and 1")
    log("Starting debarcode.py script ")
    log("Input files: \n{}".format("".join(["    " + i + "\n" for i in args.input])))
    if args.barcode != "None":
        log("Provided barcodes to demultiplex: \n{}".format(args.barcode))
    log("Output prefix: {}/".format(args.out_prefix))
    
    main(args)

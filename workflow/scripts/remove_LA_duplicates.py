import pysam
import sys
import json
import argparse

# Marek Bartosovic 13/06/2023
# Script which removes linear amplification duplicates from possorted_bam.bam file from cellranger
# Script keeps and reports PCR duplicates as these are handled fine in downstream analysis

# Warning: Input possorted_bam.bam file needs to be re-sorted by read name otherwise the script won't work
# This can be done using following command:
# samtools sort -@ 10 -n -o namesorted_bam.bam possorted_bam.bam


parser = argparse.ArgumentParser(
    description="Remove linear amplification duplicates from a name-sorted Cell Ranger BAM."
)
parser.add_argument("input_bam", help="Input BAM sorted by read name")
parser.add_argument("output_bam", help="Output BAM without linear amplification duplicates")
args = parser.parse_args()

samfile = pysam.AlignmentFile(args.input_bam, "rb")
sam_out = pysam.AlignmentFile(args.output_bam, "wb", header=samfile.header)
stats_out = args.output_bam + "_stats.txt"

reads_unique = {}
stats = {
    'not mapped': 0,
    'unique': 0,
    'not proper': 0,
    'LA duplicates': 0,
    'PCR duplicates': 0,
}

n = 0
# n_target = 1000000  # test for 1 million reads
# n_target = 10000000 # test for 10 million reads

read1 = False
read2 = False

while True:
    try:
        read1 = next(samfile)
    except StopIteration:
        break

    n += 1
    if n % 100000 == 0:
        sys.stderr.write('*** {} lines processed\n'.format(n))

    # # Only take the first n_target reads # Debugging
    # if n == n_target:
    #     # print(reads_unique)
    #     print(stats)
    #     break

    try:
        read2 = next(samfile)
    except StopIteration:
        sys.exit("*** Error: input BAM ended with an unpaired read: {}\n"
                 "*** Perhaps the BAM is not name-sorted or contains incomplete pairs?\n".format(read1.query_name))

    if read1.query_name != read2.query_name:
        sys.exit("*** Error: paired reads are not adjacent in input BAM.\n"
                 "Read 1: {}\nRead 2: {}\n"
                 "*** Perhaps the BAM is not name-sorted?\n".format(read1.query_name, read2.query_name))
    if not read1.is_read1 or not read2.is_read2:
        sys.exit("*** Error: unexpected read-pair flags for read '{}'.\n"
                 "Expected first record to be R1 and second record to be R2.\n".format(read1.query_name))

    if read1.is_unmapped or read2.is_unmapped:
        stats['not mapped'] += 1
        continue

    if not read1.is_proper_pair or not read2.is_proper_pair:
        stats['not proper'] += 1
        continue

    try:
        cell_barcode = read1.get_tag('CR')
    except KeyError:
        sys.exit("*** Error: read '{}' is missing the CR cell-barcode tag.\n"
                 "Expected Cell Ranger ATAC-style BAM tags.\n".format(read1.query_name))
    read1_position = (read1.reference_name, read1.reference_start, cell_barcode)
    read2_position = (read1.next_reference_name, read1.next_reference_start, cell_barcode)

    if read1_position not in reads_unique:
        reads_unique[read1_position] = {read2_position}
        stats['unique'] += 1
        sam_out.write(read1)
        sam_out.write(read2)
        continue

    if read2_position not in reads_unique[read1_position]:
        reads_unique[read1_position].add(read2_position)
        stats['LA duplicates'] += 1
        # LA duplicates are not reported in the final bam file
        continue

    # Case read2_position list is not empty and read2_position is in the list        == PCR duplicate
    if read2_position in reads_unique[read1_position]:
        stats['PCR duplicates'] += 1
        # PCR duplicates can be reported as they are handled in downstream analysis
        sam_out.write(read1)
        sam_out.write(read2)
        continue

    # Sanity check if something escapes
    sys.exit('found exception')

with open(stats_out,'w') as f:
    json.dump(stats,f)

samfile.close()
sam_out.close()

def get_fastq_info_from_folder(fastq_folder,sample):
    import glob
    all_fastq_files  = glob.glob(fastq_folder + "/" + sample + "/**/*.fastq.gz",recursive=True)
    all_fastq_parsed = [parse_fastq(x) for x in all_fastq_files]
    return(all_fastq_parsed)

def parse_fastq(path):
    import os
    import re
    pattern = re.compile(
        r'^(?P<id>.+)_(?P<number>S[0-9]+)_(?P<lane>L[0-9]+)_(?P<read>[RI][0-9]+)_(?P<suffix>[0-9]+)\.(?P<fastq>fastq|fq)\.(?P<compress>gz)$'
    )
    fastq = os.path.basename(path)
    match = pattern.match(fastq)
    if not match:
        raise ValueError(
            "Wrong input file name '{}'. Expected 10X-style '*_S1_L001_R1_001.fastq.gz' or '.fq.gz'.".format(fastq)
        )
    return match.groupdict()

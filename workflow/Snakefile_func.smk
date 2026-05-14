def get_fastq_info_from_folder(fastq_folder,sample):
    import glob
    all_fastq_files = []
    all_fastq_files.extend(glob.glob(fastq_folder + "/" + sample + "/**/*.fastq.gz", recursive=True))
    all_fastq_files.extend(glob.glob(fastq_folder + "/" + sample + "/**/*.fq.gz", recursive=True))
    all_fastq_files.extend(glob.glob(fastq_folder + "/" + sample + "/**/*.fastq", recursive=True))
    all_fastq_files.extend(glob.glob(fastq_folder + "/" + sample + "/**/*.fq", recursive=True))
    all_fastq_files = sorted(all_fastq_files)
    all_fastq_parsed = [parse_fastq(x) for x in all_fastq_files]
    return(all_fastq_parsed)

def parse_fastq(path):
    import os
    import re
    pattern = re.compile(
        r'^(?P<id>.+)_(?P<number>S[0-9]+)_(?P<lane>L[0-9]+)_(?P<read>[RI][0-9]+)_(?P<suffix>[0-9]+)\.(?P<ext>(?:fastq|fq)(?:\.gz)?)$'
    )
    fastq = os.path.basename(path)
    match = pattern.match(fastq)
    if not match:
        raise ValueError(
            "Wrong input file name '{}'. Expected 10X-style '*_S1_L001_R1_001.fastq.gz', '.fq.gz', '.fastq', or '.fq'.".format(fastq)
        )
    return match.groupdict()

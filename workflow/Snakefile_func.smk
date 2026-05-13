def get_fastq_info_from_folder(fastq_folder,sample):
    import glob
    all_fastq_files  = glob.glob(fastq_folder + "/" + sample + "/**/*.fastq.gz",recursive=True)
    all_fastq_parsed = [parse_fastq(x) for x in all_fastq_files]
    return(all_fastq_parsed)

def parse_fastq(path):
    import os
    import re
    result = {}
    fastq = os.path.basename(path)
    result['number'] = re.findall('_S[0-9]+_', fastq)[0].strip("_")
    result['lane']   = re.findall('_L[0-9]+_', fastq)[0].strip("_")
    result['read']   = re.findall('_[RI][0-9]+_', fastq)[0].strip("_")
    result['id']     = re.split('_S[0-9]+_',fastq)[0].strip("_")
    result['suffix']  = re.split("\.", re.split('_[RI][0-9]+_',fastq)[1].strip("_"))[0]
    result['fastq']  = re.split('\.',fastq)[1]
    result['compress']  = re.split('\.',fastq)[2]
    return(result)
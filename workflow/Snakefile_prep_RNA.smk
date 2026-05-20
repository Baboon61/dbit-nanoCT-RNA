import os
import sys

include: 'Snakefile_func.smk'

configfile: workflow.basedir + '/../config/config_RNA.yaml'

repo_root = os.path.abspath(os.path.join(workflow.basedir, '..'))

include: 'Snakefile_validation.smk'

validate_general_config_schema('RNA')

#Put RNA samples into variables
samples_list  = config['samples']
print("RNA samples to be processed : {}".format(samples_list))

barcode_r1, barcode_r2, Rfile_barcode, Rfile_sequence = get_read_layout()

#Clean path for raw and processed data to have / at the end of the paths
raw_dir = with_trailing_slash(config['general']['rawData_dir'])
proc_dir = with_trailing_slash(config['general']['processedData_dir'])
config['general']['spatial_barcodes_file'] = resolve_repo_path(config['general']['spatial_barcodes_file'])

validate_common_config("config/config_RNA.yaml")

#Check Cellranger format of fastq files
fastq_info_by_sample = {sample_id: validate_sample_fastqs(sample_id, ['R1', 'R2'], raw_dir, "config/config_RNA.yaml") for sample_id in samples_list}
res = [res for sample_id in samples_list for res in fastq_info_by_sample[sample_id]]
validate_fastq_records(res)

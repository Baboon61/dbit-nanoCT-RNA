import itertools
import os
import glob
import datetime

include: 'Snakefile_func.smk'

configfile: workflow.basedir + '/../config/config_CT.yaml'

#Put CutTag samples into variables
samples_list  = list(config['samples'].keys())
print("CutTag samples to be processed : {}".format(samples_list))

#Put barcodes into variables
barcodes_dict = {sample: config['samples'][sample]['barcodes'] for sample in samples_list}
print("List of barcodes to be processed : {}".format(barcodes_dict))

#Put antibodies into variables
antibodies_list = list(set(itertools.chain(*[barcodes_dict[sample].keys() for sample in barcodes_dict.keys()])))
print("List of antibodies to be processed : {}".format(antibodies_list))

#Register if barcodes construction is R1 or R2
if config['general']['read_barcode']=="R1":
  barcode_r1="f"
  barcode_r2="t"
  Rfile_barcode="R1"
  Rfile_sequence="R2"
elif config['general']['read_barcode']=="R2":
  barcode_r1="t"
  barcode_r2="f"
  Rfile_barcode="R2"
  Rfile_sequence="R1"
else:
  print("Error in barcode in R1/R2 in config file, please double check\n")

#Clean path for raw and processed data to have / at the end of the paths
raw_dir = config['general']['rawData_dir'] if config['general']['rawData_dir'][-1] == '/' else config['general']['rawData_dir']+"/"
proc_dir = config['general']['processedData_dir'] if config['general']['processedData_dir'][-1] == '/' else config['general']['processedData_dir']+"/"

#Check Cellranger format of fastq files
res = [res for sample_id in samples_list for res in get_fastq_info_from_folder(raw_dir, sample_id)]

#Check that the cores argument are lower than the rules cores parameters
if workflow.cores < config['general']['cellranger']['core']:
  #sys.exit(1)
  sys.stderr.write("Not enough core is given as parameters to the snakemake command line to run cellranger rule (core in Cellrnager section in the config file)\n")

if workflow.cores < config['general']['core']:
  #sys.exit(1)
  sys.stderr.write("Not enough core is given as parameters to the snakemake command line to run the general rules (core in the config file)\n")

#Remove the back up file from the cellranger barcodes to allow the copy of the new ones in rule get_barcodes_cellranger
if os.path.exists(config['general']['cellranger']['barcodes_path'] + 'BAK_' + config['general']['cellranger']['barcodes_file']):
  os.remove(config['general']['cellranger']['barcodes_path'] + 'BAK_' + config['general']['cellranger']['barcodes_file'])

#Print the 10X nomenclature error message
Error_message="*** Error: Wrong input files specified. The input must follow 10X fastq files nomenclature, see  https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/fastq-input***\n" +\
        "The files should be placed in the same folder specify as raw folder in the config.yaml file\n"

for sample in samples_list:
  file_name_list = glob.glob(raw_dir + "/" + sample + "/**/*.fastq.gz",recursive=True)
  for file_name in file_name_list:
    if not len(re.split("\.",os.path.basename(file_name))) == 3:
      sys.exit(1)
      sys.stderr.write(Error_message)

for file_struc in res:
  if not re.match("S[0-9]+",file_struc['number']):
    sys.exit(1)
    sys.stderr.write(Error_message)
  elif not re.match("L[0-9]+",file_struc['lane']):
    sys.exit(1)
    sys.stderr.write(Error_message)
  elif not re.match("[RI][0-9]+",file_struc['read']):
    sys.exit(1)
    sys.stderr.write(Error_message)
  elif not re.match("[0-9]+",file_struc['suffix']):
    sys.exit(1)
    sys.stderr.write(Error_message)
  elif not file_struc['fastq'] == 'fastq' or file_struc['fastq'] == 'fq':
    sys.exit(1)
    sys.stderr.write(Error_message)
  elif not file_struc['compress'] == 'gz':
    sys.exit(1)
    sys.stderr.write(Error_message)

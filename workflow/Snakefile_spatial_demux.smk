import glob

# Keep read pairs where the PCR primer is found within the allowed mismatch rate.
rule filter_primer:
  input:
    in1 = lambda wildcards: glob.glob(raw_fastq_glob(wildcards.sample, wildcards.number, wildcards.lane, 'R1', wildcards.suffix, wildcards.ext), recursive=True),
    in2 = lambda wildcards: glob.glob(raw_fastq_glob(wildcards.sample, wildcards.number, wildcards.lane, 'R2', wildcards.suffix, wildcards.ext), recursive=True)
  output:
    out1 = qc_fastq('{sample}', '{number}', '{lane}', 'R1', 'raw_qc_primer', '{suffix}', '{ext}'),
    out2 = qc_fastq('{sample}', '{number}', '{lane}', 'R2', 'raw_qc_primer', '{suffix}', '{ext}'),
    outfile = touch(qc_done('{sample}', '{number}', '{lane}', '{suffix}', '{ext}', 'filter_primer'))
  log:
    sample_lane_log('{sample}', '{number}', '{lane}', '{suffix}', '{ext}', 'filter_primer')
  params:
    k_length = len(config['general']['PCRprimer_sequence']),
    k_length_restrict = len(config['general']['PCRprimer_sequence'])+8,
    barcode_r1 = barcode_r1,
    barcode_r2 = barcode_r2,
    hamming_distance = config['general']['hamming_distance_linkers'],
    PCRprimer_sequence = config['general']['PCRprimer_sequence'],
    out_dir = proc_dir
  threads: config['general']['core']
  conda: '../envs/dbit-spatial-nanoct-rna-bbduk.yaml'
  shell:
    '''
    {{
    bbduk.sh \
    in1={input.in1} \
    in2={input.in2} \
    outm1={output.out1} \
    outm2={output.out2} \
    k={params.k_length} mm=f rcomp=f restrictleft={params.k_length_restrict} \
    skipr1={params.barcode_r1} skipr2={params.barcode_r2} \
    hdist={params.hamming_distance} \
    stats={params.out_dir}{wildcards.sample}/tmp_data/qc_raw_data/{wildcards.sample}_{wildcards.lane}_stats.primer.txt \
    threads={threads} \
    literal={params.PCRprimer_sequence}
    }} > {log} 2>&1
    '''

# Keep read pairs where linker 1 is found after primer filtering.
rule filter_L1:
  input:
    in1 = lambda wildcards: qc_fastq(wildcards.sample, wildcards.number, wildcards.lane, 'R1', 'raw_qc_primer', wildcards.suffix, wildcards.ext),
    in2 = lambda wildcards: qc_fastq(wildcards.sample, wildcards.number, wildcards.lane, 'R2', 'raw_qc_primer', wildcards.suffix, wildcards.ext),
    infile = lambda wildcards: qc_done(wildcards.sample, wildcards.number, wildcards.lane, wildcards.suffix, wildcards.ext, 'filter_primer')
  output:
    out1 = qc_fastq('{sample}', '{number}', '{lane}', 'R1', 'raw_qc_linker1', '{suffix}', '{ext}'),
    out2 = qc_fastq('{sample}', '{number}', '{lane}', 'R2', 'raw_qc_linker1', '{suffix}', '{ext}'),
    outfile = touch(qc_done('{sample}', '{number}', '{lane}', '{suffix}', '{ext}', 'filter_L1'))
  log:
    sample_lane_log('{sample}', '{number}', '{lane}', '{suffix}', '{ext}', 'filter_L1')
  params:
    k_length = len(config['general']['linker1_sequence']),
    k_length_restrict = len(config['general']['linker1_sequence'])+80,
    barcode_r1 = barcode_r1,
    barcode_r2 = barcode_r2,
    hamming_distance = config['general']['hamming_distance_linkers'],
    linker1_sequence = config['general']['linker1_sequence'],
    out_dir = proc_dir
  threads: config['general']['core']
  conda: '../envs/dbit-spatial-nanoct-rna-bbduk.yaml'
  shell:
    '''
    {{
    bbduk.sh \
    in1={input.in1} \
    in2={input.in2} \
    outm1={output.out1} \
    outm2={output.out2} \
    k={params.k_length} mm=f rcomp=f restrictleft={params.k_length_restrict} \
    skipr1={params.barcode_r1} skipr2={params.barcode_r2} \
    hdist={params.hamming_distance} \
    stats={params.out_dir}/{wildcards.sample}/tmp_data/qc_raw_data/{wildcards.sample}_{wildcards.lane}_stats.linker1.txt \
    threads={threads} \
    literal={params.linker1_sequence}
    }} > {log} 2>&1
    '''

# Keep read pairs where linker 2 is found after linker 1 filtering.
rule filter_L2:
  input:
    in1 = lambda wildcards: qc_fastq(wildcards.sample, wildcards.number, wildcards.lane, 'R1', 'raw_qc_linker1', wildcards.suffix, wildcards.ext),
    in2 = lambda wildcards: qc_fastq(wildcards.sample, wildcards.number, wildcards.lane, 'R2', 'raw_qc_linker1', wildcards.suffix, wildcards.ext),
    infile = lambda wildcards: qc_done(wildcards.sample, wildcards.number, wildcards.lane, wildcards.suffix, wildcards.ext, 'filter_L1')
  output:
    out1 = qc_fastq('{sample}', '{number}', '{lane}', 'R1', 'raw_qc_linker2', '{suffix}', '{ext}'),
    out2 = qc_fastq('{sample}', '{number}', '{lane}', 'R2', 'raw_qc_linker2', '{suffix}', '{ext}'),
    outfile = touch(qc_done('{sample}', '{number}', '{lane}', '{suffix}', '{ext}', 'filter_L2'))
  log:
    sample_lane_log('{sample}', '{number}', '{lane}', '{suffix}', '{ext}', 'filter_L2')
  params:
    k_length = len(config['general']['linker2_sequence']),
    k_length_restrict = len(config['general']['linker2_sequence'])+40,
    barcode_r1 = barcode_r1,
    barcode_r2 = barcode_r2,
    hamming_distance = config['general']['hamming_distance_linkers'],
    linker2_sequence = config['general']['linker2_sequence'],
    out_dir = proc_dir
  threads: config['general']['core']
  conda: '../envs/dbit-spatial-nanoct-rna-bbduk.yaml'
  shell:
    '''
    {{
    bbduk.sh \
    in1={input.in1} \
    in2={input.in2} \
    outm1={output.out1} \
    outm2={output.out2} \
    k={params.k_length} mm=f rcomp=f restrictleft={params.k_length_restrict} \
    skipr1={params.barcode_r1} skipr2={params.barcode_r2} \
    hdist={params.hamming_distance} \
    stats={params.out_dir}/{wildcards.sample}/tmp_data/qc_raw_data/{wildcards.sample}_{wildcards.lane}_stats.linker2.txt \
    threads={threads} \
    literal={params.linker2_sequence}
    }} > {log} 2>&1
    '''

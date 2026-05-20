def raw_fastq_glob(sample, number, lane, read, suffix, ext):
  return raw_dir + '{sample}/**/{sample}_{number}_{lane}_{read}_{suffix}.{ext}'.format(
    sample=sample,
    number=number,
    lane=lane,
    read=read,
    suffix=suffix,
    ext=ext
  )

def qc_fastq(sample, number, lane, read, step, suffix, ext):
  return proc_dir + '{sample}/tmp_data/qc_raw_data/{sample}_{number}_{lane}_{read}_{step}_{suffix}.{ext}'.format(
    sample=sample,
    number=number,
    lane=lane,
    read=read,
    step=step,
    suffix=suffix,
    ext=ext
  )

def qc_done(sample, number, lane, suffix, ext, step):
  return sample_lane_done(sample, number, lane, suffix, ext, step)

def tmp_fastq(sample, number, lane, read, suffix, ext):
  return proc_dir + '{sample}/tmp_data/{sample}_{number}_{lane}_{read}_{suffix}.{ext}'.format(
    sample=sample,
    number=number,
    lane=lane,
    read=read,
    suffix=suffix,
    ext=ext
  )

def tmp_done(sample, number, lane, suffix, ext, step):
  return sample_lane_done(sample, number, lane, suffix, ext, step)

def sample_done(sample, name):
  return proc_dir + '{sample}/.done/{name}.done'.format(
    sample=sample,
    name=name
  )

def sample_lane_done(sample, number, lane, suffix, ext, name):
  return proc_dir + '{sample}/.done/{name}_{number}_{lane}_{suffix}.{ext}.done'.format(
    sample=sample,
    number=number,
    lane=lane,
    suffix=suffix,
    ext=ext,
    name=name
  )

def modality_done(sample, modality, barcode, name):
  return proc_dir + '{sample}/{modality}_{barcode}/.done/{name}.done'.format(
    sample=sample,
    modality=modality,
    barcode=barcode,
    name=name
  )

def modality_lane_done(sample, modality, barcode, number, lane, suffix, ext, name):
  return proc_dir + '{sample}/{modality}_{barcode}/.done/{name}_{number}_{lane}_{suffix}.{ext}.done'.format(
    sample=sample,
    modality=modality,
    barcode=barcode,
    number=number,
    lane=lane,
    suffix=suffix,
    ext=ext,
    name=name
  )

def global_done(name):
  return proc_dir + '.done/{name}.done'.format(
    name=name
  )

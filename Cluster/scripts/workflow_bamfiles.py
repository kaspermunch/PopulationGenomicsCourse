import os, sys, re
from gwf import Workflow, AnonymousTarget
from gwf.workflow import collect
from pathlib import Path

gwf = Workflow(defaults={'account': 'populationgenomics'})

def modpath(p, parent=None, base=None, suffix=None):
    par, name = os.path.split(p)
    name_no_suffix, suf = os.path.splitext(name)
    if type(suffix) is str:
        suf = suffix
    if parent is not None:
        par = parent
    if base is not None:
        name_no_suffix = base
    new_path = os.path.join(par, name_no_suffix + suf)
    if type(suffix) is tuple:
        assert len(suffix) == 2
        new_path, nsubs = re.subn(r'{}$'.format(suffix[0]), suffix[1], new_path)
        assert nsubs == 1, nsubs
    return new_path

def bam_index(path):
    inputs = {'path': path}
    outputs = {'path': path + '.bai'}
    options = {'memory': '4g',
               'walltime': '0-02:00:00'}
    spec = f'samtools index {path}'
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def bam_subset(path, index):
    sample_name = os.path.basename(os.path.dirname(path))
    output_path = modpath(path, base=sample_name, suffix='.region.bam')
    # output_path = modpath(path, base=sample_name, suffix='.chr2.bam')
    inputs = {'path': path, 'index': index}
    outputs = {'path': output_path}
    options = {'memory': '4g',
               'walltime': '0-02:00:00'}
    spec = f'samtools view -b -h {path} -o {output_path} 2:135000000-145000000'
    # spec = f'samtools view -f 3 -F 4 -b -h {path} -o {output_path} 2'
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def bam2fastq(path):
    output_path = modpath(path, suffix='.fq')
    inputs = {'path': path}
    outputs = {'path': output_path}
    options = {'memory': '4g',
               'walltime': '0-02:00:00'}
    spec = f'bamToFastq -i {path} -fq {output_path}'
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# make sure downloaded bam files are sorted:
#samtools view -H test.bam | grep @HD # shows SO:coordinate if sorted

# bam files
bam_data_dir = Path('/home/kmt/populationgenomics/people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq')
bam_files = list(map(str, bam_data_dir.glob('**/*srt.aln.bam')))

# index bam files
bam_index_targets = gwf.map(bam_index, bam_files)

# subset bam files
bam_index_files = [d['path'] for d in bam_index_targets.outputs]
bam_subset_files = []
for i, (bam_file, bam_index_file) in enumerate(zip(bam_files, bam_index_files)):
    target = gwf.target_from_template(f'bam_subset_{i}', bam_subset(path=bam_file, index=bam_index_file))
    bam_subset_files.append(target.outputs['path'])
#bam_subset_targets = gwf.map(bam_subset, bam_files)

# make fastq files
bam2fastq_targets = gwf.map(bam2fastq, bam_subset_files)







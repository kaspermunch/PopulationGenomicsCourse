import os, sys, re
from gwf import Workflow, AnonymousTarget
from gwf.workflow import collect
from pathlib import Path

gwf = Workflow(defaults={'account': 'primatediversity'})

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

def bam_subset(path):
    sample_name = os.path.basename(os.path.dirname(path))
    output_path = modpath(path, base=sample_name, suffix='.region.bed')
    inputs = {'path': path}
    outputs = {'path': output_path}
    options = {'memory': '4g',
               'walltime': '0-02:00:00'}
    spec = f'samtools view -h {path} -o {output_path} 2:135000000-145000000'
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

bam_data_dir = Path('/home/kmt/populationgenomics/data/bamfiles')

bam_files = list(map(str, bam_data_dir.glob('**/*.bam')))

bam_index_tasks = gwf.map(bam_index, bam_files)

bam_subset_tasks = gwf.map(bam_subset, bam_files)

bam2fastq_tasks = gwf.map(bam2fastq, bam_files)

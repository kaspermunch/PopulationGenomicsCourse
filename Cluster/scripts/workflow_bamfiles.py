import os, sys
from gwf import Workflow
from gwf.workflow import collect
from templates_vcfmerge import *

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
    inputs = {'path': path}
    outputs = {'path': path + '.region'}
    options = {'memory': '4g',
               'walltime': '0-01:00:00'}
    spec = f'samtools view -h {path} -o {path}.region 2:135000000-145000000'
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def bam2fastq(path):
    inputs = {'path': path}
    outputs = {'path': path.replace('.bam', '.fq')}
    options = {'memory': '4g',
               'walltime': '0-02:00:00'}
    spec = f'bamToFastq -i filename.sorted.bam -fq filename.sorted.fq {path}'
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# Execute /path_to_bedtools/bin/bamToFastq -i filename.sorted.bam -fq filename.sorted.fq

samtools view -H test.bam | grep @HD # shows SO:coordinate if sorted

bam_files = []

bam_index_tasks = gwf.map(bam_index, bam_files)

bam_subset_tasks = gwf.map(bam_subset, bam_files)

bam2fastq_tasks = gwf.map(bam2fastq, bam_files)



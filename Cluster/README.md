
# Managing access, data, and software for students on the cluster



Data folder `data` is placed under here and them symlinked (for student access) from `populationgenomics/data`




Instructor conda environment for building data files with `workflow_bedfiles.py`:

    conda create -n popgen -c gwforg -c bioconda python=3 gwf bedtools samtools
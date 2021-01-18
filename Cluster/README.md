# Instructor zone

## Student access to the cluster

The students have temporary access to the cluster 

## Software

All command line software, except PSMC and LDhat, is installed in the `popgen` environment.  and LDhat PSMC is available in the software folder

## Running on the cluster

The students begin each exercise by running this to get an interactive session for the duration of the exercise:


    srun --mem-per-cpu=1g --time=3:00:00 --account=populationgenomics --pty bash

## Jupyter




# Managing access, data, and software for students on the cluster



Data folder `data` is placed under here and them symlinked (for student access) from `populationgenomics/data`




Instructor conda environment for building data files with `workflow_bedfiles.py`:

    conda create -n popgen -c gwforg -c bioconda python=3 gwf bedtools samtools
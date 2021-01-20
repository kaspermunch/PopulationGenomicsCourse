# Instructor zone

## Student access to the cluster

The students have temporary access to the cluster 


## the populationgenomics project folder

`instructor_data` links to a folder under the git repository. This folder holds all the intermediary files needed for building the files that the students need

`data` is the folder where the student finds the data files he/she needs for the exercises. This folder will contain symlinks to the few required files in the instructor_data folder.

`software`: builds and binaries for LDhat and PSMC

`people`: teachers and instructors folders

`students`: student folders

## Software

All command line software, except PSMC and LDhat, is installed in the `popgen` environment.  and LDhat PSMC is available in the software folder

## Running on the cluster

The students begin each exercise by running this to get an interactive session for the duration of the exercise:


    srun --mem-per-cpu=1g --time=3:00:00 --account=populationgenomics --pty bash

## Jupyter




# Managing access, data, and software for students on the cluster



Data folder `data` is placed under here and them symlinked (for student access) from `populationgenomics/data`




Instructor conda environment for building data files with `workflow_bedfiles.py`:

    conda create -n popgen_data -c gwforg -c bioconda python=3 gwf bedtools samtools



# Instructor zone

This page is for practical information related exercises on the cluster.
## The populationgenomics project folder

The project folder `populationgenomics` has the following folders:

- `instructor_data`: This links to a folder under the git repository [PopulationGenomicsCourse](https://github.com/kaspermunch/PopulationGenomicsCourse). The folder contains the original and intermediary files needed for building the files that the students need for the exercises.
- `data` is the folder where the student finds the data files he/she needs for the exercises. This folder will contain symlinks to the few required files in the instructor_data folder.
- `software`: builds and binaries for LDhat and PSMC
- `people`: teachers and instructors folders
- `students`: student folders

# Student data

To build bam and fastq files for the chr2 region run this in `data/bamfiles`:

    conda create -n popgen_data -c gwforg -c bioconda python=3 gwf bedtools samtools

    conda activate popgen_data

    sbatch download_bamfiles.sh

    gwf -f ../../scripts/workflow_bamfiles.py run

## Software

All command line software, except PSMC and LDhat, is installed in the `popgen` environment.  and LDhat PSMC is available in the software folder

## Running on the cluster

The students begin each exercise by running this to get an interactive session for the duration of the exercise:


    srun --mem-per-cpu=1g --time=3:00:00 --account=populationgenomics --pty bash

## Jupyter








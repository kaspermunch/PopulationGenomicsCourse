# Inferring the genotype phase (AKA Phasing)

During base calling, we identified the two bases at each position in each diploid individual. However, we do not know which base goes on which of the two chromosomes. That means that we do not know if the two haploid chromosomes look like the left or right example below:

    -----A-----C------         -----T-----C------ 
    -----T-----G------    or   -----A-----G------


To do that we use the program Beagle, which uses a clusering algorithm to call the genotype phase.

We put the jointly called bases for Africans, West Eurasians, and East Asians in these three files:

- Africa (9 individuals): `~/populationgenomics/data/vcf/Allvariants_africa.vcf`
- West Eurasia (10 individuals): `~/populationgenomics/data/vcf/Allvariants_westeurasia.vcf`
- Eash Asia (8 individuals): `~/populationgenomics/data/vcf/Allvariants_eastasia.vcf`

make a soft link to your own folder of the files above and this one:
~/populationgenomics/data/genetic_map/plink.chr2.GRCh37.map


In this exercise we will just use the Africans and the West Eurasians.

## Log into the cluster and request a compute node

Log into the cluster. Then request a machine for your computations. You need five gigabytes (`5g`) in this exercise so you need to run this command (see also the explanation in the previous exercise):

```bash
srun --mem-per-cpu=5g --time=3:00:00 --account=populationgenomics --pty bash
```

## Install and Activate Todays Environment 
```bash
    conda create env -f ~/populationgenomics/env/exercise_envs/phasing_wk4
    conda activate phasing_wk4
```

## Running Beagle Without a Reference Panel

For additional information see the [Beagle 4.1 manual](https://faculty.washington.edu/browning/beagle/beagle_4.1_03Oct15.pdf)

In order to obtain phased data, Beagle needs a genetic map. You can find the genetic map for chr2 (hg19 assembly) here:
~/populationgenomics/data/genetic_map/plink.chr2.GRCh37.map

Africa:

    beagle gt=Allvariants_africa.vcf map=plink.chr2.GRCh37.map out=Allvariants_africa_phased

West Eurasia

    beagle gt=Allvariants_westeurasia.vcf map=plink.chr2.GRCh37.map out=Allvariants_westeurasia_phased
    
Vcf files can both be compressed (gz) or uncompressed. IGV needs it in an uncompressed format, so decompress using

```bash
gunzip -c Allvariants_africa_phased.vcf.gz > Allvariants_africa_phased_t.vcf
```
This command outputs the decompressed to stdout, which then pipes into your file name of choice.

## Visualize the Phasing
before we can download the file and visualize it in IGV, we need to create another file, using [whatshap](https://whatshap.readthedocs.io/en/latest/guide.html#visualizing-phasing-results) 

Run this code below and replace <phased> with the file name you what to visualize.
```bash
    whatshap stats --gtf=<phased>.gtf <phased>.vcf
```



## Browsing the phased results

Download the phased VCF files and the gft files to your computer and open them in IGV (integrative genomics viewer): 
    
1. Choose Human hg19 as the reference genome.
2. Click `File > Load from File...` and select you phased VCF file.

Explore phases of haplotypes at two positions in the alignment:

Select chr2, zoom all the way in and select find the base at position 136608646. First, take a look at the WestEurasian samples. Consider these questions while zooming further out:

1. What does the haplotypes look like?
2. Do you see any long streches of homozygosity?
3. Which haplotypes agree?
4. How wide is the region where they agree?

To help derive your answers, make use of the metadata file: ~/populationgenomics/data/metadata/Sample_meta_subset.tsv

Now, compare it with the African samples.

Try to search the position chr2:136608646 in the [UCSC genome browser](https://genome-euro.ucsc.edu/cgi-bin/hgGateway?redirect=manual&source=genome.ucsc.edu). Remember we are using the Hg19 assembly version of the reference human genome. Can you find anything that explains your observations? (HINT: https://omim.org/entry/601806#0004)




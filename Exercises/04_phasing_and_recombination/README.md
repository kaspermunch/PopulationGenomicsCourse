# Inferring the genotype phase (AKA Phasing)

During base calling, we identified the two bases at each position in each diploid individual. However, we do not know which base goes on which of the two chromosomes. That means that we do not know if the two haploid chromosomes look like the left or right example below:

    -----A-----C------         -----T-----C------ 
    -----T-----G------    or   -----A-----G------


To do that we use the program Beagle, which uses a clusering algorithm to call the genotype phase.

We put the jointly called bases for Africans, West Eurasians, and East Asians in these three files:

- Africa (9 individuals): `~/populationgenomics/data/vcf/Allvariants_africa.vcf`
- West Eurasia (10 individuals): `~/populationgenomics/data/vcf/Allvariants_westeurasia.vcf`
- Eash Asia (8 individuals): `~/populationgenomics/data/vcf/Allvariants_eastasia.vcf`

In this exercise we will just use the Africans and the West Eurasians.

## Log into the cluster and request a compute node

Log into the cluster. Then request a machine for your computations. You need five gigabytes (`5g`) in this exercise so you need to run this command (see also the explanation in the previous exercise):

```bash
srun --mem-per-cpu=5g --time=3:00:00 --account=populationgenomics --pty bash
```

## Running Beagle

For additional information see the [Beagle 4.1 manual](https://faculty.washington.edu/browning/beagle/beagle_4.1_03Oct15.pdf)

In order to obtain phased data, Beagle needs a genetic map. You can find the genetic map for chr2 (hg19 assembly) here:
~/populationgenomics/data/genetic_map/plink.chr2.GRCh37.map

Africa:

    beagle gt=Allvariants_africa.vcf map=plink.chr2.GRCh37.map out=Allvariants_africa_phased

West Eurasia

    beagle gt=Allvariants_westeurasia.vcf map=plink.chr2.GRCh37.map out=Allvariants_westeurasia_phased

## Browsing the phased results

Download the phased VCF files to your computer and open them in IGV (integrative genomics viewer): 
    
1. Choose Human hg19 as the reference genome.
2. Click `File > Load from File...` and select you phased VCF file.

Explore phases of haplotypes at two positions in the alignment:

Select chr2, zoom all the way in and select find the base at position 136608646. First, take a look at the WestEurasian samples. Sort alignments by genotype (right-click on the base in the alignment tracks to get a popup menu). Consider these questions:

1. What does the haplotypes look like?
2. Do you see any long streches of homozygosity?
3. Which haplotypes agree?
4. How wide is the region where they agree?

To derive your answers, make use of the metadata file: ~/populationgenomics/data/metadata/Sample_meta_subset.tsv

Now, compare it with the African samples.

Try to search the position chr2:136608646 in the [UCSC genome browser](https://genome-euro.ucsc.edu/cgi-bin/hgGateway?redirect=manual&source=genome.ucsc.edu). Remember we are using the Hg19 assembly version of the reference human genome. Can you find anything that explains your observations? (HINT: https://omim.org/entry/601806#0004)

# Estimating a recombination map

[LDhat  manual](https://github.com/auton1/LDhat/blob/master/manual.pdf)

## Format input data for LDhat

LDhat needs its input data in a particular format. We will use vcftools to produce these input files from the phased VCF file.

First, we will need to install vcftools in our conda environment. To do so, run:

    conda install vcftools

Africa:

    vcftools --gzvcf Allvariants_africa_phased.vcf.gz --chr 2 --ldhat --out recmap_data_africa

West Eurasia:

    vcftools --gzvcf Allvariants_westeurasia_phased.vcf.gz --chr 2 --ldhat --out recmap_data_westeurasia

Have a look at the two files produced using `less`. E.g.:

    less recmap_data_africa.ldhat.sites 

> NB: press `q` to quit `less`

How do you think the information is encoded in these files?

## Running LDhat

To speed up computations you can make a lookup table first. That takes a while, so we did if for you. But it is done using the `complete` program that comes with LDhat:

    ~/populationgenomics/software/complete -n 20 -rhomax 100 -n_pts 101 -theta 0.0001

- `-n 20`:the number of haplotypes (2 * 10).
- `-rhomax 100`: maximum rho ($4N_e r$) alowed: 100 (recommended).
- `-n_pts 101`: number of points in grid: 101 (recommended).
- `-theta 0.0001`: human theta ($4N_e \mu$).

It produces a file that will serve as a look up table for the algorithm. It includes coalescent likelihoods for each pairs of SNPs using a grid of recombination rates. You can find it here:

~/populationgenomics/data/ldhat/new_lk.txt

The next step is to calculate the recombination map. This command will take some time to run (~ 6 min). We suggest you team up in pairs, so that one of you runs the command for the African samples while other runs the command for the West Eurasian samples. 

Africa:

    ~/populationgenomics/software/rhomap -seq recmap_data_africa.ldhat.sites -loc recmap_data_africa.ldhat.locs -lk new_lk.txt -its 100000 -samp 500 -burn 0

<!-- 3m9.888s -->
*or* West Eurasia:

    ~/populationgenomics/software/rhomap -seq recmap_data_westeurasia.ldhat.sites -loc recmap_data_westeurasia.ldhat.locs -lk new_lk.txt -its 100000 -samp 500 -burn 0

- `-lk`: likelihood lookup table.
- `-its`: number of iterations of the MCMC.
- `-samp`: number of iterations between sampling events, i.e how often to sample from the MCMC.
- `-burn`: how many of the initial iterations to discard. Here we set it to zero to leave keep all samples. Then we look later how much burnin to discard.

When rhomap completes it writes three files:

- `acceptance_rates.txt`: acceptance rates of the MCMC. If they are lower than 1%. The program should be run with more iterations.
- `summary.txt`: (quoting the manual) for each SNP interval, the estimated genetic map position, the estimated recombination rate, and the hotspot density (the number of hotspots per kb per iteration).
- `rates.txt`: (quoting the manual) is the output from each sample detailing the recombination rate (expressed in $4N_e r$ per kb) between each SNP. 

Rename the `rates.txt` to `rates_africa.txt` or `rates_westeurasia.txt` depending on which analysis you do. E.g.:

```bash
mv acceptance_rates.txt acceptance_rates_africa.txt
mv summary.txt summary_africa.txt
mv rates.txt rates_africa.txt
```      

## Analyze results

Open a jupyter notebook in R and paste this set of commands:

<!-- This needs to run first to allow source to pull from web -->
<!-- Sys.setenv(https_proxy = "http://in:3128", http_proxy = "http://in:3128") -->

```R
source("~/populationgenomics/data/ldhat/ldhat.r")
```

That loads a lot of R functions written by the author of LDhat.

Now run this code (if you analyze west eurasians you need to change the arguments accordingly):

```R
summarise.rhomap(rates.file = "rates_africa.txt", locs.file="recmap_data_africa.ldhat.locs")
```

The summary produces two plots:

- A graph of the recombination rate across on each polymorphic loci, along with confidence intervals.
- A plot showing how estimation of recombination rate has progressed with each MCMC sample. Notice that the initial run of MCMC samples are atypical. This is the "burn-in" of the MCMC. We want to remove that, so take notice of how many samples it corresponds to. We can produce a new set of estimates that excludes this burn-in using the `stat` program that comes with LDhat:

```bash
~/populationgenomics/software/stat -input rates_africa.txt -loc recmap_data_africa.ldhat.locs -burn 60
```

This produces a file called `res.txt` that describes the confidence in the estimated recombination rate along the sequence. Rename the `res.txt` to `res_africa.txt` or `res_westeurasia.txt` depending on which analysis you do. E.g.:

```bash
mv res.txt res_africa.txt
```

Now try to plot the final results (if you analyze West Eurasians you must read `res_westeurasia.txt`). In order to have the positions of the loci for which we have estimated mean recombination rates, we will merge the new dataset created with the summary files generated by LDhat:

```R
library(ggplot2)
library(dplyr)
library(tidyr)
library(magrittr)

summary <- read.table('summary_africa.txt', header = T)
rates <- read.table("res_africa.txt", header=T)
rates %>%
  filter(Loci > 0) %>%
  mutate(pos=summary$Position.kb.*1000) %>%
  ggplot(aes(x=pos, y=Mean_rho, ymin=L95, ymax=U95)) +  
  geom_line(color='blue') +
  geom_ribbon(alpha=0.1) +
  theme_bw()
```

Look at the plots and ponder the following questions:

- Are there any recombination hotspots?
- Are there any regions where the estimated recombination rate is really low? 
- Can you see any hotspots in Africans that are not found in West Eurasians - other the othe way around?
- What does the recombination rate look like around the lactase gene?

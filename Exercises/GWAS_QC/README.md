# GWAS Quality Control using PLINK
___

One of the most important aspects of a GWAS is to do thorough Quality control (QC) of the data before analysing it.

### Software:
We will be using plink 1.9. Plink is a comprehensive tool for handling and analyzing SNP data that can perform many different kinds of analyses. Check the documentation [here](https://www.cog-genomics.org/plink/1.9/)
If you want info about a specific command you can also use help command:
```
plink --help <name_of_command>
```

### Data:
In this practical, we will go through the steps in performing quality control (QC) of genotype data from a simulated genome-wide association study of 1000 cases and 1000 controls, typed for 317,503 autosomal and X chromosome SNPs.

We will use the following files:

```bash
~/populationgenomics/data/GWAS/GWAS_QC_data/GWA-data.bed
~/populationgenomics/data/GWAS/GWAS_QC_data/GWA-data.bim
~/populationgenomics/data/GWAS/GWAS_QC_data/GWA-data.fam
```

The practical is based on “Data quality control in genetic case-control association studies” (Anderson et al. 2010, Nature Protocols 5: 1564-73).

### Exercise contents:
We will begin by performing sample QC, including calculation of call rates, heterozygosity, and sex discordance.  We will then perform SNP QC, including calculation of call rates and deviation from Hardy-Weinberg equilibrium.  

## Sample QC

### Identification of individuals with discordant sex information
Remember to start an interactive job:

```
srun --mem-per-cpu=1g --time=3:00:00 --account=populationgenomics --pty bash
```

At the shell prompt, type:

```
plink --bfile GWA-data --check-sex --out GWA-QC
```

This command will infer the sex of the sample by looking at the mean homozygosity rate across X-chromosome markers and compare it the sex stated in the `.fam` file.

*1) Take a look at the output file “GWA-data.sexcheck”. How many problematic samples are there?*

Problematic samples can be removed by copying the family ID (FID) and individual ID (IID) of the samples into a text file (e.g. wrong_sex.txt) and using the remove command:

HINT: To filter for problematic inds, either load it into a notebook or use grep.

```
plink --bfile GWA-data --remove wrong_sex.txt --make-bed --out GWA-QC
```

The `--out` option in plink specifies the prefix of the output files that plink generates. And when we use the --make-bed command we are writing the output to the specified prefix. In this case, all our output files will have the prefix: "GWA-QC". 

*2) Each time a plink command is run it writes a summary to a log file (the file name ends with `.log`). Look at the log file after removing the problematic individuals. How many cases and controls are left in the data set?*

### Identification of individuals with elevated missing data rates or outlying heterozygosity rate
At the shell prompt, type:

```
plink --bfile GWA-QC --missing --out GWA-QC
```

This command will create the files `GWA-QC.imiss` and `GWA-QC.lmiss`.  The fourth column in the file `GWA-data.imiss` (N_MISS) denotes the number of missing SNPs and the sixth column (F_MISS) denotes the proportion of missing SNPs per individual.

At the shell prompt type:

```
plink --bfile GWA-QC --het --out GWA-QC 
```

This command will create the file `GWA-data.het`, in which the third column denotes the observed number of homozygous genotypes [O(Hom)] and the fifth column denotes the number of non-missing genotypes [N(NM)] per individual.

You can calculate the observed heterozygosity rate per individual using the formula:

Het = (N(NM) − O(Hom))/N(NM)

*3) Open a R jupyter notebook and create a plot in which the observed heterozygosity rate per individual is plotted on the x axis and the proportion of missing SNPs per individuals is plotted on the y axis. Interpret your results.*

**Hint:** You can merge the two tables by running something like this:
```R
library(dplyr)
d_miss <- read.table("GWA-QC.imiss",header=T)
d_het <- read.table("GWA-QC.het",header=T)
d <- inner_join(d_miss,d_het)
```

We will filter out outliers for either of the variables, and then save the file like this: 

```
write.table(d_m, file = "wrong_het_missing_values.txt", col.names = F, row.names = F)
```

*4) Make a file with the FID and IID of all individuals that have a genotype missing rate >=0.03 or a heterozygosity rate that is more than 3 s.d. from the mean. Then use plink to remove these individuals from the data set.*

```bash
plink --bfile GWA-QC --remove wrong_het_missing_values.txt --make-bed --out GWA-QC
```
Note that by providing the same prefix and running --make-bed, we will overwrite our bed/bim/fam files. This is OK in the context of this exercise and helps simplifying things, but can be dangerous in real life.

### Identification of duplicated or related individuals
To identify duplicated or related individuals we will calculate the identity by descent (IBD) matrix. This works best if it is done on a set of non-correlated SNPs. So first we will “prune” the data and create a list of SNPs where no pair (within a given genomic interval) has an r2 value greater than a given threshold, typically chosen to be 0.2.  This can be done by the `indep-pairwise` command, using 500kb as window size and 5 variants as step size:

```
plink --bfile GWA-QC --indep-pairwise 500kb 5 0.2 --out GWA-QC
```

It saves the list of independent SNPs as `GWA-QC.prune.in` (This data set was simulated without LD so in this case there will not be a lot of variants removed.)
To calculate IBD between each pair of individuals, type the following command at the shell prompt:

```
plink --bfile GWA-QC --extract GWA-QC.prune.in --genome --min 0.185 --out GWA-QC
```

The `--min 0.185` option means that it will only print the calculated IBD if it is above 0.185 (Mean between second-degree relatives:0.25 and third-degree relatives:0.125). The PI_HAT value in column 10 of the output file will be a number between 0 and 1 saying how much of the genome the two individuals share (1 for identical twins, 0.5 for siblings etc.). This command will produce a file called GWA-data.genome .

*5) Remove a member from each of the pairs that are too closely related from the data set. To keep it simple you can just always remove the individual mentioned first. *

**Hint:** Open the GWA-QC.genome file in your jupyter notebook and obtain unique IDs with something of the type:
```R
ibd <- read.table('GWA-QC.genome', header = TRUE)
members <- ibd$FID1
members <- unique(members)
write.table(cbind(members,members), file = 'wrong_ibd.txt', col.names = F, row.names = F)
```
Here we are exploiting the fact that FID and IID are the same in our dataset, but we would have to contition on both if that wasn't the case.

To remove these individuals, we will use again the --remove parameter and create updated bed/bim/fam files:

```
plink --bfile  GWA-QC --remove wrong_ibd.txt --make-bed --out GWA-QC
```

## SNP QC
### SNPs with an excessive missing data rate
Run the `--missing` command again to generate the `GWA-data.lmiss` with the missing data rate for each SNP.

*6) Use R to make a histogram of the missing data rates (F_MISS).*

The `--test-missing` command tests for association between missingness and case/control status, using Fisher's exact test. It produces a file with ".missing" suffix.

*7) Run the test-missing command and make a list of all the names of all SNPs where the differential missingness p-value is less than 1e-5. Save the list as `fail-diffmiss-qc.txt`.*

To remove low-quality SNPs, type the following command at the shell prompt:

```
plink --bfile GWA-QC --exclude fail-diffmiss-qc.txt --geno 0.05 --hwe 0.00001 --maf 0.01 --make-bed --out GWA-QC
```

In addition to removing SNPs identified with differential call rates between cases and controls, this command removes SNPs with call rate less than 95% with `--geno` option and deviation from HWE (p<1e-5) with the `--hwe` option. It also removes all SNPs with minor allele frequency less than a specified threshold using the `--maf` option.

*8) How many SNPs are left in the clean data set?*

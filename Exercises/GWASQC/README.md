# GWAS Quality Control using PLINK
___

One of the most important aspects of a GWAS is to do thorough Quality control (QC) of the data before analysing it.

### Software:
We will be using plink 1.9. Plink is a comprehensive tool for handling and analyzing SNP data that can perform many different kinds of analyses. Check the documentation [here](https://www.cog-genomics.org/plink/1.9/)
If you want info about a specific command you can also use help command:
```
plink --help <name_of_command>
```
We will also be using R and Rstudio to make plots and make simple calculations.

### Data:
In this practical, we will go through the steps in performing quality control (QC) of genotype data from a simulated genome-wide association study of 1000 cases and 1000 controls, typed for 317,503 autosomal and X chromosome SNPs.

The data set can be found [here](https://drive.google.com/open?id=1XC-BF1ikrhoJ6JrklaJDjCijVMqfwmYe)

The practical is based on “Data quality control in genetic case-control association studies” (Anderson et al. 2010, Nature Protocols 5: 1564-73).

### Exercise contents:
We will begin by performing sample QC, including calculation of call rates, heterozygosity, and sex discordance.  We will then perform SNP QC, including calculation of call rates and deviation from Hardy-Weinberg equilibrium.  

## Sample QC

### Identification of individuals with discordant sex information
At the shell prompt, type:
```
plink --bfile GWA-data --check-sex --out GWA-data
```
This command will infer the sex of the sample by looking at the mean homozygosity rate across X-chromosome markers and compare it the sex stated in the `.fam` file.

*1) Take a look at the output file “GWA-data.sexcheck”. How many problematic samples are there?*

Problematic samples can be removed by copying the family ID (FID) and individual ID (IID) of the samples into a text file (e.g. wrong_sex.txt) and using the remove command:
```
plink --bfile GWA-data --remove wrong_sex.txt --make-bed --out GWA-data
```
The `--out` option in plink specifies the prefix of the output files that plink generates. And when we use the --make-bed command with the same prefix as the input we are actually overwriting the input files. This is OK for these exercises, but on a “real” data set you might not want to do that.

*2) Each time a plink command is run it writes a summary to a log file (the file name ends with `.log`). Look at the log file after removing the problematic individuals. How many cases and controls are left in the data set?*

### Identification of individuals with elevated missing data rates or outlying heterozygosity rate
At the shell prompt, type:
```
plink --bfile GWA-data --missing --out GWA-data
```
This command will create the files `GWA-data.imiss` and `GWA-data.lmiss`.  The fourth column in the file `GWA-data.imiss` (N_MISS) denotes the number of missing SNPs and the sixth column (F_MISS) denotes the proportion of missing SNPs per individual.

At the shell prompt type:
```
plink --bfile GWA-data --het --out GWA-data
```
This command will create the file `GWA-data.het`, in which the third column denotes the observed number of homozygous genotypes [O(Hom)] and the fifth column denotes the number of non-missing genotypes [N(NM)] per individual.

You can calculate the observed heterozygosity rate per individual using the formula:

Het = (N(NM) − O(Hom))/N(NM)

*3) Use R to create a plot in which the observed heterozygosity rate per individual is plotted on the x axis and the proportion of missing SNPs per individuals is plotted on the y axis.*

*4) Use R to make a file with the FID and IID of all individuals that have a genotype missing rate >=0.03 or a heterozygosity rate that is more than 3 s.d. from the mean. Then use plink to remove these individuals from the data set.*

### Identification of duplicated or related individuals
To identify duplicated or related individuals we will calculate the identity by descent (IBD) matrix. This works best if it is done on a set of non-correlated SNPs. So first we will “prune” the data and create a list of SNPs where no pair (within a given genomic interval) has an r2 value greater than a given threshold, typically chosen to be 0.2.  This can be done by the `indep-pairwise` command:
```
plink --bfile GWA-data --indep-pairwise 500kb 5 0.2 --out GWA-data
```
It saves the list of independent SNPs as `GWA-data.prune.in` (This data set was simulated without LD so in this case there will not be a lot of variants removed.)
To calculate IBD between each pair of individuals, type the following command at the shell prompt:
```
plink --bfile GWA-data --extract GWA-data.prune.in --genome --min 0.185 --out GWA-data
```
The `--min 0.185` option means that it will only print the calculated IBD if it is above 0.185. The PI_HAT value in column 10 of the output file will be a number between 0 and 1 saying how much of the genome the two individuals share (1 for identical twins, 0.5 for siblings etc.).

*5) Remove a member from each of the pairs that are too closely related from the data set. To keep it simple you can just always remove the individual mentioned first.*

## SNP QC
### SNPs with an excessive missing data rate
Run the `--missing` command again to generate the `GWA-data.lmiss` with the missing data rate for each SNP.

*6) Use R to make a histogram of the missing data rates.*

The `--test-missing` command tests all markers for differences in the call rate between cases and controls.

*7) Run the test-missing command and make a list of all the names of all SNPs where the differential missing-ness p-value is less than 10e-5. Save the list as `fail-diffmiss-qc.txt`.*

To remove low-quality SNPs, type the following command at the shell prompt:
```
plink --bfile GWA-data --exclude fail-diffmiss-qc.txt --geno 0.05 --hwe 0.00001 --maf 0.01 --make-bed --out GWA-data
```
In addition to removing SNPs identified with differential call rates between cases and controls, this command removes SNPs with call rate less than 5% with `--geno` option and deviation from HWE (p<10-5) with the `--hwe` option. It also removes all SNPs with minor allele frequency less than a specified threshold using the `--maf` option.

*8) How many SNPs are left in the clean data set?*

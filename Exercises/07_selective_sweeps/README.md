Background
----------

In this exercise section you will analyse positive selection in human populations. You will be looking at the X chromosome of individuals from the Simons genome diversity project. The advantage of X chromosomes in males is that they are haploid, and therefore, fully phased. We also have other reasons to belive that they often experience natural selection. Finally, they have typically not been investigated to the same extent as the autosomes.

You will have data from the following populations:

| Population  | Individuals |
|-------------|:-----------:|
| Africa      |   26 males  |
| West Europe |   44 males  |
| South Asia  |   30 males  |
| East Asia   |   24 males  |

![](img/unnamed-chunk-1-1.png)


The data consists of **24198** SNPs from the region 73-81 Mb on the X chromosome and there is no missing data. The haplotype data for each population is found in separate files (**genotypes360\_400\_.**), whereas they use a common SNP identity file **snps360\_400\_filtered.snp**. You can find all the necessary files here:

~/populationgenomics/data/haplotypes_chrX

Package
-------

The package requires 2 inputs:

1.  The snp file (`snp_ID`, `Chromosome`, `position`, `reference allele` and `derived allele`).

2.  The haplotype file: haplotype of each individual encoded as reference allele, alternative allele and missing data N.

Analysis
--------

You will perform a genome wide scan and then focus on candidate SNPs. The package that you will be using on these analysis is `rehh`. The manual of the package can be found [here](https://cran.r-project.org/web/packages/rehh/rehh.pdf). 

### Reading in data in REHH format

``` r
install.packages("rehh")
library(rehh)

> hap360_400_AF <-data2haplohh(hap_file="genotypes360_400_AF",map_file="snps360_400_filtered",
                               allele_coding="map", 
                               min_perc_geno.mrk=100,
                               min_perc_geno.hap=100,
                               haplotype.in.columns=TRUE,
                               chr.name=1)
```

#### Q1. How many haplotypes and snps are found in each population?

Scan the region using iHS, Rsb and XP-EHH
-----------------------------------------

You should first perform a scan for each of the regions for extreme values of integrated haplotype score (iHS). This calculation is done by standardizing to differences in allele frequency. You can do that by using the functions `scan_hh` and followed by `ihh2ihs`. Then you can visualize the results with ggplot.

#### Q2. Try also to plot a histogram of the allele frequencies of the SNPs in each population. Do you observe population differences?

Hint: Allele frequencies are calculated and stored as part of the dataframe resulted from scan\_hh. Use (par mfrow) function to combine the 4 different population plots. 

``` r
res.scanAF<-scan_hh(hap360_400_AF)
head(res.scanAF) # Inspect the dataframe
```
``` r
library(ggplot2)
library(dplyr)
res.scanAF %>% ggplot() +
  geom_histogram(aes(x=FREQ_A))
res.scanAF %>% ggplot() +
  geom_histogram(aes(x=FREQ_D)) 
```

Once you have obtained the data\_frame produced by `scan_hh` you can compute the standardized iHS (iHH), as described in [Voight et al. (2006)](http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0040072).

#### Q3. For what reason do they standardize iHS measure?

``` r
wg.ihsAF<-ihh2ihs(res.scanAF, freqbin = 0.05) 
manhattanplot(wg.ihsAF)
```

#### Q4. Do you find outliers?

If so, record the SNP positions of the most significant SNPs for later analysis, using e.g. which.max() and which.min().

Repeat the analysis for the other two populations.

Perform pairwise population tests
---------------------------------

There are two possible functions to be used, `ies2rsb` and `ies2xpehh`, both require the dataframes of `scan_hh` class.

#### Q5. Can you see the differences between the different methods?

``` r
# Have a look at the differences between the functions:
#?ies2rsb()
#?ies2xpehh()

# Computing ies2rsb between African population over West European population:
wg.rsbAFWE <- ines2rsb(res.scanAF,res.scanWE, popname1 = "Africa", popname2 = "West_Europe")
manhattanplot(wg.rsbAFWE)
```

``` r
wg.XPEHHAFWE <- ies2xpehh(res.scanAF,res.scanWE, popname1 = "Africa", popname2 = "West_Europe")
manhattanplot(wg.XPEHHAFWE)
```

Zooming in interesting markers
------------------------------

From the scan you can find SNPs that give extreme values of iEHS or of XPEHH for a set of populations. You can then analyse the haplotype structure around them. This is done by including the index position of the interested marker in the functions `calc_ehhs` and `bifurcation.diagram`.

Try to plot markers that show outlier values in the above statistics and compare populations. Hint: use which.max() and which.min() (especially when using XPEHH or Rsb). Select a SNP that shows some interesting results.

``` r
marker = which.min(wg.XPEHHAFWE$XPEHH_Africa_West_Europe)
snp = row.names(wg.XPEHHAFWE)[marker]
a = calc_furcation(hap360_400_AF, mrk=snp)
plot(a))
```

#### Q6. What is the biological function of the region around this snp?

Have a look at UCSC Genome Browser and remember that this dataset belongs to chromosome X and Hg19 as reference assembly.

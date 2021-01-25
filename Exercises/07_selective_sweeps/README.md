Background
----------

In this exercise section you will analyse positive selection in human populations. You will be looking at the X chromosome of individuals from the Simons genome diversity project. The advantage of X chromosomes in males is that they are haploid, and therefore, fully phased. We also have other reasons to belive that they often experience natural selection. Finally, they have typically not been investigated previously to the same extent as the autosomes.

You will have data from the following populations:

| Population  | Individuals |
|-------------|:-----------:|
| Africa      |   26 males  |
| West Europe |   44 males  |
| South Asia  |   30 males  |
| East Asia   |   24 males  |

![](img/unnamed-chunk-1-1.png)

<!-- TODO: Update file names -->

The data consists of **24198** SNPs from the region 73-81 Mb on the X chromosome and there is no missing data. The haplotype data for each population is found in separate files (**genotypes360\_400\_.**), whereas they use a common SNP identity file **snps360\_400\_filtered.snp**.

Package
-------

The package requires 2 inputs:

1.  The snp file (`snp_ID`, `Chromosome`, `position`, `reference allele` and `derived allele`).

2.  The haplotype file: haplotype of each individual encoded as reference allele, alternative allele and missing data N.

Analysis
--------

You will perform a genome wide scan and then focus on candidate SNPs. The package that you will be using on these analysis is `rehh`. The manual of the package can be found [here](https://cran.r-project.org/web/packages/rehh/rehh.pdf). 

### Reading in data in REHH format


<!-- TODO: Update file names -->

``` r
install.packages("rehh")
library(rehh)

> hap360_400_AF <-data2haplohh(hap_file="/home/Data/genotypes360_400_AF",map_file="/home/Data/snps360_400_filtered",
+                              allele_coding="map", 
+                              min_perc_geno.mrk=100,
+                              min_perc_geno.hap=100,
+                              haplotype.in.columns=TRUE,
+                              chr.name=1)
* Reading input file(s) *
Map info: 24198 markers declared for chromosome 1 .
Haplotype input file in transposed format assumed.
Alleles are being recoded according to fourth and fifth column of map file.
* Filtering data *
Discard haplotypes with less than 100 % of genotyped markers.
No haplotype discarded.
Discard markers genotyped on less than 100 % of haplotypes.
No marker discarded.
Data consists of 26 haplotypes and 24198 markers.
Number of mono-, bi-, multi-allelic markers:
1 2 
10018 14180 
```

#### Q1. How many haplotypes and snps are found in each population?

Scan the region using iHS, Rsb and XP-EHH
-----------------------------------------

You should first perform a scan for each of the regions for extreme values of **i**ntegrated **h**aplotype **s**core (iHS). This calculation is done by standardizing to differences in allele frequency. This implies using the functions `scan_hh`, followed by `ihh2ihs` and then plotting the results using either ihsplot or plotting by yourself.

#### Q2. Try also to plot a histogram of the allele frequencies of the SNPs in each population. Do you observe population differences?

Hint: Allele frequencies are calculated and stored as part of the dataframe resulted from scan\_hh. Use (par mfrow) function to combibe the 4 different population plots. Once you have obtained the data\_frame produced by `scan_hh` you can compute the standardized iHS (iHH), as described in [Voight et al. (2006)](http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0040072).

``` r
# Producing the required input dataframe:
> res.scanAF<-scan_hh(hap360_400_AF)
> head(res.scanAF)
             CHR POSITION     FREQ_A    FREQ_D NHAPLO_A NHAPLO_D IHH_A IHH_D IES INES
X:X_73263128   1 73263128 0.00000000 1.0000000        0       26     0    NA  NA   NA
X:X_73264012   1 73264012 0.03846154 0.9615385        1       25     0    NA  NA   NA
X:X_73264487   1 73264487 0.03846154 0.9615385        1       25     0    NA  NA   NA
X:X_73264603   1 73264603 0.03846154 0.9615385        1       25     0    NA  NA   NA
X:X_73265256   1 73265256 0.00000000 1.0000000        0       26     0    NA  NA   NA
X:X_73266095   1 73266095 1.00000000 0.0000000       26        0    NA     0  NA   NA

library(ggplot2)
res.scanAF %>% ggplot() +
  geom_histogram(aes(x=FREQ_A))
res.scanAF %>% ggplot() +
  geom_histogram(aes(x=FREQ_D)) 
```

#### Q3. For what reason do they standardize iHS measure?

``` r
# Scanning each population at time:
wg.ihsAF<-ihh2ihs(res.scanAF, freqbin = 0.05) 

# Plotting the results:
manhattanplot(wg.ihsAF)
```

#### Q4. Do you find outliers?

If so, record the SNP positions of the most significant SNPs for later analysis, using e.g. which.max() or which.min().

Perform pairwise population tests
---------------------------------

There are two possible functions to be used, `ies2rsb` and `ies2xpehh`, both require the dataframes of `scan_hh` class.

#### Q5. Can you see the differences between the different methods?

``` r
# Have a look at the differences between the functions:
#?ies2rsb()
#?ies2xpehh()

# Computing ies2rsb between African population over West European population:
wg.rsbAFWE <- ines2rsb(res.scanAF,res.scanWE, popname1 = "Africa", popname2 = "W Europe")
manhattanplot(wg.rsbAFWE)
```

Use the function rsbplot() and xpehhplot() to explore and plot your results:

``` r
wg.XPEHHAFWE <- ies2xpehh(res.scanAF,res.scanWE, popname1 = "Africa", popname2 = "W Europe")
manhattanplot(wg.XPEHHAFWE)
```

Zooming in interesting markers
------------------------------

From the scan you can find SNPs that give extreme values of iEHS or of XPEHH for a set of populations. You can then analyse the haplotype structure around them. This is done by including the index position of the interested marker in the functions `calc_ehhs` and `bifurcation.diagram`.

Try to plot markers that show outlier values in the above statistics and compare populations. Hint: use which.max() and which.min() (especially when using XPEHH or Rsb). Select a SNP that shows some interesting results.

``` r
marker = which.min(wg.XPEHHAFWE["XPEHH_Africa_W Europe"])
snp = row.names(wg.XPEHHAFWE)[marker]
a = calc_furcation(hap360_400_AF, mrk=snp)
plot(a)
```

#### Q6. What is the biological function of the region around this snp?

Have a look at UCSC Genome Browser and remember that this dataset belongs to chromosome X and Hg19 as reference assembly.

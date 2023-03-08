# vcfR exercises

## Data
We will devote this exercise session to analyze variants called for our dataset of 28 individuals from across the world. For that we will make use of a VCF file, such as the one we generated last week, but including all 28 samples, thus making possible the study of population genetic's summary statistics from the variants called. Additionally, we will need the metadata file for annotating the samples in the VCF file.

- vcf file: ~/populationgenomics/data/vcf/chr2_135_145.vcf.gz
- annotation file: ~/populationgenomics/data/metadata/sample_infos_accessionnb.csv

We will use a jupyter notebook in R to do the analysis. Run the command to open jupyter.

## Useful links
During the exercise, we will use two very useful packages: [dplyr](https://cran.r-project.org/web/packages/dplyr/vignettes/dplyr.html) and [ggplot](https://monashbioinformaticsplatform.github.io/r-more/topics/tidyverse.html).

## vcfR
First, to read vcf-files into R, we'll need the R package `vcfR` . We can read the vcf-file like this:

```r
library(vcfR)
vcf <- read.vcfR("~/populationgenomics/data/vcf/chr2_135_145.vcf.gz")
```

We will use the `dplyr` package to analyse the data so we want the data to be in "tidy" format. We can get that by using the vcfR2tidy function.

```r
library(dplyr)
tvcf <- vcfR2tidy(vcf, 
          single_frame = TRUE,
          info_fields = c("TR"),
          format_fields = c("GT","GQ","NR"),
          info_types = TRUE,
          format_types = c(NR="i",GQ="i"))
```

One of the measures provided for each variant call is the Phred score, i.e the quality associated with that given base pair. It is computed by the following equation:

<a href="https://www.codecogs.com/eqnedit.php?latex=Q&space;=&space;-10&space;log10(P)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?Q&space;=&space;-10&space;log10(P)" title="Q = -10 log10(P)" /></a>

Being P the probability of error of a given the base call.

We also want to add some info about the samples and have all heterozygous sites in the same format:

```r
info <- read.csv2("~/populationgenomics/data/metadata/sample_infos_accessionnb.csv")
d <- inner_join(tvcf$dat,info, by= c("Indiv" = "ENA.RUN")) %>%
  mutate(gt_GT=replace(gt_GT, gt_GT=="1/0", "0/1"))
```

Now `d` contains the relevant data in "tidy" format. We can count the fraction of missing genotypes for each individual:
```r
d %>% 
  group_by(Indiv) %>%
  summarise(missing = mean(is.na(gt_GT))) %>%
  arrange(desc(missing))
```

And we can plot it using ggplot:
```r
library(ggplot2)
d %>% 
  group_by(Indiv) %>%
  summarise(missing = mean(is.na(gt_GT))) %>%
  ggplot(aes(x=Indiv,y=missing)) + geom_col() + coord_flip()
```

As you can see, individual `ERR1025639` is an outlier, so we want to remove that individual from the analyses.
```r
d <- d %>% filter(Indiv!="ERR1025639")
```

1. *Which individual has the largest amount of missing data now?*
2. *Some variants have values different from "PASS" in the FILTER column. These variants should be removed. Does that change the fraction of missing data?*
3. *The column gt_NR contains the number of reads covering the position. What is the average depth in the data set?*
4. *If a variant is heterozygous the gt_GT variable will have the value "0/1". Make a plot of the number of variants that are heterozygous for each individual. Which population has the highest fraction of heterozygous variants in this genomic region?*

All the variants we look at are bi-allelic so we can calculate the allele frequences in each sub-population (region) like this:
```r
d2 <- d %>%
  group_by(POS,region) %>%
  summarise(na=sum(gt_GT=="0/1",na.rm=T)+2*sum(gt_GT=="0/0",na.rm=T),
            nA=sum(gt_GT=="0/1",na.rm=T)+2*sum(gt_GT=="1/1",na.rm=T))  %>%
  mutate(pS=na/(na+nA), qS= nA/(na+nA))
```

We can then plot the number of polymorphic sites for each region:
```r
d2 %>% 
  filter(na!=0, nA!=0) %>%
  ggplot(aes(x=region, fill=region)) + geom_bar() + coord_flip()
```

## Calculating F-statistics

`FST` can be calculated as `FST = 1- (H_S/H_T)`. Where `H_T` is the expected heterozygosity of the entire population and `H_S` is the mean expected heterozygosity across subpopulations.
Using `d2` from above you can calculate `H_S` as `H_S=mean(2*pS*qS)` and `H_T=2*mean(pS)*mean(qS)`.

5. *Use `d2` from above to calculate `FST` for each position. What is the median FST value?*
6. *Make a histogram of the `FST` values. (Hint: use geom_histogram() instead of geom_bar())*
7. *Calculate `FST` using only the European and African samples and make a histogram of the values.*
8. *Make a plot with the genomic position on the x axis and the `FST` value on the y axis. (Hint: use geom_point() or geom_line())*

We can also look at bins of a given size along the genome (to keep it simple we will just plot non-overlapping bins instead of sliding windows). We can fx. plot the fraction of sites in each bin that are polymorphic for each subpopulation:
```r
bin_width = 500000
d2 %>% 
  mutate(binmid=((POS %/% bin_width)*bin_width)) %>%
  filter(na!=0, nA!=0) %>%
  group_by(binmid, region) %>%
  summarise(frac_polymorph=n()/bin_width) %>%
  ggplot(aes(x=binmid, y=frac_polymorph, color=region)) + geom_line()
```

9. *Make a plot with average FST in bins along the genome.*


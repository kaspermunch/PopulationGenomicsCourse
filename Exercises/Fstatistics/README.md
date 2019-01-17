# vcfR exercises

## Data
In this exercise we will analyse a vcf-file like the one you created yourselves last week. We will only use R so you can do the analyses on your own machine this week. The vcf-file is called `chr2_135_145.vcf.gz` and you can download it from Materials/Week3 on blackboard.

## Useful links
In this exercise we will be using two very useful packages: `dplyr` and `ggplot`. The learning curve can be quite steep, but once you learn it, it becomes very natural and useful! :)
You can find nice links to go deeper into your learning: [ggplot](https://monashbioinformaticsplatform.github.io/r-more/topics/tidyverse.html) and [dplyr](https://cran.r-project.org/web/packages/dplyr/vignettes/dplyr.html).

## vcfR
In this exercise will use the R package `vcfR` to read vcf-files into R. We can read the vcf-file like this:

```r
install.packages('vcfR')
library(vcfR)
vcf <- read.vcfR("chr2_135_145.vcf.gz")
```

We will use the `dplyr` package to analyse the data so we want the data to be in "tidy" format. We can get that by using the vcfR2tidy function.

```r
install.packages('dplyr')
library(dplyr)
tvcf <- vcfR2tidy(vcf, 
          single_frame = TRUE,
          info_fields = c("TR"),
          format_fields = c("GT","GQ","NR"),
          info_types = TRUE,
          format_types = c(NR="i",GQ="i"))
```

We also want ta add some info about the samples:
```r
library(dplyr)
info <- read.csv2("sample_infos_accessionnb.csv")
d <- inner_join(tvcf$dat,info, by= c("Indiv" = "ENA.RUN")) %>%
  mutate(gt_GT=replace(gt_GT, gt_GT=="1/0", "0/1"))
```

Now `d` contains the relevant data in "tidy" format. We can fx. count the fraction of missing genotypes for each individual:
```r
d %>% 
  group_by(Indiv) %>%
  summarise(missing = mean(is.na(gt_GT))) 
```

And we can plot it using ggplot:
```r
library(ggplot2)
d %>% 
  group_by(Indiv) %>%
  summarise(missing = mean(is.na(gt_GT))) %>%
  ggplot(aes(x=Indiv,y=missing)) + geom_col() + coord_flip()
```

As you can see individual `ERR1025639` is an outlier so we want to remove that individual from the analyses.
```r
d <- d %>% filter(Indiv!="ERR1025639")
```

1. *Which individual has the largest amount of missing data now?*
2. *Some variants have values different from "PASS" in the FILTER column. These variants should be removed. Does that change the fraction of missing data?*
3. *The column gt_NR contains the number of reads covering the position. What is the average depth in the data set?*
4. *If a variant is heterozygous the gt_GT variable will have the value "0/1". Make a plot of the number of variants that are heterozygous for each individual. Which population has the highest fraction of heterozygous variants in this genomic region?*

All the variants we look at are bi-allelic so we can calculate the allele frequences in each sub-population (region) like this (it's a bit slow):
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

A description of `F_ST` can be found in Box 5.2. on page 147 of HEG.

`F_ST` can be calculated as `F_ST = 1- (H_S/H_T)`. Where `H_T` is the expected heterozygocity of the entire population and `H_S` is the mean expected heterozygocity across subpopulations.
Using `d2` from above you can calculate `H_S` as `H_S=mean(2*pS*qS)` and `H_T=2*mean(pS)*mean(qS)`.

5. *Use `d2` from above to calculate `F_ST` for each position. What is the median F_ST value?*
6. *Make a histogram of the `F_ST` values.*
7. *Calculate `F_ST` using only the European and African samples and make a histogram of the values.*
8. *Make a plot with the genomic position on the x axis and the `F_ST` value on the y axis.*

We can also look at bins of a given size along the genome (to keep it simple we will just plot non-overlapping bins instead of sliding windows). We can fx. plot the fraction of sites in each bin that are polymorphic for each subpopulation:
```r
bin_width = 500000
d2 %>% 
  mutate(binmid=((as.integer(POS/bin_width)+0.5)*bin_width)) %>%
  filter(na!=0, nA!=0) %>%
  group_by(binmid, region) %>%
  summarise(frac_polymorph=n()/bin_width) %>%
  ggplot(aes(x=binmid, y=frac_polymorph, color=region)) + geom_line()
```

9. *Make a plot with average F_ST in bins along the genome.*


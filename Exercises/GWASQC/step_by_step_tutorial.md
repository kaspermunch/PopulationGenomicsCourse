# GWAS QC: step-by-step tutorial

First I create a new directory in my home folder in the cluster named GWAS_QC and copy the input files there:

```bash
mkdir GWAS_QC
cp /home/Data/GWAS_QC_data/* ./GWAS_QC/
```
I check that plink is installed by running plink like this:

```bash
plink
```
This will also show me the version of the program (1.90) and all the tools inside the program.

Works? Good, we have everything we need to start. The filtering has two parts: Sample filtering and SNP filtering.
I'll use the commands provided in the README file but I'll change the output prefix each time we filter out samples or variants, thus avoiding the dreadful "~" files, that complicated your lives earlier. That way we'll also be able to trace back any problematic step.

## Sample filtering

#### Wrong sex

I check for samples with wrong sex aannotation with the following:

```bash
ls -lth
plink --bfile GWA-data --check-sex --out GWA-data
ls -lth
head GWA-data.sexcheck
grep -v "OK" GWA-data.sexcheck > wrong_sex.txt
cat wrong_sex.txt
```

I included ls -lh steps before and after running plink --check-sex, so that I could check that the original files had the same size, i.e we it did not modify them so you shoudn't get "~" files. Plus, we can see the only outputs of the command. The -l option allows us to see all the information of the different files in the directory, the -h option makes the memory easier to read and the -t option orders the files by time, so we can see the ones added more recently first.
If you read the standard output (std out), i.e what was written on the screen when you ran the command, or the log file, you will see that there are 3 identified samples with wrong sex information. 
I inspected the contents and structure of the file "GWA-data.sexcheck" using head. Then, I used grep to susbset the lines that did not have an "OK", since I saw that these was the field to define good variants from the head command -longer than 3 lines. I could also have run something like:

```bash
awk '{if (NR > 1) print $5}' GWA-data.sexcheck|sort -k1,1 -u
```

This would have omitted the header and printed unique values for column 5, which is the one to store the status of the sample, i.e what we are looking for.

By doing cat "wrong_sex.txt" I can see the wrong samples and I'll filter them out with:

```bash
plink --bfile GWA-data --remove wrong_sex.txt --make-bed --out GWA-data_sflt
```

By changing the output prefix, you'll avoid again "~" files. You can see from std out or the log file that we have the same number of variants and 3 samples less. You can also inspect the new files created by running:

```bash
ls -lth
```

#### Missing data

To get an analysis of missing data and heterozygosity, I run the following:

```
plink --bfile GWA-data_sflt --missing --out GWA-data_sflt
plink --bfile GWA-data_sflt --het --out GWA-data_sflt
ls -lth
```

Note that the -bfile is now changed, since I modified the prefix in the last step and we want to keep using the filtered data. I keep the same out prefix because I'm not changing the files, but just generating new files from these, as you can see from the output of "ls -lht".

I'll go to Rstudio now and analyze these results:

```R
Imiss = read.table('GWA-data_sflt.imiss', header = TRUE)
mean(Imiss$N_MISS)
mean(Imiss$F_MISS)

het = read.table('GWA-data_sflt.het', header = T)
het$het_cal = (het$N.NM. - het$O.HOM.)/ het$N.NM.
mean_het = mean(het$het_cal)
sd_het = sd(het$het_cal)

plot(x = het$het_cal,y = Imiss$F_MISS, lwd =4, col = 'light blue', xlab = 'Heterozygosity rate per individual', ylab = 'Proportion of missing SNPs per individual')
abline(h = 0.03, col = 'red', lty = 3)
abline(v = mean_het + 3*sd_het, col = 'red', lty = 3)
abline(v = mean_het - 3*sd_het, col = 'red', lty = 3)
```

By running this I get the mean of number of missing genotypes and missing frequencies for each sample. I also compute the observed heterozygosity on each individual. Then, I plot these results, outlying extreme values on both variables, i.e missing frequency greater than 3% and heterozygosity lower or greaater than 3 standard deviations from the mean.
A high proportion of missing genotypes is indicative of low DNA quality or concentration, while higher or lower proportion of heterozygotes is indicative of contamination or inbreeding, respectively.

Thus, I proceed to filter samples based on these criteria, with:

```R
right_tail = mean_het + 3*sd(het$het_cal)
left_tail =  mean_het - 3*sd(het$het_cal)

filtering = cbind(Imiss, het)
outlier_ind = subset(filtering, filtering$F_MISS >= 0.03 | filtering$het_cal > right_tail | filtering$het_cal < left_tail)

nrow(outlier_ind)

write.table(outlier_ind[,c(1,2)], 'wrong_het_missing_values.txt', col.names = FALSE, row.names = FALSE)
```

I save the outlier 43 samples in a file named 'wrong_het_missing_values.txt'. I only need the IDs for filtering, so I just save the first two columns.

Back to the command line, I filter the dataset based on the file I've just created:

```bash
plink --bfile GWA-data_sflt --remove wrong_het_missing_values.txt --make-bed --out GWA-data_sflt_mflt_hflt
ls -lht
```

After that, you should have the same number of variants and 1954 samples left. You should also see that we have created bed, bim and fam files for GWA-data_sflt_mflt_hflt, containing our filtered dataset.

#### Relatedness

We want to filter out samples with high identity by descent (IBD), indicative of relatedness between individuals or repeaated samples, that could bias our GWAS.
Before that, we have to filter out variants in LD that could bias the estimates of IBD.

On the command line, I run:

```bash
plink --bfile GWA-data_sflt_mflt_hflt --indep-pairwise 500kb 5 0.2 --out GWA-data_sflt_mflt_hflt
plink --bfile GWA-data_sflt_mflt_hflt --extract GWA-data_sflt_mflt_hflt.prune.in --genome --min 0.185 --out GWA-data_sflt_mflt_hflt
ls -lht
```

From the first file you will get a set of variants that show r2 lower than 0.2, using windows of 500Kb and a step size of 5, named "GWA-data_sflt_mflt_hflt.prune.in". Using the second command, I compute IDB with these variants in linkage equilibrium and output the pairs of individuals with IBD greater than 0.185 (mean bewteen second degree relatives and third degree relatives). Again, this is not going to change our files so we don't need to change the prefix on the output files.
In R, I proceed in analyzing the output file, named "GWA-data_sflt_mflt_hflt.genome":

```R
ibd = read.table('GWA-data_sflt_mflt_hflt.genome', header = TRUE)
members = ibd[,1]
members = unique(members)
length(members)
write.table(cbind(members,members), file = 'wrong_ibd.txt', col.names = F, row.names = F)
```

With these lines of code, I obtain a unique set of 14 samples that show outlying IDB with another sample. Back to the command line, I filter out samples in our dataset using this file by running:

```bash
plink --bfile GWA-data_sflt_mflt_hflt --remove wrong_ibd.txt --make-bed --out GWA-data_sflt_mflt_hflt_iflt
ls -lth
```

I get a new set of files, containing the same number of variables and the samples that survived the wrong sex, the missing data and ibd filtering. At this point, you should have 1940 samples, 988 cases and 952 controls.

## Variant filtering

We need to re-compute the missing data for each variant, since we filtered some samples. I do that by running:

```bash
plink --bfile GWA-data_sflt_mflt_hflt_iflt --missing --out GWA-data_sflt_mflt_hflt_iflt
ls -lht
```

With "ls -lht" I can see again the new files produced. I open "GWA-data_sflt_mflt_hflt_iflt.lmiss" on R to plot an histogram of the variant's missing frequency:

```R
miss_data = read.table('GWA-data_sflt_mflt_hflt_iflt.lmiss', header = T)
hist(miss_data$F_MISS, breaks = 50, col = 'light blue', main = 'Missing Data Distribution', xlab = 'Missing Data rate')
```

A lot of care must be taken when filtering variants, since we can lose potential variants associaated with the phenotypic trait. Therefore, we only filter these variants in which the missigness is associated with the phenotype, i.e case or control, so they could be a source of bias in our study. I compute the Fisher exact test between these two variables by running:

```bash
plink --bfile GWA-data_sflt_mflt_hflt_iflt --test-missing --out GWA-data_sflt_mflt_hflt_iflt
```
This outputs a file named "GWA-data_sflt_mflt_hflt_iflt.missing". Using R and applying a significance criterion of 10e-5, I get significant variants by:

```R
test_missing = read.table('GWA-data_sflt_mflt_hflt_iflt.missing', header = TRUE)
fail_diffmiss_qc = test_missing[test_missing$P < 10e-5, 2]
write.table(fail_diffmiss_qc, file = 'fail-diffmiss-qc.txt', row.names = F, col.names = F)
```

I filter out these variants and also variants that show deviations from HWE and variants with low allele frequency, by running:

```bash
plink --bfile GWA-data_sflt_mflt_hflt_iflt --exclude fail-diffmiss-qc.txt --geno 0.05 --hwe 0.00001 --maf 0.01 --make-bed --out GWA-data_sflt_mflt_hflt_iflt_vflt
```

Deviations from HWE might indicate genotype error, but we choose a very low siginifcance value, since deviations from HWE might also be a product of selection implying that we could be filtering out variants potentially associated with the phenotypic trait. 
Samples with low allele frequency are also more likely to be miss-calls. Additionally, these variants have low power to be significant associations in the GWAS study due to their low small sample size.

After all filtering, you should be left with:
- 313896 variants.
- 1940 people (988 cases and 952 controls).








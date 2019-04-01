# Association testing using PLINK

### Software:
We will be using plink 1.9, you cna see the documentation [here](https://www.cog-genomics.org/plink/1.9/).
We will also be using R and Rstudio to make plots and make simple calculations.

### Data:
We will use a simulated date set of 47 cases and 41 controls. It is uploaded to Blackboard (Materials/week10).

### Exercise contents:
In this practical, we will go through the steps of performing association tests in plink and adjusting for principle components.

### Test for association with disease status using a Fisher’s exact test
To test for association between SNPs and disease status using an allelic Fisher’s exact test, type the following command at the shell prompt:
```
plink --bfile gwa --fisher --out gwa
```
*1) Take a look at the output file “gwa.assoc.fisher”. What is the p-value and location of the most significant variant?*

*2) Is the most significant variant significant if you do Bonferroni correction?*

### Make plots
We will use the R package “qqman” to make Manhattan plots and qq-plots. The package is available in CRAN so it can be installed using the “install.packages” commando. You can read in the association results and make a Manhattan plot by typing:
```
d <- read.table('gwa.assoc.fisher', head=T)
manhattan(d)
```

*3) Are there other variants close to the most significant variant that are also associated with the disease?*

To make a QQ-plot you should use the “qq” function:
```
qq(d$P)
```
*4) Is there a general inflation of the test statistic?*

### Genomic Control.
The inflation factor, 𝝺 , can be calculated as the median chi-square value divided by 0.456. Given a p-value, p, the corresponding Chi-square can be calculated as:
```
qchisq(p, df=1, lower.tail = F)
```
*5) What is the inflation factor?*

To do genomic control (to adjust for inflated test statistic) you divide the chi-square values with the inflation factor. To turn a chi-square value, q, into a p-value you use the “pchisq” function:
```
pchisq(q, df=1, lower.tail = F)
```
*6) What is the p-value of the most significant marker after genomic control?*

### PCA
It is best to perform the PCA on a LD-pruned set of SNPs:
```
plink --bfile gwa --indep-pairwise 500kb 5 0.2 --out gwa
```
To use the pruned set of SNPs to calculate the relationship matrix and calculate the first 20 principle components (PCs) type:
```
plink --bfile gwa --extract gwa.prune.in --pca 20 --out gwa
```
This calculates the eigenvalues and the eigenvectors, and stores them in two files (gwa.eigenval, gwa.eigenvec).

*7) Load gwa.eigenvec into R and make a plot with the first PC on the x-axis and the second PC on the y-axis. Does it look like there is population structure in the data? How many populations?*

The eigenvalues divided by the number of individuals should correspond approximately to the variance explained by the eigenvectors.
 
*8) How large a percentage of the variance does the first PC approximately explain?*

### Adjusting for PCs
We can use a logistic regression test to perform an association test while correcting for covariates. To include the first PC as a covariate type:
```
plink --bfile gwa --logistic --covar gwa.eigenvec --covar-number 1
```
The resulting file “plink.assoc.logistic” contains p-values for both the SNPs and the covariates. To get the p-values for the SNPs should look at the rows with the value “ADD” in the “TEST” column. (It is possible to include more PCs. To include the first x covariates you can write “--covar-number 1-x”.)

*9) Create Manhattan plot and QQ-plot for the new results. Does the QQ-plot look better?*

*10) What is the inflation factor now?*

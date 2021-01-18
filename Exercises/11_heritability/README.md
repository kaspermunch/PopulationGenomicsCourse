# Estimating heritability using GCTA

### Software:

In this exercise we will be using GCTA. You can see the documentation [here](
http://cnsgenomics.com/software/gcta/#Download). We will be also using plink 1.9, you can see the documentation [here](https://www.cog-genomics.org/plink/1.9/).

### Exercise contents:

We will estimate the amount of variance explained by the SNPs in a GWAS dataset. You can find the data here:

<!-- TODO: Update file names -->
```bash
/home/Data/GWAS_heritability
```
Copy the content of the directory to your home.

### Calculating the genetic relationship matrix

We will use plink to calculate the genetic relationship matrix (GRM) since it is faster than gcta. At the shell prompt, type:

<!-- TODO: Update file names -->
```
plink --make-grm-gz --bfile gwa --out gwa
```

 This will save the genetic relationship matrix in the zipped file gwa.grm.gz. Try to read it into R:

<!-- TODO: Update file names -->
```
d <- read.table(gzfile('gwa.grm.gz'))
```

*1) If you exclude the lines where an individual is compared to itself (column 1 is equal to column 2) what is the highest value in the GRM then?*

### Estimating variance components

We can use gcta to estimate how much of the variance in the phenotype in gwa.phen is explained by the SNPs:

<!-- TODO: Update file names -->
```
gcta64 --grm-gz gwa --pheno gwa.phen --reml --out test
```

*2) How much of the phenotypic variance (Vp) is explained by the genetic variance (V(G))?*

*3) Is this number larger or smaller than the narrow-sense heritability (h^2) of the trait?*

### Estimating variance components for groups of SNPs

The estimation of variance components can be used to answer questions about how much of the heritability is explained by different parts of the genome (for example different chromosomes or different functional annotations).

 *4) How much of the phenotypic variance can be explained by the genetic variants on chromosome 6? (You can use the “--chr” flag in plink to build a GRM only using variants from a particular chromosome)*

*5) Does chromosome 6 contribute more to the heritability than would be expected? How many of the genetic variants in the data set are located on chr 6? (you can use the genetic map in gwa.bim to see the location of the variants).*

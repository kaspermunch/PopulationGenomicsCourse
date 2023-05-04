
# Project description and requirements

The last part of this course focus on applying the knowledge acquired during these 10 weeks into real genetic datasets. We encourage you to work on the project in class Monday and Wednesday where one (or several of us) will be there to guide you and answer questions. 

The project is **mandatory** and needs to be handed in as a report. 

The requirements of the report are:

 * It should be at most 10 pages (including reference and figures)
 
 * It must to be divided into sections: abstract, introduction,
 results/discussion, conclusion and references. 
 
 * It must to be cohesive and coerent.
 
 * Source code must be provided (appended or linked/github repository). 
 
 * It must be in a PDF format.

You can choose between three different projects:

1. **Selective sweeps on chromosome 3**

2. **GWAS of eye color or height**

3. **Exploring non-African archaic segments**

Each project lists two relevant papers. The two papers, which go with the project you choose, are included in your curriculum for the oral exam.

## Deadline
The report must be handed in on the 27th of May at 9 AM.

## Submission

The project needs to be submitted through Brightspace. 
The name of your report must states your *name* and the chosen *project*:

MYNAME_archaic.pdf, MYNAME_xchromosome.pdf or MYNAME_GWAS.pdf

-----------------------------------------------------------------------------------------------

# Positive selection in the chromosome 3 region 3p21.31

The 3p21.31 region on the human chromosome 3 spans about five megabases where positive selection seems to act recurrently. Previously published papers suggest that genes in the region have been under selection on multiple occasions in both African humans, in the ancestors of humans and chimpanzees, and more generally across primates. Why strong selection so often affects this region and which genes this selection affects is not not really known. With your newly aquired skills, you can apply the the most advanced population genomic methods and produce an updated inference of selection in Africans. For this project you have phased genotypes for chr3:46000000-54000000 individuals from the following populations:

```
YRI 	Yoruba      Yoruba in Ibadan, Nigeria
LWK 	Luhya       Luhya in Webuye, Kenya
GWD 	Gambian     Gambian in Western Division, The Gambia 
MSL 	Mende       Mende in Sierra Leone
ESN 	Esan        Esan in Nigeria
```

Make yourself familar with the study populations. Where in Africa are they? How are they related?

## Investigate the following

A. Perform an Fst scan between sets of populations in a sliding window of 100 SNP positions, comparing at least five pairs of populations. Identify the Fst outlier regions in each case.

B. Use Relate on all the individuals and visualize trees (using Relate or tskit) to get an impression of the relationship between the populations. How does this relate to your Fst retults?

C. Use Relate on each population seperately to infer positive selection.

D. Run one or more additional methods for selecion inference. If possible this should be  another tree sequence based method such as CLUES. Compare the results to those obtained using Relate.

E. Identify genes potentially under selection and any known function of these genes. Consider what may drive recurrent selection in this region.

## Papers

[Patterns of Ancestry, Signatures of Natural Selection, and Genetic Association with Stature in Western African Pygmies](https://doi.org/10.1371/journal.pgen.1002641)

[https://academic.oup.com/mbe/article/33/12/3065/2450101](https://academic.oup.com/mbe/article/33/12/3065/2450101)

Perhaps: [An approximate full-likelihood method for inferring selection and allele frequency trajectories from DNA sequence data](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1008384)



## Data

Data for the project can be found in this folder on the cluster:

```
~/populationgenomics/project_data/chr3region/chr3_all_460_540_phased.vcf.gz
~/populationgenomics/project_data/chr3region/chr3_ESN_460_540_phased.vcf.gz
~/populationgenomics/project_data/chr3region/chr3_GWD_460_540_phased.vcf.gz
~/populationgenomics/project_data/chr3region/chr3_LWK_460_540_phased.vcf.gz
~/populationgenomics/project_data/chr3region/chr3_MSL_460_540_phased.vcf.gz
~/populationgenomics/project_data/chr3region/chr3_YRI_460_540_phased.vcf.gz
~/populationgenomics/project_data/chr3region/all_inds.txt
~/populationgenomics/project_data/chr3region/ESN_inds.txt
~/populationgenomics/project_data/chr3region/GWD_inds.txt
~/populationgenomics/project_data/chr3region/LWK_inds.txt
~/populationgenomics/project_data/chr3region/MSL_inds.txt
~/populationgenomics/project_data/chr3region/YRI_inds.txt
~/populationgenomics/project_data/chr3region/20140520.chr3.strict_mask.fasta.gz
~/populationgenomics/project_data/chr3region/human_ancestor_3.fa
```

The files are VCF files for all individuals and for each African population seperately. Each file has a corresponding file with the individuals included. The last two files are a mask file and an ancestor sequence file used by Relate. The files are named as in the Relate exercises.

<!-- 
## Population Genetics on X-chromosome 

-----------------------------------------------------------------------------------------------

In this project you will perform an extension of the analysis you did in the exercise on selective sweeps. You will identify regions affected by positive selection on the X chromosome by comparing patterns of genetic variation within and between human populations. 

### Investigate the following

A. Perform an Fst scan between sets of populations in a sliding window of 100 SNP positions, including at least the contrast between Africa and Europe, between Europe and East Asia, and between East Asia and Africa. Identify the 10 strongest Fst outlier regions in each case. Identify their genomic position and the genes covered by these Fst peaks. Discuss potential adaptive explanations. 

You can obtain the allele frequencies at each position from the haplotype data. Then, you can estimate FST by comparing the expected heterozygosity within population (Hs) and across populations (Ht). 

Fst = (Ht - Hs)/Ht. 

B. Perform an iHS scan of the whole X chromosome for at least three populations. Identify the 10 most significant regions and associated with genes as in A.

C. Perform an XP-EHH scan of the whole X chromosome for at least three populations. Identify the 10 most significant regions and associated with genes as in A. 

D. Intersect the analysis of Fst and XP-EHH.

E. Perform any additional analysis of your own choice, such as (diversity along the X chromosome).

### Papers:

[Extreme selective sweeps displaced archaic admixture across the human X chromosome around 50,000 years ago](https://www.biorxiv.org/content/10.1101/503995v2)

### Data

Data for the project can be found in this folder on the cluster:

```
~/populationgenomics/project_data/Xchromosome
```

The files are:


The dataset is composed of 150 male individuals from the Simons Diversity Project. African individuals with high probability of a non-african component (Masai, Somali, Mozabite and Saharawi) have been excluded. We will use a total of 411892 SNPs with no missing data.

Files description:

snps_filtered.txt: Text file containing id, chromosome, position, derived allele and ancestral allele.

genotypes_* : File containing genotypes for each individual and position. Rows correspond to SNP positions and columns to individuals.
WE = WestEurasia
AF = Africa
EA = EastAsia
SA = South Asia
Am = America
CAS = CentralAsiaSiberia
O = Oceania

metadata.txt: Text file containing the population and region of all individuals used.

gencode.v30lift37.annotation.gtf.gz: Gtf file containing the gene annotation for human (assembly version Hg19; the same used in Simons Diversity Project). 

----------------------------------------------------------------------------------------------- 

-->

-----------------------------------------------------------------------------------------------

# GWAS of eye color or height

In this project you will be looking at GWAS data from [openSNP](https://opensnp.org), which is a web site where users of direct-to-customer genetic tests can share their personal data with other users. The phenotypes we will be looking at is self-reported eye color and height. 
When looking at the data you should be aware that:
- The data comes from different companies that use different chips so there are markers that are missing from some individuals because they were not present on the chip used by their company.
- The gender information is missing from the file and by default plink will ignore the phenotype of individuals without gender information. So you have to use “--allow-no-sex” option in plink.

## Investigate the following

A. Are there any closely related individuals in the sample?

B. Do a PCA plot. What does it tell you about the samples?

C. The files eye_color.txt and height.txt contains the self-reported eye color and height for the individuals in the data. Do a GWAS on one these traits. There are 12 eye color categories and you can group some of them together to create a binary phenotype. How many significant loci do you find? 

D. Do additional analyses using plink, GCTA, R or any other tool you might find relevant. The list below contains some suggestions for possible analyses, but you can also come up with your ideas

Suggestions for further analyses:
- Use mixed model for GWAS.	
- Do imputation (either of the whole genome or the region around the most significant SNP) and see if you can then find variants with lower p-values.
- If you use half of the data set to calculate a polygenic score, how well does that score predict height on the other half?
- Find a trained height PRS on the internet. How well does it predict the height in this data set?
- Test for epistasis.
- What are the distribution of phenotypes for each of the genotypes at th most significant SNP? If you want to analyse it in R you can use the "--recode A" together with the "--snp" and "--window" option in plink to get the variants around a specific SNP written to a text file that it is easy to load in R. 
- How many of the significant variants found in the largest published GWAS study can you replicate those hits in this data set?
- Make association tests where you condition on most significant variant (you can use the --condition option in plink)

## Papers

* [Genome-wide association study in almost 195,000 individuals identifies 50 previously unidentified genetic loci for eye color](https://advances.sciencemag.org/content/7/11/eabd1239)
* [A saturated map of common genetic variants associated with human height | Nature](https://www.nature.com/articles/s41586-022-05275-y)

## Data

Data for the project can be found in this folder on the cluster:

```
~/populationgenomics/project_data/GWAS
```


-----------------------------------------------------------------------------------------------

# Exploring non-African archaic segments

In this project you will be looking at segments of archaic genomes idendentified individual modern humans. You will investigate how alike  non-African genomes are in terms of how archaic segments are distributed. You will be working with an extended version of the data set that you worked with in the admixture exercise. In this version you also have an addtional file with the positions of candidate archaic SNPs.

## Investigate the following

In your project you must address the questions below, but you are also expected to expand the project to answer your own questions. How you do this is up to you. You do not need to answer them in the order they are listed. Make a project plan with a set of analyses that  will allow you to answer the questions. 

A. To what extent do individuals share SNPs contributed by archaic human introgression? In other words, how correlated are the archaic contents in two individuals?

B. How does this correlation chagne when you compare individuals from different populations?

C. How does it chagne when comparing individuals from different geographical regions?

D. Does the region containing the EPAS1 gene stand out in any way? (Redo the analysis above for a 1Mb window surrounding this gene).

E. What is the total amount of admixture (archiac genomic sequence) in each non-African individual genome?

F. What is the total amount of admixture (archiac genomic sequence) in the region around EPAS1 in each individual?

G. Do individuals with large admixture totals have more correlated admixture patterns? Do individuals with large admixture totals *in the EPAS1 region* have more correlated admixture patterns *in the EPAS1 region*? Can you find any evidence of adaptive introgression?

H. Perform any additional analyses of your own choice.

## Papers

- [Analysis of Human Sequence Data Reveals Two
Pulses of Archaic Denisovan Admixture](https://doi.org/10.1016/j.cell.2018.02.031)
- [Altitude adaptation in Tibetans caused by introgression of Denisovan-like DNA](https://doi.org/10.1038/nature13408)

## Data

Data for the project can be found in this folder on the cluster:

```
~/populationgenomics/project_data/ArchaicAdmixture
```

The files are:

- ArchaicSegments.txt: This file is formatted the same was as the one you used for the archiac admixture exercise.
- SNP.txt: This file has all the non-African SNPs that remain after removing all SNPs found in (Subsaharan) Africa. When each SNP is found in any of the high coverage archaic genomes (Altai or Vindija Neanderthals or the Denisova) it is labelled as such. Otherwise it is labeled "human".


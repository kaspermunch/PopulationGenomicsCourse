
# Project description and requirements

The last part of this course focus on applying the knowledge acquired during these 10 weeks into real genetic datasets. We encourage you to work on the project in class Monday and Thursday where one (or several of us) will be there to guide you and answer questions. 

The project is **mandatory** and needs to be handed in as a report. 

The requirements of the report are:

 * It should be at most 10 pages (including references and figures)
 
 * It must be divided into sections: abstract, introduction,
 results/discussion, conclusion, and references. 
 
 * It must be cohesive and coherent.
 
 * Source code must be provided (appended or linked/GitHub repository). 
 
 * It must be in a PDF format.

You can choose between three different projects:

1. **Selective sweeps on chromosome 3**

2. **GWAS of eye color or height**

3. **Exploring non-African archaic segments**

Each project lists two relevant papers. The two papers that go with the project you choose are included in your curriculum for the oral exam.

## Deadline
To be eligible for the exam, the report must be handed in before **midnight on June 5.**.

## Submission

The project needs to be submitted through Brightspace. 
The name of your report must state your *name* and the chosen *project*:

MYNAME_archaic.pdf, MYNAME_xchromosome.pdf or MYNAME_GWAS.pdf

-----------------------------------------------------------------------------------------------

# Positive selection in the chromosome 3 region 3p21.31

The 3p21.31 region on the human chromosome 3 spans about five megabases where positive selection seems to act recurrently. Previously published papers suggest that genes in the region have been under selection on multiple occasions in both African humans, in the ancestors of humans and chimpanzees, and more generally across primates. Why strong selection so often affects this region and which genes this selection affects is not really known. With your newly acquired skills, you can apply the most advanced population genomic methods and produce an updated inference of selection in Africans. For this project, you have phased genotypes for chr3:46000000-54000000 individuals from the following populations:

```
YRI 	Yoruba      Yoruba in Ibadan, Nigeria
LWK 	Luhya       Luhya in Webuye, Kenya
GWD 	Gambian     Gambian in Western Division, The Gambia 
MSL 	Mende       Mende in Sierra Leone
ESN 	Esan        Esan in Nigeria
```

Make yourself familiar with the study populations. Where in Africa are they? How are they related?

## Investigate the following

A. Perform an Fst scan between sets of populations in a sliding window of 100 SNP positions, comparing at least five pairs of populations. 
<!-- You can also use [tskit](https://tskit.dev/tskit/docs/stable/introduction.html) to compute Fst from the tree sequences you produce with Relate.  -->
Identify the Fst outlier regions in each case.

B. Use Relate on all the individuals and visualize trees (using Relate or [tskit](https://tskit.dev/tskit/docs/stable/introduction.html)) to get an impression of the relationship between the populations. How does this relate to your Fst results?

C. Run Relate on each population separately to infer positive selection. You can run Relate the 46-54Mb region only but also on the entire chr3 (see below). Even if you only investigate the 46-54Mb region you should consider the latter option, as it will produce a more reliable estimation of population demography, which is used in the selection inference. 

D. Do one or more additional analyses to address selection and compare the results to those obtained using Relate. If possible, these should also based on tree sequences and might involve:

- Alternative statistics computed based on the trees already obtained from Relate. [tskit](https://tskit.dev/tskit/docs/stable/introduction.html) has ready-made ones, or you can devise one yourself.
- Use another method of inference, such as Clues.
- or do something else...

E. Identify genes potentially under selection and any known function of these genes. Consider what may drive recurrent selection in this region.

## Papers

[Patterns of Ancestry, Signatures of Natural Selection, and Genetic Association with Stature in Western African Pygmies](https://doi.org/10.1371/journal.pgen.1002641)

[Selective Sweeps across Twenty Millions Years of Primate Evolution](https://academic.oup.com/mbe/article/33/12/3065/2450101)

Perhaps: [An approximate full-likelihood method for inferring selection and allele frequency trajectories from DNA sequence data](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1008384)



## Data

Data for the project can be found in this folder on the cluster:

```
~/populationgenomics/project_data/chr3region
```

<!-- 
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
-->

In your project directory, run these commands to create links to the main data files:

```
ln -s ~/populationgenomics/project_data/chr3region/1000g_chr3_AFR.vcf
ln -s ~/populationgenomics/project_data/chr3region/1000g_chr3_46_54_AFR.vcf
ln -s ~/populationgenomics/project_data/chr3region/1000g_chr3_46_54_YRI.vcf
ln -s ~/populationgenomics/project_data/chr3region/1000g_chr3_46_54_LWK.vcf
ln -s ~/populationgenomics/project_data/chr3region/1000g_chr3_46_54_GWD.vcf
ln -s ~/populationgenomics/project_data/chr3region/1000g_chr3_46_54_MSL.vcf
ln -s ~/populationgenomics/project_data/chr3region/1000g_chr3_46_54_ESN.vcf
```

The first one is a VCF with all Africans for the entire chr3. The second one has all Africans but only for the region 46Mb-54Mb. The last four are also only the region 46Mb-54Mb but these contain only data for individuals from each of the five African populations.

## Running Relate

Run these commands to create links to the additional files you need for the Relate analyses:

```
ln -s ~/populationgenomics/project_data/chr3region/AFR.poplabels
ln -s ~/populationgenomics/project_data/chr3region/YRI.poplabels
ln -s ~/populationgenomics/project_data/chr3region/LWK.poplabels
ln -s ~/populationgenomics/project_data/chr3region/GWD.poplabels
ln -s ~/populationgenomics/project_data/chr3region/MSL.poplabels
ln -s ~/populationgenomics/project_data/chr3region/ESN.poplabels
ln -s ~/populationgenomics/project_data/chr3region/all_except_YRI.txt
ln -s ~/populationgenomics/project_data/chr3region/all_except_LWK.txt
ln -s ~/populationgenomics/project_data/chr3region/all_except_GWD.txt
ln -s ~/populationgenomics/project_data/chr3region/all_except_MSL.txt
ln -s ~/populationgenomics/project_data/chr3region/all_except_ESN.txt
ln -s ~/populationgenomics/project_data/chr3region/human_ancestor_3.fa
ln -s ~/populationgenomics/project_data/chr3region/20140520.chr3.strict_mask.fasta.gz
ln -s ~/populationgenomics/project_data/chr3region/genetic_map_chr3_combined_b37.txt
```

The first ten are meta data. The last three are human ancestor, quality mask, and recombination map.


### Relate on region 46-54Mb only

Run this command to create files in the Relate input file format for all the African individuals:

```
~/populationgenomics/software/relate/bin/RelateFileFormats --mode ConvertFromVcf --haps 1000g_chr3_46_54_AFR.haps --sample 1000g_chr3_46_54_AFR.sample -i 1000g_chr3_46_54_AFR --poplabels AFR.poplabels
```

To run the remaining steps on all Africans you do:

```
~/populationgenomics/software/relate/scripts/PrepareInputFiles/PrepareInputFiles.sh --haps 1000g_chr3_46_54_AFR.haps --sample 1000g_chr3_46_54_AFR.sample --ancestor human_ancestor_3.fa --mask 20140520.chr3.strict_mask.fasta.gz -o 1000g_chr3_46_54_AFR_ALL --poplabels AFR.poplabels 
~/populationgenomics/software/relate/bin/Relate --mode All -m 1.25e-8 -N 20000 --sample 1000g_chr3_46_54_AFR_ALL.sample.gz --haps 1000g_chr3_46_54_AFR_ALL.haps.gz --map genetic_map_chr3_combined_b37.txt --annot 1000g_chr3_46_54_AFR_ALL.annot --dist 1000g_chr3_46_54_AFR_ALL.dist.gz --memory 20 -o 1000g_chr3_46_54_AFR_ALL
~/populationgenomics/software/relate/scripts/EstimatePopulationSize/EstimatePopulationSize.sh -m 1.25e-8 -N 20000 -i 1000g_chr3_46_54_AFR_ALL --poplabels AFR.poplabels -o 1000g_chr3_46_54_AFR_ALL_popsize --threshold 0 --num_iter 5 --years_per_gen 29 --threads 14 --threshhold 0
~/populationgenomics/software/relate/scripts/DetectSelection/DetectSelection.sh -i 1000g_chr3_46_54_AFR_ALL -m 1.25e-8 --poplabels AFR.poplabels -o 1000g_chr3_46_54_AFR_ALL_selection
~/populationgenomics/software/relate/bin/RelateFileFormats --mode ConvertToTreeSequence -i 1000g_chr3_46_54_AFR_ALL -o 1000g_chr3_46_54_AFR_ALL
```

To run separate analyses for each of the populations populations: Yoruba in Ibadan, Nigeria (YRI), Luhya in Webuye, Kenya (LWK), Gambian in Western Division – Mandinka (GWD), Mende in Sierra Leone (MSL), and Esan in Nigeria (ESN). The commands below run Relate on the individuals from the Luhya population (Notice the LWK-part of file names):

```
~/populationgenomics/software/relate/scripts/PrepareInputFiles/PrepareInputFiles.sh --haps 1000g_chr3_46_54_AFR.haps --sample 1000g_chr3_46_54_AFR.sample --ancestor human_ancestor_3.fa --mask 20140520.chr3.strict_mask.fasta.gz --remove_ids all_except_LWK.txt -o 1000g_chr3_46_54_AFR_LWK --poplabels LWK.poplabels
~/populationgenomics/software/relate/bin/Relate --mode All -m 1.25e-8 -N 20000 --sample 1000g_chr3_46_54_AFR_LWK.sample.gz --haps 1000g_chr3_46_54_AFR_LWK.haps.gz --map genetic_map_chr3_combined_b37.txt --annot 1000g_chr3_46_54_AFR_LWK.annot --dist 1000g_chr3_46_54_AFR_LWK.dist.gz --memory 20 -o 1000g_chr3_46_54_AFR_LWK
~/populationgenomics/software/relate/scripts/EstimatePopulationSize/EstimatePopulationSize.sh -m 1.25e-8 -N 20000 -i 1000g_chr3_46_54_AFR_LWK --poplabels LWK.poplabels -o 1000g_chr3_46_54_AFR_LWK_popsize --threshold 0 --num_iter 5 --years_per_gen 29 --threads 14 --threshhold 0
~/populationgenomics/software/relate/scripts/DetectSelection/DetectSelection.sh -i 1000g_chr3_46_54_AFR_LWK -m 1.25e-8 --poplabels LWK.poplabels -o 1000g_chr3_46_54_AFR_LWK_selection
~/populationgenomics/software/relate/bin/RelateFileFormats --mode ConvertToTreeSequence -i 1000g_chr3_46_54_AFR_LWK -o 1000g_chr3_46_54_AFR_LWK
```

### Relate on entire chr3

Run this command to create files in the Relate input file format for all the African individuals:

```
~/populationgenomics/software/relate/bin/RelateFileFormats --mode ConvertFromVcf --haps 1000g_chr3_AFR.haps --sample 1000g_chr3_AFR.sample -i 1000g_chr3_AFR --poplabels AFR.poplabels
```

To run the remaining steps on all Africans (will take a *long* time) you do:

```
~/populationgenomics/software/relate/scripts/PrepareInputFiles/PrepareInputFiles.sh --haps 1000g_chr3_AFR.haps --sample 1000g_chr3_AFR.sample --ancestor human_ancestor_3.fa --mask 20140520.chr3.strict_mask.fasta.gz -o 1000g_chr3_AFR_ALL --poplabels AFR.poplabels 
~/populationgenomics/software/relate/bin/Relate --mode All -m 1.25e-8 -N 20000 --sample 1000g_chr3_AFR_ALL.sample.gz --haps 1000g_chr3_AFR_ALL.haps.gz --map genetic_map_chr3_combined_b37.txt --annot 1000g_chr3_AFR_ALL.annot --dist 1000g_chr3_AFR_ALL.dist.gz --memory 20 -o 1000g_chr3_AFR_ALL
~/populationgenomics/software/relate/scripts/EstimatePopulationSize/EstimatePopulationSize.sh -m 1.25e-8 -N 20000 -i 1000g_chr3_AFR_ALL --poplabels AFR.poplabels -o 1000g_chr3_AFR_ALL_popsize --threshold 0 --num_iter 5 --years_per_gen 29 --threads 14 --threshhold 0
~/populationgenomics/software/relate/scripts/DetectSelection/DetectSelection.sh -i 1000g_chr3_AFR_ALL -m 1.25e-8 --poplabels AFR.poplabels -o 1000g_chr3_AFR_ALL_selection
~/populationgenomics/software/relate/bin/RelateFileFormats --mode ConvertToTreeSequence -i 1000g_chr3_AFR_ALL -o 1000g_chr3_AFR_ALL
```

To run separate analyses for each of the populations populations: Yoruba in Ibadan, Nigeria (YRI), Luhya in Webuye, Kenya (LWK), Gambian in Western Division – Mandinka (GWD), Mende in Sierra Leone (MSL), and Esan in Nigeria (ESN). The commands below run Relate on the individuals from the Luhya population (Notice the LWK-part of file names):

```
~/populationgenomics/software/relate/scripts/PrepareInputFiles/PrepareInputFiles.sh --haps 1000g_chr3_AFR.haps --sample 1000g_chr3_AFR.sample --ancestor human_ancestor_3.fa --mask 20140520.chr3.strict_mask.fasta.gz --remove_ids all_except_LWK.txt -o 1000g_chr3_AFR_LWK
~/populationgenomics/software/relate/bin/Relate --mode All -m 1.25e-8 -N 20000 --sample 1000g_chr3_AFR_LWK.sample.gz --haps 1000g_chr3_AFR_LWK.haps.gz --map genetic_map_chr3_combined_b37.txt --annot 1000g_chr3_AFR_LWK.annot --dist 1000g_chr3_AFR_LWK.dist.gz --memory 20 -o 1000g_chr3_AFR_LWK
~/populationgenomics/software/relate/scripts/EstimatePopulationSize/EstimatePopulationSize.sh -m 1.25e-8 -N 20000 -i 1000g_chr3_AFR_LWK --poplabels LWK.poplabels -o 1000g_chr3_AFR_LWK_popsize --threshold 0 --num_iter 5 --years_per_gen 29 --threads 14 --threshhold 0
~/populationgenomics/software/relate/scripts/DetectSelection/DetectSelection.sh -i 1000g_chr3_AFR_LWK -m 1.25e-8 --poplabels LWK.poplabels -o 1000g_chr3_AFR_LWK_selection
~/populationgenomics/software/relate/bin/RelateFileFormats --mode ConvertToTreeSequence -i 1000g_chr3_AFR_LWK -o 1000g_chr3_AFR_LWK
```








<!-- 
## Population Genetics on X-chromosome 

-----------------------------------------------------------------------------------------------

In this project you will perform an extension of the analysis you did in the exercise on selective sweeps. You will identify regions affected by positive selection on the X chromosome by comparing patterns of genetic variation within and between human populations. 

### Investigate the following

A. Perform an Fst scan between sets of populations in a sliding window of 100 SNP positions, including at least the contrast between Africa and Europe, between Europe and East Asia, and between East Asia and Africa. Identify the 10 strongest Fst outlier regions in each case. Identify their genomic position and the genes covered by these Fst peaks. Discuss potential adaptive explanations. 

You can obtain the allele frequencies at each position from the haplotype data. Then, you can estimate FST by comparing the expected heterozygosity within populations (Hs) and across populations (Ht). 

Fst = (Ht - Hs)/Ht. 

B. Perform an iHS scan of the whole X chromosome for at least three populations. Identify the 10 most significant regions associated with genes as in A.

C. Perform an XP-EHH scan of the whole X chromosome for at least three populations. Identify the 10 most significant regions associated with genes as in A. 

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


The dataset is composed of 150 male individuals from the Simons Diversity Project. African individuals with a high probability of a non-African component (Masai, Somali, Mozabite, and Saharawi) have been excluded. We will use a total of 411892 SNPs with no missing data.

File descriptions:

snps_filtered.txt: A text file containing id, chromosome, position, derived allele, and ancestral allele.

genotypes_* : File containing genotypes for each individual and position. Rows correspond to SNP positions and columns to individuals.
WE = WestEurasia
AF = Africa
EA = EastAsia
SA = South Asia
Am = America
CAS = CentralAsiaSiberia
O = Oceania

metadata.txt: A text file containing the population and region of all individuals used.

gencode.v30lift37.annotation.gtf.gz: Gtf file containing the gene annotation for humans (assembly version Hg19; the same used in the Simons Diversity Project). 

----------------------------------------------------------------------------------------------- 

-->

-----------------------------------------------------------------------------------------------

# GWAS of eye color or height

In this project, you will be looking at GWAS data from [openSNP](https://opensnp.org), which is a website where users of direct-to-customer genetic tests can share their personal data with other users. The phenotypes we will be looking at are self-reported eye color and height. 
When looking at the data, you should be aware that:
- The data comes from different companies that use different chips so there are markers that are missing from some individuals because they were not present on the chip used by their company.
- The gender information is missing from the file and by default plink will ignore the phenotype of individuals without gender information. So you have to use “--allow-no-sex” option in plink.

## Investigate the following

A. Do QC. Are there any closely related individuals in the sample?

B. Do a PCA plot. What does it tell you about the samples?

C. The files eye_color.txt and height.txt contain the self-reported eye color and height of the individuals in the data. Do a GWAS on one of these traits. There are 12 eye color categories, and you can group some of them together to create a binary phenotype. How many significant loci do you find? 

D. Do additional analyses using Plink, GCTA, R, or any other tool you might find relevant. The list below contains some suggestions for possible analyses, but you can also come up with your ideas

Suggestions for further analyses:
- Use mixed model for GWAS.	
- Do imputation (either of the whole genome or the region around the most significant SNP) and see if you can then find variants with lower p-values.
- If you use half of the data set to calculate a polygenic score, how well does that score predict height on the other half?
- Find a trained height PRS on the internet. How well does it predict the height in this data set?
- Test for epistasis.
- What are the distribution of phenotypes for each of the genotypes at th most significant SNP? If you want to analyse it in R you can use the "--recode A" together with the "--snp" and "--window" option in plink to get the variants around a specific SNP written to a text file that it is easy to load in R. 
- How many of the significant variants found in the largest published GWAS study can you replicate those hits in this data set?
- Make association tests where you condition on the most significant variant (you can use the --condition option in plink)

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

In this project, you will be looking at segments of archaic genomes identified in individual modern humans. You will investigate how alike non-African genomes are in terms of how archaic segments are distributed. You will be working with an extended version of the data set that you worked with in the admixture exercise. In this version you also have an addtional file with the positions of candidate archaic SNPs.

## Investigate the following

In your project, you must address the questions below, but you are also expected to expand the project to answer your own questions. How you do this is up to you. You do not need to answer them in the order they are listed. Make a project plan with a set of analyses that  will allow you to answer the questions. 

A. To what extent do individuals share SNPs contributed by archaic human introgression? In other words, how correlated are the archaic contents in two individuals?

B. How does this correlation change when you compare individuals from different populations?

C. How does it change when comparing individuals from different geographical regions?

D. Does the region containing the EPAS1 gene stand out in any way? (Redo the analysis above for a 1Mb window surrounding this gene).

E. What is the total amount of admixture (archaic genomic sequence) in each non-African individual genome?

F. What is the total amount of admixture (archaic genomic sequence) in the region around EPAS1 in each individual?

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

- ArchaicSegments.txt: This file is formatted the same way as the one you used for the archaic admixture exercise.
- SNP.txt: This file has all the non-African SNPs that remain after removing all SNPs found in (Subsaharan) Africa. When each SNP is found in any of the high-coverage archaic genomes (Altai or Vindija Neanderthals or the Denisova), it is labeled as such. Otherwise, it is labeled "human".


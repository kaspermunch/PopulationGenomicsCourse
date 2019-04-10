
### Project description and requirements

The last part of this course focus on applying the knowledge acquired during these 10 weeks into real genetic datasets. The project is **mandatory** and needs to be handed-in as a report. The requirements of the report are:

 * Report of at maximum 10 pages (including reference and figures)
 
 * It needs to be divided by sections: abstract, introduction,
 results/discussion, conclusion and references. 
 
 * It needs to be cohesive and coerent.
 
 * Source code must be provided (appended or linked/github repository).
 
 * It must be in a pdf format.
 
Three choices are given and you must choose one and produce the required:

1. **Population Genetics on X-chromosome**

2. **Analysis of GWAS summary statistics**

3. **Exploring non-African archaic segments**

#### Deadline
The report must be handed in on the 22th of May at 9 AM.

#### Submission

The project needs to be submitted through Blackboard. 
The name of your report must states your *name* and the chosen *project*:

MICA_archaic.pdf, MICA_xchromosome.pdf or MICA_GWAS.pdf

#### Q&A session

We encourage you to help each other during the process. Feel free to approach the TA when necessary. We will also give you the opportunity to ask questions in a QA session that will happen on the 9th of May. 

-----------------------------------------------------------------------------------------------

## Population Genetics on X-chromosome 

-----------------------------------------------------------------------------------------------

The data consists of a vcf file of 150 male full X chromosomes, a bed file with callable regions, a gif gene annotation file, a metafile with information about the samples and a set of files for use with REHH.

Gene annotation:

Gene annotation (gtf format) for Hg19 can be found in the following website
https://www.gencodegenes.org/releases/17.html
It was also uploaded to the dropbox as **gencode.v17.annotation.gtf**

Fst Calculation:

You will do the analysis from scratch by reading the genotype file of each population into different tables (rememeber rows are SNP positions and columns are individuals), the information about the snps are in the .snp file (ancestral and derived alleles). The data is haploid (n) therefore calculating Fst consists of estimating the allele frequencies for each position and calculating the expected heterozygosity within population Hs and contrasting Expected Heterozygosity across populations Ht. 

Fst = (Ht - Hs)/Ht. 

This can be done by averaging Fst values for a set of consecutive markers in a given window size (100 SNPs).


### Investigate the following

A. Perform an Fst scan between sets of populations in a sliding window of 100 SNP positions, including at least the contrast between Africa and Europe, between Europe and East Asia, and between East Asia and Africa. Identify the 10 strongest Fst outlier regions in each case. Identify their genomic position and the genes covered by these Fst peaks. Discuss potential adaptive explanations. 

B. Perform an iHS scan of the whole X chromosome for at least three populations. Identify the 10 most significant regions and associated with genes as in A.

C. Perform an XP-EHH scan of the whole X chromosome for at least three populations. Identify the 10 most significant regions and associated with genes as in A. 

D. Intersect the analysis of Fst and XP-EHH

E. Perform any additional analysis of your own choice, such as (diversity along the C X chromosome)

### Data

A DropBox link uploaded on BlackBoard (./Materials/Week 12: Projects Materials)

-----------------------------------------------------------------------------------------------

## Analysis of GWAS summary statistics

-----------------------------------------------------------------------------------------------

### Summary statistics

To protect the privacy of the participants in GWA studies the raw genotype data is usually not made public. It has, however, become the norm that studies will publish the summary statistics (p-value, effect size, frequency etc.) for all the variants they have tested. For this reason there are many new methods designed to do further analyses based on summary statistics. In this project you will analyse summary statistcs from a large meta-analysis of BMI and:
- Do a per gene test and look for enriched gene sets.
- Estimate the heritability of BMI.

### Data
You can download summary statistics from a BMI study with ~700,000 individuals [here](http://portals.broadinstitute.org/collaboration/giant/images/6/63/Meta-analysis_Wood_et_al%2BUKBiobank_2018.txt.gz).

### Per gene test.
The gcta tool that you used in the last GWAS exercise can calculate a per gene test based on summary statistics (see [here](http://gcta.freeforums.net/thread/309/gcta-fastbat-based-association-analysis)). In order to see if the presence of multiple significant variants in the same gene is due to LD or multiple independent signals it is necessary to provide a data set in plink-format that the program can use to estimate the LD between variants. For this purpose you can use data from the 1000 genomes project that can be downloaded in plink format [here](https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_plinkfiles.tgz). The BMI study and LD data is based on hg19 so you should use the gene list file called "glist-hg19.txt". The genotype data consists of separate files for each chromosome so you can do a test per chromosome and then concatenate the results in the end.

#### Gene set enrichment
Try to see if the genes affecting BMI are enriched in specific Gene Ontologies or Pathways or Tissue types. You can for example use the tool enrichR: http://amp.pharm.mssm.edu/Enrichr/

### Estimating heritability
The method called LD-score regression can be used to estimate heritability using summary statistics. The method is described in [this article](https://www.nature.com/articles/ng.3211). A description of how to use the software to estimate heritability can be found [here](https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation). The software can be downloaded [here](https://github.com/bulik/ldsc).

-----------------------------------------------------------------------------------------------


## Exploring non-African archaic segments

-----------------------------------------------------------------------------------------------

### How alike are non-African genomes in terms of how archaic segments are distributed?

In this project you will be working with an extended version of the data set that you worked with in the admixture exercise. In this version you also have an addtional file with the positions of candidate archaic SNPs.

In your project you must address the questions below, but you are also expected to expand the project to answer your own questions. How you do this is up to you. You do not need to answer them in the order they are listed. Make a project plan with a set of analyses that  will allow you to answer the questions. 

A. To what extent do individuals share SNPs contributed by archaic human introgression? In other words, how correlated are the archaic contents in two individuals?

B. Is this correlation stronger when you compare individuals from different populations?

C. Is it even stronger when comparing individuals from different geographical regions?

D. Does the region containing the EPAS1 gene stand out in any way? (Redo the analysis above for a 1Mb window surrounding this gene).

E. What is the total amount of admixture in each non-African individual?

F. What is the total amount of admixture in the region around EPAS1 in each individual?

G Do individuals with large admixture totals have more correlated admixture patterns? Do individuals with large admixture totals *in the EPAS1 region* have more correlated admixture patterns *in the EPAS1 region*? Can you find any evidence of adaptive introgression?

### Data

[Google Drive folder with data](https://drive.google.com/open?id=1lrRfFcoxpyVpXgOi4RYP2_vauM_-rRz_)

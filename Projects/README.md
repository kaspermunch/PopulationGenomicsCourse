
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

2. **GWAS of eye color**

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

In this project you will perform an extension of the analysis done during the exercise session on selective sweeps. Thereby, you will try to infer signs of positive selection on the X chromosome, by comparing intraspecific variation among populations. 

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

[Google drive folder with data](https://drive.google.com/open?id=16TKX5WJ0vPDVttb5bHAwcVla2o4ARlay)

The dataset is composed of 150 male individuals from the Simons Diversity Project. African individuals with high probability of a non-african component (Masai, Somali, Mozabite and Saharawi) have been excluded. We will use a total of 411892 SNPs with no missing data.

Files description:

snps_filtered.txt: Text file containing id, chromosome, position, ancestral allele and derived allele.

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

## GWAS of eye color

-----------------------------------------------------------------------------------------------

In this project you will be looking at GWAS data from [openSNP](https://opensnp.org), which is a web site where users of direct-to-customer genetic tests can share their personal data with other users. The phenotype we will be looking at is self-reported eye color. 
When looking at the data you should have in mind that:
- The data comes from different companies that use different chips so there are markers that are missing from some individuals because they were not present on the chip used by their company.
- Many users have not reported their gender. By default plink will ignore the phenotype of these individuals unless you tell it not to.

### Investigate the following

A. Are there any closely related individuals in the sample?

B. Do a PCA plot. What does it tell you about the samples?

C. The file eye_color.txt contains the self-reported eye colors for the individuals in the data set. First try to do a GWAS for Brown vs. Blue eyes. How many loci do you find? 

D. Try to look at the SNPs at the most significant locus. If you want to analyse it in R you can use the "--recode A" together wither the "--snp" and "--window" option in plink to get a the variants around a specific SNP written to a text file that it is easy to load in R. How is the distribution of eye colors for each genotype of the most significant SNP? Is the effect additive, dominant or recessive?

E. Are there any significant variants after you condition on the most significant hit? (you can use the "--condition" option in plink)

F. Can you find genome-wide significant hits that helps distinguish between other colors (fx. green or hazel)? What if you restrict the analysis to a region around the most significant SNP in the blue vs brown comparison?

G. Do any additional analysis. Using plink, GCTA or any other tool you might find relevant.

### Papers:

[A GWAS in Latin Americans highlights the convergent evolution of lighter skin pigmentation in Eurasia | Nature Communications](https://www.nature.com/articles/s41467-018-08147-0)

### Data:

[Google drive folder with data](https://drive.google.com/open?id=1gOzJAh-2lJsZSMqN92M4xLeTCMCGolVn)

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

### Papers

- [Analysis of Human Sequence Data Reveals Two
Pulses of Archaic Denisovan Admixture](https://doi.org/10.1016/j.cell.2018.02.031)
- [Detecting archaic introgression using an unadmixed outgroup](https://doi.org/10.1371/journal.pgen.1007641
)
- [Altitude adaptation in Tibetans caused by introgression of Denisovan-like DNA](https://doi.org/10.1038/nature13408)

### Data

[Google Drive folder with data](https://drive.google.com/open?id=1lrRfFcoxpyVpXgOi4RYP2_vauM_-rRz_)

- ArchaicSegments.txt: This files is formatted the same was as the one you used for the archiac admixture exercise.
- SNP.txt: This file has all the non-African SNPs that remain after removing all SNPs found in (Subsaharan) Africa. Each SNP is found in any of the high coverage archaic genomes (Altai or Vindija Neanderthals or the Denisova) it is labelled as such. Otherwise it is labeled "human".


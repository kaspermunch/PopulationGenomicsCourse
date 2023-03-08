# Population Genomics 2021

## Description
The participants will after the course have detailed knowledge of the methods and applications required to perform a typical population genomic study.

The participants must at the end of the course be able to:

* Identify an experimental platform relevant to a population genomic analysis.
* Apply commonly used population genomic methods.
* Explain the theory behind common population genomic methods.
* Reflect on strengths and limitations of population genomic methods.
* Interpret and analyze results of population genomic inference.
* Formulate population genetics hypotheses based on data

The course introduces key concepts in population genomics from generation of population genetic data sets to the most common population genetic analyses and association studies. The first part of the course focuses on generation of population genetic data sets. The second part introduces the most common population genetic analyses and their theoretical background. Here topics include analysis of demography, population structure, recombination and selection. The last part of the course focus on applications of population genetic data sets for association studies in relation to human health.

## Curriculum
The curriculum for each week is listed below. "Coop" refers to a set of [lecture notes by Graham Coop](https://github.com/cooplab/popgen-notes/releases/download/v1.2/minicoop.pdf) that we will use throughout the course.

## Padlet
We will use a padlet for shared communication about the curriculum. There I may post questions to guide your studies and you can comment on which parts of the curriculum you find most challenging so we can focus on that.

link to padlet

## Access to computing cluster
You will do the exercises on the GenomeDK computing cluster. So before the course begins you must request a user acount by applying [here](https://console.genome.au.dk/user-requests/create/). You need to fill in some information. Most of it is selv-explanatory. For "Organization" fill in "Aarhus University", for "Department" fill in "BiRC", for "Zone" choose "Open", for "Reason" fill in "Population genomics course", for "Username" fill in a short username that you think might be unique.

## Student presentations
You will each do two student presentations together with a fellow student. If possible, you should sign up for one presentation in the first half of the course and one in the last half. In this [Google Sheet](https://docs.google.com/spreadsheets/d/1XuTLhy8Kx14y9XGm_fK9hz6CqvQy79-Mh8IONdPs_PE/edit?usp=sharing), you can see the available dates and the topics to choose from on each date. Fill in your name as "student one" or "student two" for two dates.

## Lectures
Lectures/discussions are on Tuesdays: 12-14. You can see the curriculum for each lecture in the weekly plan below. Each lecture session will be structured roughly as follows:

- 15 min student presentation on a topic related to the past week's week's curriculum.
- 30 min lecture based on the current week's curriculum.
- 15 break
- 45 min student discussion on the current week's curriculum.

## Exercises
Exercises are on Thursdays: 12-15. You can find the exercises on the GitHub page for the course.

## Week plan	

1. Course intro and overview: 
   - Lecture (Kasper): Coop chapters 1, 2, 3, [Paper: Genome Diversity Project](https://www.nature.com/articles/nature18964)
    - Exercise (Jilong): [Cluster practicals](https://github.com/kaspermunch/PopulationGenomicsCourse/tree/master/Exercises/01_cluster_practicals)
2. Drift and the coalescent:
    - Lecture (Kasper): Coop chapter 4; [Paper: Platypus](https://www.nature.com/articles/ng.3036)
    - Exercise (Jilong): [Read mapping and base calling](https://github.com/kaspermunch/PopulationGenomicsCourse/tree/master/Exercises/02_mapping_and_calling)
3. Recombination:
    - Lecture (Juraj): [Review: Recombination in eukaryotes](https://royalsocietypublishing.org/doi/10.1098/rstb.2016.0455), [Review: Recombination rate estimation](https://www.nature.com/articles/s41576-020-0240-1)
    - Exercise (Jilong): [Phasing and recombination rate](https://github.com/kaspermunch/PopulationGenomicsCourse/tree/master/Exercises/04_phasing_and_recombination)
4. Population strucure and incomplete lineage sorting:
    - Lecture (Kasper): Coop chapter 6, [Review: Incomplete lineage sorting](https://doi.org/10.1146/annurev-genet-120213-092532)
    - Exercise (Jilong): [Working with VCF files](https://github.com/kaspermunch/PopulationGenomicsCourse/tree/master/Exercises/03_f_statistics)
5. Hidden Markov models:
    - Lecture (Kasper): Durbin chapter 3, [Paper: population structure](https://www.nature.com/articles/nature07331)
    - Exercise (Jilong): [Inference of population structure and admixture](https://github.com/kaspermunch/PopulationGenomicsCourse/tree/master/Exercises/05_population_structure)
6. Ancestral recombination graphs:
    - Lecture (Kasper): [Paper: Approximating the ARG](https://bmcgenomdata.biomedcentral.com/articles/10.1186/1471-2156-7-16), [Paper: Tree inference](https://www.nature.com/articles/s41588-019-0484-x)
    - Exercise (Jilong): ARG dashboard exercises + Inference of trees along sequence
7. Past population demography:
    - Lecture (Juraj): Coop chapter 4, [Paper: PSMC](https://www.nature.com/articles/nature10231), revisit [Paper: Tree inference](https://www.nature.com/articles/s41588-019-0484-x)
    - Exercise (Jilong): [Inferring historical populations](https://github.com/kaspermunch/PopulationGenomicsCourse/tree/master/Exercises/06_historical_population_size)
8. Direct and linked selection:
    - Lecture (Kasper): Coop chapters 12, 13, revisit [Paper: Tree inference](https://www.nature.com/articles/s41588-019-0484-x)
    - Exercise (Jilong): TBA
8. Admixture:
   - Lecture (Mikkel): [Review: Admixture](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007349), [Paper: Admixture inference](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007641)
   - Exercise (Jilong): [Detecting archaic ancestry in modern humans](https://github.com/kaspermunch/PopulationGenomicsCourse/tree/master/Exercises/08_archaic_humans)
10. Genome-wide association study (GWAS):
    - Lecture (Søren): GWAS review, Population stratification review, [Coop lecture notes](https://github.com/cooplab/popgen-notes/releases/download/v1.2/release_popgen_notes.pdf) 99-120
    - Exercise (Jilong): [GWAS quality control](https://github.com/kaspermunch/PopulationGenomicsCourse/tree/master/Exercises/09_GWAS_QC)
11. Heritability:
    - Lecture (Søren): [Missing heritability and mixed models review]() ; Coop Lecture notes Sec. 2.2 (p23-36) + Chap. 7 (p119-142)		
    - Exercise (Jilong): [Association testing](https://github.com/kaspermunch/PopulationGenomicsCourse/tree/master/Exercises/10_GWAS_association)
12. Evolution and disease:
    - Lecture (Søren): [Genetic architecture review]() ; [Article about "omnigenic" model]() ; Coop Lecture notes Sec. 11.0.1 (p217-221)	
    - Exercise (Jilong): [Estimating heritability](https://github.com/kaspermunch/PopulationGenomicsCourse/tree/master/Exercises/11_heritability)
13. Project presentations:	
    - Lecture (Kasper, Mikkel, Søren): None
    - Exercise (Jilong):  None, Focus on projects
14. Project guidance:	
    - Lecture (Kasper, Mikkel, Søren): None
    - Exercise (Jilong):  None, Focus on projects
15. Project guidance:	
    - Lecture (Kasper, Mikkel, Søren): None
    - Exercise (Jilong):  None, Focus on projects

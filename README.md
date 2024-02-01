# Population genomics course

## Description
After the course, the participants will have detailed knowledge of the methods and applications required to perform a typical population genomic study.

The participants must, at the end of the course, be able to:

* Identify an experimental platform relevant to a population genomic analysis.
* Apply commonly used population genomic methods.
* Explain the theory behind standard population genomic methods.
* Reflect on the strengths and limitations of population genomic methods.
* Interpret and analyze results of population genomic inference.
* Formulate population genetics hypotheses based on data

The course introduces fundamental concepts in population genomics, from generating population genetic data sets to the most common population genetic analyses and association studies. The course's first part focuses on generating population genetic data sets. The second part introduces the most common population genetic analyses and their theoretical background. Here, topics include analysis of demography, population structure, recombination, and selection. The last part of the course focuses on applications of population genetic data sets for association studies relating to human health.

## Curriculum
The curriculum for each week is listed below. "Coop" refers to a set of [lecture notes by Graham Coop](https://github.com/cooplab/popgen-notes/releases/download/v1.2/minicoop.pdf) that we will use throughout the course.

## Padlet
We will use a padlet for shared communication about the curriculum. There, I may post questions to guide your studies, and you can comment on which parts of the curriculum you find most challenging so we can focus on that.

[link to padlet](https://padlet.com/kaspermunch/population-genomics-2024-404m16k90hk54n97)

## Access to the computing cluster
You will do the exercises on the GenomeDK computing cluster. So before the course begins, you must request a user account by applying [here](https://console.genome.au.dk/user-requests/create/). You need to fill in some information. Most of it is self-explanatory. For "Organization" fill in "Aarhus University", for "Department" fill in "BiRC", for "Zone" choose "Open", for "Reason" fill in "Population genomics course", for "Username" fill in a short username that you think might be unique.

## Student presentations
You will each do two student presentations together with a fellow student. You should sign up for one presentation in the first half of the course and one in the last half. In this [Google Sheet](https://docs.google.com/spreadsheets/d/1XuTLhy8Kx14y9XGm_fK9hz6CqvQy79-Mh8IONdPs_PE/edit?usp=sharing), you can see the available dates and the topics to choose from on each date. Fill in your name as "student one" or "student two" for two dates.

## Lectures
Lectures/discussions are on Mondays from 12:15 to 15:00. You can see the curriculum for each lecture in the weekly plan below. Each lecture session will be structured like this:
* 10-minute student presentation on a topic related to the past week's curriculum.
* 10-minute student presentation on a topic connected to the past week's week's curriculum.
* 30-minute lecture based on the current week's curriculum.
* 15 break
* 45-minute lecture/discussion on the current week's curriculum.
* 15 break
* 45-minute lecture/discussion on the current week's curriculum.



## Exercises
Exercises are on Thursdays from 12:15 to 14:00. The week plan below has links to the exercises hosted on the GitHub page for the course.

## Week plan    

<!-- 
TODO:
Change lectures to the 3-hour slot
Maybe drop either calling or phasing
Maybe change from LD hat to pyro
Drop PCA and Admixture and do MOSAIC instead
Take over the admixture lecture from Mikkel
Maybe start the projects earlier so they work on them on the side for longer
-->

<!-- 5. Hidden Markov models:
    - Lecture : Durbin chapter 3, [Paper: population structure](https://www.nature.com/articles/nature07331)
    - Exercise (Bjarke): [Inference of population structure and admixture](https://github.com/kaspermunch/PopulationGenomicsCourse/tree/master/Exercises/05_population_structure)
 -->


<!-- - Exercise (Bjarke): [Working with VCF files](https://github.com/kaspermunch/PopulationGenomicsCourse/tree/master/Exercises/03_f_statistics) -->

<!-- 5. Ancestral recombination graphs:
    - Lecture : [Paper: Approximating the ARG](https://bmcgenomdata.biomedcentral.com/articles/10.1186/1471-2156-7-16), [Paper: Tree inference](https://www.nature.com/articles/s41588-019-0484-x)
    - Exercise (Bjarke): ARG dashboard exercises + Inference of trees along sequence -->

1. Course intro and overview: 
   - Lecture (Kasper): Coop chapters 1, 2, 3, [Paper: Simons Genome Diversity Project](https://www.nature.com/articles/nature18964)
    - Exercise (Bjarke): [Cluster practicals](https://github.com/kaspermunch/PopulationGenomicsCourse/tree/master/Exercises/cluster_practicals)
2. Drift and the coalescent:
    - Lecture (Kasper): Coop chapter 4; 
    - Exercise (Bjarke): [Read mapping and base calling](https://github.com/kaspermunch/PopulationGenomicsCourse/tree/master/Exercises/mapping_and_calling)
3. Recombination and Ancestral recombination graphs:
    - Lecture (Kasper): [Paper: Platypus](https://www.nature.com/articles/ng.3036), [Review: Recombination in eukaryotes](https://royalsocietypublishing.org/doi/10.1098/rstb.2016.0455), [Review: Recombination rate estimation](https://www.nature.com/articles/s41576-020-0240-1)
    - Exercise (Bjarke): [Phasing and recombination rate](https://github.com/kaspermunch/PopulationGenomicsCourse/tree/master/Exercises/phasing_and_recombination)
4. Past population demography, Ancestral recombination graph:
    - Lecture (Kasper): Revisit Coop chapter 4, [Paper: Approximating the ARG](https://bmcgenomdata.biomedcentral.com/articles/10.1186/1471-2156-7-16), [Paper: Tree inference](https://www.nature.com/articles/s41588-019-0484-x)
    - Exercise (Bjarke): [Tree sequences and inference of demography](https://github.com/kaspermunch/PopulationGenomicsCourse/tree/master/Exercises/tree_sequences)
5. Population structure and incomplete lineage sorting:
    - Lecture: Coop chapter 6, [Review: Incomplete lineage sorting](https://doi.org/10.1146/annurev-genet-120213-092532)
    - Exercise: [Inference of population structure and admixture](https://github.com/kaspermunch/PopulationGenomicsCourse/tree/master/Exercises/05_population_structure)
6. Direct and linked selection:
    - Lecture: Coop chapters 12, 13, revisit [Paper: Tree inference](https://www.nature.com/articles/s41588-019-0484-x)
    - Exercise (Bjarke): [Inference of positive selection](https://github.com/kaspermunch/PopulationGenomicsCourse/tree/master/Exercises/selection)
7. Admixture:
   - Lecture (Kasper): [Review: Admixture](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007349), [Paper: Admixture inference](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007641)
   - Exercise (Bjarke): [Detecting archaic ancestry in modern humans](https://github.com/kaspermunch/PopulationGenomicsCourse/tree/master/Exercises/archaic_humans)
8. Genome-wide association study (GWAS):
    - Lecture (Søren): GWAS review, Population stratification review, [Coop lecture notes](https://github.com/cooplab/popgen-notes/releases/download/v1.2/release_popgen_notes.pdf) 99-120
    - Exercise (Bjarke): [GWAS quality control](https://github.com/kaspermunch/PopulationGenomicsCourse/tree/master/Exercises/GWAS_QC)
9.  Heritability:
    - Lecture (Søren): [Missing heritability and mixed models review]() ; Coop Lecture notes Sec. 2.2 (p23-36) + Chap. 7 (p119-142)     
    - Exercise (Bjarke): [Association testing](https://github.com/kaspermunch/PopulationGenomicsCourse/tree/master/Exercises/GWAS_association)
10. Trees and GWAS:
    - Lecture (Kasper): [TBA]()
    - Exercise (Bjarke): [TBA]()
11.  Evolution and disease:
    - Lecture (Søren): [Genetic architecture review]() ; [Article about "omnigenic" model]() ; Coop Lecture notes Sec. 11.0.1 (p217-221)    
    - Exercise (Bjarke): [Estimating heritability](https://github.com/kaspermunch/PopulationGenomicsCourse/tree/master/Exercises/heritability)
12.  Project presentations:  
     - Lecture (Kasper, Søren): None
     - Exercise (Bjarke):  None. Focus on projects
13.  Project guidance:   
     - Lecture (Kasper, Søren): None
     - Exercise (Bjarke):  None. Focus on projects
14. Project guidance:   
    - Lecture (Kasper, Søren): None
    - Exercise (Bjarke):  None. Focus on projects

# Introduction to Population Genomics
**A course of the danish health data science sandbox**

This course is based on the material developed for the Population Genomics course at Aarhus university. The material is organized in four separated jupyter notebooks in both `R`, `bash`, and `python`, where you will benefit of an interactive coding setup. 

If you use any of this material for your research, please cite this course with the DOI below, and acknowledge the Health Data Science Sandbox project of the Novo Nordisk Foundation (grant number NNF20OC0063268). It is of great help to support the project and the creation of new courses.
[![DOI](https://zenodo.org/badge/468293635.svg)](https://zenodo.org/badge/latestdoi/468293635)


## Course description

The course introduces key concepts in population genomics from generation of population genetic data sets to the most common population genetic analyses and association studies. The first part of the course focuses on generation of population genetic data sets. The second part introduces the most common population genetic analyses and their theoretical background. Here topics include analysis of demography, population structure, recombination and selection. The last part of the course focus on applications of population genetic data sets for association studies in relation to human health.

----------------------

### Prerequisites

This is an introductory course that needs a basic understanding of genomics, and not necessarily programming experience (thought that helps).
  
### Learning Outcomes

After the course, you will have detailed knowledge of the methods and applications required to perform a typical population genomic study.
You will be able to:

* Identify an experimental platform relevant to a population genomic analysis.
* Apply commonly used population genomic methods.
* Explain the theory behind common population genomic methods.
* Reflect on strengths and limitations of population genomic methods.
* Interpret and analyze results of population genomic inference.
* Formulate population genetics hypotheses based on data


### Supporting material

* jupyter notebooks for interactive coding
* Structure of the course with lecture list

The curriculum for each week of the course is listed below. "Coop" refers to a set of [lecture notes by Graham Coop](https://github.com/cooplab/popgen-notes/releases/download/v1.2/minicoop.pdf) that are used throughout the course.

### Course duration and structure

This course is one-semester long.


1. Course intro and overview: 
   - Lecture (Kasper): Coop chapters 1, 2, 3, [Paper: Genome Diversity Project](https://www.nature.com/articles/nature18964)
    - Exercise: Cluster practicals
2. Drift and the coalescent:
    - Lecture: Coop chapter 4; [Paper: Platypus](https://www.nature.com/articles/ng.3036)
    - Exercise: [Read mapping and base calling]()
3. Recombination:
    - Lecture: [Review: Recombination in eukaryotes](https://royalsocietypublishing.org/doi/10.1098/rstb.2016.0455), [Review: Recombination rate estimation](https://www.nature.com/articles/s41576-020-0240-1)
    - Exercise: [Phasing and recombination rate](https://github.com/kaspermunch/PopulationGenomicsCourse/tree/master/Exercises/04_phasing_and_recombination)
4. Population strucure and incomplete lineage sorting:
    - Lecture: Coop chapter 6, [Review: Incomplete lineage sorting](https://zh.booksc.eu/book/32923932/889942)
    - Exercise: [Working with VCF files](https://github.com/kaspermunch/PopulationGenomicsCourse/tree/master/Exercises/03_f_statistics)
5. Hidden Markov models:
    - Lecture : Durbin chapter 3, [Paper: population structure](https://www.nature.com/articles/nature07331)
    - Exercise: [Inference of population structure and admixture](https://github.com/kaspermunch/PopulationGenomicsCourse/tree/master/Exercises/05_population_structure)
6. Ancestral recombination graphs:
    - Lecture: [Paper: Approximating the ARG](https://bmcgenomdata.biomedcentral.com/articles/10.1186/1471-2156-7-16), [Paper: Tree inference](https://www.nature.com/articles/s41588-019-0484-x)
    - Exercise: ARG dashboard exercises + Inference of trees along sequence
7. Past population demography:
    - Lecture (Juraj): Coop chapter 4, [Paper: PSMC](https://www.nature.com/articles/nature10231), revisit [Paper: Tree inference](https://www.nature.com/articles/s41588-019-0484-x)
    - Exercise: [Inferring historical populations](https://github.com/kaspermunch/PopulationGenomicsCourse/tree/master/Exercises/06_historical_population_size)
8. Direct and linked selection:
    - Lecture: Coop chapters 12, 13, revisit [Paper: Tree inference](https://www.nature.com/articles/s41588-019-0484-x)
9. Admixture:
   - Lecture: [Review: Admixture](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007349), [Paper: Admixture inference](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007641)
   - Exercise: [Detecting archaic ancestry in modern humans](https://github.com/kaspermunch/PopulationGenomicsCourse/tree/master/Exercises/08_archaic_humans)
10. Genome-wide association study (GWAS):
    - Lecture: Coop lecture notes 99-120
    - Exercise: [GWAS quality control](https://github.com/kaspermunch/PopulationGenomicsCourse/tree/master/Exercises/09_GWAS_QC)
11. Heritability:
    - Lecture: Coop Lecture notes Sec. 2.2 (p23-36) + Chap. 7 (p119-142)		
    - Exercise: [Association testing](https://github.com/kaspermunch/PopulationGenomicsCourse/tree/master/Exercises/10_GWAS_association)
12. Evolution and disease:
    - Lecture : Coop Lecture notes Sec. 11.0.1 (p217-221)	
    - Exercise: [Estimating heritability](https://github.com/kaspermunch/PopulationGenomicsCourse/tree/master/Exercises/11_heritability)


### Course authors

Head of the course: Kasper Munch.

Contact: Samuele Soraggi (samuele at birc.au.dk).

-----------------------------



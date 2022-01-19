# Instructor zone

This page is for practical information related exercises on the cluster.
## The populationgenomics project folder

The project folder `populationgenomics` has the following dirs, many of which are symlinks to dirs under the [PopulationGenomicsCourse](https://github.com/kaspermunch/PopulationGenomicsCourse) git project checked out to my directory in `people/kmt/`:

```
populationgenomics
├── data
├── fasta -> people/kmt/PopulationGenomicsCourse/Cluster/data/fasta/
├── instruktor_data -> people/kmt/PopulationGenomicsCourse/Cluster/data
├── metadata -> people/kmt/PopulationGenomicsCourse/Cluster/metadata
├── people
├── README.md -> people/kmt/PopulationGenomicsCourse/Cluster/README.md
├── shared_results
├── software
└── students
```

- `software`: builds and binaries for LDhat and PSMC
- `people`: teachers and instructors folders
- `students`: student folders
- `instructor_data` links to a folder that contains all data files including the original and intermediary files needed for building the files that the students need for the exercises.
- `data` is the folder where the student finds the data files he/she needs for the exercises. This folder contains symlinks to only the relevant data folders and files:
- 
```
data
├── archaic -> ../people/kmt/PopulationGenomicsCourse/Cluster/data/archaic
├── assignment -> ../people/kmt/PopulationGenomicsCourse/Cluster/data/assignment
├── bam
│   ├── README
│   ├── S_Ami-1.chr2.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Ami-1/S_Ami-1.chr2.bam
│   ├── S_Ami-1.chr2.bam.bai -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Ami-1/S_Ami-1.chr2.bam.bai
│   ├── S_Ami-1.region.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Ami-1/S_Ami-1.region.bam
│   ├── S_Atayal-1.chr2.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Atayal-1/S_Atayal-1.chr2.bam
│   ├── S_Atayal-1.chr2.bam.bai -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Atayal-1/S_Atayal-1.chr2.bam.bai
│   ├── S_Atayal-1.region.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Atayal-1/S_Atayal-1.region.bam
│   ├── S_Bulgarian-1.chr2.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Bulgarian-1/S_Bulgarian-1.chr2.bam
│   ├── S_Bulgarian-1.chr2.bam.bai -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Bulgarian-1/S_Bulgarian-1.chr2.bam.bai
│   ├── S_Bulgarian-1.region.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Bulgarian-1/S_Bulgarian-1.region.bam
│   ├── S_Bulgarian-2.chr2.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Bulgarian-2/S_Bulgarian-2.chr2.bam
│   ├── S_Bulgarian-2.chr2.bam.bai -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Bulgarian-2/S_Bulgarian-2.chr2.bam.bai
│   ├── S_Bulgarian-2.region.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Bulgarian-2/S_Bulgarian-2.region.bam
│   ├── S_Cambodian-1.chr2.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Cambodian-1/S_Cambodian-1.chr2.bam
│   ├── S_Cambodian-1.chr2.bam.bai -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Cambodian-1/S_Cambodian-1.chr2.bam.bai
│   ├── S_Cambodian-1.region.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Cambodian-1/S_Cambodian-1.region.bam
│   ├── S_Cambodian-2.chr2.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Cambodian-2/S_Cambodian-2.chr2.bam
│   ├── S_Cambodian-2.chr2.bam.bai -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Cambodian-2/S_Cambodian-2.chr2.bam.bai
│   ├── S_Cambodian-2.region.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Cambodian-2/S_Cambodian-2.region.bam
│   ├── S_Druze-2.chr2.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Druze-2/S_Druze-2.chr2.bam
│   ├── S_Druze-2.chr2.bam.bai -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Druze-2/S_Druze-2.chr2.bam.bai
│   ├── S_Druze-2.region.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Druze-2/S_Druze-2.region.bam
│   ├── S_English-1.chr2.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_English-1/S_English-1.chr2.bam
│   ├── S_English-1.chr2.bam.bai -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_English-1/S_English-1.chr2.bam.bai
│   ├── S_English-1.region.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_English-1/S_English-1.region.bam
│   ├── S_Esan-2.chr2.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Esan-2/S_Esan-2.chr2.bam
│   ├── S_Esan-2.chr2.bam.bai -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Esan-2/S_Esan-2.chr2.bam.bai
│   ├── S_Esan-2.region.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Esan-2/S_Esan-2.region.bam
│   ├── S_Georgian-1.chr2.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Georgian-1/S_Georgian-1.chr2.bam
│   ├── S_Georgian-1.chr2.bam.bai -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Georgian-1/S_Georgian-1.chr2.bam.bai
│   ├── S_Georgian-1.region.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Georgian-1/S_Georgian-1.region.bam
│   ├── S_Georgian-2.chr2.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Georgian-2/S_Georgian-2.chr2.bam
│   ├── S_Georgian-2.chr2.bam.bai -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Georgian-2/S_Georgian-2.chr2.bam.bai
│   ├── S_Georgian-2.region.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Georgian-2/S_Georgian-2.region.bam
│   ├── S_Hungarian-1.chr2.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Hungarian-1/S_Hungarian-1.chr2.bam
│   ├── S_Hungarian-1.chr2.bam.bai -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Hungarian-1/S_Hungarian-1.chr2.bam.bai
│   ├── S_Hungarian-1.region.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Hungarian-1/S_Hungarian-1.region.bam
│   ├── S_Hungarian-2.chr2.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Hungarian-2/S_Hungarian-2.chr2.bam
│   ├── S_Hungarian-2.chr2.bam.bai -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Hungarian-2/S_Hungarian-2.chr2.bam.bai
│   ├── S_Hungarian-2.region.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Hungarian-2/S_Hungarian-2.region.bam
│   ├── S_Icelandic-1.chr2.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Icelandic-1/S_Icelandic-1.chr2.bam
│   ├── S_Icelandic-1.chr2.bam.bai -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Icelandic-1/S_Icelandic-1.chr2.bam.bai
│   ├── S_Icelandic-1.region.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Icelandic-1/S_Icelandic-1.region.bam
│   ├── S_Iranian-1.chr2.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Iranian-1/S_Iranian-1.chr2.bam
│   ├── S_Iranian-1.chr2.bam.bai -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Iranian-1/S_Iranian-1.chr2.bam.bai
│   ├── S_Iranian-1.region.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Iranian-1/S_Iranian-1.region.bam
│   ├── S_Japanese-2.chr2.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Japanese-2/S_Japanese-2.chr2.bam
│   ├── S_Japanese-2.chr2.bam.bai -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Japanese-2/S_Japanese-2.chr2.bam.bai
│   ├── S_Japanese-2.region.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Japanese-2/S_Japanese-2.region.bam
│   ├── S_Ju_hoan_North-2.region.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Ju_hoan_North-2/S_Ju_hoan_North-2.region.bam
│   ├── S_Ju_hoan_North-3.chr2.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Ju_hoan_North-3/S_Ju_hoan_North-3.chr2.bam
│   ├── S_Ju_hoan_North-3.chr2.bam.bai -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Ju_hoan_North-3/S_Ju_hoan_North-3.chr2.bam.bai
│   ├── S_Ju_hoan_North-3.region.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Ju_hoan_North-3/S_Ju_hoan_North-3.region.bam
│   ├── S_Kinh-2.chr2.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Kinh-2/S_Kinh-2.chr2.bam
│   ├── S_Kinh-2.chr2.bam.bai -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Kinh-2/S_Kinh-2.chr2.bam.bai
│   ├── S_Kinh-2.region.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Kinh-2/S_Kinh-2.region.bam
│   ├── S_Korean-2.chr2.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Korean-2/S_Korean-2.chr2.bam
│   ├── S_Korean-2.chr2.bam.bai -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Korean-2/S_Korean-2.chr2.bam.bai
│   ├── S_Korean-2.region.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Korean-2/S_Korean-2.region.bam
│   ├── S_Luhya-2.chr2.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Luhya-2/S_Luhya-2.chr2.bam
│   ├── S_Luhya-2.chr2.bam.bai -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Luhya-2/S_Luhya-2.chr2.bam.bai
│   ├── S_Luhya-2.region.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Luhya-2/S_Luhya-2.region.bam
│   ├── S_Luo-2.chr2.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Luo-2/S_Luo-2.chr2.bam
│   ├── S_Luo-2.chr2.bam.bai -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Luo-2/S_Luo-2.chr2.bam.bai
│   ├── S_Luo-2.region.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Luo-2/S_Luo-2.region.bam
│   ├── S_Mandenka-1.chr2.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Mandenka-1/S_Mandenka-1.chr2.bam
│   ├── S_Mandenka-1.chr2.bam.bai -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Mandenka-1/S_Mandenka-1.chr2.bam.bai
│   ├── S_Mandenka-1.region.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Mandenka-1/S_Mandenka-1.region.bam
│   ├── S_Mandenka-2.chr2.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Mandenka-2/S_Mandenka-2.chr2.bam
│   ├── S_Mandenka-2.chr2.bam.bai -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Mandenka-2/S_Mandenka-2.chr2.bam.bai
│   ├── S_Mandenka-2.region.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Mandenka-2/S_Mandenka-2.region.bam
│   ├── S_Miao-1.chr2.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Miao-1/S_Miao-1.chr2.bam
│   ├── S_Miao-1.chr2.bam.bai -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Miao-1/S_Miao-1.chr2.bam.bai
│   ├── S_Miao-1.region.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Miao-1/S_Miao-1.region.bam
│   ├── S_Naxi-2.chr2.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Naxi-2/S_Naxi-2.chr2.bam
│   ├── S_Naxi-2.chr2.bam.bai -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Naxi-2/S_Naxi-2.chr2.bam.bai
│   ├── S_Naxi-2.region.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Naxi-2/S_Naxi-2.region.bam
│   ├── S_Yoruba-1.chr2.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Yoruba-1/S_Yoruba-1.chr2.bam
│   ├── S_Yoruba-1.chr2.bam.bai -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Yoruba-1/S_Yoruba-1.chr2.bam.bai
│   ├── S_Yoruba-1.region.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Yoruba-1/S_Yoruba-1.region.bam
│   ├── S_Yoruba-2.chr2.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Yoruba-2/S_Yoruba-2.chr2.bam
│   ├── S_Yoruba-2.chr2.bam.bai -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Yoruba-2/S_Yoruba-2.chr2.bam.bai
│   └── S_Yoruba-2.region.bam -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Yoruba-2/S_Yoruba-2.region.bam
├── consensus_fastq -> ../people/kmt/PopulationGenomicsCourse/Cluster/data/consensus_fastq
├── fasta
│   ├── chr2.fa
│   └── chr2.fa.fai
├── fastq
│   ├── README
│   ├── S_Ami-1.region.fq -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Ami-1/S_Ami-1.region.fq
│   ├── S_Atayal-1.region.fq -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Atayal-1/S_Atayal-1.region.fq
│   ├── S_Bulgarian-1.region.fq -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Bulgarian-1/S_Bulgarian-1.region.fq
│   ├── S_Bulgarian-2.region.fq -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Bulgarian-2/S_Bulgarian-2.region.fq
│   ├── S_Cambodian-1.region.fq -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Cambodian-1/S_Cambodian-1.region.fq
│   ├── S_Cambodian-2.region.fq -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Cambodian-2/S_Cambodian-2.region.fq
│   ├── S_Druze-2.region.fq -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Druze-2/S_Druze-2.region.fq
│   ├── S_English-1.region.fq -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_English-1/S_English-1.region.fq
│   ├── S_Esan-2.region.fq -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Esan-2/S_Esan-2.region.fq
│   ├── S_Georgian-1.region.fq -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Georgian-1/S_Georgian-1.region.fq
│   ├── S_Georgian-2.region.fq -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Georgian-2/S_Georgian-2.region.fq
│   ├── S_Hungarian-1.region.fq -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Hungarian-1/S_Hungarian-1.region.fq
│   ├── S_Hungarian-2.region.fq -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Hungarian-2/S_Hungarian-2.region.fq
│   ├── S_Icelandic-1.region.fq -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Icelandic-1/S_Icelandic-1.region.fq
│   ├── S_Iranian-1.region.fq -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Iranian-1/S_Iranian-1.region.fq
│   ├── S_Japanese-2.region.fq -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Japanese-2/S_Japanese-2.region.fq
│   ├── S_Ju_hoan_North-3.region.fq -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Ju_hoan_North-3/S_Ju_hoan_North-3.region.fq
│   ├── S_Kinh-2.region.fq -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Kinh-2/S_Kinh-2.region.fq
│   ├── S_Korean-2.region.fq -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Korean-2/S_Korean-2.region.fq
│   ├── S_Luhya-2.region.fq -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Luhya-2/S_Luhya-2.region.fq
│   ├── S_Luo-2.region.fq -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Luo-2/S_Luo-2.region.fq
│   ├── S_Mandenka-1.region.fq -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Mandenka-1/S_Mandenka-1.region.fq
│   ├── S_Mandenka-2.region.fq -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Mandenka-2/S_Mandenka-2.region.fq
│   ├── S_Miao-1.region.fq -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Miao-1/S_Miao-1.region.fq
│   ├── S_Naxi-2.region.fq -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Naxi-2/S_Naxi-2.region.fq
│   ├── S_Yoruba-1.region.fq -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Yoruba-1/S_Yoruba-1.region.fq
│   └── S_Yoruba-2.region.fq -> ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq/S_Yoruba-2/S_Yoruba-2.region.fq
├── genetic_map
│   ├── plink.chr2.GRCh37.map
│   ├── plink.GRCh37.map.zip
│   ├── plink.README.txt
│   ├── README
│   └── README.txt
├── GWAS -> ../people/kmt/PopulationGenomicsCourse/Cluster/data/GWAS
├── haplotypes_chrX
│   ├── genotypes360_400_AF
│   ├── genotypes360_400_EA
│   ├── genotypes360_400_SA
│   ├── genotypes360_400_WE
│   └── snps360_400_filtered
├── ldhat
│   ├── ldhat.r
│   └── new_lk.txt
├── metadata -> ../people/kmt/PopulationGenomicsCourse/Cluster/metadata
└── vcf -> ../people/kmt/PopulationGenomicsCourse/Cluster/data/vcf
```

The symlinks in the `bam` folder are made like this:

```
find ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq -name '*.region.bam' -exec ln -s {} \;
find ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq -name '*.chr2.bam' -exec ln -s {} \;
find ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq -name '*.chr2.bam.bai' -exec ln -s {} \;
```

and the ones in the `fastq` folder are made like this:

```
find ../../people/kmt/PopulationGenomicsCourse/Cluster/data/bam_and_fastq -name '*region.fq' -exec ln -s {} \;
```




# Data

Large data files for exercises and projects are placed under the `PopulationGenomicsCourse`. They are NOT not under version control but are backed up on the cluster. Only the huge bam files are not backed up. To re-build data files for chr2 not under version control (but backed up), run this in `data/bam_and_fastq`:

    conda create -n popgen_data -c gwforg -c bioconda python=3 gwf bedtools samtools bcftools

    conda activate popgen_data

    sbatch download_bamfiles.sh

    gwf -f ../../scripts/workflow_bamfiles.py run

## Software

All command line software, except PSMC and LDhat, is installed in the `popgen` environment.  and LDhat PSMC is available in the software folder

## Running on the cluster

The students begin each exercise by running this to get an interactive session for the duration of the exercise:

    srun --mem-per-cpu=1g --time=3:00:00 --account=populationgenomics --pty bash

## Jupyter








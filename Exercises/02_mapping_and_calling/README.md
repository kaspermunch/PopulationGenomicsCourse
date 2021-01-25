

# Mapping and SNP calling exercise

As we learned last week, high-throughput sequencing technologies have in the past few years been producing millions of reads of human genome and other species. To be useful, this genetic information has to be 'put together' in a smart way, in the same way as the pieces of a puzzle (reads) need to be mounted according to a picture (reference genome). In this exercise section you will be exposed to different softwares used for mapping reads to a reference sequence and calling variants from the produced alignments. We will use a dataset composed of 30 individuals from 3 different regions: Africa, EastAsia and WestEurasia.

    ##       X...ID    ENA.RUN    population region country latitude longitude
    ## 1 ERS1042176 ERR1019075 Ju_hoan_North Africa Namibia    -18.9      21.5
    ## 2 ERS1042177 ERR1019076 Ju_hoan_North Africa Namibia    -18.9      21.5
    ## 3 ERS1042248 ERR1025622          Esan Africa Nigeria      6.5       6.0
    ## 4 ERS1042265 ERR1025639         Luhya Africa   Kenya      1.3      36.8
    ## 5 ERS1042266 ERR1025640      Mandenka Africa Senegal     12.0     -12.0
    ## 6 ERS1042267 ERR1025641      Mandenka Africa Senegal     12.0     -12.0
    ##      Sex       Illumina.ID
    ## 1   male LP6005441-DNA_B11
    ## 2   male LP6005441-DNA_A11
    ## 3 female LP6005442-DNA_B10
    ## 4   male LP6005442-DNA_E11
    ## 5   male LP6005441-DNA_E07
    ## 6 female LP6005441-DNA_F07

![](img/unnamed-chunk-1-1.png)

This dataset is a subset of the Simons Diversity Project (discussed last week).

## Log into the cluster and request a compute node on the cluster

Log into the cluster:

```bash
    ssh username@login.genome.au.dk
```

Now request a machine for your computations. You need five gigabytes (`5g`) in this exercise so you need to run this command (see also the explanation in the previous exercise):

```bash
srun --mem-per-cpu=5g --time=3:00:00 --account=populationgenomics --pty bash
```

## Data source

You will be separated in pairs so you can help each other out with the commands. Each of you will be responsible for 2 individuals and at the end of this exercise we will estimate the mean SNP heterozygosity per individual of the 10 MB region in chromosome 2. The data is placed in a folder called **Data** in the same directory as users folder. You should introduce you results [here](https://docs.google.com/spreadsheets/d/1OEHI1tNiwHrwKkl9L5rPtbVKCHQDpCZtKzpnZ1sWKJY/edit?usp=sharing) 

The following tutorial is based on **ERR1019076**.

## Mapping reads against the reference

The first step when dealing with raw reads is mapping (aligning) them to a reference sequence. For this, we will be using the BWA mapper. BWA stands for Burrows-Wheeler aligner, which allows for fast and accurate alignment of short reads to an indexed reference sequence. [Here](http://bio-bwa.sourceforge.net/bwa.shtml) is the manual. We decided to focus on a 10 MB region of chromosome 2 (from 135MB to 145MB), which contains the lactase gene.

Two input files are needed to do genome mapping:

- Fasta file containing your reference genome Which was downloaded from:
    ([GRCh37](ftp://ftp.ensembl.org/pub/release-75//fasta/homo_sapiens/dna/))
- The reads in fastq format

We have prepared them for you and can be found in the shared Data folder: /home/Data/

We will create a soft-link of fasta reference to your folder, so that we don't need to type in the full path to the reference everytime we want to use it and avoid copying it to out own directory.

<!-- TODO: Update file names -->
```bash
    ln -s /home/Data/Homo_sapiens.GRCh37.75.dna.chromosome.2.fa /home/[user]/
```

First we need to index the reference file for later use. This creates index files used by bwa mem to perform the alignment. To produce these files, run the following command:

<!-- TODO: Update file names -->
```bash
    bwa index -p Homo_sapiens.GRCh37.75.dna.chromosome.2 -a bwtsw Homo_sapiens.GRCh37.75.dna.chromosome.2.fa
```

where -a bwtsw specifies that we want to use the indexing algorithm that is capable of handling the human genome.

You also need to generate a fasta file index. This can be done using **samtools**:

<!-- TODO: Update file names -->
```bash
    samtools faidx Homo_sapiens.GRCh37.75.dna.chromosome.2.fa
```

Now you can map the reads to the reference. This will take around 10 minutes. You can start installing the software that will be used later in this tutorial (IGV) while you wait for it.

<!-- TODO: Update file names -->
```bash
    bwa mem -t 16 -p Homo_sapiens.GRCh37.75.dna.chromosome.2 /home/Data/sorted_ERR1019076_reads_135_145.fq | \
    samtools sort -O BAM -o ERR1019076.bam
```

<!-- TODO: Update file names -->
This command is composed of two sub-commands where the output of the "bwa mem" command is piped ("|" is the pipe symbol) into the "samtools sort" command. The output of the "bwa mem" command is an unsorted bam file, which is then used as input into the "samtools sort" command to produce a sorted bam file, which is necessary for further analysis. We could also run the two commands separately, but by using piping we save disc space, as we do not have to save the intermediate unsorted bam file, and altogether speed up the analysis.

You can have a look at the bam file generated:

<!-- TODO: Update file names -->
```bash
    samtools view ERR1019076.bam | head
```

Bam files follow this structure:

![Alt text](https://us.v-cdn.net/5019796/uploads/editor/f4/uuzmf2cbau1y.png)

- Read name: ID for the given read.
- Flags: Combination of bitwise FLAGs that provide information on how the read is mapped  [Extra information](https://broadinstitute.github.io/picard/explain-flags.html).
- Position: Chromosome and position of the first base in the alignment. 
- MAPQ: Probability of wrong mapping of the read. It's in Phred scale, so higher numbers mean lower probabilities:
![Alt_text](https://genome.sph.umich.edu/w/images/math/e/9/d/e9dc88d1834c4579de12153a67ac3afa.png)
- CIGAR: summary of the alignment, including start position on the reference sequence, matches, mismatches, deletions and insertions. It may also include information on soft/hard clipping, i.e bases in the 3' and 5' ends of teh read that are not part of the alignment.
[Extra information](https://wiki.bits.vib.be/index.php/CIGAR)
- Mate information: chromosome and start position of teh read pair, and inferred insert size.
- Quality scores: base qualities of the read.
- Metadata: optional extra information. 

For more information, read [this](https://samtools.github.io/hts-specs/SAMv1.pdf)

Get some useful stats of your mapping:

<!-- TODO: Update file names -->
```bash
    samtools flagstat ERR1019076.bam
```

Once the map is generated, you can index the bam file to visualize it using IGV. Indexing a genome sorted BAM file allows one to quickly extract alignments overlapping particular genomic regions. Moreover, indexing is required by genome viewers such as IGV so that the viewers can quickly display alignments in each genomic region to which you navigate.

<!-- TODO: Update file names -->
```bash
    samtools index ERR1019076.bam
```

## Downloading via terminal

You can download the data via terminal by the following:

<!-- TODO: Update file names -->
```bash
    scp username@login.genome.au.dk:populationgenomics/students/username/ERR1019076.bam folder_on_your_computer
```

## IGV software

IGV is an Integrative Genomics viewer and can be very useful to look at the results of Mapping and SNP calling. We have not installed it in the cluster, so you can download it to your machine you can go to its [website](http://software.broadinstitute.org/software/igv/). Three files are necessary to look at this dataset: a reference sequence and the
**.bam** and **.bai** files, download it from the cluster in a specific directory. Since we are using a human reference, the sequence is already available in the software:

Go to Genomes ----> Load Genome from server... ----> Filter by human and choose the Human hg19 reference (which is the GRCh37).

After it you will the chromosomes and genes. Now you can download the mapping results by typing: File ----> Load from File... ----> ERR1019076.bam.

When you zoom in to the lactase region on chromosome 2, you will see something like this: ![](img/IGV_example.png)

Try to understand what are the different attributes present in the viewer. If you zoom in very much you will find single nucleotide polymorphisms (SNPs), where the reference sequence does not have the same nucleotide as the data mapped to.

## Plotting results

<!-- TODO: This should be run on the cluster -->

One of the attributes one could learn from mapping reads back to the reference is the coverage of reads across the genome. In order to calculate the coverage depth you can use the command **samtools depth**.

<!-- TODO: Update file names -->
```bash
    samtools depth ERR1019076.bam > ERR1019076.coverage
```

You can have a look at the resulted file. What do you find in the three different columns?

<!-- TODO: Update file names -->
```bash
    less ERR1019076.coverage
```

<!-- This should be done in a slurm-jupyter R notebook. Add instructions on how -->

Using Rstudio or R, run something like this:

```R
    library(ggplot2)
    library(dplyr)

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

    scaf <- read.table("./ERR1019076.coverage",header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE, col.names = c("Scaffold", "locus", "depth"))
      
    head(scaf)

    ##   Scaffold  locus depth
    ## 1        2 833855     1
    ## 2        2 833856     1
    ## 3        2 833857     1
    ## 4        2 833858     1
    ## 5        2 833859     1
    ## 6        2 833860     1

    # Compressing the dataframe in windows
    scaf %>% 
    mutate(rounded_position = round(locus, -2)) %>%
        group_by(rounded_position) %>% 
            summarize(mean_cov = mean(depth)) -> compressed

    # Plotting the data
    p <- ggplot(data =  compressed, aes(x=rounded_position, y=mean_cov)) + geom_area() + theme_classic() + ylim(0, 400)

    #p

    # Saving your coverage plot
    ggsave("ERR1019076.coverage.pdf",p)

    ## Saving 7 x 5 in image

```

You can do it by opening Rsturio on the terminal or transfer the file to your local machine using scp and then use your locally installed Rstudio. Alternatively, you could also open a text editor in the terminal (e.g vi) and run the script with Rscript [path to the script].R .

What are the conclusions you can extract from these analysis? Does the coverage match with what you observed with IGV? Does it match with what you would expect, i.e what you know from the data? 

## SNP calling

Even though just a tiny portion (around 2%) of our genomes are based of protein coding regions, this partition contains most of the disease causal variants (mutations), and that is why variant calling is so important from a medical point of view. From the population genetics side of view it is also possible to use these variants to establish differences between individuals, populations and species. It can also be used to clarify the genetic basis of adaptation. These topics will come back to your mind during the following weeks.

Once we have mapped our reads we can now start with variant detection. For now we will be using the software **Platypus**: a tool designed for efficient and accurate variant-detection in high-throughput sequencing data. You can access their website [here](http://www.well.ox.ac.uk/platypus). 

To run platypus, we can use this line of code:

<!-- TODO: Update file names -->
```bash
    platypus callVariants --bamFile=ERR1019076.bam --refFile=Homo_sapiens.GRCh37.75.dna.chromosome.2.fa --output=AllVariants.vcf
```

The output will be a single [VCF](http://samtools.github.io/hts-specs/VCFv4.2.pdf) file containing all the variants that Platypus identified, and a 'log.txt' file, containing log information.

Look at the output vcf file. What does the format look like? Does that match with what you observed in the IGV? Download the VCF file to the IGV browser.

<!-- TODO: Update file names -->
```bash
    less -S AllVariants.vcf
```

You will be using this format further in the course, for now let's just count the number of heterozygous SNPs in each individual:

<!-- TODO: Update file names -->
```bash
    grep -o '0/1\|1/0' AllVariants.vcf  | wc -l
```

-   0/0 - the sample is homozygous to the reference (note that these sites usually won't be listed in single sample vcf files as they are not variants)
-   0/1 OR 0/1 - the sample is heterozygous, carrying 1 copy of each of the REF
    and ALT alleles
-   1/1 - the sample is homozygous for the alternate allele

Given this information you are now able to estimate the mean SNP heterozygosity for your individual of the 10 MB region in chromosome 2.

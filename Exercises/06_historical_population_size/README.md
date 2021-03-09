# Inferring historical populations sizes using PSMC

The Pairwise Sequentially Markovian Coalescent (PSMC) model uses information in the complete diploid sequence of a single individual to infer the history of population size changes. The method was published in 2011 ([Li and Durbin 2011](https://www.nature.com/articles/nature10231)) in the paper that you discussed in class. It has become a very popular tool in the world of genomics. In this exercise, we first walk through the steps to generate the necessary input data for PSMC. Then we run PSMC on chromosome 2 of an individual from the Simons Diversity Panel and plot the results.

For additional detail on how to run PSMC see the [GitHub page](https://github.com/lh3/psmc) for PSMC source code.

The bam files and reference genome necessary to run the following scripts can be found at: `/home/Data'.

The method you used for base calling in an earlier exercise is state of the art. Unfortunately, to produce the input data for PSMC we cannot just use the base calls or VCF files that we already produced. The first reason is that PSMC required more data than the 10Mb of chromosome 2 that you called bases on. The second reason is that the way you did your base calls do not let us easily produce input data for PSMC, which is a consensus sequence that represents the diploid genome. 

The example individual used below is a Hungarian individual with id ERR1025630. You can replace that to run the same analysis on another individual.

## Calling consensus sequence
Starting from mapped reads, the first step is to produce a consensus sequence in FASTQ format, which stores both the sequence and its corresponding quality scores, that will be used for QC filtering. The consensus sequence has A, T, C or G at homozygous sites, and other letters [IUPAC codes](https://www.bioinformatics.org/sms/iupac.html) to represent heterozygotes. To make the consensus calls, we use the samtools/bcftools suite. We first use `samtools mpileup` to take the mapped reads and produce a VCF file. We then generate a consensus sequence with `bcftools`, which we convert to FASTQ (with some additional filtering) by `vcfutils.pl`. We take advantage of Unix pipes and the ability of `samtools` to work with streaming input and output to run the whole pipeline (`samtools` -> `bcftools` -> `vcfutils.pl`) as one command. We run our consensus calling pipeline, consisting of a linked set of `samtools`, `bcftools`, and `vcfutils.pl` commands:

```bash
samtools mpileup -Q 30 -q 30 -u -v -f ~/populationgenomics/data/fasta/chr2.fa -r 2 ~/populationgenomics/data/bam/S_Ami-1.chr2.bam | ~/populationgenomics/software/bcftools call -c | ~/populationgenomics/software/vcfutils.pl vcf2fq -d 5 -D 100 -Q 30 > S_Ami-1.chr2.fq
```

The command takes as input an aligned bam file and a reference genome, generates a summary of the coverage of mapped reads on a reference sequence at a single base pair resolution using `samtools mpileup`, then calls the consensus sequence with `bcftools`, and then filters and converts the consensus to FASTQ format. Some parameter explanations:

1. `samtools`:
    - `-Q` and `-q` in mpileup determine the cutoffs for baseQ and mapQ, respectively
    - `-v` tells mpileup to produce vcf output, and `-u` says that should be uncompressed
    - `-f` is the reference fasta used 
    - `-r` is the region to call the mpileup for (in this case, a particular chromosome)
    - `S_Ami-1.chr2.bam` is the bam file to use
2. `bcftools`:
    - call `-c` calls a consensus sequence from the mpileup using the original calling method
3. `vcfutils.pl`:
    - `-d 5` and `-D 100` determine the minimum and maximum coverage to allow for `vcf2fq`, anything outside that range is filtered
    - `-Q 30` sets the root mean squared mapping quality minimum to 30

This takes a long to run (about 5-6 hours) so if you get tired of waiting you can get it here:

    ~/populationgenomics/data/consensus_fastq/S_Ami-1.chr2.fq

There you can also find FASTQ files for all the other individuals we have been working with.

## Creating a PSMC input file
PSMC takes the consensus FASTQ file, and infers the history of population sizes, but first we need to convert this FASTQ file to the input format for PSMC:

```bash
    fq2psmcfa -q20 S_Ami-1.chr2.fq > S_Ami-1.chr2.psmcfa
```

This transforms the consensus sequence into a fasta-like format where the i-th character in the output sequence indicates whether there is at least one heterozygote in the bin [100i, 100i+100). Have a look at the file using `less`.

## Running PSMC

Now we are finally ready to run PSMC. You do that like this:

```bash
    psmc -N50 -t15 -r5 -p "4+25*2+4+6" -o S_Ami-1.chr2.psmc S_Ami-1.chr2.psmcfa
```

The comman d line in the example above has been shown to be suitable for modern humans, inappropiate settings might lead to under/over-fitting. The `-p` and `-t` options are used to specify the length and number of time intervals. The `-r` option is used to specify the initial theta/rho ratio. The `-N` option sets the maximum number of EM iterations in the fitting of model parameters.

This PSMC analysis takes about 25 minutes to complete. 

## Plot your results

When the PSMC completes you can make the PSMC plot. You have to specify the per-generation mutation rate using `-u` and the generation time in years using `-g`. To make the plotting script work must first run the following command so the plotting routine knows where to find a file it needs:

	export GNUPLOT_PS_DIR=~/anaconda3/envs/popgen/share/gnuplot/5.0/PostScript

Then you can generate the plot like this:

```bash
    psmc_plot.pl -R -u 1.2e-08 -g 25 -p S_Ami-1_plot S_Ami-1.chr2.psmc
```

The `-u` option specifies the per year mutation rate and the `-g` the generation time. The `-p` option makes the script produce pdf format. The `-R` option preserves the intermediate files the script produces. The latter is handy if you want to make plots yourself combining several PSMC analyses. The two arguments are the base name for the output files and the input file with the psmc data to plot.

Does the plot resemble the ones in Li and Durbin. We used a different (more correct) mutation rate than Li et al.. What do you think that does to the plot?

## Compare individuals from different regions of the world

Now compare individuals from different regions. You can find all individuals at:

<!-- TODO: Update file names -->
```bash
/home/Data/consensus_files/PSMC
```

Do they look different? If so, what do you think could explain this?

Try to make a plot that show all your results together. Try out the code below:

<!-- TODO: Updata file names -->
```R
psmc_data1 <- read.table("ERR1025630_sort_dedup_consensus_plot.0.txt", header=F)
# you can read data from more than one psmc run:
# psmc_data2 <- read.table( ...
# psmc_data3 <- read.table( ...

plot(c(log(1e04), log(1e7)), c(0,3), type="n",  xaxt="n", bty="n", main="Results of PSMC", xlab="Years", ylab='Effective population size', las=1)
with(psmc_data1, lines(log(V1), V2, type="S", lw=1.5, col='red'))
# with(psmc_data2, lines(log(V1), V2, type="S", lw=1.5, col='green'))
# with(psmc_data3, lines(log(V1), V2, type="S", lw=1.5, col='blue'))
axis(side=1, at=log((c(1:9*1e4, 1:9*1e5, 1:9*1e6, 1:9*1e7))), labels=F) 
labs <- c(expression("10"^"4"), expression("10"^"5"), expression("10"^"6"), expression("10"^"7"))
text(x=log(c(1e4, 1e5, 1e6, 1e7)), y=1.5*par("usr")[3], pos=1, adj=1, labels = labs, xpd = TRUE)

## ggplot version
library(ggplot2)
psmc_data1 <- read.table("x_sort_dedup_consensus_plot.0.txt", header=F, col.names = c('Years', 'Effective_pop_size', 'X', 'Y', 'C'))
psmc_data2 <- read.table("y_sort_dedup_consensus_plot.0.txt", header =F, col.names = c('Years', 'Effective_pop_size', 'X', 'Y', 'C'))


# If data 1 is African and data2 is European you can type: 

psmc_data1$type = 'African'
psmc_data2$type = 'European'

df3 = data.frame(type=c(psmc_data1$type, psmc_data2$type), Years = c(psmc_data1$Years, psmc_data2$Years), 'Effective_pop_size'=c(psmc_data1$Effective_pop_size,psmc_data2$Effective_pop_size))
# Ggplot version
library(ggplot2)

g <- ggplot(df3, aes(x=Years, y=Effective_pop_size, color="NCL-08")) + geom_line(aes(color=type), size=1.5) + 
  theme_bw() + 
  labs(x= expression(paste("Years (g=25, ", mu, "=2,5*", 10^-8,")")), y="Effective population size", title='Results of PSMC') +
  scale_x_log10(breaks=c(1000, 10000, 100000, 1000000), minor_breaks=c(500, 5000, 50000, 500000)) +
  scale_y_continuous(limits = c(0,3))
g

```

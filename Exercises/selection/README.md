
# Inference of positive selection

In this exercise we will work on inference of selection using the Relate tool that we also used in the [exercise about tree sequences](https://github.com/kaspermunch/PopulationGenomicsCourse/tree/master/Exercises/tree_sequences). In that exercise we analyzed individuals from three populations and inferred all the trees along the chr2 region we have been working on. In this exercise we are going to use those trees to take the analysis a step further and look for positive selection in our chr2 region. So we pick up the [tree sequences exercise](https://github.com/kaspermunch/PopulationGenomicsCourse/tree/master/Exercises/tree_sequences) where we left it and you should continue this exericse in the same folder. That way you have access to the files you you produced the other exercise.

## Request an interactive session on a compute node:

Begin by requesting an interactive job:

```
srun --mem-per-cpu=5g --time=3:00:00 --account=populationgenomics --pty bash
```

I will take a while to finish, so just leave the terminal for now and go on with the exercise.

## Build the environment for the exercise:

For this you do not need to wait for `srun` to return. Just open another terminal on the cluster and run this command: 

```
conda env create -f ~/populationgenomics/env/pg-relate.yml
```

I will take a while to finish, so just leave the terminal for now and go on with the exercise.

## Know what you are doing

Since we already have the trees along our genomic alignment we can now use Relate to analyze each tree and compute the likelihood that it was shaped by positive selection. But before running the analysis you need to make sure you understand how Relate quantifies evidence of selelction and how that relates to what you know about positive selection. It just more fun when you know what is happening. So team up with one or more fellow students and make sure you understand the relevant Methods section in the [Relate paper](https://www.nature.com/articles/s41588-019-0484-x) as well as the section "A tree-based statistic for detecting positive selection" in the [supplementary note for the paper](https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-019-0484-x/MediaObjects/41588_2019_484_MOESM1_ESM.pdf).

## Running Relate on thee populations

When you are ready, you should have a command prompt from srun. Use that terminal to activate the `pg-relate` environment.

```
conda activate pg-relate
```

The Relate command below detects positive selection. At this point, we are detecting selection based on three distinct populations. 

```
~/populationgenomics/software/relate/scripts/DetectSelection/DetectSelection.sh -i popsize -m 1.25e-8 --poplabels 60_inds.txt -o selection_relate
```

Have a look at the `selection_relate.sele` produced, and see what is in there. Column 35 is the log10 p-value for selection. Run slurm-jupyter with this command on your own computer:

```
slurm-jupyter -e pg-relate -A populationgenomics -m 8g -t 3h -u username
```

and open an R notebook (click big blue `+` button top left and click R notebook). Then this R code to to plot the p-values across the chr2 region:

```R
library(ggplot2)
options(repr.plot.width=14, repr.plot.height=5)
relate_data <- read.table("selection_relate.sele", header=1)#, col.names = c('Years', 'Effective_pop_size', 'X', 'Y', 'C'))
ggplot(relate_data, aes(x=pos, y=-when_mutation_has_freq2, color=when_mutation_has_freq2)) + 
  geom_point(size=1.5) + 
  scale_x_continuous(limits = c(130000000,145000000)) +
  scale_y_continuous(limits = c(0,7)) +
  theme_bw()
  ```

Now try to zoom in on the region 135000000-136000000 by changing the arguments to scale_x_continuous. What do you see? It is impossible to read the actual coordinates of each SNP, but you can quickly list the sites with the strongest signal of selection in the terminal: You reverse-sort the result file on column 35 (which is the one with the log10 p-value for selection) using `sort -k35r` and then pipe the output into `head -n 10` which prints the first 10 lines of output:

```
sort -k35r selection_relate.sele | head -n 10
```

If you only want columns 1, 2, and 35 you can cut those out like this:

```
sort -k35r selection_relate.sele | cut -f 1,2,35 -d ' ' | head -n 10
```

**Do you see any of the top SNPs in the interval you just looked at?**

Try plotting some of the sites with the following command. Use -o to determine name, and bp_of_interest for the position (remember to chagne both each time you run the command). 

```
~/populationgenomics/software/relate/scripts/TreeView/TreeView.sh --haps chr2.haps --sample chr2.sample --anc popsize.anc --mut popsize.mut --poplabels 60_inds.txt --years_per_gen 28 -o tree_135472847 --bp_of_interest 135472847
```

Each command produces a pdf file (`tree_135472847.pdf` in the above case) that you can view using slurm-jupyter.

**Are the tree populations affected by the sweep in the same way?**

**Do the tree populations for seperate clades? Are they expected to, sweep or not?**

**Find the name of the SNP at 135472847 using the command extracting the top SNPs. The name is in column 2 and starts with "rs".**

**Use the UCSC genome browser to look for genes near this SNP. Can you find any good candidates?**

## Analysis of only the British population (GBR)

**Consider if it using samples from different populations could cause problems for the way Relate detects selection**

In general, selection scans of this kind will perform better if it is for a single population. Try rerunning the analysis with pure brits. First link these two files to your working directory: 

```
ln -s ~/populationgenomics/data/relate_data/50_GBR_inds.txt
ln -s ~/populationgenomics/data/relate_data/chr2_GBR_phased.vcf.gz 
```

Now you can rerun all the relate commands with these two new files, including the ones you did in the [exercise about tree sequences](https://github.com/kaspermunch/PopulationGenomicsCourse/tree/master/Exercises/tree_sequences). All the other files you need are still the same.


**Do you get different or similar results?**

**Try plotting trees at the same positions in the two different analyses, and see whether the trees in the analysis of only GBR are similar to the ones including all three populations.**



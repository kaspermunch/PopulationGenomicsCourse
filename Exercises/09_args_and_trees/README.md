Background
----------

As we saw in Week 7, male X chromosomes are easy to use for haplotype-based analysis, as they are haploid and therefore does not require phasing. However, autosomes are diploid, and we would like to use haplotype-based analysis here as well. This can be done through phasing.

Multiple phasing methods are known, such as Beagle (used in Week 4) and Shapeit. If you are lucky, the dataset you are using has already been phased, which is where we start this week. We are back to looking at chr2.

The data is taken from 60 individuals in the 1000Genomes project. In this dataset, there are 20 individuals from the following populations: GBR (Brits/Scots), JPT (Japanese), YRI (Yoruban).

Relate
------

![billede](https://user-images.githubusercontent.com/47324240/158781125-b0d4af85-69dd-4d4a-b722-30da62e8c18f.png)

The documentation for Relate can be found [here](https://myersgroup.github.io/relate/).

To start off, request an interactive job.

```
srun --mem-per-cpu=5g --time=3:00:00 --account=populationgenomics --pty bash
```

Relate will produce plots with its population size and marginal tree scripts, and this requires some specific r packages. To keep the popgen environment light and snappy, let us install this in a separate environment.

```
conda create --name relate_r -c conda-forge r-base r-tidyverse r-ggplot2 r-cowplot r-gridextra
```

All the scripts can be run in this environment, so stay in relate_r for this exercise.

A specific path needs to be manually set, due to the underlying cluster structure. If you at some point get an error referring to GLIBC, run this command again.

```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/miniconda3/envs/popgen/lib
```

Create some symlinks that lead to the following files:

```
ln -s ~/populationgenomics/data/relate_data/20140520.chr2.strict_mask.fasta
ln -s ~/populationgenomics/data/relate_data/genetic_map_chr2_combined_b37.txt
ln -s ~/populationgenomics/data/relate_data/human_ancestor_2.fa
ln -s ~/populationgenomics/data/relate_data/60_inds.txt
ln -s~/populationgenomics/data/relate_data/chr2_130_145_phased.vcf.gz
```

The first file is a mask containing areas that either have abnormal read depth or has been identified to contain repetitive elements.
The second file is a recombination map.
The third file is the ancestral state of every site, based on an alignment with Gorilla, Chimp and Human genomes.
The fourth file is a metadata file detailing the population and region for each sample.
The last file is the phased genotype vcf.

Relate does not accept the standard vcf file format, but instead uses a haps/sample format. You can read up on in the Relate documentation. The authors have been so kind as to supply a script to transform it.

### Note that input files in Relate always have to be supplied without file extensions.

First, the vcf is converted to another file format (haplotype file format). If you want to know how it is structured, you can read about it [here](https://www.cog-genomics.org/plink/2.0/formats#haps).

```
~/populationgenomics/software/relate/bin/RelateFileFormats --mode ConvertFromVcf --haps chr2.haps --sample chr2.sample -i chr2_130_145_phased
```

Then, repetitive regions are masked, as well as state flipped, based on what is the inferred ancestral state.

```
~/populationgenomics/software/relate/scripts/PrepareInputFiles/PrepareInputFiles.sh --haps chr2.haps --sample chr2.sample --ancestor human_ancestor_2.fa --mask 20140520.chr2.strict_mask.fasta -o prep.chr2
```

### Q1: How many SNPs were removed due to non-matching nucleotides?
### Q2: How many were removed due to the mask?

Now, the input is fully prepared, and Relate can be run. This should take less than a minute.

```
~/populationgenomics/software/relate/bin/Relate --mode All -m 1.25e-8 -N 30000 --haps prep.chr2.haps.gz --sample prep.chr2.sample.gz --map genetic_map_chr2_combined_b37.txt -o chr2_relate
```

### Q3: How many SNPs are left per haplotype?

The lengthiest process is this step, in which population size is estimated, and the population size is re-estimates branch lengths. This takes around 20 minutes. While you are waiting, look at the ARG notebook at the bottom of this page.

```
~/populationgenomics/software/relate/scripts/EstimatePopulationSize/EstimatePopulationSize.sh -i chr2_relate -m 1.25e-8 --poplabels 60_inds.txt -o popsize --threshold 0
```

### Q4: At the end, Relate outputs estimated mutation rate and coalescence times along the region - can this tell us anything?

Selection
---------

With the coalescence rates estimated, it is possible to detect selection.

```
~/populationgenomics/software/relate/scripts/DetectSelection/DetectSelection.sh -i popsize_relate -m 1.25e-8 --poplabels 60_inds.txt -o selection_relate
```

At this point, we are detecting selection based on three distinct populations. 

### Q5: Consider whether the underlying assumptions for this make sense.

To find the sites with the strongest selection, we will sort based on col 35 (which is the one with the log10 p-value for selection).

```
sort -k35r selection_relate.sele | head -n 10
```

Try plotting some of the sites with the following command. Use -o to determine name, and bp_of_interest for the position.

```
~/populationgenomics/software/relate/scripts/TreeView/TreeView.sh --haps chr2.haps --sample chr2.sample --anc popsize.anc --mut popsize.mut --poplabels 60_inds.txt --years_per_gen 28 -o tree --bp_of_interest 14000000
```

### Q6: Is it the same populations in which the sweep seem to occur?

### Q7: The various populations are generally not fully split. Consider why Relate can end up inferring this.

### Q8: For a site which seems promising, figure out what the nearby genes code for.

In general, selection scans of this kind will perform better if it is for a single population. Try rerunning the analysis with pure brits.

Necessary files for this are:

```
ln -s ~/populationgenomics/data/relate_data/50_GBR_inds.txt
ln -s ~/populationgenomics/data/relate_data/chr2_GBR_phased.vcf.gz 
```

The rest of the files are still the same.

### Q9: Compared to the results in Q6, do you get different or similar results?

### Q10: Try plotting some of the same sites, and see whether the trees (when considering GBR) are similar.


Ancestral Recombination Graphs and Sequence Genealogies
-------------------------------------------------------


```
conda install -c conda-forge -c plotly -c kaspermunch popgen-dashboards
```

Download the notebook by right-clicking <a href="https://raw.githubusercontent.com/kaspermunch/PopulationGenomicsCourse/master/Notebooks/arg-dashboard.ipynb" download="arg-dashboard.ipynb">
this link
</a> and "choose save link as".
To run the dashboard, also download this [file:](https://github.com/kaspermunch/popgen-dashboards/blob/main/popgen_dashboards/arg_dashboard.py) and place it in the same folder.

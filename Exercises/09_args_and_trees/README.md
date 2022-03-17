# Background
------------

As we saw in Week 7, male X chromosomes are easy to use for haplotype-based analysis, as they are haploid and therefore does not require phasing. However, autosomes are diploid, and we would like to use haplotype-based analysis here as well. This can be done through phasing.

Multiple phasing methods are known, such as Beagle (used in Week 4) and Shapeit. If you are lucky, the dataset you are using has already been phased, which is where we start this week. We are back to looking at chr2.

# Relate
--------

![billede](https://user-images.githubusercontent.com/47324240/158781125-b0d4af85-69dd-4d4a-b722-30da62e8c18f.png)

The documentation for Relate can be found [here](https://myersgroup.github.io/relate/).

To start off, request an interactive job.

```
srun --mem-per-cpu=3g --time=3:00:00 --account=populationgenomics --pty bash
```

A specific path needs to be manually set, due to the underlying cluster structure.

```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/miniconda3/envs/popgen/lib
```

Create some symlinks that lead to the following files:

```
~/populationgenomics/data/relate_data/20140520.chr2.strict_mask.fasta
~/populationgenomics/data/relate_data/genetic_map_chr2_combined_b37.txt
~/populationgenomics/data/relate_data/human_ancestor_2.fa
```

The first file is a mask containing areas that either have abnormal read depth or has been identified to contain repetitive elements.
The second file is a recombination map.
The last map is the ancestral state of every site, based on an alignment with Gorilla, Chimp and Human genomes.

Relate does not accept the standard vcf file format, but instead uses a haps/sample format. You can read up on in the Relate documentation. The authors have been so kind as to supply a script to transform it.

# Note that input files in Relate always have to be supplied without file extensions.

```
~/populationgenomics/software/relate/bin/RelateFileFormats --mode ConvertFromVcf --haps chr2.haps --sample chr2.sample -i ~/populationgenomics/data/vcf/chr2_135_145_flt
```

# Ancestral Recombination Graphs and Sequence Genealogies
---------------------------------------------------------


```
conda install -c conda-forge -c plotly -c kaspermunch popgen-dashboards
```

Download the notebook by right-clicking <a href="https://raw.githubusercontent.com/kaspermunch/PopulationGenomicsCourse/master/Notebooks/arg-dashboard.ipynb" download="arg-dashboard.ipynb">
this link
</a> and "choose save link as".

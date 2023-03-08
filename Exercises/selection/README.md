
---

## Selection

With the coalescence rates estimated, it is possible to detect selection.

```
~/populationgenomics/software/relate/scripts/DetectSelection/DetectSelection.sh -i popsize -m 1.25e-8 --poplabels 60_inds.txt -o selection_relate
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



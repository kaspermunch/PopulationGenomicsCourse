# Meta data files

ls consensus_files_for_pg2018/ERR1* | awk 'match($0, /(ERR[^_]+)/, a) {print a[1]}' > samples.txt

- `sample.txt`: contains samples used in this course.
- `Simons_meta_ENArun.txt` maps between sample ids.
- `ena.ftp.pointers.txt` has ftp urls for download of bam files.

E.g. to get sample ids and corresonding urls for bam files do:

    grep -f ../metadata/samples.txt ../metadata/Simons_meta_ENArun.txt | cut -f 5 | grep -f - ../metadata/ena.ftp.pointers.txt | cut -f 4,5


# Generate a simplified metadata file for the exercises.

grep -f <(ls ~/populationgenomics/data/fastq/|awk -F "." '{print $1}') ~/populationgenomics/data/metadata/Simons_meta_ENArun.txt |awk BEGIN'{OFS="\t"; print "ID","ENA-RUN","Sample ID (SGDP)","Genetic sex assignment","Population ID","Region","Country"}{OFS="\t";print $1,$2,$7,$8,$10,$11,$12}' > ~/populationgenomics/data/metadata/Sample_meta_subset.tsv
# Meta data files

- `sample.txt`: contains samples used in this course.
- `Simons_meta_ENArun.txt` maps between sample ids.
- `ena.ftp.pointers.txt` has ftp urls for download of bam files.

E.g. to get sample ids and corresonding urls for bam files do:

    grep -f ../metadata/samples.txt ../metadata/Simons_meta_ENArun.txt | cut -f 5 | grep -f - ../metadata/ena.ftp.pointers.txt | cut -f 4,5

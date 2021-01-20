
grep -f samples.txt Simons_meta_ENArun.txt | cut -f 5 | grep -f - ena.ftp.pointers.txt | cut -f 4,5 | column -t > ftp_urls.txt


Get bam subset for chr2 regon: http://grch37.ensembl.org/Homo_sapiens/Tools/DataSlicer. (this should also work but does not: https://www.internationalgenome.org/faq/how-do-i-get-sub-section-bam-file/)

Convert bam to fastq: https://seqome.com/convert-bam-file-fastq/



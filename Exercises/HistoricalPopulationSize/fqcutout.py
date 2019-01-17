
import sys
from Bio import SeqIO

_, start, end, input_file_name, output_file_name = sys.argv
fastq_parser = SeqIO.parse(input_file_name, "fastq")
record = next(fastq_parser)
part = record[int(start):int(end)]
SeqIO.write(part, output_file_name, "fastq")


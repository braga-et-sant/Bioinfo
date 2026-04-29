from Bio import motifs
from Bio.SeqIO import parse

with open("example2.fasta") as file:
    records = list(parse(file, "fasta"))

m = motifs.create(records)
print(m)
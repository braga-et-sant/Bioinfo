from Bio.SeqIO import parse
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction, molecular_weight, six_frame_translations

file = open("example.fasta")

records = parse(file, "fasta")


for record in records:
    my_seq = record.seq
    print(f"ID: {record.id}")
    print(f"Description: {record.description}")
    print(f"Sequence Data: {record.seq}")
    print(gc_fraction(my_seq))
    print(molecular_weight(my_seq))
    print(six_frame_translations(my_seq))

coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
coding_dna.transcribe()

print(coding_dna.transcribe())
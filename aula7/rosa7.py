from Bio import SeqIO

q = 28
count = 0

for record in SeqIO.parse("example.fastq", "fastq"):
    print("PHRED scores:", record.letter_annotations["phred_quality"])
    for i in record.letter_annotations["phred_quality"]:
        if i >= q:
            count += 1

print(count)
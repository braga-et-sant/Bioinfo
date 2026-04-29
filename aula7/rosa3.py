from Bio import Entrez
from Bio import SeqIO

sampl = "FJ817486 JX069768 JX469983"

Entrez.email = "your_name@your_mail_server.com"
handle = Entrez.efetch(db="nucleotide", id=sampl, rettype="fasta")
records = list (SeqIO.parse(handle, "fasta"))
print(records[0].id)
print(len(records[-1].seq))


#Given: A genus name, followed by two dates in YYYY/M/D format.

#Return: The number of Nucleotide GenBank entries for the given genus that were published between the dates specified.

from Bio import Entrez

gename = "Anthoxanthum"
mindate = "2003/7/25"
maxdate =  "2005/12/27"

Entrez.email = "your_name@your_mail_server.com"
handle = Entrez.esearch(
    db="nucleotide",
    term=f'"{gename}"[Organism]',
    mindate=mindate,
    maxdate=maxdate,
    datetype="pdat"
)
record = Entrez.read(handle)

print(record["Count"])
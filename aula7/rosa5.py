from Bio import Entrez, SeqIO, Align
from Bio.Align import substitution_matrices


Entrez.email = "your_email@example.com"

ID1 = "JX205496.1"
ID2 = "JX469991.1"


def fetch_sequence(genbank_id):
    handle = Entrez.efetch(
        db="nucleotide",
        id=genbank_id,
        rettype="fasta",
        retmode="text"
    )

    record = SeqIO.read(handle, "fasta")
    handle.close()

    return str(record.seq)


seq1 = fetch_sequence(ID1)
seq2 = fetch_sequence(ID2)

aligner = Align.PairwiseAligner()
aligner.mode = "global"

# EMBOSS Needle settings:
# Gap opening penalty = 10
# Gap extension penalty = 1
aligner.open_gap_score = -10
aligner.extend_gap_score = -1

# Closest Biopython equivalent to EMBOSS DNAfull
aligner.substitution_matrix = substitution_matrices.load("NUC.4.4")

score = aligner.score(seq1, seq2)

print(round(score))
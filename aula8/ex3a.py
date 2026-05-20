from Bio import Entrez, SeqIO
import requests

Entrez.email = "teu_email@example.com"

distant_ncbi = [
    "NP_000509.1",
    "XP_508242",
    "NP_001257813",
    "NP_058652",
    "NP_990820",
]

close_ncbi = [
    "NP_000509",
    "NP_005359",
    "NP_067080",
]

def download_ncbi_proteins(accessions, output_file):
    with Entrez.efetch(
        db="protein",
        id=",".join(accessions),
        rettype="fasta",
        retmode="text"
    ) as handle:
        records = list(SeqIO.parse(handle, "fasta"))

    SeqIO.write(records, output_file, "fasta")
    print(f"Guardadas {len(records)} sequências em {output_file}")


def append_fasta_from_url(url, output_file):
    response = requests.get(url)
    response.raise_for_status()

    with open(output_file, "a") as f:
        f.write("\n")
        f.write(response.text)

    print(f"Sequência adicionada a {output_file}: {url}")


# 1. Homólogos distantes
download_ncbi_proteins(distant_ncbi, "distant_globins.fasta")

# 2. Homólogos próximos, parte NCBI
download_ncbi_proteins(close_ncbi, "close_globins.fasta")

# 3. Adicionar UniProt P02237
append_fasta_from_url(
    "https://rest.uniprot.org/uniprotkb/P02237.fasta",
    "close_globins.fasta"
)

# 4. Adicionar sequência da estrutura PDB 1D8U
append_fasta_from_url(
    "https://www.rcsb.org/fasta/entry/1D8U",
    "close_globins.fasta"
)
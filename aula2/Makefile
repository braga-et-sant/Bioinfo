PY = python3
CONDA_ENV = bioinfo_lab2
HUMAN_CDNA_URL = https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
HUMAN_CHRY_URL = https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.Y.fa.gz

DMEL_CDNA_URL = https://ftp.ensembl.org/pub/release-115/fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP6.54.cdna.all.fa.gz
DMEL_CHRY_URL = https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
PACKAGES= python=3.11 emboss seqkit

all: results/gc_comparison.tsv env init

init:
	mkdir -p data
	mkdir -p results

env:
	conda create -c bioconda -n $(CONDA_ENV) $(PACKAGES) -y

download_data: data/human_cdna.fa.gz data/dmel_cdna.fa.gz data/dmel_chry.fa.gz data/human_chry.fa.gz
# Local files
HUMAN_CDNA=data/human_cdna.fa.gz
HUMAN_DNA=data/human_chry.fa.gz
DMEL_CDNA=data/dmel_cdna.fa.gz
DMEL_DNA=data/dmel_chry.fa.gz

$(HUMAN_CDNA):
	mkdir -p $(dir $@)
	curl -L $(HUMAN_CDNA_URL) -o $@

# COMPLETE
$(DMEL_CDNA):
	mkdir -p $(dir $@)
	curl -L $(DMEL_CDNA_URL) -o $@

# COMPLETE
$(HUMAN_DNA):
	mkdir -p $(dir $@)
	curl -L $(HUMAN_CHRY_URL) -o $@

# COMPLETE
$(DMEL_DNA):
	mkdir -p $(dir $@)
	curl -L $(DMEL_CHRY_URL) -o $@


# =========================================================
# Exercise 5: Exploratory Analysis
# TODO: Add a target that outputs total base count and transcript count to a text file
exploratory_stats: num_nucleotides.csv num_seqs.csv num_Ns.csv

# Exercise 5.1
num_nucleotides.csv: $(HUMAN_CDNA).nnuc.csv $(DMEL_CDNA).nnuc.csv
	cat $^ > $@

# COMPLETE
%.fa.gz.nnuc.csv: %.fa.gz
	cat $^ > grep -v "^>" dmel_cdna.fa | tr -d "\n" | wc -c


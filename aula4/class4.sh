#!/usr/bin/env bash
#conda init
#conda activate bioinfo_lab4

samtools view -h \
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/alignment/HG00096.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam \
20:1000000-1000500 > subset_region.sam

#2.2
#Software: bwa
#Version: 0.5.9-r16

#2.3
#@CO	$known_indels_file(s) = ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_mapping_resources/ALL.wgs.indels_mills_devine_hg19_leftAligned_collapsed_double_hit.indels.sites.vcf.gz
#@CO	$known_indels_file(s) .= ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_mapping_resources/ALL.wgs.low_coverage_vqsr.20101123.indels.sites.vcf.gz
#@CO	$known_sites_file(s) = ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_mapping_resources/ALL.wgs.dbsnp.build135.snps.sites.vcf.gz

#2.4.1
#19
#read is paired = 1
#read is mapped in a proper pair = 2
#query read maps to reverse strand = 16

#2.4.2
samtools view -h -f 19 subset_region.sam > reverse_proper_pairs.sam
#count them
samtools view -c -f 19 subset_region.sam

#2.4.3
samtools view -b -f 8 subset_region.sam > mate_unmapped.bam

#2.5
samtools view subset_region.sam | awk '$6 ~ /[ID]/ {print length($6), $6}' | sort -k1,1nr -k2,2 | head -1

#!/bin/bash

# ================================
# Exercise 2.6
# Convert, sort and index SAM file
# ================================

# 2.6.1 Convert the SAM file to BAM format
# -b tells samtools to output BAM instead of SAM
samtools view -b subset_region.sam > subset.bam


# 2.6.2 Sort the BAM file by genomic coordinates
# Sorting is required before indexing
samtools sort subset.bam -o subset.sorted.bam


# 2.6.3 Generate an index for the sorted BAM file
# This creates subset.sorted.bam.bai
samtools index subset.sorted.bam


# ================================
# Exercise 2.7
# Generate pileup for region
# ================================

# The SAM header shows the reference genome used is hs37d5
# Download the reference genome if not already present
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz


# Uncompress the reference genome
gunzip hs37d5.fa.gz


# Index the reference genome so samtools can query it efficiently
samtools faidx hs37d5.fa


# 2.7.1 Generate pileup for chromosome 20 positions 1,000,000–1,000,100
# -f specifies the reference genome
# -r specifies the genomic region
samtools mpileup -f hs37d5.fa -r 20:1000000-1000100 subset.sorted.bam > pileup.txt


# ================================
# Optional: display the pileup
# ================================

# Show the pileup results
cat pileup.txt

# ================================
# Exercise 3.1
# Data preparation
# ================================

# 3.1.1 Download the reference genome
wget https://raw.githubusercontent.com/BenLangmead/bowtie2/master/example/reference/lambda_virus.fa -O reference.fa


# 3.1.2 Download the FASTQ reads
wget https://raw.githubusercontent.com/BenLangmead/bowtie2/master/example/reads/reads_1.fq -O reads.fq


# ================================
# Exercise 3.2
# Align reads using minimap2
# ================================

# Align reads to the reference and pipe the SAM output to samtools
# samtools sort produces a sorted BAM file
minimap2 -a reference.fa reads.fq | samtools sort -o minimap2_alignment.bam


# Index the BAM file (optional but useful)
samtools index minimap2_alignment.bam


# Count how many reads aligned (exclude unmapped reads)
samtools view -c -F 4 minimap2_alignment.bam


# ================================
# Exercise 3.3
# Align reads using bowtie2
# ================================

# First build the bowtie2 index for the reference genome
bowtie2-build reference.fa reference_index


# Align reads using bowtie2 and pipe to samtools sort
bowtie2 -x reference_index -U reads.fq | samtools sort -o bowtie2_alignment.bam


# Index the BAM file
samtools index bowtie2_alignment.bam


# Count aligned reads
samtools view -c -F 4 bowtie2_alignment.bam


# ================================
# Exercise 3.4
# Count secondary alignments
# ================================

# Secondary alignments have flag 256
# Expected result: 0 for both aligners

samtools view -c -f 256 minimap2_alignment.bam
samtools view -c -f 256 bowtie2_alignment.bam



# ============================================
# Exercise 6 - Local BLAST+ workflow
# ============================================
# This script:
# 1) installs/checks BLAST+
# 2) creates query.fasta with the two mystery reads
# 3) downloads an E. coli reference genome
# 4) builds a local BLAST database
# 5) runs blastn locally
# 6) shows the top 5 hits with the lowest E-values
# ============================================

set -euo pipefail

# --------------------------------------------
# 6.1 Check BLAST+ installation
# --------------------------------------------
# If blastn is not found, install BLAST+ into the current conda environment.
if ! command -v blastn >/dev/null 2>&1; then
    echo "blastn not found. Installing BLAST+ with conda..."
    conda install -c conda-forge -c bioconda blast -y
fi

echo "BLAST+ version:"
blastn -version

# --------------------------------------------
# 6.2 Create query FASTA file
# --------------------------------------------
cat > query.fasta << 'EOF'
>mystery_read_1
GGTAAAGCGGCAAAAAACTACCGTGAAAAGTCGGTGGATGTGGCGGGTTATGATGAACT
>mystery_read_2
GTCTCTCTGACTTCACACAGCGACACCCACTCCTCCACCTTTGACGCTGGGGCTGGCA
EOF

echo "Created query.fasta"

# --------------------------------------------
# 6.2 Download E. coli reference genome
# --------------------------------------------
# Using NCBI RefSeq assembly for E. coli K-12 MG1655
if [ ! -f ecoli.fa ]; then
    echo "Downloading E. coli reference genome..."
    wget -O ecoli.fa.gz "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz"
    gunzip -f ecoli.fa.gz
fi

echo "Reference genome ready: ecoli.fa"

# --------------------------------------------
# 6.3 Build local BLAST nucleotide database
# --------------------------------------------
# This generates database files such as:
# - ecoli_db.nhr : header data
# - ecoli_db.nin : index data
# - ecoli_db.nsq : sequence data
echo "Building local BLAST database..."
makeblastdb -in ecoli.fa -dbtype nucl -out ecoli_db

# --------------------------------------------
# 6.4 Run local BLAST search
# --------------------------------------------
# outfmt 6 = tab-delimited table:
# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
echo "Running blastn..."
blastn -query query.fasta -db ecoli_db -out blast_results.txt -outfmt 6

echo "BLAST search complete: blast_results.txt"

# --------------------------------------------
# Show the top 5 hits with the lowest E-values
# --------------------------------------------
echo
echo "Top 5 hits with lowest E-values:"
sort -k11,11g blast_results.txt | head -5

echo
echo "Top 5 hits with query id, subject id, and E-value only:"
cut -f1,2,11 blast_results.txt | sort -k3,3g | head -5

# --------------------------------------------
# Optional: quick explanation of BLAST DB files
# --------------------------------------------
echo
echo "BLAST database files created:"
ls -1 ecoli_db.n*

echo
echo ".nhr = header information"
echo ".nin = index information"
echo ".nsq = nucleotide sequence data"

# ============================================
# Exercise 6.5 - Downloading a prebuilt BLAST database
# ============================================

# 6.5.1 Locate the update_blastdb.pl script
# This script is included with BLAST+ and is used to download
# preformatted BLAST databases from NCBI.

which update_blastdb.pl


# --------------------------------------------
# 6.5.2 Download the RefSeq Select RNA database
# --------------------------------------------
# --decompress extracts the database automatically after download
# This will download about ~250MB

update_blastdb.pl --decompress refseq_select_rna


# After extraction you should see files like:
ls refseq_select_rna.*


# These files differ slightly from the manual makeblastdb output:
# Instead of .nin .nhr .nsq you will see modern BLASTv5 files like:
#   .ndb
#   .nhr
#   .nin
#   .not
#   .nsq
# They contain additional metadata and indexing information.


# --------------------------------------------
# 6.5.3 Run BLAST against the new database
# --------------------------------------------
blastn -query query.fasta \
-db refseq_select_rna \
-out blast_refseq_results.txt \
-outfmt 6


# Show the top 5 hits with the lowest E-values
sort -k11,11g blast_refseq_results.txt | head -5


# --------------------------------------------
# 6.5.4 Explanation (not code)
# --------------------------------------------
# Institutions often keep local BLAST databases because:
# - Searches run much faster locally than through the NCBI website
# - Sensitive data does not need to be uploaded externally
# - BLAST can be integrated into automated pipelines and scripts
# - Databases can be updated on a fixed schedule (e.g., weekly)
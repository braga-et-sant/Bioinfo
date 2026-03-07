#!/usr/bin/env bash
set -euo pipefail

# 1) Generate FASTA files from scratch
cat > seqA.fasta <<'EOF'
>seqA
GATCTA
EOF

cat > seqB.fasta <<'EOF'
>seqB
GAACGTA
EOF

# 2) Execute EMBOSS needle (Needleman–Wunsch global alignment) with default parameters
# -outfile writes results to the required file name.
needle -asequence seqA.fasta -bsequence seqB.fasta -auto -outfile global_alignment.txt

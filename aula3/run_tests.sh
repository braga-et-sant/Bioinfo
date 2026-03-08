#!/usr/bin/env bash

run_pair() {
    local name="$1"
    local seq1="$2"
    local seq2="$3"
    local result

    if result=$(python hamming_distance.py --seq1 "$seq1" --seq2 "$seq2" 2>/dev/null); then
        echo "$name,$result"
    else
        echo "$name,-1"
    fi
}

cat > hamming.csv << EOF
Pair A,$(run_pair "Pair A" "ATGCGTACGTAGCTA" "ATGCCTACGTAGCTA" | cut -d, -f2)
Pair B,$(run_pair "Pair B" "GAGCCTACTAACGGGAT" "CATCGTAATGACGGCCT" | cut -d, -f2)
Pair C,$(run_pair "Pair C" "ATGC" "ATG" | cut -d, -f2)
EOF
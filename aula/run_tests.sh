

result1=$(python hamming_distance.py "--seq1" "ATGCGTACGTAGCTA" "--seq2" "ATGCCTACGTAGCTA")
result2=$(python hamming_distance.py "--seq1" "GAGCCTACTAACGGGAT" "--seq2" "CATCGTAATGACGGCCT")
result3=$(python hamming_distance.py "--seq1" "ATGC" "--seq2" "ATG")
python hamming_distance.py "--seq1" "ATGCGTACGTAGCTA" "--seq2" "ATGCCTACGTAGCTA"
exit_code1=$?
python hamming_distance.py "--seq1" "GAGCCTACTAACGGGAT" "--seq2" "CATCGTAATGACGGCCT"
exit_code2=$?
python hamming_distance.py "--seq1" "ATGC" "--seq2" "ATG"
exit_code3=$?
cat << EOF >> result.csv
Pair A, $exit_code1
Pair B, $exit_code2
Pair C, $exit_code3
EOF
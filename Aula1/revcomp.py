#!/usr/bin/env python3

import sys


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    complement = {
        "A": "T", "T": "A",
        "C": "G", "G": "C",
        "N": "N",
        "R": "Y", "Y": "R",
        "S": "S", "W": "W",
        "K": "M", "M": "K",
        "B": "V", "V": "B",
        "D": "H", "H": "D"
    }
    seq = seq.upper()
    return "".join(complement.get(base, "N") for base in reversed(seq))


def test_reverse_complement():
    """Basic tests."""
    assert reverse_complement("ATGC") == "GCAT"
    assert reverse_complement("atgc") == "GCAT"
    assert reverse_complement("AANN") == "NNTT"
    print("All tests passed.")


def read_fasta(filepath: str):
    """Parse FASTA file → (header, sequence)."""
    header = None
    seq_lines = []

    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                header = line
            else:
                seq_lines.append(line)

    sequence = "".join(seq_lines)
    return header, sequence


def main():
    if len(sys.argv) != 2:
        print("Usage: python revcomp.py <fasta_file>")
        sys.exit(1)

    filepath = sys.argv[1]

    header, sequence = read_fasta(filepath)
    revcomp_seq = reverse_complement(sequence)

    print(header)
    print(revcomp_seq)


if __name__ == "__main__":
    # run tests before execution (optional but nice for submission)
    test_reverse_complement()
    main()
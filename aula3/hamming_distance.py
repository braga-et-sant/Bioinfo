import argparse
import sys

def calculate_hamming_distance(seq1, seq2):
    if seq1 is None or seq2 is None:
        raise ValueError("Error null sequence")

    if seq1 == "" or seq2 == "":
        raise ValueError("Error empty sequence")

    if len(seq1) != len(seq2):
        raise ValueError("Error sequences are not the same length")

    hd = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            hd += 1

    return hd


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--seq1", required=True)
    parser.add_argument("--seq2", required=True)
    args = parser.parse_args()

    try:
        result = calculate_hamming_distance(args.seq1, args.seq2)
        print(result)
        sys.exit(0)
    except ValueError as e:
        print(e, file=sys.stderr)
        sys.exit(1)
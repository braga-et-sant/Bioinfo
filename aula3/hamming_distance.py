import argparse

def calculate_hamming_distance(seq1, seq2):
    try:
        seq1 != ""
    except:
        raise ValueError("Error null sequence")
        exit(-1)
    if len(seq1) != len(seq2):
        raise ValueError("Error sequences are not the same length")
        exit(-1)
    if seq1 == "" and seq2 == "":
        raise ValueError("Error empty sequence")
        exit(-1)

    hd=0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            hd += 1

    print(hd)
    return(hd)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--seq1")
    parser.add_argument("--seq2")
    args = parser.parse_args()

    exit(calculate_hamming_distance(args.seq1, args.seq2))
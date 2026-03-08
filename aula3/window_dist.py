def get_kmers(seq, k):
    kmers = []

    for i in range(len(seq) - k + 1):
        kmers.append(seq[i:i+k])

    return kmers

if __name__ == "__main__":
    seq = "GATCGATC"
    k = 3

    print(get_kmers(seq, k))
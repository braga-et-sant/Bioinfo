#!/usr/bin/env python3

import argparse

def score(a, b, match, mismatch):
    return match if a == b else mismatch

def needleman_wunsch(seq1, seq2, match, mismatch, gap):
    n = len(seq1)
    m = len(seq2)

    # DP matrix
    F = [[0]*(m+1) for _ in range(n+1)]

    # Traceback matrix
    traceback = [[""]*(m+1) for _ in range(n+1)]

    # Initialization
    for i in range(1, n+1):
        F[i][0] = i * gap
        traceback[i][0] = "↑"

    for j in range(1, m+1):
        F[0][j] = j * gap
        traceback[0][j] = "←"

    # Fill matrix
    for i in range(1, n+1):
        for j in range(1, m+1):
            diag = F[i-1][j-1] + score(seq1[i-1], seq2[j-1], match, mismatch)
            up   = F[i-1][j] + gap
            left = F[i][j-1] + gap

            best = max(diag, up, left)
            F[i][j] = best

            # Store traceback direction (priority: diag > up > left)
            if best == diag:
                traceback[i][j] = "↖"
            elif best == up:
                traceback[i][j] = "↑"
            else:
                traceback[i][j] = "←"

    return F, traceback


def traceback_alignment(seq1, seq2, traceback):
    i = len(seq1)
    j = len(seq2)

    align1 = []
    align2 = []

    while i > 0 or j > 0:
        direction = traceback[i][j]

        if direction == "↖":
            align1.append(seq1[i-1])
            align2.append(seq2[j-1])
            i -= 1
            j -= 1

        elif direction == "↑":
            align1.append(seq1[i-1])
            align2.append("-")
            i -= 1

        elif direction == "←":
            align1.append("-")
            align2.append(seq2[j-1])
            j -= 1

        else:
            break

    return "".join(reversed(align1)), "".join(reversed(align2))


def print_matrix(F, seq1, seq2):
    print("\nDP Matrix:\n")

    print("    ", "  ".join(["-"] + list(seq2)))
    for i, row in enumerate(F):
        if i == 0:
            label = "-"
        else:
            label = seq1[i-1]
        print(label, row)


def print_traceback(traceback, seq1, seq2):
    print("\nTraceback Matrix:\n")

    print("    ", "  ".join(["-"] + list(seq2)))
    for i, row in enumerate(traceback):
        if i == 0:
            label = "-"
        else:
            label = seq1[i-1]
        print(label, row)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seq1", required=True)
    parser.add_argument("--seq2", required=True)
    parser.add_argument("--match", type=int, default=1)
    parser.add_argument("--mismatch", type=int, default=-1)
    parser.add_argument("--gap", type=int, default=-1)

    args = parser.parse_args()

    F, tb = needleman_wunsch(
        args.seq1, args.seq2,
        args.match, args.mismatch, args.gap
    )

    print_matrix(F, args.seq1, args.seq2)
    print_traceback(tb, args.seq1, args.seq2)

    align1, align2 = traceback_alignment(args.seq1, args.seq2, tb)

    print("\nFinal Alignment:\n")
    print(align1)
    print(align2)

    print("\nScore:", F[len(args.seq1)][len(args.seq2)])


if __name__ == "__main__":
    main()
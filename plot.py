import matplotlib.pyplot as plt
import numpy as np
import sys
from main import get_kmer_count


def histo_kmer_distrib(filename: str, k: int):
    """Plot and save a histogram of k-mer distribution for a given sequence file and k-mer length."""
    kmer_counts = get_kmer_count(filename, k)
    count_occurrences = list(kmer_counts.values())

    plt.figure(figsize=(10, 6))
    plt.hist(
        count_occurrences,
        bins=np.arange(min(count_occurrences), max(count_occurrences) + 1, 1),
        color="skyblue",
        edgecolor="black",
    )
    plt.yscale("log")
    plt.xlabel("k-mer Count")
    plt.ylabel("Frequency (log scale)")
    plt.title(f"Distribution of k-mer Frequencies (k={k})")

    plt.savefig(f"hist_k={k}.png")
    plt.show()
    plt.close()


def main(filename: str, k: int):
    histo_kmer_distrib(
        filename,
        k,
    )


if __name__ == "__main__":

    if len(sys.argv) < 3:
        print("Usage: python3 plot.py <filename> <k>")
        sys.exit(1)

    main(sys.argv[1], int(sys.argv[2]))

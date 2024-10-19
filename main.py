"""
Python tool for finding mutations in a genome sequence using k-mer analysis and BLAST search.
BLAST hits are saved to XML files for further analysis.

Usage: python main.py <wild_filename> <resistant_filename>
"""

from collections import defaultdict
from typing import Iterable
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Blast
import sys
import time


def get_kmer_count(filename: str, k: int, stop: int = -1) -> dict[str, int]:
    """
    Counts k-mers in a given FASTA file.
    Returns a dictionary with k-mer sequences as keys and their counts as values.
    """

    def generate_kmers(sequence: str, k: int) -> Iterable[str]:
        for i in range(len(sequence) - k + 1):
            yield sequence[i : i + k]

    kmer_counts = defaultdict(int)
    with open(filename, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            for kmer in generate_kmers(str(record.seq), k):
                kmer_counts[kmer] += 1
    return kmer_counts


def remove_sequencing_errors(kmers: dict[str, int], threshold: int = 10) -> set[str]:
    """Remove k-mers that appear less than the given threshold, and the ones whose reverse complement is not in the set."""
    return set(
        kmer
        for kmer, count in kmers.items()
        if count > threshold and str(Seq(kmer).reverse_complement() in kmers)
    )


def score_kmer_similarity(kmer1: str, kmer2: str):
    """Compare two k-mers and return the similarity score and the differing bases."""
    k = len(kmer1)
    score = 0
    diff1 = ""
    diff2 = ""
    for i in range(k):
        if kmer1[i] == kmer2[i]:
            score += 1
        else:
            diff1 += kmer1[i]
            diff2 += kmer2[i]

    return (score / k, diff1, diff2)


def kmers_containing_mutation(
    wild_kmers: set[str], resistant_kmers: set[str], expected_score: float = 0.8
):
    """
    Find k-mers that are unique to each set and have a given similarity score.
    Returns a list of tuples containing the differing bases and the two k-mers:
    (differing bases, differing bases, resistant k-mer, wild k-mer)
    """
    kmers = []
    for resistant_kmer in resistant_kmers:
        for wild_kmer in wild_kmers:
            cmp = score_kmer_similarity(resistant_kmer, wild_kmer)
            if cmp[0] == expected_score:
                kmers.append((cmp[1], cmp[2], resistant_kmer, wild_kmer))

    return kmers


def find_sequence_containing_kmer(kmer: str, filename: str, padding: int = 40):
    """
    Find a sequence containing the given k-mer in the FASTA file.
    Returns the sequence with padding around the k-mer.
    """
    kmer_length = len(kmer)

    for record in SeqIO.parse(filename, "fasta"):
        # Find the position of the k-mer in the sequence
        start_position = record.seq.find(kmer)

        # Check if kmer is found and if it is within the padding
        if (
            start_position != -1
            and start_position > padding
            and (start_position + kmer_length + padding) <= len(record.seq)
        ):
            return record.seq[
                start_position - padding : start_position + kmer_length + padding
            ]

    return ""


def BLAST_search(seq: str, outfile: str = "hits.xml"):
    """Blast the given sequence and save the results to a file."""

    # Look up the sequence
    result_stream = Blast.qblast("blastx", "nr", seq)

    # Save
    with open(outfile, "wb") as out_stream:
        out_stream.write(result_stream.read())
    result_stream.close()


def main(wild_filename, resistant_filename):

    k = 15  # Length of k-mers
    threshold = 10  # Threshold for k-mer counts
    expected_similarity = 0.8  # Expected similarity score for k-mers

    # Count k-mers in the wild and resistant genomes
    wild_kmers = get_kmer_count(wild_filename, k)
    resistant_kmers = get_kmer_count(resistant_filename, k)

    # Remove k-mers that appear less than the threshold or if their reverse complements are not in the set
    wild_kmers = remove_sequencing_errors(wild_kmers, threshold)
    resistant_kmers = remove_sequencing_errors(resistant_kmers, threshold)

    # Find k-mers that are unique to each set
    only_in_wild = wild_kmers - resistant_kmers
    only_in_resistant = resistant_kmers - wild_kmers

    # Find k-mers that contain a mutation
    kmer_differences = kmers_containing_mutation(
        only_in_wild, only_in_resistant, expected_similarity
    )

    # Find a sequence containing the k-mer in both the wild and resistant genomes
    seq_resistant = ""
    seq_wild = ""
    padding = 40

    for tuple in kmer_differences:
        if seq_resistant == "":
            seq_resistant = find_sequence_containing_kmer(
                tuple[2], resistant_filename, padding
            )
        if seq_wild == "":
            seq_wild = find_sequence_containing_kmer(tuple[3], wild_filename, padding)
        if seq_resistant != "" and seq_wild != "":
            break

    print(f"Mutation found! {kmer_differences[0][1]} became {kmer_differences[0][0]}")

    return seq_wild, seq_resistant


if __name__ == "__main__":

    if len(sys.argv) != 3:
        print("Usage: python main.py <wild_filename> <resistant_filename>")
        sys.exit(1)

    # Filenames for the wild and resistant genome reads
    wild_filename = sys.argv[1]
    resistant_filename = sys.argv[2]

    # Find the mutation
    mutation_exec_time = time.time()

    seq_wild, seq_resistant = main(wild_filename, resistant_filename)

    duration = (time.time() - mutation_exec_time) / 60  # convert seconds to minutes
    print(f"Mutation found in: {duration:.2f} minutes")

    # BLAST the sequences to find the gene
    blast_exec_time = time.time()

    BLAST_search(seq_wild, "wild_hits.xml")
    BLAST_search(seq_resistant, "resistant_hits.xml")

    duration = (time.time() - blast_exec_time) / 60  # convert seconds to minutes
    print(f"BLAST hit in: {duration:.2f} minutes")

from frequent_kmers_mismatches import hamming_ball
from count_kmers import get_kmers
from approximate_pattern_matching import HammingDistance
from utils import strip_list

def dna_contains_motif(motif, dna, d):
    k = len(motif)
    for sequence in dna:
        exists = False

        for i in range(len(sequence)-k+1):
            kmer = sequence[i:i+k]
            if HammingDistance(kmer, motif) <= d:
                exists = True
                break
        if not exists:
            return False
    return True

def motif_enumeration(dna, k, d):
    patterns = set()
    dna = dna.split(" ")

    for sequence in dna:
        for i in range(len(sequence)-k+1):
            kmer = sequence[i:i+k]
            neighborhood = hamming_ball(kmer, d, ['A','G','C','T'])
            for neighbor in list(neighborhood):
                if dna_contains_motif(neighbor, dna, d):
                    patterns.add(neighbor)
    return list(patterns)

if __name__ == "__main__":
    seq = "GAATTCCTCACCTCTAGCGCGAGGG ACCGATTTTTACCAGTGCGTGGATT CCTAGACTTTGAATTCTCTACTAAA GTGAGAGAGCTCTATGGCAAGTATT CCCTTGTCACGGATTTCGCCTTGGT CTACCGAGTAGTATTGATAGATCCC"
    k = 5
    d = 1
    print(strip_list(str(motif_enumeration(seq, k, d))))
import math
from approximate_pattern_matching import HammingDistance
from itertools import product
from count_kmers import get_kmers
from tqdm import tqdm

def count_motif_percentages(motifs):
    count = {}
    columns = []

    for i in range(len(motifs[0])):
        columns.append([motif[i] for motif in motifs])
    for i in range(len(columns)):
        count[i] = {
            'A':columns[i].count('A')/len(columns[i]),
            'C':columns[i].count('C')/len(columns[i]),
            'G':columns[i].count('G')/len(columns[i]),
            'T':columns[i].count('T')/len(columns[i])
        }
    return count

def motif_entropy(motifs):
    entropy = 0
    percents = count_motif_percentages(motifs)
    for i in range(len(percents)):
        for nucleotide in percents[i]:
            pin = percents[i][nucleotide]
            if pin != 0:
                entropy += pin * math.log2(pin)
    return -entropy

def d_pm(pattern, motifs):
    pattern = pattern.split(" ")
    score = 0

    for motif in motifs:
        _motif = motif.split(" ")
        score += HammingDistance(pattern, _motif)
    
    return score

def d_p_t(pattern, text):
    distances = []
    pattern = list(pattern)
    kmers = get_kmers(list(text), len(pattern))


    for kmer in kmers:
        distances.append(HammingDistance(pattern, kmer))
    return min(distances)

def sum_dpt(pattern, dna):
    count = 0
    for sequence in dna:
        count += d_p_t(pattern, sequence)
    return count

def possible_kmers(k):
    return list(map(''.join, product('ACGT', repeat = k)))
def median_string(dna, k):
    distance = 100000000
    for kmer in tqdm(possible_kmers(k)):
        dpd = sum_dpt(kmer, dna)
        if distance > dpd:
            distance = dpd
            median = kmer
    return median
    

if __name__ == "__main__":
    #motifs = ["TCGGGGGTTTTT", "CCGGTGACTTAC","ACGGGGATTTTC","TTGGGGACTTTT","AAGGGGACTTCC","TTGGGGACTTCC","TCGGGGATTCAT","TCGGGGATTCCT","TAGGGGAACTAC","TCGGGTATAACC"]
    #print(motif_entropy(motifs))

    dnas = "AGCGGCTCATAGTCCCTAACCGAGCGAGTAACTTCGCCTGAC ACATCTTGTGTGAGCTTACGAGACCGAGTACGCCGATACGAG GTGCCAACTACAAGTGCGAGTTAGGGAGTACGACCTAAGTGG TTGAGACGTGTTTCTTAGAGCGATTGAGTATCTGAGATGTAT TGAGTAAAAGGATCACCGCAAATCCTGCTTGAGTAGTCGTAG CGGGTTCGAGTACACCTCAGATTCCGAGCGTGTAACCTTTAC CTTGATGCGGCCCAGCTCCTTATCCGAGTAGAGTTACGGGCG ATCCGTTCGGCGGATGAGCGAGGCGCAACTAAGCACAGAGTA GAGAAAAAAACGTGTGCACCCGGTCATGCTCAGCATAGAGTA ATGGAGGCGCATCCAACGGAGTATGTCATTCGAGTATTAGAA".split(" ")
    k = 6
    print(median_string(dnas, k))

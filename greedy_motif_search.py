import numpy as np
from count_kmers import get_kmers

def most_probable_kmer(sequence, k, profile_matrix):
    kmers = get_kmers(sequence, k)
    permuted = list(np.transpose(np.array(profile_matrix), (1, 0)).tolist())
    INDEX_NUCLEOTIDE = {0:'A', 1:'C', 2:"G", 3:"T"}

    kmer_probabilities = []

    #best_index = 0
    for kmer in kmers:
        prob = get_kmer_probability(kmer, profile_matrix)
        kmer_probabilities.append(prob)
    
    index = kmer_probabilities.index(max(kmer_probabilities))
    
    return kmers[index]

def get_kmer_probability(kmer, profile_matrix):
    probabilities = []
    NUCLEOTIDE_INDEX = {"A":0, "C":1, "G":2, "T":3}

    for pos, nucleotide in enumerate(kmer):
        probability = profile_matrix[NUCLEOTIDE_INDEX[nucleotide]][pos]
        probabilities.append(probability)
    
    count = 1.0

    for probability in probabilities:
        #print(probability)c
        probability=float(probability)
        #print(type(probability), type(count))
        count = count*probability
    
    return count
    
def string_to_matrix(matrix_str, k):
    matrix_list = np.array(matrix_str.split(" "))
    return np.reshape(matrix_list, (4, k))


matrix_str = "0.229 0.229 0.313 0.205 0.241 0.289 0.301 0.193 0.229 0.205 0.253 0.253 0.253 0.301 0.241 0.325 0.349 0.193 0.193 0.217 0.301 0.229 0.229 0.217 0.229 0.229 0.169 0.205 0.181 0.277 0.313 0.373 0.229 0.277 0.301 0.241 0.289 0.241 0.277 0.265 0.229 0.241 0.193 0.217 0.241 0.289 0.217 0.289"
seq = "TTAATGAAAATTCGAAATCTAGCAAGTAGCCTTCTACTCTGAAACGCTCTCCTGATCAACGACGATTTAAGAGGACTATGTGCGGATTTGATAGGACGCTCCCTCTCCCCGGGCATGGTATGGCATGCAGTAGCGTCTCTGGGCTAATGTATCCCTGTCAGCCGGGGCCTCCTTATTCCGGTCCGGCTAGCATCGGTATATGCTTGTTAGACGGGGAGACAATCCACTAGATCATGTCGTTGATCCTGGGTAGTATACGGTTTTACCATAATGATGATACCCCAGAAGGCTCGGGATATTGCCCGGCCAGCGATTTTTTAGCCTTTGTATTTGCACAAACAACGTTGCGAGTAAATGGGGAAGCAGACTTAGCGAATTCGTGCTGCCCTGTCGTTCGGCGCTGACCAATCGTTACTGACTGCTGTGGAATACAGTCTTGTTGCGTAAGGAACATGACACATCACTGCTAGAGTCCCCTTCCTAGTGTTGTCGGGCGGTCACGCAAGATCGCCACCAGTTACTTGCCGAAAGCTGCCGCAAGCACCCTTCTTTCGCGTTGTCTGGTGATGTCCGGGGCAGGACCTGTACACTCCATAGAACGTAACTAAGGATAGCTTGCCCATTGCATCAAAGACTCATTACGGAACTAGACACGTCATCTAGCCGACAGCCCTCTTCACGACGTGGATTATTTGCATTTACGTGTGCTAACGGGGTCGACGGTGTTCATAAAGGTATGGCAATCACTTTATGGGTAAGTAAGCTCCCGCCCAGCTAATTTCCTCTATGGTACCTATAGACCCTTAAATTATTTAACCCTGCGCAATTACTATGCATACGCCTAGCAGTCACCATGAACGTCATTTATTGCGTATACAGAAATTGGTGGAACAGAGCACTAACTACATCCCACGCCGCTGTATCGTGTGTAGATGAGGACCCTCCTACCACGTCATTGGTTGCATGAGCGCGGTTCGGCGGAAATAGAACGCTTAAAGAC"
k = 12
print(most_probable_kmer(seq, k, string_to_matrix(matrix_str, k)))
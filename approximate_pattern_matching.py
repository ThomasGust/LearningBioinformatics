from utils import strip_list

def HammingDistance(p, q):
    d = 0
    for p, q in zip(p, q):
        if p!= q: 
            d += 1
    return d

def ApproximatePatternMatching(Pattern, Text, d):
    positions = [] # initializing list of positions
    for i in range(len(Text) - len(Pattern)+1):
        # and using distance < d, rather than exact matching
        if HammingDistance(Pattern, Text[i:i+len(Pattern)]) <= d:
            positions.append(i)
    return positions

def approximate_pattern_count(pattern, text, d):
    indicies = ApproximatePatternMatching(pattern, text, d)
    return len(indicies)

def first_symbol(sequence):
    sequence = list(sequence)
    return sequence[0]

def suffix(sequence):
    return "".join(sequence[1:])
    
def neighbors(pattern, d):
    NUCLEOTIDES = ['A', 'C', 'G', 'T']
    if d == 0:
        return pattern
    lp = len(pattern)
    if lp == 1:
        return NUCLEOTIDES
    neighborhood = []
    suffix_neighbors = neighbors(suffix(pattern),d)

    for sequence in suffix_neighbors:
        if HammingDistance(suffix(pattern), sequence) < d:
            for n in NUCLEOTIDES:
                neighborhood.append(n+sequence)
        else:
            neighborhood.append(first_symbol(pattern)+sequence)
    return list(set(neighborhood))

if __name__ == "__main__":
    seq1 = "CTTGAAGTGGACCTCTAGTTCCTCTACAAAGAACAGGTTGACCTGTCGCGAAG"
    seq2 = "ATGCCTTACCTAGATGCAATGACGGACGTATTCCTTTTGCCTCAACGGCTCCT"
    print(HammingDistance(seq1, seq2))
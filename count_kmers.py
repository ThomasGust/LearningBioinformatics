
def get_kmers(sequence, k):
    return [sequence[i:i+k] for i in range(len(sequence)-k+1)]

def count_kmers(sequence, pattern):
    k = len(list(pattern))
    kmers = get_kmers(sequence, k)
    return kmers.count(pattern)

if __name__ == "__main__":
    c = count_kmers("ACACGAT", "AC")
    print(c)

    sequence = "GACCATCAAAACTGATAAACTACTTAAAAATCAGT"
    pattern = "AAA"
    c = count_kmers(sequence, pattern)
    print(c)
from count_kmers import get_kmers
def get_kmers(sequence, k):
    return [sequence[i:i+k] for i in range(len(sequence)-k+1)]

def frequency_table(sequence, k):
    freq_map = {}
    n = len(sequence)
    kmers = get_kmers(sequence, k)
    for i in range(n-k):
        pattern = kmers[i]
        if pattern not in list(freq_map.keys()):
            freq_map[pattern] = 1
        else:
            freq_map[pattern] = freq_map[pattern]+1
    
    return freq_map

def max_freq_table(freq_map):
    return max(freq_map.values())

def most_frequent_kmers(sequence, k):
    frequent_patterns = []
    freq_map = frequency_table(sequence, k)
    m = max_freq_table(freq_map)

    for pattern in list(freq_map.keys()):
        if freq_map[pattern] == m:
            frequent_patterns.append(pattern)
    
    return frequent_patterns

if __name__ == "__main__":
    seq = "CGGAGGACTCTAGGTAACGCTTATCAGGTCCATAGGACATTCA"
    k = 3
    mfk = str(most_frequent_kmers(seq, k))
    print(mfk.replace('[', "").replace("]","").replace("'","").replace(",",""))
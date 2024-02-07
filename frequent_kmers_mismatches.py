from itertools import chain, combinations, product
from utils import strip_list
from tqdm import tqdm
from collections import defaultdict
from reverse_complement import reverse_complement

def get_kmers(sequence, k):
    return [sequence[i:i+k] for i in range(len(sequence)-k+1)]

def hamming_circle(s, n, alphabet):
    """Generate strings over alphabet whose Hamming distance from s is
    exactly n.

    >>> sorted(hamming_circle('abc', 0, 'abc'))
    ['abc']
    >>> sorted(hamming_circle('abc', 1, 'abc'))
    ['aac', 'aba', 'abb', 'acc', 'bbc', 'cbc']
    >>> sorted(hamming_circle('aaa', 2, 'ab'))
    ['abb', 'bab', 'bba']

    """
    for positions in combinations(range(len(s)), n):
        for replacements in product(range(len(alphabet) - 1), repeat=n):
            cousin = list(s)
            for p, r in zip(positions, replacements):
                if cousin[p] == alphabet[r]:
                    cousin[p] = alphabet[-1]
                else:
                    cousin[p] = alphabet[r]
            yield ''.join(cousin)

def hamming_ball(s, n, alphabet):
    """Generate strings over alphabet whose Hamming distance from s is
    less than or equal to n.

    >>> sorted(hamming_ball('abc', 0, 'abc'))
    ['abc']
    >>> sorted(hamming_ball('abc', 1, 'abc'))
    ['aac', 'aba', 'abb', 'abc', 'acc', 'bbc', 'cbc']
    >>> sorted(hamming_ball('aaa', 2, 'ab'))
    ['aaa', 'aab', 'aba', 'abb', 'baa', 'bab', 'bba']

    """
    return chain.from_iterable(hamming_circle(s, i, alphabet)
                               for i in range(n + 1))

def frequency_table_reverse_complement(sequence, hamming_threshold, k):
    NUCLEOTIDES = ['A', 'C', 'G', "T"]
    freq_map = defaultdict(int)

    # Use a sliding window to get kmers
    for i in tqdm(range(len(sequence) - k + 1)):
        pattern = sequence[i:i+k]
        for neighbor in hamming_ball(pattern, hamming_threshold, NUCLEOTIDES):
            freq_map[neighbor] += 1
        #for neighbor in hamming_ball(reverse_complement(pattern), hamming_threshold, NUCLEOTIDES):
        #    freq_map[neighbor] += 1

    max_count = max(freq_map.values())
    patterns = [kmer for kmer, count in freq_map.items() if count == max_count]

    return patterns, freq_map, max_count

if __name__ == "__main__":
    #sequence = "CATGCCATTCGCATTGTCCCAGTGA"
    #hamming_threshold = 2
    #k = 3
    #_, table, _ = frequency_table_reverse_complement(sequence, hamming_threshold, k)
    #print(table["CCC"])
    print(len(list(hamming_ball("TGCAT", 2, ['A', 'G',' C', 'T']))))
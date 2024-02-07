
def strip_list(l):
    return l.replace("[","").replace("]","").replace(",","").replace("'","")

def get_substring_indices(seq, pattern):
    seq = list(seq)
    pattern = list(pattern)
    lp = len(pattern)

    indices = []

    for i in range(len(seq)-lp):
        if seq[i:i+lp] == pattern:
            indices.append(i)
    
    return indices

if __name__ == "__main__":
    pattern = "CGC"
    seq = "ATGACTTCGCTGTTACGCGC"
    indxs = get_substring_indices(seq, pattern)
    indxs = str(indxs)
    print(strip_list(indxs))
    


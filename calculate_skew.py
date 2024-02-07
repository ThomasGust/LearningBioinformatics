from utils import strip_list

def get_skew(seq):
    seq = list(seq)
    counts = [0]

    for base in seq:
        if base == "C":
            counts.append(counts[-1]-1)
        elif base == "G":
            counts.append(counts[-1]+1)
        else:
            counts.append(counts[-1])
    
    return counts

def minimum_skew(seq):
    skew = get_skew(seq)
    indxs = []
    min_skew = min(skew)

    for i, c in enumerate(skew):
        if c == min_skew:
            indxs.append(i)
    return indxs

def max_skew(seq):
    skew = get_skew(seq)
    indxs = []
    min_skew = max(skew)

    for i, c in enumerate(skew):
        if c == min_skew:
            indxs.append(i)
    return indxs

if __name__ == "__main__":
    seq = "CATTCCAGTACTTCATGATGGCGTGAAGA"
    #print(strip_list(str(get_skew(seq))))
    print(max_skew(seq))

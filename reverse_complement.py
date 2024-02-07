
def reverse_complement(seq):
    seq = list(seq)
    nucleotide_mapping = {"A":"T","T":"A","C":"G","G":"C"}
    rev_comp = [nucleotide_mapping[seq[-i-1]] for i in range(len(seq))]
    return "".join(rev_comp)

if __name__ == "__main__":
    seq = "GATTACA"
    print(reverse_complement(seq))
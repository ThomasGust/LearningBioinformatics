from frequent_kmers_mismatches import frequency_table_reverse_complement
from utils import strip_list
from calculate_skew import get_skew, minimum_skew, max_skew
import matplotlib.pyplot as plt


with open("C:\\Users\Thomas\OneDrive\Apps\Documents\Visual studio code projects\CourseraUCSDBioinformaticsSpecialization\course_1\salmonella_enterica.txt", "r") as f:
    data = "".join(f.readlines()).replace("\n", "")

#most_frequent = frequency_table_reverse_complement(data, hamming_threshold=2, k=9)
#print(strip_list(str(most_frequent)))
skew = get_skew(data)
min_skew = minimum_skew(data)
m_skew = max_skew(data)
# min_skew for this genome is only 2 and they are right next to eachother, so we can center our search around there
print(min_skew)
min_skew = int(sum(min_skew)/len(min_skew))
print(min_skew)
print()

lower_bound = min_skew - 2000
upper_bound = min_skew + 2000
lseq = list(data)

window = lseq[lower_bound:upper_bound]
#This window is were we can look for the most frequent kmers

mfk, freq_map, count = frequency_table_reverse_complement(window, hamming_threshold=1, k=9)
print(count, mfk)
"""
plot = plt.plot(range(len(skew)), skew)
plt.title("Salmonella Enterica Skew Plot")
plt.xlabel("Genome Position")
plt.ylabel("Guanine - Cytosine")
plt.show()
"""
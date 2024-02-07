from find_clumps import find_clumps
from substring_matching import strip_list

with open("C:\\Users\Thomas\OneDrive\Apps\Documents\Visual studio code projects\CourseraUCSDBioinformaticsSpecialization\course_1\ecoli.txt", "r") as f:
    genome = f.read()

k=30
l=500
t=3

clumps = find_clumps(genome, k, l, t)

print(strip_list(str(clumps)))
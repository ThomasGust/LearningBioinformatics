from substring_matching import get_substring_indices, strip_list

with open("C:\\Users\Thomas\OneDrive\Apps\Documents\Visual studio code projects\CourseraUCSDBioinformaticsSpecialization\course_1\\vibrio_cholareae.txt", "r") as f:
    genome = f.read()

pattern = "CTTGATCAT"

print(strip_list(str(get_substring_indices(genome, pattern))))
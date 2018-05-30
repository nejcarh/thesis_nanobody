import os


candidate_ids = [16, 18, 20, 22, 68, 70, 72, 73, 74, 75, 77, 79, 81, 83, 84]  # change to whatever the main script used

mut_aa_list = ["D", "E", "H"]
dir_coord = ("O", "T")
dir_num = ("1", "2", "3", "4", "5", "6")

mutations = {}  # [num][mut_aa]

for i in candidate_ids:
    mutations[i] = {}
    for aa in mut_aa_list:
        mutations[i][aa] = 0

for d1 in dir_coord:
    for d2 in dir_num:
        path = ".//output/NUM/" + d1 + "/" + d2 + "/"

        if os.path.exists(path):
            for folder in next(os.walk(path))[1]:

                split_name = folder.split()
                for i in range(0, len(split_name) - 1):
                    num = int(split_name[i][1:3])
                    mut_aa = split_name[i][3]
                    count = int(split_name[len(split_name) - 1].replace(")", "").replace("(", ""))
                    mutations[num][mut_aa] += count

f = open("output_distribution.txt", "w")

f.write("  |")
for num in candidate_ids:
    f.write(str(num) + "|")
f.write("\n")

for mut_aa in mut_aa_list:
    f.write(mut_aa + " |")
    for num in candidate_ids:
        f.write(str(mutations[num][mut_aa]) + "|")
    f.write("\n")

f.close()

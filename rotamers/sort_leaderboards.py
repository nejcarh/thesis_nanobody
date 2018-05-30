import os
import glob
from operator import itemgetter

dir_coord = ("O", "T")
dir_num_oct = ("3", "4", "5")
dir_num_tet = ("3", "4")

mutations = []

#  Tuple structure: ( NUM(REV), VARIANTS(REV), AVG_CLASHES, MIN_CLASHES, AVG_NUM_H2O_CLASHES, MIN_H2O_CLASHES, NUM_HIS(REV), NAME)

for d in dir_num_oct:
    path = ".//output/NUM/O/" + d + "/"  # make sure it matched whatever the main script did

    if os.path.exists(path):
        for folder in next(os.walk(path))[1]:

            split = folder.split()

            num = d    # DONE
            variants = int(split[len(split) - 1].replace("(", "").replace(")", "")) # DONE
            avg_num_clashes = 0     # DONE
            min_num_clashes = 999   # DONE
            avg_num_h2o_clashes = 0
            min_num_h2o_clashes = 999
            num_his = 0
            name = folder[:-4]

            fullpath = os.path.join(path, folder)

            for log in glob.glob1(fullpath, "*.txt"):  # CHECK ALL VARIANT LOGS

                f = open(os.path.join(fullpath, log), "r")
                log_lines = f.readlines()

                firstline_split = log_lines[0].split()
                num_clash_variant = int(firstline_split[len(firstline_split) - 1])
                avg_num_clashes += num_clash_variant / variants
                if num_clash_variant < min_num_clashes:
                    min_num_clashes = num_clash_variant

                variant_num_h2o_clashes = 0

                for line in log_lines:
                    if "HOH" in line:
                        variant_num_h2o_clashes += 1

                avg_num_h2o_clashes += variant_num_h2o_clashes / variants
                if variant_num_h2o_clashes < min_num_h2o_clashes:
                    min_num_h2o_clashes = variant_num_h2o_clashes

            for i in range(0, len(split) - 1):
                if split[i][3] == "H":
                    num_his += 1

            avg_num_clashes = round(avg_num_clashes, 1)
            avg_num_h2o_clashes = round(avg_num_h2o_clashes, 1)

            mutations.append((num, variants, avg_num_clashes, min_num_clashes, avg_num_h2o_clashes, min_num_h2o_clashes, num_his, name))


mutations = sorted(mutations, key=itemgetter(6), reverse=True)
mutations = sorted(mutations, key=itemgetter(5), reverse=False)
mutations = sorted(mutations, key=itemgetter(4), reverse=False)
mutations = sorted(mutations, key=itemgetter(3), reverse=False)
mutations = sorted(mutations, key=itemgetter(2), reverse=False)
mutations = sorted(mutations, key=itemgetter(1), reverse=True)
mutations = sorted(mutations, key=itemgetter(0), reverse=True)

f = open("output_O.txt", "w")

for i in mutations:
    f.write(str(i[0]) + "|")
    f.write(str(i[7]) + "|")
    f.write(str(i[1]) + "|")
    f.write(str(i[2]) + "|")
    f.write(str(i[3]) + "|")
    f.write(str(i[4]) + "|")
    f.write(str(i[5]) + "|")
    f.write(str(i[6]) + "\n")

f.close()


mutations = []

#  Tuple structure: ( NUM(REV), VARIANTS(REV), AVG_CLASHES, MIN_CLASHES, AVG_NUM_H2O_CLASHES, MIN_H2O_CLASHES, NUM_HIS(REV), NAME)

for d in dir_num_tet:
    path = ".//output/NUM/T/" + d + "/"

    if os.path.exists(path):
        for folder in next(os.walk(path))[1]:

            split = folder.split()

            num = d    # DONE
            variants = int(split[len(split) - 1].replace("(", "").replace(")", "")) # DONE
            avg_num_clashes = 0     # DONE
            min_num_clashes = 999   # DONE
            avg_num_h2o_clashes = 0
            min_num_h2o_clashes = 999
            num_his = 0
            name = folder[:-4]

            fullpath = os.path.join(path, folder)

            for log in glob.glob1(fullpath, "*.txt"):   # CHECK ALL VARIANT LOGS

                f = open(os.path.join(fullpath, log), "r")
                log_lines = f.readlines()

                firstline_split = log_lines[0].split()
                num_clash_variant = int(firstline_split[len(firstline_split) - 1])
                avg_num_clashes += num_clash_variant / variants
                if num_clash_variant < min_num_clashes:
                    min_num_clashes = num_clash_variant

                variant_num_h2o_clashes = 0

                for line in log_lines:
                    if "HOH" in line:
                        variant_num_h2o_clashes += 1

                avg_num_h2o_clashes += variant_num_h2o_clashes / variants
                if variant_num_h2o_clashes < min_num_h2o_clashes:
                    min_num_h2o_clashes = variant_num_h2o_clashes

            for i in range(0, len(split) - 1):
                if split[i][3] == "H":
                    num_his += 1

            avg_num_clashes = round(avg_num_clashes, 1)
            avg_num_h2o_clashes = round(avg_num_h2o_clashes, 1)

            mutations.append((num, variants, avg_num_clashes, min_num_clashes, avg_num_h2o_clashes, min_num_h2o_clashes, num_his, name))


mutations = sorted(mutations, key=itemgetter(6), reverse=True)
mutations = sorted(mutations, key=itemgetter(5), reverse=False)
mutations = sorted(mutations, key=itemgetter(4), reverse=False)
mutations = sorted(mutations, key=itemgetter(3), reverse=False)
mutations = sorted(mutations, key=itemgetter(2), reverse=False)
mutations = sorted(mutations, key=itemgetter(1), reverse=True)
mutations = sorted(mutations, key=itemgetter(0), reverse=True)

f = open("output_T.txt", "w")

for i in mutations:
    f.write(str(i[0]) + "|")
    f.write(str(i[7]) + "|")
    f.write(str(i[1]) + "|")
    f.write(str(i[2]) + "|")
    f.write(str(i[3]) + "|")
    f.write(str(i[4]) + "|")
    f.write(str(i[5]) + "|")
    f.write(str(i[6]) + "\n")

f.close()

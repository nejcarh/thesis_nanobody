import numpy as np
import os
import sys
import time
from Bio.PDB.DSSP import *
from os import listdir
from os.path import isfile, join


def to_key(pair):
    return (pair[1].get_id(), pair[0].get_parent().get_id())

parser = PDBParser(PERMISSIVE=1)
files = [f for f in listdir(os.path.dirname(os.path.realpath(__file__))) if isfile(join(os.path.dirname(os.path.realpath(__file__)), f))]

log_full = []
log_hits_details = []
log_hits = []

now = time.time()
searched = 0
found = 0
found_contacts = 0
skipped = 0

maxdist = 4   # in Ångstroms

# iterate trough all .pdb files in same directory

for i in files:

    searched += 1
    dssp_success = True

    if not '.pdb' in i:
        skipped += 1
        print('{} (models={})   SKIPPED'.format(structure.get_id(), len(structure)))
        log_full.append('{} (models={})   SKIPPED'.format(structure.get_id(), len(structure)))
        continue

    structure = PDBParser(QUIET=True).get_structure(i.replace('.pdb', ''), i)

    print('{} (models={})'.format(structure.get_id(), len(structure)))
    log_full.append('{} (models={})'.format(structure.get_id(), len(structure)))

# Handle DSSP failure

    try:
        dssp = DSSP(structure[0], os.path.normpath(sys.path[0] + '/' + i))
    except:
        dssp_success = False

# get matrices for current structure

    lines_symm = []

    with open(i) as f:
        content = f.readlines()
        for line in content:
            if 'SMTRY' in line:
                lines_symm.append(line)

    lines_mat = []

# extract matrices

    for l in lines_symm:
        lines_mat.append(l.split()[4:8])

    list_matrix = []

    for i in range(0, int(len(lines_symm) - 1), 3):
        list_matrix.append(np.array(lines_mat[i:(i + 3)], dtype="float64"))

    list_transformed = []

# create transformed copies of all chains

    for i3 in range(1, len(list_matrix)-1):

        t = []

        N = 0

        for chain in structure[0]:
            for res in chain:
                N += 1
                for a in res:
                    a_transformed = a.copy()
                    a_transformed.set_parent(res)

                    coord = np.array([a_transformed.get_vector()[0], a_transformed.get_vector()[1], a_transformed.get_vector()[2]], dtype="float64")
                    rotated_coord = np.dot(list_matrix[i3][0:3, 0:3], coord) + list_matrix[i3][:3, 3]

                    a_transformed.set_coord(rotated_coord)
                    t.append((a_transformed, chain))

        list_transformed.append(t)

    first_found = True

# compare transformed atoms with original structure

    for a_original in structure[0].get_atoms():

        if a_original.get_id() == 'CO':

            for i_matrix in range(1, len(list_transformed)-1):

                for tuple_transformed in list_transformed[i_matrix]:

                    if tuple_transformed[0].get_parent().get_resname() != 'HOH':
                        distance = round(a_original - tuple_transformed[0], 3)

                        if distance < maxdist:
                            found_contacts += 1

                            if first_found:
                                found += 1
                                log_hits_details.append(structure.get_id())
                                log_hits.append(structure.get_id())
                                first_found = False

                            note = ''
                            if distance < 0.5:
                                note += '   ...CLASH?'

                            if dssp_success:
                                try:
                                    struct = dssp[to_key(tuple_transformed)][2]

                                    message = 'res = {},   atom = {:>3},   dist = {:5.3f},   struct: {}  {}'.format(tuple_transformed[0].get_parent().get_resname(), tuple_transformed[0].get_id(), distance, struct, note)
                                    print(message)
                                    log_full.append(message)
                                    log_hits_details.append(message)

                                except:
                                    message = 'res = {},   atom = {:>3},   dist = {:5.3f},   struct: ?{}'.format(tuple_transformed[0].get_parent().get_resname(), tuple_transformed[0].get_id(), distance, note)
                                    print(message)
                                    log_full.append(message)
                                    log_hits_details.append(message)

                            else:
                                message = 'res = {},   atom = {:>3},   dist = {:5.3f},   struct = ?{}   ...file error?'.format(tuple_transformed[0].get_parent().get_resname(), tuple_transformed[0].get_id(), distance, note)
                                print(message)
                                log_full.append(message)
                                log_hits_details.append(message)

print('Searched all.')
elapsed = time.time() - now

# Write results to all output files

f = open('log_full.txt', 'w')
f.write('searched {} file(s), found {} contact(s) in {} candidate(s), skipped {} file(s), maxdist = {}A, t = {}s\n'.format(searched, found_contacts, found, skipped, maxdist, round(elapsed, 0)))
f.write('\n')

for line in log_full:
    f.write(line + '\n')

f.close()

f = open('log_hits_details.txt', 'w')
f.write('searched {} file(s), found {} contact(s) in {} candidate(s), skipped {} file(s), maxdist = {}A, t = {}s\n'.format(searched, found_contacts, found, skipped, maxdist, round(elapsed, 0)))
f.write('\n')

for line in log_hits_details:
    f.write(line + '\n')

f.close()

f = open('log_hits.txt', 'w')
f.write('searched {} file(s), found {} contact(s) in {} candidate(s), skipped {} file(s), maxdist = {}A, t = {}s\n'.format(searched, found_contacts, found, skipped, maxdist, round(elapsed, 0)))
f.write('\n')

for line in log_hits:
    f.write(line + '\n')

f.close()

print('Saved all.')

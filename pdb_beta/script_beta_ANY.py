import numpy as np
import os
import sys
import time
from Bio.PDB.DSSP import *
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from os import listdir
from os.path import isfile, join


def check_for(first, second):
    set1 = set([first])
    set2 = set(second)
    diff = set2 - set1
    return [first] + list(diff)


def return_atom(res, name):
    for atom in res:
        if atom.get_id() == name:
            return atom
    return 'failed'


class Defective( Exception ):
    pass

parser = PDBParser(PERMISSIVE=1)
files = [f for f in listdir(os.path.dirname(os.path.realpath(__file__))) if isfile(join(os.path.dirname(os.path.realpath(__file__)), f))]

AA = ["GLY", "ALA", "LEU", "MET", "PHE", "TRP", "LYS", "GLN", "GLU", "SER", "PRO", "VAL", "ILE", "CYS", "TYR", "HIS", "ARG", "ASN", "ASP", "THR"]

log_full = []
log_hits_details = []
log_hits = []

start_time = time.time()

found_pdb = 0
found_res = 0
found_contacts = 0
found_defective = 0
dssp_error_count = 0

# Change for different beta sheet runs
RES_NAME = 'ANY'
DISTANT_ATOM_CONDITION = (8, 13)
DSSP_CONDITION = 'E'

# iterate trough all .pdb files in same directory

for file in files:

    if '.pdb' not in file:
        print('{}   SKIPPED'.format(file))
        log_full.append('{}   SKIPPED'.format(file))
        continue

    found_pdb += 1

    structure = PDBParser(QUIET=True).get_structure(file.replace('.pdb', ''), file)

    print('{} (m={})'.format(structure.get_id(), len(structure)))
    log_full.append('{} (m={})'.format(structure.get_id(), len(structure)))

# Handle DSSP failure

    try:
        dssp = DSSP(structure[0], os.path.normpath(sys.path[0] + '/' + file))
    except:
        dssp_error_count += 1
        message = '---dssp exception! (initial)'
        print(message)
        log_full.append(message)
        continue


# get matrices for current structure

    lines_symm = []

    with open(file) as f:
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

    list_original = []
    list_transformed = []

# collect relevant residues

    for chain in structure[0]:
        for res in chain:
            try:
                struct = dssp[(chain.get_id(), res.get_id())][2]

                if struct == DSSP_CONDITION:

                    # for atom_id in check_for("CB", "CA"):
                    #     if all(atom.get_id() != atom_id for atom in res):
                    #         raise Defective

                    found_res += 1
                    list_original.append(res)

            except Defective:
                found_defective += 1
                res_info = '{}.{}.{}'.format(res.get_resname(), res.get_full_id()[3][1], res.get_full_id()[2])
                log_full.append('---defective ' + res_info)
                print('---defective ' + res_info)
                continue

            except:
                dssp_error_count += 1
                # message = '---dssp exception!' + str((chain.get_id(), res.get_id()))
                # print(message)
                # log_full.append(message)
                # log_hits_details.append(message)

# create transformed copies of the relevant residues

    for i in range(1, len(list_matrix) - 1):
        list_transformed.append([])

        for res in list_original:
            res_transformed = res.copy()
            res_transformed.set_parent(res.get_parent())
            res_transformed.transform(list_matrix[i][0:3, 0:3], list_matrix[i][:3, 3])

            list_transformed[i - 1].append(res_transformed)

# compare transformed atoms with original structure

    added_hit = False

    local_contacts = 0

    for res_original in list_original:

        for i in range(len(list_transformed) - 1):

            for res_transformed in list_transformed[i]:

                if res_transformed.get_resname() == "GLY" or res_original.get_resname() == "GLY":

                    distance_d = return_atom(res_original, "CA") - return_atom(res_transformed, "CA")

                    if distance_d < 13:

                        found_contacts += 1
                        local_contacts += 1

                        if not added_hit:
                            log_hits.append(structure.get_id())
                            log_hits_details.append(structure.get_id())
                            added_hit = True

                        res_o_info = '{}.{}.{}'.format(res_original.get_resname(), res_original.get_full_id()[3][1], res_original.get_full_id()[2])
                        res_t_info = '{}.{}.{}({})'.format(res_transformed.get_resname(), res_transformed.get_full_id()[3][1], res_transformed.get_full_id()[2], i)

                        message = '{}/{}  s: {}  d: {:5.2f}'.format(res_o_info, res_t_info, DSSP_CONDITION, distance_d)

                        print(message)
                        log_full.append(message)
                        log_hits_details.append(message)

                else:

                    distance_d = return_atom(res_original, "CA") - return_atom(res_transformed, "CA")

                    try:

                        distance_d = return_atom(res_original, "CA") - return_atom(res_transformed, "CA")
                        distance_c = return_atom(res_original, "CB") - return_atom(res_transformed, "CB")

                        if distance_c < distance_d < 13:

                            res_o_vect = return_atom(res_original, "CB").get_vector() - return_atom(res_original, "CA").get_vector()
                            res_t_vect = return_atom(res_transformed, "CB").get_vector() - return_atom(res_transformed, "CA").get_vector()

                            if res_o_vect.angle(res_t_vect) / np.pi * 180 > 60:
                                found_contacts += 1
                                local_contacts += 1

                            if not added_hit:
                                log_hits.append(structure.get_id())
                                log_hits_details.append(structure.get_id())
                                added_hit = True

                            res_o_info = '{}.{}.{}'.format(res_original.get_resname(), res_original.get_full_id()[3][1], res_original.get_full_id()[2])
                            res_t_info = '{}.{}.{}({})'.format(res_transformed.get_resname(), res_transformed.get_full_id()[3][1], res_transformed.get_full_id()[2], i)

                            message = '{}/{}  s: {}  d: {:5.2f}  a: {:.1f}°'.format(res_o_info, res_t_info, DSSP_CONDITION, distance_d, res_o_vect.angle(res_t_vect) / np.pi * 180)

                            print(message)
                            log_full.append(message)
                            log_hits_details.append(message)

                    except:

                        if distance_d < 13:

                            try:

                                res_o_vect = return_atom(res_original, "CB").get_vector() - return_atom(res_original, "CA").get_vector()
                                res_t_vect = return_atom(res_transformed, "CB").get_vector() - return_atom(res_transformed, "CA").get_vector()

                            except:
                                message = "broken  residues: " + str(res_original.get_full_id()) + str(res_transformed.get_full_id())

                                print(message)
                                log_full.append(message)
                                continue

                            if res_o_vect.angle(res_t_vect) / np.pi * 180 > 60:
                                found_contacts += 1
                                local_contacts += 1

                            if not added_hit:
                                log_hits.append(structure.get_id())
                                log_hits_details.append(structure.get_id())
                                added_hit = True

                            res_o_info = '{}.{}.{}'.format(res_original.get_resname(), res_original.get_full_id()[3][1], res_original.get_full_id()[2])
                            res_t_info = '{}.{}.{}({})'.format(res_transformed.get_resname(), res_transformed.get_full_id()[3][1], res_transformed.get_full_id()[2], i)

                            # message = '{}/{}  s: {}  d: {:5.2f}  a: {:.1f}°'.format(res_o_info, res_t_info, DSSP_CONDITION, distance_d, res_o_vect.angle(res_t_vect) / np.pi * 180)
                            message = 's: {}  d: {:5.2f}'.format(DSSP_CONDITION, distance_d)

                            print(message)
                            log_full.append(message)
                            log_hits_details.append(message)

    ln = 0

    for cn in structure[0]:
        for rs in cn:
            if rs.get_resname() in AA:
                ln += 1

    if local_contacts > 0:
        try:
            message = "N = " + str(local_contacts) + "/" + str(ln) + "   " + str(round(local_contacts / ln * 100, 1)) + "%"
        except:
            message = "N = " + str(local_contacts) + "/" + str(ln) + "   "
        print(message)
        log_hits_details.append(message)
        log_hits.append(message)
        print("")
        log_hits_details.append("")
        log_hits.append("")
        continue

print('Searched all.')
elapsed = time.time() - start_time

# Write results to all output files

run_name = RES_NAME + '-' + "beta"

print('searched {} pdb(s), found {} contact(s) among {} {}/{} residues; {} dssp failures, {} defective residues, atoms=({} / {}), t = {}s\n'.format(found_pdb, found_contacts, found_res, RES_NAME, DSSP_CONDITION, dssp_error_count, found_defective, "CA", "CB", round(elapsed, 1)))

f = open('log_full ' + run_name + '.txt', 'w')
f.write('searched {} pdb(s), found {} contact(s) among {} {}/{} residues; {} dssp failures, {} defective residues, atoms=({} / {}), t = {}s\n'.format(found_pdb, found_contacts, found_res, RES_NAME, DSSP_CONDITION, dssp_error_count, found_defective, "CA", "CB", round(elapsed, 1)))
f.write('\n')

for line in log_full:
    f.write(line + '\n')

f.close()

f = open('log_hits_details ' + run_name + '.txt', 'w')
f.write('searched {} pdb(s), found {} contact(s) among {} {}/{} residues; {} dssp failures, {} defective residues, atoms=({} / {}), t = {}s\n'.format(found_pdb, found_contacts, found_res, RES_NAME, DSSP_CONDITION, dssp_error_count, found_defective, "CA", "CB", round(elapsed, 1)))
f.write('\n')

for line in log_hits_details:
    f.write(line + '\n')

f.close()

f = open('log_hits ' + run_name + '.txt', 'w')
f.write('searched {} pdb(s), found {} contact(s) among {} {}/{} residues; {} dssp failures, {} defective residues, atoms=({} / {}), t = {}s\n'.format(found_pdb, found_contacts, found_res, RES_NAME, DSSP_CONDITION, dssp_error_count, found_defective, "CA", "CB", round(elapsed, 1)))
f.write('\n')

for line in log_hits:
    f.write(line + '\n')

f.close()

print('Saved all.')

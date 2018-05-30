import numpy as np
import time
import os
import shutil
import glob
from Bio.PDB.DSSP import *
from Bio.PDB.Atom import *
from Bio.PDB.Residue import *
from Bio.PDB.Chain import *
from Bio.PDB.Model import *
from Bio.PDB.Vector import *
from Bio.PDB.PDBIO import *

#   Important parameters:

show_only_best_set_residues = True
check_for_all_steric_clashes = True
search_for_HIS_by_keeping_all = True

search_volume_x = 40 #40
search_volume_y = 28 #28
search_volume_z = 2.5
search_spacing = 0.5
coord_radius = 2.1      # Tetrahedral ~ 2.06 (2.12/1.94) [2.33/1.85]    Octahedral ~ 2.17 (2.27/2.07) [2.73/1.74]
coord_tolerance = 0.5   # About 0.2 from crystal structure deviation and 0.2 from rotamer alignment deviation   ---> CHANGE LATER, AND/OR INTRODUCE DIMER SHIFTS
angular_tolerance_ratio = 0.15
clash_dist = 2.5
CA_CO_max = 8   # Approximate maximum for GLU-Co (~6.8 recorded)


#   Helpful lists

candidate_ids = [17, 19, 21, 23, 68, 70, 72, 73, 74, 75, 77, 79, 81, 83, 84]  # Initial: [17, 19, 21, 23, 68, 70, 72, 73, 74, 75, 77, 79, 81]
candidate_chains = ["A", "B"]
residue_names_1 = ["H", "D", "E"]  # Full list: ["C", "D", "E", "H", "K", "N", "Q", "R", "S", "T", "Y"], perhaps try HIS only?
common_atoms = ["N", "CA", "C", "O", "OXT"]
common_atoms_steric = ["N", "CA", "C", "O"]

#   Path variables

path_dict = "./input/"
path_rotamers = "./input/rotamers/PDB/"
path_structure = "./input/structure/"
path_output = "./output/"

#   Initiating global variables

now = time.time()
res_name_dict = {}      # 'ID' : 'ID_3'
res_rev_name_dict = {}  # 'ID_3' : 'ID'
coord_name_dict = {}    # 'ID' : [ 'coord_atom_ID' , 'coord_atom_ID' ]
ref_dict = {}           # 'ID' : 'atom_ID'
angle_dict = {}         # 'ID' : [ angle , tolerance ]
all_moving = {}
all_rotamers = {}       # 'FIXED_RES_CHAIN' : {'FIXED_RES_NUM' : { CANDIDATE_RES_1 : [ transformed_residues ] }}
rotamer_dict = {}       # 'ID' : [ alt_loc_1, alt_loc_2,... ] - Gets initiated in build_residues
messages = []           # [max_oct, opt_oct, max_tet, opt_tet, num_his, num_clash, [messages]]  # DEPRECATED ?
mutations = {}          # [sorted IDs] : [cobalt], [residues], num_HIS, clash_info_list, num, type, coords      # CHECK USAGE!!!

#   Initiating classes

parser = PDBParser(PERMISSIVE=1)
structure = PDBParser(QUIET=False).get_structure("dimer", path_structure + "dimer.pdb")

#   Defining functions


def to_key(pair):
    return (pair[1].get_id(), pair[0].get_parent().get_id())


def to_vector(t):
    return Vector(t[0], t[1], t[2])


def to_tuple(v):
    return (v[0], v[1], v[2])


def cross(v1, v2):
    return Vector(v1[1] * v2[2] - v1[2] * v2[1], v1[2] * v2[0] - v1[0] * v2[2], v1[0] * v2[1] - v1[1] * v2[0])


def interpret_rotamer(i):
    file = open(path_rotamers + i + ".pdb")
    lines = file.readlines()

    atom_data = []  # name, coord, bfactor, occupancy, altloc, fullname, serial, element
    # PDB: ATOM, SERIAL, NAME, ALTLOC IND., RESNAME, CHAIN, RESIDUE SEQ., INSERT?, X, Y, Z, OCCUPANCY, TEMP, SEGMENT, ELEMENT

    for line in lines:
        line = line.split()

        if line[0] == "ATOM":
            if len(line[2]) > 3:
                line_atom = [line[2][0:3], line[2][3], line[2][4:7], line[3], line[4], " ", (float(line[5]), float(line[6]), float(line[7])), float(line[8]), 0.0, " ", line[10]]
                atom_data.append(line_atom)
            else:
                if len(line[3]) == 3:
                    line_atom = [line[2], " ", line[3], line[4], line[5], " ", (float(line[6]), float(line[7]), float(line[8])), float(line[9]), 0.0, " ", line[11]]
                    atom_data.append(line_atom)
                else:
                    line_atom = [line[2], line[3][0], line[3][1:4], line[4], line[5], " ", (float(line[6]), float(line[7]), float(line[8])), float(line[9]), 0.0, " ", line[11]]
                    atom_data.append(line_atom)

                    # PDB NORMAL:                       PDB SHORTENED:
                    # 2 - ATOM NAME                     2 - ATOM NAME
                    # 3 - ALT. LOC. INDICATOR           2 - ALT. LOC. INDICATOR
                    # 3 - RESNAME                       2 - RESNAME
                    # 4 - CHAIN                         3 - CHAIN
                    # 5 - RES. SEQ. = 1                 4 - RES. SEQ. = 1
                    #   - INSERT INDICATOR = " "          - INSERT INDICATOR = " "
                    # 6 - X                             5 - X
                    # 7 - Y                             6 - Y
                    # 8 - Z                             7 - Z
                    # 9 - OCCUPANCY                     8 - OCCUPANCY
                    # 10 - TEMP = 0.0                   9 - TEMP = 0.0
                    #    - SEG. ID = "    "               - SEG. ID = "    "
                    # 11 - ELEMENT                      10 - ELEMENT

                    # ATOM_DATA_NORMAL:
                    # 0 - ATOM NAME
                    # 1 - ALT. LOC. INDICATOR
                    # 2 - RESNAME
                    # 3 - CHAIN
                    # 4 - RES. SEQ. = 1
                    # 5 - INSERT INDICATOR = " "
                    # 6 - (X, Y, Z)
                    # 7 - OCCUPANCY
                    # 8 - TEMP = 0.0
                    # 9 - SEG. ID = "    "
                    # 10 - ELEMENT

    return atom_data


def build_residues(list_blueprints, list_names):
    list_residues = []
    resname = res_rev_name_dict[list_blueprints[0][0][2]]

    rotamer_dict[resname] = []

    for list_rotamer in list_blueprints:

        new_res = Residue((" ", 279, " "), resname, " ")  # RESET TOTAL_RES_ID TO 1 ON SUBSEQUENT RUNS?
        rotamer_dict[resname].append(list_rotamer[5][1])

        for blueprint_atom in list_rotamer:
            new_atom = Atom(blueprint_atom[0], blueprint_atom[6], blueprint_atom[8], blueprint_atom[7], blueprint_atom[1], blueprint_atom[0], 1, blueprint_atom[10])
            # INIT:         NAME,              COORDINATES,       BFACTOR,           OCCUPANCY,         ALTLOC,            FULL NAME,         SERIAL NUMBER, ELEMENT
            new_atom.set_parent(new_res)
            new_res.add(new_atom)

        list_residues.append(new_res)

    return list_residues


def find_residues(data):
    list_residues = []

    for atom_data in data:
        if atom_data[1] != " ":
            if not atom_data[1] in list_residues:
                list_residues.append(atom_data[1])

    return list_residues


def assemble_residues(data, list_residues):
    res_blueprints = []

    for current_name in list_residues:
        current = []

        for i in common_atoms:
            for item in data:
                if item[0] == i:
                    current.append(item)

        for item in data:
            if item[1] == current_name:
                current.append(item)
        res_blueprints.append(current)

    return res_blueprints


def init_res_dict():
    f = open(path_dict + "res.dict")
    lines = f.readlines()

    for line in lines:
        line = line.rstrip().split(":")
        res_name_dict[line[0]] = line[1]


def init_res_rev_dict():
    f = open(path_dict + "res_rev.dict")
    lines = f.readlines()

    for line in lines:
        line = line.rstrip().split(":")
        res_rev_name_dict[line[0]] = line[1]


def init_coord_dict():
    f = open(path_dict + "coord.dict")
    lines = f.readlines()

    for line in lines:
        line = line.rstrip().split(":")
        coord_name_dict[line[0]] = []
        for i in line[1].split(","):
            coord_name_dict[line[0]].append(i)


def init_ref_dict():
    f = open(path_dict + "ref.dict")
    lines = f.readlines()

    for line in lines:
        line = line.rstrip().split(":")
        ref_dict[line[0]] = line[1]


def init_angle_dict():
    f = open(path_dict + "angle.dict")
    lines = f.readlines()

    for line in lines:
        line = line.rstrip().split(":")
        angle_dict[line[0]] = []
        angle_dict[line[0]].append(line[1])
        angle_dict[line[0]].append(line[2])


def return_fixed(chain_id, res_id):
    for model in structure:
        for chain in model:
            for res in chain:
                info = res.get_full_id()

                if info[2] == chain_id and info[3][1] == res_id:
                    return res


def return_all_moving_res(res_id):

    atom_data = interpret_rotamer(res_id)
    list_names = find_residues(atom_data)
    blueprints = assemble_residues(atom_data, list_names)
    res_list = build_residues(blueprints, list_names)

    return res_list


def init_all_moving():
    for res_id in residue_names_1:
        res_list = return_all_moving_res(res_id)
        res_name = res_name_dict[res_id]
        all_moving[res_name] = res_list


def refresh_vectors(res_moving, res_fixed):
    ca_m = to_vector(res_moving["CA"].get_coord())
    cb_m = to_vector(res_moving["CB"].get_coord())
    n_m = to_vector(res_moving["N"].get_coord())
    ca_f = to_vector(res_fixed["CA"].get_coord())
    cb_f = to_vector(res_fixed["CB"].get_coord())
    n_f = to_vector(res_fixed["N"].get_coord())

    n_ca_m = n_m - ca_m
    n_ca_f = n_f - ca_f
    cb_ca_m = cb_m - ca_m
    cb_ca_f = cb_f - ca_f

    return ca_m, cb_m, n_m, ca_f, cb_f, n_f, n_ca_m, n_ca_f, cb_ca_m, cb_ca_f


def res_translate(res, translation):
    for a in res:
        v = to_vector(a.get_coord())
        a.set_coord(v + translation)


def res_rotate(res, rot):
    for a in res:
        v = to_vector(a.get_coord())
        v = v.left_multiply(rot)
        a.set_coord(v)


def rotmat_n_ca(n_ca_m, n_ca_f):
    angle = n_ca_m.angle(n_ca_f)
    perp = to_vector(np.cross(to_tuple(n_ca_m), to_tuple(n_ca_f)))

    return rotaxis(angle, perp)


def rotmat_cb_ca(cb_ca_m, cb_ca_f, n_ca_m, n_ca_f):
    #   Correction for mobile residue
    # theta1 = cb_ca_m.angle(n_ca_m)
    # x1 = cb_ca_m.norm() * np.cos(theta1 - np.pi)
    # cb_ca_m_cor = cb_ca_m + n_ca_m.normalized().right_multiply(x1)
    # #   Correction for fixed residue
    # theta2 = cb_ca_f.angle(n_ca_f)
    # x2 = cb_ca_f.norm() * np.cos(theta2 - np.pi)
    # cb_ca_f_cor = cb_ca_f + n_ca_f.normalized().right_multiply(x2)
    #   Angle between corrected vectors

    cross1_m = cross(n_ca_m, cb_ca_m)
    cross2_m = cross(cross1_m, n_ca_m)
    cross1_f = cross(n_ca_f, cb_ca_f)
    cross2_f = cross(cross1_f, n_ca_f)

    if n_ca_m.angle(cross(cross2_m, cross2_f)) < 0.1:
        angle = cross2_m.angle(cross2_f)
    else:
        angle = - cross2_m.angle(cross2_f)

    return rotaxis(angle, n_ca_m)


# def rotmat_cb_ca_cor(cb_ca_m, cb_ca_f, n_ca_m, n_ca_f):  # NOT NEEDED SINCE ROTATION FIX
#     # #   Correction for mobile residue
#     # theta1 = cb_ca_m.angle(n_ca_m)
#     # x1 = cb_ca_m.norm() * np.cos(theta1 - np.pi)
#     # cb_ca_m_cor = cb_ca_m + n_ca_m.normalized().right_multiply(x1)
#     # #   Correction for fixed residue
#     # theta2 = cb_ca_f.angle(n_ca_f)
#     # x2 = cb_ca_f.norm() * np.cos(theta2 - np.pi)
#     # cb_ca_f_cor = cb_ca_f + n_ca_f.normalized().right_multiply(x2)
#     # #   Angle between corrected vectors
#     # angle = - cb_ca_m_cor.angle(cb_ca_f_cor)
#
#     cross1_m = cross(n_ca_m, cb_ca_m)
#     cross2_m = cross(cross1_m, n_ca_m)
#     cross1_f = cross(n_ca_f, cb_ca_f)
#     cross2_f = cross(cross1_f, n_ca_f)
#
#     angle = - 2 * cross2_m.angle(cross2_f)
#
#     print(n_ca_m.angle(cross(cross2_m, cross2_f)))
#
#     print("fixing ", np.round(angle * 180 / np.pi, 2))
#
#     return rotaxis(angle, n_ca_m)


def return_aligned(res_moving, res_fixed):
    ca_m, cb_m, n_m, ca_f, cb_f, n_f, n_ca_m, n_ca_f, cb_ca_m, cb_ca_f = refresh_vectors(res_moving, res_fixed)

    res_translate(res_moving, -ca_m)

    r1 = rotmat_n_ca(n_ca_m, n_ca_f)
    res_rotate(res_moving, r1)  # OK, ALIGNED

    ca_m, cb_m, n_m, ca_f, cb_f, n_f, n_ca_m, n_ca_f, cb_ca_m, cb_ca_f = refresh_vectors(res_moving, res_fixed)  # OK

    r2 = rotmat_cb_ca(cb_ca_m, cb_ca_f, n_ca_m, n_ca_f)
    res_rotate(res_moving, r2)

    ca_m, cb_m, n_m, ca_f, cb_f, n_f, n_ca_m, n_ca_f, cb_ca_m, cb_ca_f = refresh_vectors(res_moving, res_fixed)

    # if cb_ca_m.angle(cb_ca_f) > 0.1:  # NOT NEEDED SINCE ROTATION FIX
    #     r2 = rotmat_cb_ca_cor(cb_ca_m, cb_ca_f, n_ca_m, n_ca_f)
    #     res_rotate(res_moving, r2)
    #     ca_m, cb_m, n_m, ca_f, cb_f, n_f, n_ca_m, n_ca_f, cb_ca_m, cb_ca_f = refresh_vectors(res_moving, res_fixed)

    res_translate(res_moving, ca_f)

    return res_moving


def scale(v, a):
    return (v[0]*a, v[1]*a, v[2]*a)


def return_middle(res1, res2):
    v1 = res1["CA"].get_coord()
    v2 = res2["CA"].get_coord()
    return (v1 + v2)/2


def place_atom(element, position, resname, chain):
    global runtime_mark_id
    res_temp = Residue((" ", runtime_mark_id, " "), resname, " ")
    runtime_mark_id += 1

    atom_temp = Atom("X", position, 0.0, 1, " ", "X", runtime_mark_id, element)
    res_temp.add(atom_temp)
    chain.add(res_temp)


def place_atoms(atoms, element, resname, chain):
    global runtime_mark_id
    res_temp = Residue((" ", runtime_mark_id, " "), resname, " ")
    runtime_mark_id += 3

    n = 0
    for a in atoms:
        n += 1
        atom_temp = Atom("X" + str(n), a.get_coord(), 0.0, 1, " ", "X", runtime_mark_id, element)
        res_temp.add(atom_temp)

    chain.add(res_temp)


def place_line(element, position1, position2, resname, chain):
    global runtime_mark_id
    res_temp = Residue((" ", runtime_mark_id, " "), resname, " ")
    runtime_mark_id += 1

    atom_temp1 = Atom("X1", position1, 0.0, 1, " ", "X", runtime_mark_id, element)
    atom_temp2 = Atom("X2", position2, 0.0, 1, " ", "X", runtime_mark_id, element)
    res_temp.add(atom_temp1)
    res_temp.add(atom_temp2)
    chain.add(res_temp)


def get_base_vectors():

    x1 = return_middle(structure[0]["A"][73], structure[1]["B"][70])
    x2 = return_middle(structure[0]["A"][70], structure[1]["B"][73])
    y1 = return_middle(structure[0]["A"][75], structure[1]["B"][75])
    y2 = return_middle(structure[0]["A"][70], structure[1]["B"][70])

    base_x = x1 - x2
    base_x = to_vector(base_x)
    base_x.normalize()

    base_y = y2 - y1
    base_y = to_vector(base_y)
    base_y.normalize()

    base_z = to_vector(np.cross(to_tuple(base_y), to_tuple(base_x)))
    origin = y1 - scale(base_x, 20) - scale(base_y, 5)
    axis = y1 - scale(base_y, 5)

    return base_x, base_y, base_z, origin, axis


def idm(chain):
    if chain == "A":
        return 0
    elif chain == "B":
        return 1
    else:
        return 2


def init_all_rotamers():  # 'FIXED_RES_CHAIN' : {'FIXED_RES_NUM' : { CANDIDATE_RES_1 : [ transformed_residues ] }}
    for chain in candidate_chains:
        all_rotamers[chain] = {}

        for res_fixed_num in candidate_ids:
            all_rotamers[chain][res_fixed_num] = {}

            for rot_name_1 in residue_names_1:
                all_rotamers[chain][res_fixed_num][rot_name_1] = []

                for moving in all_moving[res_name_dict[rot_name_1]]:

                    aligned = return_aligned(moving, structure[idm(chain)][chain][res_fixed_num])
                    aligned_copy = Residue((" ", res_fixed_num, " "), res_name_dict[rot_name_1], " ")

                    #   Residue((" ", 279, " "), resname, " ")

                    for i in aligned:
                        aligned_copy.add(i.copy())

                    all_rotamers[chain][res_fixed_num][rot_name_1].append(aligned_copy)


def dist(a, b):
    return np.sqrt((a[0] - b[0])*(a[0] - b[0]) + (a[1] - b[1])*(a[1] - b[1]) + (a[2] - b[2])*(a[2] - b[2]))


def add_markers(chain):
    #   Marking the symmetry axis
    place_line("HE", axis, axis + scale(by, 1), "A", chain)
    place_line("HE", axis + scale(by, search_volume_y), axis + scale(by, search_volume_y - 1), "A", chain)
    #   Marking the origin and base vectors
    place_line("O", o, o + scale(bx, 1), "BX", chain)
    place_line("CL", o, o + scale(by, 1), "BY", chain)
    place_line("N", o, o + scale(bz, 1), "BZ", chain)


def create_cobalt(v):
    global runtime_mark_id
    #res_co = Residue((" ", runtime_mark_id, " "), "Co3", " ")
    co = Atom("CO", v, 0.0, 1, " ", "X", runtime_mark_id, "CO")
    #runtime_mark_id += 1
    #res_co.add(co)
    #new_chain_mark.add(res_co)
    return co


def evaluate_rotamer(co, aligned, coord):                                 # Improve?
    aligned_resname = res_rev_name_dict[aligned.get_resname()]
    # Check for clashes with cobalt
    for a in aligned:
        if a.get_id() not in coord_name_dict[aligned_resname]:
            temp = a - co

            if temp < clash_dist:
                return False

    # Check for angles
    v_a = to_vector(aligned[coord].get_coord())
    v_ref = to_vector(aligned[ref_dict[aligned_resname]].get_coord())
    v_co = to_vector(co.get_coord())
    co_a = v_co - v_a
    ref_a = v_ref - v_a
    angle_deg = co_a.angle(ref_a) * 180 / np.pi

    if np.abs(angle_deg - float(angle_dict[aligned_resname][0])) > float(angle_dict[aligned_resname][1]):
        return False

    if aligned_resname == "H":
        if coord == "NE2":
            v_r1 = v_a - to_vector(aligned["CD2"].get_coord())
            v_r2 = v_a - to_vector(aligned["CE1"].get_coord())
        elif coord == "ND1":
            v_r1 = v_a - to_vector(aligned["CG"].get_coord())
            v_r2 = v_a - to_vector(aligned["CE1"].get_coord())
        v_sum = v_r1 + v_r2
        angle_deg = v_sum.angle(co_a) * 180 / np.pi

        if angle_deg > 25:
            return False

    # for mod_id in [0, 1]:
    #     for res_fixed in structure[mod_id].get_residues():
    #         for a_fixed in common_atoms_steric:
    #             for a_aligned in coord_name_dict[aligned_resname]:
    #                 if dist(res_fixed[a_fixed].get_coord(), aligned[a_aligned].get_coord()) < clash_dist:
    #                     return False

    return True


def evaluate_result(list_res, co, coords):
    build_octahedron(list_res, co, coords)
    build_tetrahedron(list_res, co, coords)


def add_edge(angle_param, list_res, set_edges):
    for res_other in list_res:
        for atom_other_name in coord_name_dict[res_rev_name_dict[res_other.get_resname()]]:
            atom_other = res_other[atom_other_name]

            co_dist = atom_other - set_edges[0]
            if np.abs(co_dist - coord_radius) > coord_tolerance:
                continue

            if atom_other not in set_edges:
                candidate = True
                identifier_other = (atom_other.get_parent().get_id()[0], atom_other.get_parent().get_id()[2])

                for i in range(1, len(set_edges)):
                    a_co_edge = to_vector(set_edges[i].get_coord()) - to_vector(set_edges[0].get_coord())
                    a_co_other = to_vector(atom_other.get_coord()) - to_vector(set_edges[0].get_coord())
                    angle = a_co_edge.angle(a_co_other) * 180 / np.pi
                    identifier_set = (set_edges[i].get_parent().get_id()[0], set_edges[i].get_parent().get_id()[2])

                    if identifier_set == identifier_other:
                        candidate = False
                        break

                    if angle_param == 90:
                        if np.abs(angle - angle_param) > angle_param * angular_tolerance_ratio and np.abs(angle - angle_param * 2) > angle_param * angular_tolerance_ratio:
                            candidate = False
                            break
                    else:
                        if np.abs(angle - angle_param) > angle_param * angular_tolerance_ratio:
                            candidate = False
                            break

                if candidate:
                    set_edges.append(atom_other)
                    return True

    return False


def build_octahedron(list_res, co, coords):
    len_max = 0
    temp_candidates = []

    for res_origin in list_res:
        for atom_origin_name in coord_name_dict[res_rev_name_dict[res_origin.get_resname()]]:
            atom_origin = res_origin[atom_origin_name]

            co_dist = atom_origin - co
            if np.abs(co_dist - coord_radius) > coord_tolerance:
                continue

            set_edges = [co, atom_origin]

            while add_edge(90, list_res, set_edges):
                pass

            set_edges.remove(co)
            invalid = False

            for set_other in temp_candidates:
                if invalid:
                    break
                same = True

                for atom in set_edges:
                    for atom_other in set_other[1]:
                        if atom.get_id() == atom_other.get_id():
                            if atom.get_parent().get_resname() == atom_other.get_parent().get_resname():
                                if atom.get_parent().get_id()[0] == atom_other.get_parent().get_id()[0]:
                                    if atom - atom_other != 0:
                                        same = False
                                        break
                if same:
                    invalid = True

            if invalid:
                continue

            histidines = []
            for a in set_edges:
                a_parent = a.get_parent()
                if a_parent.get_resname() == "HIS":
                    if a_parent not in histidines:
                        histidines.append(a_parent)

            if len(set_edges) > len_max:
                len_max = len(set_edges)
                if not search_for_HIS_by_keeping_all or len(histidines) == 0:
                    del temp_candidates[:]
                temp_candidates.append([co, set_edges, len_max, len(histidines)])
            elif len(set_edges) == len_max:
                temp_candidates.append([co, set_edges, len_max, len(histidines)])

    for candidate in temp_candidates:

        if len(candidate[1]) > 1:
            invalid = True  # Check for structures coming from the same chain
            for a1 in candidate[1]:
                for a2 in candidate[1]:
                    if a1 != a2:
                        if a1.get_parent().get_full_id()[0][2] != a2.get_parent().get_full_id()[0][2]:
                            invalid = False
            if invalid:
                print("Skipping candidate set (O, {}):".format(len(candidate[1])))
                for a in candidate[1]:
                    print(a.get_parent().get_full_id())
                continue

        if len(candidate[1]) > 2:
            water_id = 3000

            positions = []
            v_co = to_vector(co.get_coord())

            for a in candidate[1]:
                positions.append(to_vector(a.get_coord()) - v_co)

            a1 = positions[0]
            a2 = 0

            for i in range(1, len(positions)):
                if np.abs(a1.angle(positions[i]) * 180 / np.pi - 90) < 15:
                    a2 = positions[i]
                    break

            a3 = v_co - a1
            a4 = v_co - a2
            a5 = v_co + cross(a1, a2).normalized().right_multiply(coord_radius)
            a6 = v_co - cross(a1, a2).normalized().right_multiply(coord_radius)

            water_positions = [a3, a4, a5, a6]

            for wt_p in water_positions:
                empty = True

                for p in positions:
                    if dist(wt_p, v_co + p) < 1:
                        empty = False

                if empty:
                    water_residue = Residue((" ", water_id, " "), "HOH", " ")
                    water = Atom("CA", wt_p, 0.0, 1, " ", "CA", water_id, "SR")
                    water.set_parent(water_residue)
                    water_residue.add(water)
                    water_id += 1
                    candidate[1].append(water)

            save_candidate(candidate[0], candidate[1], candidate[2], candidate[3], "O", coords)

        elif len(candidate[1]) == 2:
            water_id = 3500
            positions = []
            v_co = to_vector(co.get_coord())

            for a in candidate[1]:
                positions.append(to_vector(a.get_coord()) - v_co)

            if np.abs(positions[0].angle(positions[1]) * 180 / np.pi - 90) < 15:

                a1 = positions[0]
                a2 = positions[1]
                a3 = v_co - a1
                a4 = v_co - a2
                a5 = v_co + cross(a1, a2).normalized().right_multiply(coord_radius)
                a6 = v_co - cross(a1, a2).normalized().right_multiply(coord_radius)

                water_residue = Residue((" ", water_id, " "), "HOH", " ")
                water = Atom("CA", a3, 0.0, 1, " ", "CA", water_id, "SR")
                water.set_parent(water_residue)
                water_residue.add(water)
                water_id += 1
                candidate[1].append(water)

                water_residue = Residue((" ", water_id, " "), "HOH", " ")
                water = Atom("CA", a4, 0.0, 1, " ", "CA", water_id, "SR")
                water.set_parent(water_residue)
                water_residue.add(water)
                water_id += 1
                candidate[1].append(water)

                water_residue = Residue((" ", water_id, " "), "HOH", " ")
                water = Atom("CA", a5, 0.0, 1, " ", "CA", water_id, "SR")
                water.set_parent(water_residue)
                water_residue.add(water)
                water_id += 1
                candidate[1].append(water)

                water_residue = Residue((" ", water_id, " "), "HOH", " ")
                water = Atom("CA", a6, 0.0, 1, " ", "CA", water_id, "SR")
                water.set_parent(water_residue)
                water_residue.add(water)
                water_id += 1
                candidate[1].append(water)

                save_candidate(candidate[0], candidate[1], candidate[2], candidate[3], "O", coords)


def build_tetrahedron(list_res, co, coords):  # Optimize
    len_max = 0
    temp_candidates = []

    for res_origin in list_res:
        for atom_origin_name in coord_name_dict[res_rev_name_dict[res_origin.get_resname()]]:
            atom_origin = res_origin[atom_origin_name]

            co_dist = atom_origin - co
            if np.abs(co_dist - coord_radius) > coord_tolerance:
                continue

            set_edges = [co, atom_origin]

            while add_edge(120, list_res, set_edges):
                pass

            set_edges.remove(co)
            invalid = False

            for set_other in temp_candidates:
                if invalid:
                    break
                same = True

                for atom in set_edges:
                    for atom_other in set_other[1]:
                        if atom.get_id() == atom_other.get_id():
                            if atom.get_parent().get_resname() == atom_other.get_parent().get_resname():
                                if atom.get_parent().get_id()[0] == atom_other.get_parent().get_id()[0]:
                                    if atom - atom_other != 0:
                                        same = False
                                        break
                if same:
                    invalid = True

            if invalid:
                continue

            histidines = []
            for a in set_edges:
                a_parent = a.get_parent()
                if a_parent.get_resname() == "HIS":
                    if a_parent not in histidines:
                        histidines.append(a_parent)

            if len(set_edges) > len_max:
                len_max = len(set_edges)
                if not search_for_HIS_by_keeping_all or len(histidines) == 0:
                    del temp_candidates[:]
                temp_candidates.append([co, set_edges, len_max, len(histidines)])
            elif len(set_edges) == len_max:
                temp_candidates.append([co, set_edges, len_max, len(histidines)])  # SHOULD MAYBE ALSO BE UNDER CONTROL OF "HIS" SETTING? OR IMPLEMENT A FAST SEARCH? OR CLASH TOGGLE

    for candidate in temp_candidates:

        if len(candidate[1]) > 1:

            invalid = True  # Check for structures coming from the same chain
            for a1 in candidate[1]:
                for a2 in candidate[1]:
                    if a1 != a2:
                        if a1.get_parent().get_full_id()[0][2] != a2.get_parent().get_full_id()[0][2]:
                            invalid = False
            if invalid:
                print("Skipping candidate set (T):")
                for a in candidate[1]:
                    print(a.get_parent().get_full_id())
                continue

            water_id = 4000
            positions = []
            v_co = to_vector(co.get_coord())

            for a in candidate[1]:
                positions.append(to_vector(a.get_coord()) - v_co)

            a1 = positions[0]
            a2 = positions[1]
            r = rotaxis(np.pi * 2 / 3, a1)
            a3 = v_co + a2.left_multiply(r)
            r = rotaxis(-np.pi * 2 / 3, a1)
            a4 = v_co + a2.left_multiply(r)

            water_positions = [a3, a4]

            for wt_p in water_positions:
                empty = True

                for p in positions:
                    if dist(wt_p, v_co + p) < 1:
                        empty = False

                if empty:
                    water_residue = Residue((" ", water_id, " "), "HOH", " ")
                    water = Atom("CA", wt_p, 0.0, 1, " ", "CA", water_id, "SR")
                    water.set_parent(water_residue)
                    water_residue.add(water)
                    water_id += 1
                    candidate[1].append(water)

            save_candidate(candidate[0], candidate[1], candidate[2], candidate[3], "T", coords)


def return_rmsd(res1, res2):
    d = 0
    n = 0
    if res1.get_resname() == res2.get_resname():
        for a1 in res1:
            n += 1
            d += np.square(a1 - res2[a1.get_id()])
        return np.sqrt(d / n)
    else:
        return 7   # The residues are different


# mutations = {}     [sorted IDs] : [cobalt], [residues], num_HIS, clash_info_list, num, type, coords


def save_candidate(co, set, num, num_HIS, coord_type, coords):
    mut_ids = []
    res_ids = []
    conflicting = []

    for a in set:
        parent_id = a.get_parent().get_segid()
        if a.get_parent().get_resname() != "HOH":
            parent_name = res_rev_name_dict[a.get_parent().get_resname()]
        temp_id = str(parent_id) + ":" + str(parent_name)
        if parent_id not in res_ids:
            res_ids.append(parent_id)
        if temp_id not in mut_ids:
            mut_ids.append(temp_id)

    for id1 in mut_ids:
        for id2 in mut_ids:
            if id1 != id2:
                temp_id1 = str(id1).split(":")
                temp_id2 = str(id2).split(":")
                if temp_id1[0] == temp_id2[0]:
                    if temp_id1[1] != temp_id2[1]:
                        if str(id1) not in conflicting:
                            conflicting.append(str(id1))
                        if str(id2) not in conflicting:
                            conflicting.append(str(id2))

    mut_ids.sort()
    # if len(conflicting) == 0:
    ids = tuple(mut_ids)
    num_clashes = -1
    rmsd_clashes = -1

    residues = []
    for edge in set:
        if edge.get_parent() not in residues:
            residues.append(edge.get_parent())

    clashes = []
    if check_for_all_steric_clashes:
        clashes = evaluate_result_steric(residues, res_ids)  # No clash details?

    if ids in mutations.keys():
        mutations[ids].append([co, residues, num_HIS, clashes, num, coord_type, coords])
    else:
        mutations[ids] = []
        mutations[ids].append([co, residues, num_HIS, clashes, num, coord_type, coords])

    # else:
    #     print("CONFLICTING MUTATIONS - SPLITTING MODELS")
    #     print("Split to:")
    #     for conflict in conflicting:
    #         ids = copy.deepcopy(mut_ids)
    #         num_clashes = -1
    #         rmsd_clashes = -1
    #
    #         residues = []
    #         for edge in set:
    #             parent_id = edge.get_parent().get_segid()
    #             parent_name = res_rev_name_dict[edge.get_parent().get_resname()]
    #             temp_id = str(parent_id) + ":" + str(parent_name)
    #
    #             if conflict == temp_id:
    #
    #                 ids.remove(temp_id)
    #             else:
    #                 residues.append(edge.get_parent())
    #
    #         ids = tuple(ids)
    #
    #         if check_for_all_steric_clashes:
    #             num_clashes, rmsd_clashes = evaluate_result_steric(residues, res_ids)  # No clash details?
    #
    #         print(ids)
    #
    #         if ids in mutations.keys():
    #             mutations[ids].append([co, residues, num_HIS, num_clashes, rmsd_clashes, num, coord_type, coords])
    #         else:
    #             mutations[ids] = []
    #             mutations[ids].append([co, residues, num_HIS, num_clashes, rmsd_clashes, num, coord_type, coords])


def evaluate_result_steric(set_res, parent_ids):
    clashes = []

    for chain_name in ["A", "B"]:
        chain = structure[idm(chain_name)][chain_name]

        for old_res in chain:
            if old_res.get_full_id()[3][1] not in parent_ids:

                for new_res in set_res:
                    if dist(old_res["CA"].get_coord(), new_res["CA"].get_coord()) < 10:

                        for old_atom in old_res:
                            for new_atom in new_res:

                                # if old_atom.get_id() in common_atoms or new_atom.get_id() in common_atoms and new_res.get_resname() != "HOH":
                                #     continue

                                if new_atom.get_id() in common_atoms and new_res.get_resname() != "HOH":
                                    continue

                                dst = new_atom - old_atom  # Details probably not needed?
                                if dst < clash_dist:

                                    # info tuple format: (new_res, new_atom, old_res, old_atom, distance)
                                    clash_info = (new_res, new_atom.get_id(), old_res, old_atom.get_id(), np.round(dst, 2))
                                    clashes.append(clash_info)

    return clashes


def print_clash_data(fullpath, clashes):
    f = open(fullpath, "w")
    f.write("Number of clashes: {}\n".format(len(clashes)))

    rmsd = 0
    for i in clashes:
        rmsd += np.square(i[4]) / len(clashes)
    rmsd = np.round(np.sqrt(rmsd), 2)
    f.write("RMSD: {}\n".format(rmsd))
    f.write("\n")
    f.write("Clashes: (A)\n")

    for i in clashes:
        f.write("{}\tbetween: new {}{}({})\tand old {}{}({})\n".format(i[4], i[0].get_resname(), i[0].get_segid(), i[1], i[2].get_resname(), i[2].get_full_id()[3][1], i[3]))
    f.close()


#   ------------------------------------------------------------------------------------------ Main ------------------------------------------------------------------------------------------


print("Script started.")

#   Initiating larger lists
init_res_dict()
init_res_rev_dict()
init_coord_dict()
init_ref_dict()
init_angle_dict()


print("Interpreting possible rotamers.")

init_all_moving()

print("Creating all possible rotamers for all candidate residues.")

init_all_rotamers()


now = time.time()
print("Script initiated in " + str(np.round(time.time() - now, 2)) + "s")


#   Establishing a search space coordinate system
bx, by, bz, o, axis = get_base_vectors()


folder = './/output'
if os.path.exists(folder):
    print("Cleaning output path.")

    for the_file in os.listdir(folder):
        file_path = os.path.join(folder, the_file)
        try:
            if os.path.isfile(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path): shutil.rmtree(file_path)
        except Exception as e:
            print(e)

print("Starting search.")

#  Main loop: Iterating trough the search space
for x in np.arange(0, search_volume_x + search_spacing / 2, search_spacing):
    for y in np.arange(0, search_volume_y + search_spacing / 2, search_spacing):
        for z in np.arange(-search_volume_z, search_volume_z + search_spacing / 2, search_spacing):

            structure = PDBParser(QUIET=False).get_structure("dimer", path_structure + "dimer.pdb")
            message = [0, 0, 0, 0, 0, 0, []]
            messages.append(message)

            print("")
            print("Current position: " + str((x, y, z)))

            runtime_residue_id = 1000
            runtime_mark_id = 300
            candidates = []
            hit_chains = []

            cobalt = create_cobalt(o + scale(bx, x) + scale(by, y) + scale(bz, z))

            for current_chain in candidate_chains:
                for current_res_id in candidate_ids:
                    reference = structure[idm(current_chain)][current_chain][current_res_id]

                    if reference["CA"] - cobalt < CA_CO_max:    # Check only for close residues
                        for aligned_rotamer_group in all_rotamers[current_chain][current_res_id]:
                            for aligned_rotamer_resname in residue_names_1:
                                for aligned_rotamer in all_rotamers[current_chain][current_res_id][aligned_rotamer_resname]:
                                    for coord_atom in coord_name_dict[aligned_rotamer_resname]:

                                        distance = aligned_rotamer[coord_atom] - cobalt
                                        if np.abs(distance - coord_radius) < coord_tolerance:

                                            # Check for steric constraints middles of larger chains
                                            if not evaluate_rotamer(cobalt, aligned_rotamer, coord_atom):
                                                break

                                            aligned_rotamer_copy = Residue((current_res_id, runtime_residue_id, current_chain), res_name_dict[aligned_rotamer_resname], current_res_id)
                                            runtime_residue_id += 2

                                            for atom in aligned_rotamer:
                                                temp_atom = atom.copy()
                                                temp_atom.set_parent(aligned_rotamer_copy)
                                                aligned_rotamer_copy.add(temp_atom)

                                            identical = False

                                            for res in candidates:   # Check for duplication source?
                                                if return_rmsd(res, aligned_rotamer_copy) == 0:
                                                    identical = True
                                            if identical:
                                                break

                                            candidates.append(aligned_rotamer_copy)

                                            if current_chain not in hit_chains:
                                                hit_chains.append(current_chain)
                                            break

            if len(hit_chains) < 2:
                print("Hits in chain " + str(hit_chains))
                continue

            print("Hits in chains", hit_chains, "- evaluating.")

            # RATE CURRENT HITS
            evaluate_result(candidates, cobalt, (x, y, z))   # Maybe rename function?

            list_mut_names = []

            for key in mutations.keys():

                item = mutations[key]
                # print(key, item)
                mut_name_full = ""
                for i in key:
                    mut_id = str(i).split(":")[0]
                    if mut_id != " ":
                        original_name = res_rev_name_dict[structure[0]["A"][int(mut_id)].get_resname()]
                        mut_name = str(i).split(":")[1]
                        mut_name_full += original_name + mut_id + mut_name + " "
                list_mut_names.append(mut_name_full.rstrip())

            for key in mutations.keys():

                invalid = False
                item = mutations[key]
                mut_name_full = ""
                for i in key:
                    mut_id = str(i).split(":")[0]
                    if mut_id != " ":
                        original_name = res_rev_name_dict[structure[0]["A"][int(mut_id)].get_resname()]
                        mut_name = str(i).split(":")[1]
                        mut_name_full += original_name + mut_id + mut_name + " "

                mut_name_full = mut_name_full.rstrip()

                for name in list_mut_names:
                    if name != mut_name_full:
                        if mut_name_full in name:
                            invalid = True

                if invalid:
                    continue

                num_items = {"O": 0, "T": 0}
                for hit in item:
                    num_items[hit[5]] += 1

                for hit in item:

                    if len(hit[1]) > 1:

                        structure = PDBParser(QUIET=False).get_structure("dimer", path_structure + "dimer.pdb")
                        new_model = Model(2)
                        new_chain = Chain("C")
                        new_model.add(new_chain)
                        new_chain_mark = Chain("M")
                        add_markers(new_chain_mark)
                        new_model.add(new_chain_mark)
                        structure.add(new_model)

                        # mutations = {}     [sorted IDs] : [cobalt], [residues], num_HIS, clash_info_list, num, type, coords

                        res_co = Residue((" ", runtime_mark_id, " "), "Co3", " ")
                        res_co.add(hit[0])
                        new_chain.add(res_co)  # Add cobalt in its own residue

                        for res in hit[1]:  # Add rotamers
                            if res not in new_chain:
                                new_chain.add(res)

                        if hit[2] > 0:

                            directory = path_output + "/HIS/" + str(hit[2]) + "/" + mut_name_full + "/"
                            if not os.path.exists(directory):
                                os.makedirs(directory)

                            count = len(glob.glob1(directory, "*.pdb"))

                            io = PDBIO()
                            io.set_structure(structure)
                            io.save(directory + str(hit[6]) + " " + str(count + 1) + ".pdb")

                            if check_for_all_steric_clashes:
                                print_clash_data(directory + str(hit[6]) + " " + str(count + 1) + ".txt", hit[3])

                        directory = path_output + "/NUM/" + hit[5] + "/" + str(hit[4]) + "/" + mut_name_full + "/"

                        if not os.path.exists(directory):
                            os.makedirs(directory)

                        count = len(glob.glob1(directory, "*.pdb"))

                        io = PDBIO()
                        io.set_structure(structure)
                        io.save(directory + str(hit[6]) + " " + str(count + 1) + ".pdb")

                        if check_for_all_steric_clashes:
                            print_clash_data(directory + str(hit[6]) + " " + str(count + 1) + ".txt", hit[3])

            mutations = {}


#   --------------------------------------------------------------------------------------- Renaming folders ---------------------------------------------------------------------------------------


# Filesystem:

#   output/type/NUM/num_hits/mutations (num_struct)/coords.pdb
# example:
#   output/octahedral/4/R19E A74E K75D (6)/(15, 6.5, -2).pdb

#   output/type/HIS/num_HIS/mutations (num_struct)/coords.pdb
# example:
#   output/tetrahedral/2/R19E A74E K75D (3)/(15, 6.5, -2).pdb


dir_coord = ("O", "T")
dir_num = ("1", "2", "3", "4", "5", "6")

for d1 in dir_coord:
    for d2 in dir_num:
        path = ".//output/NUM/" + d1 + "/" + d2 + "/"

        if os.path.exists(path):
            for folder in next(os.walk(path))[1]:
                fullpath = os.path.join(path, folder)
                num = str(len(glob.glob1(fullpath, "*.pdb")))
                os.rename(fullpath, fullpath + " (" + num + ")")

for d1 in dir_num:
    path = ".//output/HIS/" + d1 + "/"

    if os.path.exists(path):
        for folder in next(os.walk(path))[1]:
            fullpath = os.path.join(path, folder)
            num = str(len(glob.glob1(fullpath, "*.pdb")))
            os.rename(fullpath, fullpath + " (" + num + ")")


now = time.time()
print("Script Finished in " + str(np.round(time.time() - now, 2)) + "s")

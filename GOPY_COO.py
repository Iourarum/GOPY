# -----------------------------------------------------------
# GOPY - A tool for building 2D graphene-based computational models
#
# Released under GNU Public License (GPL)
# email sebmuraru@gmail.com
# -----------------------------------------------------------

import os
from scipy import spatial
import random
import math 
import numpy as np
import sys

class Typical_Bond:
    def __init__(self, length, specific_id):
        self.length = length
        self.identity = specific_id
        
CX_CX_sg = Typical_Bond(1.418, "CXCXGGGGGG")
CX_CY_sg_C1A = Typical_Bond(1.418, "CXCYGGGC1A")
CX_CY_sg_H1A = Typical_Bond(1.418, "CXCYGGGH1A")
CX_CY_sg_E1A = Typical_Bond(1.418, "CXCYGGGE1A")
CX_COOH = Typical_Bond(1.520, "CXC4GGGC1A")
CY_COOH = Typical_Bond(1.520, "CYC4C1AC1A")
CX_OH = Typical_Bond(1.49, "CXOLGGGH1A")
CY_OH = Typical_Bond(1.49, "CYOLH1AH1A")
CX_OS = Typical_Bond(1.460, "CXOEGGGE1A")
CY_OS = Typical_Bond(1.460, "CYOEE1AE1A")
C4_OOH = Typical_Bond(1.210, "C4OJC1AC1A")
C4_OH = Typical_Bond(1.32 ,"C4OKC1AC1A")
O_H = Typical_Bond(0.967, "OLHKH1AH1A")
O_H_2 = Typical_Bond(0.967, "OKHKC1AC1A")

bond_list_1 = [CX_CX_sg, CX_CY_sg_C1A, CX_CY_sg_E1A, CX_CY_sg_H1A, CX_COOH, CY_COOH, CX_OH, CY_OH, CX_OS, CY_OS, CX_OH, CY_OH, C4_OOH, C4_OH, O_H, O_H_2] 

CX_N2 = Typical_Bond(1.475, "CXN2GGGP1A")
CY_N2 = Typical_Bond(1.475, "CYN2P1AP1A")
N2_C6 = Typical_Bond(1.475, "N2C6P1AP1A")
C6_C5 = Typical_Bond(1.540, "C6C5P1AP1A")
C5_O2 = Typical_Bond(1.430, "C5O2P1AP1A")
O2_C4 = Typical_Bond(1.430, "O2C4P1AP1A")
C4_C3 = Typical_Bond(1.540, "C4C3P1AP1A")
C3_O1 = Typical_Bond(1.430, "C3O1P1AP1A")
O1_C2 = Typical_Bond(1.430, "O1C2P1AP1A")
C2_C1 = Typical_Bond(1.540, "C2C1P1AP1A")
C1_N1 = Typical_Bond(1.475, "C1N1P1AP1A")
N2_H15 = Typical_Bond(1.01, "N2H15P1AP1A")
C6_H14 = Typical_Bond(1.09, "C6H14P1AP1A")
C6_H13 = Typical_Bond(1.09, "C6H13P1AP1A")
C5_H12 = Typical_Bond(1.09, "C5H12P1AP1A")
C5_H11 = Typical_Bond(1.09, "C5H11P1AP1A")
C4_H10 = Typical_Bond(1.09, "C4H10P1AP1A")
C4_H9 = Typical_Bond(1.09, "C4H9P1AP1A")
C3_H8 = Typical_Bond(1.09, "C3H8P1AP1A")
C3_H7 = Typical_Bond(1.09, "C3H7P1AP1A")
C2_H6 = Typical_Bond(1.09, "C2H6P1AP1A")
C2_H5 = Typical_Bond(1.09, "C2H5P1AP1A")
C1_H4 = Typical_Bond(1.09, "C1H4P1AP1A")
C1_H3 = Typical_Bond(1.09, "C1H3P1AP1A")
N1_H2 = Typical_Bond(1.01, "N1H2P1AP1A")
N1_H1 = Typical_Bond(1.01, "N1H1P1AP1A")

bond_list_2 = [CX_CX_sg, CX_CY_sg_C1A, CX_CY_sg_E1A, CX_CY_sg_H1A, CX_COOH, CY_COOH, CX_OH, CY_OH, CX_OS, CY_OS, CX_OH, CY_OH, C4_OOH, C4_OH, O_H, O_H_2, CX_N2, CY_N2, N2_C6, C6_C5, C5_O2, O2_C4, C4_C3, C3_O1, O1_C2, C2_C1, C1_N1, N2_H15, C6_H14, C6_H13, C5_H12, C5_H11, C4_H10, C4_H9, C3_H8, C3_H7, C2_H6, C2_H5, C1_H4, C1_H3, N1_H2, N1_H1]

gen_NH = Typical_Bond(1.010, "NH")
gen_NO = Typical_Bond(1.060, "NO")
gen_NC = Typical_Bond(1.475, "NC")
gen_NN = Typical_Bond(1.450, "NN")
gen_OH = Typical_Bond(0.970, "OH")
gen_OC = Typical_Bond(1.160, "OH")
gen_OO = Typical_Bond(1.490, "OO")
gen_CH = Typical_Bond(1.090, "CH")
gen_CC = Typical_Bond(1.540, "CC")
gen_HH = Typical_Bond(0.740, "HH")

bond_list_3 = [gen_NH, gen_NO, gen_NC, gen_NN, gen_OH, gen_OC, gen_OO, gen_CH, gen_CC, gen_HH]

class Atom:
    """
    This is the Atom class. All atom objects are expected to contain the usual parameters
    written in a PDB file: atom number, atom name, residue name, residue number and XYZ
    coordinates.
    """
    def __init__(self, atom_number, atom_name, residue_name, residue_number, x, y, z):
        self.atom_number = atom_number
        self.atom_name = atom_name
        self.residue_name = residue_name
        self.residue_number = residue_number
        self.x = x
        self.y = y
        self.z = z
        
### PRSTINE GRAPHENE ###   

def add_atom(atom_list, atom_name, residue_name, residue_number, x, y, z, atom_number):
    """
    Adds an atom to the atom list. The atom comes with all the information usually found in a PDB file.
    """
    atom_list.append(Atom(atom_number, atom_name, residue_name, residue_number, x, y, z))
    return atom_list

def remove_atom(atom_list, atom):
    """
    Removes an atom from the atom list. 
    """
    del atom_list[atom.atom_number - 1]
    del atom
    return atom_list
    
def generate_pristine_graphene(x_dim, y_dim, filename1): 
    """
    This is the function used to generate a pristine graphene sheet. To fit PDB standards, all 
    distances are expressed in Angstroms. 
    """
    y_number = round(y_dim / 1.228)
    x_number = int(x_dim / 2.127)
    x_addition = (x_dim / 2.127 ) % 1
    list_of_coords = []
    a = 0
    b = 0
    c = 0
    list_of_coords = fill_row(list_of_coords, y_number, a,b,c, [], 5, prev = False)
    for i in range(1,x_number):
        if (i == x_number-1):
            if (i % 2 == 1):
                a += 1.228
                b += 2.127
                list_of_coords = fill_row(list_of_coords, y_number, a, b, c, [], 6, prev = True)
                fill_hexagon(list_of_coords, -1.228, b, c, [0, 1, 3, 4, 5], full=6, prev=False)
            if (i % 2 == 0):
                a -= 1.228
                b += 2.127
                list_of_coords = fill_row(list_of_coords, y_number, a, b, c, [], 6, prev = False)
                fill_hexagon(list_of_coords, y_number*1.228, b, c, [0, 1, 3, 4, 5], full=6, prev=False)
        elif (i % 2 == 1):
            a += 1.228
            b += 2.127
            list_of_coords = fill_row(list_of_coords, y_number, a, b, c, [], 6, prev = True)
        elif (i % 2 == 0):
            a -= 1.228
            b += 2.127
            list_of_coords = fill_row(list_of_coords, y_number, a, b, c, [], 6, prev = False)
    list_x_steps = [0, 0.33, 0.66, 1]
    x_step = min(list_x_steps, key=lambda x:abs(x-x_addition))
    if (x_step == 0.33):
        list_of_coords = fill_row(list_of_coords, y_number, 0, 0, 0, [], 6, prev = False)
        fill_hexagon(list_of_coords, y_number*1.228, 0, 0, [0, 1, 2, 3, 4], full=6, prev=False)
    elif (x_step == 0.66):
        if (x_number % 2 == 1):
            a += 1.228
            b += 2.127
            list_of_coords = fill_row(list_of_coords, y_number, a, b, c, [2], 6, prev = True)
        elif (x_number % 2 == 0):
            a -= 1.228
            b += 2.127
            list_of_coords = fill_row(list_of_coords, y_number, a, b, c, [2], 6, prev = False)
    elif (x_step == 1):
        if (x_number % 2 == 1):
            a += 1.228
            b += 2.127
            list_of_coords = fill_row(list_of_coords, y_number, a, b, c, [], 6, prev = True)
        elif (x_number % 2 == 0):
            a -= 1.228
            b += 2.127
            list_of_coords = fill_row(list_of_coords, y_number, a, b, c, [], 6, prev = False)
    writepdb3(list_of_coords, filename1)
    print('done.')
    return list_of_coords

def writepdb3(list_of_coords, name):
    """
    This is the function used to write a PDB file using the list of coordinates from generate_pristine_graphene. 
    """
    list_of_coords2 = []
    for elem in range(len(list_of_coords)):
        if (list_of_coords[elem] not in list_of_coords2):
            list_of_coords2.append(list_of_coords[elem])
    os.chdir(os.getcwd())
    if ((".pdb" not in name) and (".PDB" not in name)):
        string = str(name) + ".pdb"
    else:
        string = str(name)
    with open(string, 'a') as le_file:
        for element in range(len(list_of_coords)):
            temp_atom = Atom(element, "CX", "GGG", element, list_of_coords[element][0], list_of_coords[element][1], list_of_coords[element][2])
            line = "ATOM" + lw2(7, str(temp_atom.atom_number)) + str(temp_atom.atom_number) + lw(4, str(temp_atom.atom_name)) + str(temp_atom.atom_name) + "  " + str(temp_atom.residue_name) + lw(6, str(temp_atom.residue_number)) + str(temp_atom.residue_number) + lw(12, str(temp_atom.x)) + str(temp_atom.x) + lw(8, str(temp_atom.y)) + str(temp_atom.y) + lw(8, str(temp_atom.z)) + str(temp_atom.z) + "  1.00  0.00             "
            del temp_atom
            le_file.write(line + '\n')
    return list_of_coords

def lw2(max_no, str_obj):
    """Function used to add white spaces when writing the PDB file from writepdb3."""
    x = max_no - len(str_obj)
    y = 0
    string = ''
    for y in range(x):
        string = string + ' '
    return string

def fill_row(list_of_coords, halves, a,b,c, takeout, full=6, prev = False):
    """Used to fill a row of hexagons given a starting point XYZ coords. (here a,b,c) and the number of "halves" of a hexagon to add."""
    for i in range(int(halves / 2 )):
        list_of_coords = fill_hexagon(list_of_coords, a, b, c, takeout, full, prev)
        a += 2.456
    if ((halves % 2 == 1) and (prev == False)):
        takeout.append(3)
        takeout.append(4)
        list_of_coords = fill_hexagon(list_of_coords, a, b, c, takeout, full, prev)
    return list_of_coords

def toString(triplet):
    return([str(triplet[0]), str(triplet[1]), str(triplet[2])])
    
def check_me(triplet, list_of_coords):
    """Used to make sure atoms are within range of where they are supposed to be - and to provide a small error tolerance if necessary."""
    c = True
    for element in list_of_coords:
        if (float(triplet[0])*0.99 <= float(element[0]) <= float(triplet[0])*1.01):
            if (float(triplet[1])*0.99 <= float(element[1]) <= float(triplet[1])*1.01):
                if (float(triplet[2])*0.99 <= float(element[2]) <= float(triplet[2])*1.01):
                    c = False
    return c
    
def fill_hexagon(list_of_coords, a, b, c, takeout, full = 6, prev = False):
    """Given the XYZ coordinates (here a,b,c) of a carbon atom, this function is used to find the coordinates of the other atoms making up the hexagon.
    Takeout takes out one of the 6 atoms making up the hexagon. "prev" is used to know whether the row is an odd number row or not."""
    list_of_vertices = [[a, b, c], [a, b + 1.418, c], [a + 1.228, b + 2.127, c], [a + 2.456, b + 1.418, c], [a + 2.456, b, c], [a + 1.228, b - 0.705, c]]
    for element in range(full):
        if (check_me(list_of_vertices[element], list_of_coords) and (int(element) not in takeout)):
            list_of_coords.append(list_of_vertices[element])
    if ((prev == True) and (check_me([a - 1.228, b - 0.705, 0], list_of_coords))):
        if ((2 not in takeout) and (check_me([a - 1.228, b + 2.127, 0], list_of_coords))):
            list_of_coords.append([a - 1.228, b + 2.127, 0])
        list_of_coords.append([a - 1.228, b - 0.705, 0])
    for elem in range(len(list_of_coords)):
        list_of_coords[elem] = ["%.3f" % float(elem2) for elem2 in list_of_coords[elem]]
    return list_of_coords

def pristine_coords_to_objects(list_of_coords):
    """Function used to take the list of coordinates and add an Atom object instance for each coordinate triplet"""
    list_of_objects = []
    for element in range(len(list_of_coords)):
        list_of_objects.append(Atom(element, "CX", "GGG", element, list_of_coords[element][0], list_of_coords[element][1], list_of_coords[element][2]))
    return list_of_objects

def read_in_graphene(pdbfile):
    """Reads in a PG layer where C atom name is 'CX', residue name is 'GGG'"""
    with open(pdbfile, "r") as f:
        filedata = f.read()
        filedata = filedata.replace("C   GRA X", "CX  GGG  ")
        content = filedata.splitlines()
        atom_lines = [x.split() for x in content if (('ATOM' in str(x)) and ('GGG' in str(x)) and ('CX' in str(x)))]
        atoms = [Atom(int(str(atom_lines[x][1])), str(atom_lines[x][2]), str(atom_lines[x][3]), int(str(atom_lines[x][1])), float(str(atom_lines[x][5])), float(str(atom_lines[x][6])), float(str(atom_lines[x][7]))) for x in range(len(atom_lines))] 
    return atoms

def read_in_GO(pdbfile):
    """Reads in a GO layer where C atom name is "CX", residue name is 'GGG'
    Expected carboyl residue name: C1A; Expected epoxy residue name: E1A; Expected hydroxyl residue name: H1A;"""
    with open(pdbfile, "r") as f:
        filedata = f.read()
        filedata = filedata.replace("C   GRA X", "CX  GGG  ")
        content = filedata.splitlines()
        atom_lines = [x.split() for x in content if (('ATOM' in str(x)) and (('C1A' in str(x)) or ('E1A' in str(x)) or ('H1A' in str(x)) or ('GGG' in str(x))))]
        atoms = [Atom(int(str(atom_lines[x][1])), str(atom_lines[x][2]), str(atom_lines[x][3]), int(str(atom_lines[x][4])), float(str(atom_lines[x][5])), float(str(atom_lines[x][6])), float(str(atom_lines[x][7]))) for x in range(len(atom_lines))] 
    return atoms
 
def calculate_3D_distance_2_atoms(atom1, atom2):
    return spatial.distance.euclidean((atom1.x, atom1.y, atom1.z), (atom2.x, atom2.y, atom2.z))

def calculate_3D_distance_2_centers(x1, y1, z1, x2, y2, z2):
    return spatial.distance.euclidean((x1, y1, z1), (x2, y2, z2))

def get_bond_id(atom1, atom2): 
    """Creates bond_id of two atoms, taking into account their residue names and atom_names."""
    id1 = [str(atom1.atom_name) + str(atom2.atom_name) + str(atom1.residue_name) + str(atom2.residue_name), str(atom2.atom_name) + str(atom1.atom_name) + str(atom2.residue_name) + str(atom1.residue_name)]
    return id1

def check_bond(atom1, atom2):
    """A potential bond between two atoms is verified in the bond_list taking into account its id resulting from get_bond_id and distance between the two atoms."""
    check = False
    for bond in bond_list:
        if (((bond.identity == get_bond_id(atom1, atom2)[0]) or (bond.identity == get_bond_id(atom1, atom2)[1])) and 0.975 * bond.length <= calculate_3D_distance_2_atoms(atom1, atom2) <= 1.025 * bond.length):
            check = True
            break
    return check

def check_if_no_bond(atom1, atom2, bond_list, bond_generic):
    """Similar to check_bond, but using two specified bond_lists from the three available. Used when placing hydrogen atoms in add_rGO_PEG_NH2"""
    check = False
    for bond in bond_list:
        if ((bond.identity == get_bond_id(atom1, atom2)[0]) or (bond.identity == get_bond_id(atom1, atom2)[1]) and calculate_3D_distance_2_atoms(atom1, atom2) > 1.05 * bond.length):
            check = True
    for bond in bond_generic:
        if (((atom1.atom_name[0] + atom2.atom_name[0]) == bond.identity) or (atom2.atom_name[0] + atom1.atom_name[0] == bond.identity) and (calculate_3D_distance_2_atoms(atom1, atom2) > 1.05 * bond.length)):
            check = True            
    return check
    
def check_connected(chosen_atom, identified_bonds):
    """Based on the get_bond_id output, determine whether an atom is connected to a functional group. Used when removing functional groups."""
    check = False
    for bond in identified_bonds:
        if (("E1AE1A" in str(get_bond_id(chosen_atom, bond[0])[0])) or ("C1AC1A" in str(get_bond_id(chosen_atom, bond[0])[0])) or ("H1AH1A" in str(get_bond_id(chosen_atom, bond[0])[0])) or ("P1AP1A" in str(get_bond_id(chosen_atom, bond[0])[0]))):
            check = True
    return check
            
def identify_bonds(chosen_atom, atom_list):
    """The identify_bonds function finds the neighbours of the given atom in a distance based manner. 
    This distance is different for when placing GO groups to when placing PEG-NH2 chains and also when placing hydrogens vs the skeleton of the PEG-NH2 chains.
    The function then uses the check_bond function to identify valid bonds."""
    list_of_hydrogens = ['H15', 'H14', 'H13', 'H12', 'H11', 'H10', 'H9', 'H8', 'H7', 'H6', 'H5', 'H4', 'H3', 'H2', 'H1']    
    if ((chosen_atom.atom_name not in list_of_hydrogens) and (chosen_atom.residue_name != "P1A")):
        nearby_atoms_crude = [atom for atom in atom_list if ((abs(chosen_atom.x - atom.x) <= 2) and (abs(chosen_atom.y - atom.y) <= 2) and (abs(chosen_atom.z - atom.z) <= 2))]
        nearby_atoms = [atom for atom in nearby_atoms_crude if (0 < calculate_3D_distance_2_atoms(chosen_atom,atom) <= 2)]
        identified_bonds = [[atom, calculate_3D_distance_2_atoms(chosen_atom, atom)] for atom in nearby_atoms if (check_bond(chosen_atom, atom) == True)]              
    elif ((chosen_atom.atom_name not in list_of_hydrogens) and (chosen_atom.residue_name == "P1A")):
        nearby_atoms_crude = [atom for atom in atom_list if ((abs(chosen_atom.x - atom.x) <= 2) and (abs(chosen_atom.y - atom.y) <= 2) and (abs(chosen_atom.z - atom.z) <= 2))]
        nearby_atoms = [atom for atom in nearby_atoms_crude if (0 < calculate_3D_distance_2_atoms(chosen_atom,atom) <= 1.8)]
        identified_bonds = [[atom, calculate_3D_distance_2_atoms(chosen_atom, atom)] for atom in nearby_atoms if (check_bond(chosen_atom, atom) == True)]               
    else:
        nearby_atoms_crude = [atom for atom in atom_list if ((abs(chosen_atom.x - atom.x) <= 1.6) and (abs(chosen_atom.y - atom.y) <= 1.6) and (abs(chosen_atom.z - atom.z) <= 1.6))]
        nearby_atoms = [atom for atom in nearby_atoms_crude if (0 < calculate_3D_distance_2_atoms(chosen_atom,atom) <= 1.6)]
        identified_bonds = [[atom, calculate_3D_distance_2_atoms(chosen_atom, atom)] for atom in nearby_atoms if (check_bond(chosen_atom, atom) == True)]               
        for elements in nearby_atoms:
            if (check_if_no_bond(chosen_atom, elements, bond_list, bond_list_3) == True):
                nearby_atoms.remove(elements)
    if (len(nearby_atoms) == len(identified_bonds)):
        return identified_bonds
    else:
        return []

def top_or_down():
    """Random top or below draw."""
    ct = random.randint(1,2)
    if (ct == 1):
        return 1
    else:
        return -1
    
def get_map_anywhere(atom_list):
    """Creates a list of all available atoms which are not connected to functional groups."""
    anywhere_map = [atom for atom in atom_list if (check_connected(atom, identify_bonds(atom, atom_list)) == False)]
    return anywhere_map
    
def get_map_central(atom_list):
    """Creates a list of all available atoms which are not connected to a functional group and do not represent an edge of any kind (thus 3 bonds required)."""
    central_map = [atom for atom in atom_list if ((len(identify_bonds(atom, atom_list)) == 3) and (check_connected(atom, identify_bonds(atom, atom_list)) == False))]
    return central_map

def get_map_edge(atom_list):
    """Creates a list of all edge atoms (thus 1 or 2 bonds) which are not connected to functional groups."""
    edge_map = [atom for atom in atom_list if ((0 < len(identify_bonds(atom, atom_list)) < 3) and (check_connected(atom, identify_bonds(atom, atom_list)) == False))]
    return edge_map

def pick_to_add(no_cooh, no_epoxy, no_hydroxyl):
    """Random pick."""
    make_list = ["carboxyl", "epoxy", "hydroxyl"]
    if (no_cooh == 0):
        make_list.remove("carboxyl")
    if (no_epoxy == 0):
        make_list.remove("epoxy")
    if (no_hydroxyl == 0):
        make_list.remove("hydroxyl")
    if (make_list == []):
        return "done"
    else:
        return make_list
    
def create_GO(init_file, no_COOH, no_epoxy, no_OH, filename1):
    """Function used to create the GO morphology. The three functional groups are placed in no specific order, picking the next to be added at random.
    One may want to place COOH groups first due to limited space, however this is rarely an issue. 50 attempts are made for each placement. After 50 attempts
    the placement is considered done (despite there being no successful placement) and the function moves to the next functional group to place.
    After placement, there is a rewriting of the PDB file where each CX atom connected to a functional group becomes a CY atom and new atom numbers are given
    to take into account the CX atoms replaced by CY. """
    global atoms
    global bond_list
    bond_list = bond_list_1
    atoms = read_in_graphene(init_file)
    global anywhere_map
    anywhere_map = get_map_anywhere(atoms)
    global edge_map
    edge_map = get_map_edge(atoms)
    
    list_residue_numbers = [x.residue_number for x in atoms]
    added_functional_groups = max(list_residue_numbers)
    
    must_add = no_COOH + no_epoxy + no_OH
    while (must_add > 0):
        print("Left to add: ", "cooh: ", no_COOH, "epoxy: ", no_epoxy, "hydroxyl: ", no_OH)
        chosen = random.choice(pick_to_add(no_COOH, no_epoxy, no_OH))
        if (chosen == "carboxyl"):
            attempt = 0
            while (attempt < 50):
                old_length = len(atoms)
                new_atoms = add_carboxyl(random_pick_spot("carboxyl", edge_map, anywhere_map), atoms, added_functional_groups, top_or_down())
                if (old_length != len(new_atoms)):
                    atoms = new_atoms
                    added_functional_groups += 1
                    must_add -= 1
                    no_COOH -= 1
                    attempt = 1888
                else:
                    attempt += 1
            if (attempt == 50):
                must_add = -1
        elif (chosen == "epoxy"):     
            attempt = 0
            while (attempt < 50):
                old_length = len(atoms)
                new_atoms = add_epoxy(random_pick_spot("epoxy", edge_map, anywhere_map), atoms, added_functional_groups, top_or_down())
                if (old_length != len(new_atoms)):
                    atoms = new_atoms
                    added_functional_groups += 1
                    must_add -= 1
                    no_epoxy -= 1
                    attempt = 1888
                else:
                    attempt += 1
            if (attempt == 50):
                must_add = -1
        elif (chosen == "hydroxyl"):
            attempt = 0
            while (attempt < 50):
                old_length = len(atoms)
                new_atoms = add_hydroxyl(random_pick_spot("hydroxyl", edge_map, anywhere_map), atoms, added_functional_groups, top_or_down())
                if (old_length != len(new_atoms)):
                    atoms = new_atoms
                    added_functional_groups += 1
                    must_add -= 1
                    no_OH -=1
                    attempt = 1888                    
                else:
                    attempt += 1
            if (attempt == 50):
                must_add = -1
    atno = 1
    new_list = []
    for atom in atoms:
        if (atom.atom_name == "CX"):
            New_CX = Atom(atno, "CX", "GGG", atno, atom.x, atom.y, atom.z)
            new_list.append(New_CX)
            atno += 1 
    
    for atom in atoms:
        if (atom.atom_name == "C4"):
            check = False
            for atom_CY in atoms:
                if ((atom_CY.atom_name == "CY") and (atom_CY.residue_name == "C1A") and (atom_CY.residue_number == atom.residue_number)):
                    for atom_OJ in atoms:
                        if ((atom_OJ.atom_name == "OJ") and (atom_OJ.residue_name == "C1A") and (atom_OJ.residue_number == atom.residue_number)):
                            for atom_OK in atoms:
                                if ((atom_OK.atom_name == "OK") and (atom_OK.residue_name == "C1A") and (atom_OK.residue_number == atom.residue_number)):
                                    # for atom_HK in atoms:
                                        # if ((atom_HK.atom_name == "HK") and (atom_HK.residue_name == "C1A") and (atom_HK.residue_number == atom.residue_number)):
                                    New_CY = Atom(atno + 0, "CY", "C1A", atom.residue_number, atom_CY.x, atom_CY.y, atom_CY.z )
                                    New_C4 = Atom(atno + 1, "C4", "C1A", atom.residue_number, atom.x, atom.y, atom.z)
                                    New_OJ = Atom(atno + 2, "OJ", "C1A", atom.residue_number, atom_OJ.x, atom_OJ.y, atom_OJ.z)
                                    New_OK = Atom(atno + 3, "OK", "C1A", atom.residue_number, atom_OK.x, atom_OK.y, atom_OK.z)
                                        # New_HK = Atom(atno + 4, "HK", "C1A", atom.residue_number, atom_HK.x, atom_HK.y, atom_HK.z)
                                    atno += 4 #5
                                    new_list.append(New_CY); new_list.append(New_C4); new_list.append(New_OJ); new_list.append(New_OK); # new_list.append(New_HK);
                                    check = True
                                    break
                                if (check == True):
                                    break
                        if (check == True):
                            break
                if (check == True):
                    break                
                            
        elif (atom.atom_name == "OE"):  
            check = False
            for atom_CY in atoms:
                if ((atom_CY.atom_name == "CY") and (atom_CY.residue_name == "E1A") and (atom_CY.residue_number == atom.residue_number)):
                    for atom_CY2 in atoms: 
                        if ((atom_CY2.atom_name == "CZ") and (atom_CY2.residue_name == "E1A") and (atom_CY2.residue_number == atom.residue_number) and (atom_CY2 != atom_CY)):
                            New_CY = Atom( atno + 0, "CY", "E1A", atom.residue_number, atom_CY.x,  atom_CY.y,  atom_CY.z)
                            New_CY2 = Atom(atno + 1, "CZ", "E1A", atom.residue_number, atom_CY2.x, atom_CY2.y, atom_CY2.z)
                            New_OE = Atom( atno + 2, "OE", "E1A", atom.residue_number, atom.x,     atom.y,     atom.z)
                            atno += 3
                            new_list.append(New_CY); new_list.append(New_CY2); new_list.append(New_OE);
                            check = True
                            break
                if (check == True):
                    break
        elif (atom.atom_name == "OL"):
            check = False
            for atom_CY in atoms:
                if ((atom_CY.atom_name == "CY") and (atom_CY.residue_name == "H1A") and (atom_CY.residue_number == atom.residue_number)):
                    for atom_HK in atoms:
                        if ((atom_HK.atom_name == "HK") and (atom_HK.residue_name == "H1A") and (atom_HK.residue_number == atom.residue_number)):
                            New_CY = Atom(atno + 0, "CY", "H1A", atom.residue_number, atom_CY.x, atom_CY.y, atom_CY.z)
                            New_OL = Atom(atno + 1, "OL", "H1A", atom.residue_number, atom.x,    atom.y,    atom.z)
                            New_HK = Atom(atno + 2, "HK", "H1A", atom.residue_number, atom_HK.x, atom_HK.y, atom_HK.z)
                            atno += 3
                            new_list.append(New_CY); new_list.append(New_OL); new_list.append(New_HK);
                            check = True
                            break
                if (check == True):
                    break
            
    atoms = new_list.copy()
    writepdb(atoms, filename1)
    sum_c1a = 0; sum_e1a = 0; sum_h1a = 0; sum_ggg = 0
    for atom in atoms:
        if (atom.residue_name == "C1A"):
            sum_c1a += 1
        elif (atom.residue_name == "E1A"):
            sum_e1a += 1
        elif (atom.residue_name == "H1A"):
            sum_h1a += 1
        elif (atom.residue_name == "GGG"):
            sum_ggg += 1
    print("Placed:")
    print("carboxyl: ", sum_c1a/4)
    print("epoxy: ", sum_e1a/3)
    print("hydroxyl: ", sum_h1a/3)
    print("graphene atoms (CX - GGG) left: ", sum_ggg)
    return 'done.'

def random_pick_spot(spot, map_edge, map_anywhere):
    if (spot == "carboxyl"):
        atom = random.choice(map_edge)
    elif (spot == "epoxy"):
        atom = random.choice(map_anywhere)
    elif (spot == "hydroxyl"):
        atom = random.choice(map_anywhere)
    return atom
        
def add_carboxyl(atom, atom_list, added_functional_groups, ct):
    """Function used to add carboxyl atoms. Atom placement follows a geometrical approach.
    Atom - chosen CX atom. Atom_list - the atom list so far. Added_functional_groups - the number of functional groups added so far. Ct - the random top or below coordinate.
    Adds 4 atoms and replaces 1 (5 in total).
    """
    
    global anywhere_map
    global edge_map
    current_size = len(atom_list)
    placed = 0
    alpha = random.randint(0,359)
    x = atom.x
    y = atom.y
    z = atom.z
    while (placed <= 359):
        alpha += 1
        carbon_atom = Atom(current_size + 1, 'C4', 'C1A', str(added_functional_groups + 1), float("{0:.3f}".format(x)), float("{0:.3f}".format(y)), float("{0:.3f}".format(ct * 1.52 + z)))
        atom_list.append(carbon_atom)
        if ((len(identify_bonds(carbon_atom, atom_list)) == 1) and (identify_bonds(carbon_atom, atom_list)[0][0].atom_number == atom.atom_number)):
            h = math.sin(math.radians(60)) * 1.20
            oxygen_atom_1 = Atom(current_size + 2, 'OJ', 'C1A', str(added_functional_groups + 1), float("{0:.3f}".format(carbon_atom.x - math.cos(math.radians(alpha)) * h)), float("{0:.3f}".format(carbon_atom.y - math.sin(math.radians(alpha)) * h)), float("{0:.3f}".format(carbon_atom.z + ct * math.cos(math.radians(60)) * 1.20)))
            atom_list.append(oxygen_atom_1)
            if ((len(identify_bonds(oxygen_atom_1, atom_list)) == 1) and (identify_bonds(oxygen_atom_1, atom_list)[0][0].atom_number == carbon_atom.atom_number)):
                h = math.sin(math.radians(60)) * 1.34
                oxygen_atom_2 = Atom(current_size + 3, 'OK', 'C1A', str(added_functional_groups + 1), float("{0:.3f}".format(carbon_atom.x - math.cos(math.radians(alpha + 180)) * h)), float("{0:.3f}".format(carbon_atom.y - math.sin(math.radians(alpha+180)) * h)), float("{0:.3f}".format(carbon_atom.z + ct * math.cos(math.radians(60)) * 1.34)) )
                atom_list.append(oxygen_atom_2)
                if ((len(identify_bonds(oxygen_atom_2, atom_list)) == 1) and (identify_bonds(oxygen_atom_2, atom_list)[0][0].atom_number == carbon_atom.atom_number)):
                    # hydrogen_atom = Atom(current_size + 4, 'HK', 'C1A', str(added_functional_groups + 1), float("{0:.3f}".format(oxygen_atom_2.x)),  float("{0:.3f}".format(oxygen_atom_2.y)),  float("{0:.3f}".format(oxygen_atom_2.z+ct*0.98)))
                    # atom_list.append(hydrogen_atom)
                    # if ((len(identify_bonds(hydrogen_atom, atom_list)) == 1) and (identify_bonds(hydrogen_atom, atom_list)[0][0].atom_number == oxygen_atom_2.atom_number)):
                    placed = 888
                    if atom in edge_map: edge_map.remove(atom)
                    if atom in anywhere_map: anywhere_map.remove(atom)
                    CY = Atom(atom.atom_number, 'CY', 'C1A', carbon_atom.residue_number, atom.x, atom.y, atom.z)
                    atom_list.append(CY)
                    atom_list.remove(atom)
                    del atom
                    # else:
                        # placed += 5
                        # del hydrogen_atom
                        # del atom_list[current_size + 3]
                        # del oxygen_atom_2
                        # del atom_list[current_size + 2]
                        # del oxygen_atom_1
                        # del atom_list[current_size + 1]
                        # del carbon_atom
                        # del atom_list[current_size + 0]
                else:
                    placed += 5
                    del oxygen_atom_2
                    del atom_list[current_size + 2]
                    del oxygen_atom_1
                    del atom_list[current_size + 1]
                    del carbon_atom
                    del atom_list[current_size + 0]
            else:
                placed += 5
                del oxygen_atom_1
                del atom_list[current_size + 1]
                del carbon_atom
                del atom_list[current_size + 0]
        else:
            placed += 5
            del carbon_atom 
            del atom_list[current_size + 0]
    return atom_list

def add_epoxy(atom, atom_list, added_functional_groups, ct):
    """Function used to add epoxy atoms. Atom placement follows a geometrical approach.
    Atom - chosen CX atom. Atom_list - the atom list so far. Added_functional_groups - the number of functional groups added so far. Ct - the random top or below coordinate.
    Ads 1 atom and replaces 2 (3 in total)."""
    global anywhere_map
    global edge_map
    current_size = len(atom_list)
    list_of_n = [x[0] for x in identify_bonds(atom, atom_list) if (x[0].atom_name != "CY")]
    if (len(list_of_n) != 0):
        atom2 = random.choice(list_of_n)
        epoxy_atom = Atom( current_size + 1, 'OE', 'E1A', str(added_functional_groups + 1), float("{0:.3f}".format(abs(atom.x - atom2.x) / 2 + min(atom.x, atom2.x))), float("{0:.3f}".format(abs(atom.y - atom2.y) / 2 + min(atom.y, atom2.y))), float("{0:.3f}".format(ct * 1.46 * math.sin(math.radians(60)))) )
        atom_list.append(epoxy_atom)
        if ((len(identify_bonds(epoxy_atom, atom_list)) == 2) and ((identify_bonds(epoxy_atom, atom_list)[0][0].atom_number == atom.atom_number) or (identify_bonds(epoxy_atom, atom_list)[1][0].atom_number == atom.atom_number)) and ((identify_bonds(epoxy_atom, atom_list)[0][0].atom_number == atom2.atom_number) or ((identify_bonds(epoxy_atom, atom_list)[1][0].atom_number == atom2.atom_number)))):
            CY = Atom(atom.atom_number, 'CY', 'E1A', epoxy_atom.residue_number, atom.x, atom.y, atom.z)
            CY2 = Atom(atom2.atom_number, 'CZ', 'E1A', epoxy_atom.residue_number, atom2.x, atom2.y, atom2.z)
            atom_list.remove(atom); atom_list.remove(atom2)
            atom_list.append(CY); atom_list.append(CY2) 
            if atom in edge_map: edge_map.remove(atom)
            if atom in anywhere_map: anywhere_map.remove(atom)
            if atom2 in edge_map: edge_map.remove(atom2)
            if atom2 in anywhere_map: anywhere_map.remove(atom2)
        else:
            atom_list.remove(epoxy_atom)
            del epoxy_atom
    return atom_list

def add_hydroxyl(atom, atom_list, added_functional_groups, ct): 
    """Function used to add hydroxyl atoms. Atom placement follows a geometrical approach.
    Atom - chosen CX atom. Atom_list - the atom list so far. Added_functional_groups - the number of functional groups added so far. Ct - the random top or below coordinate.
    Ads 2 atom and replaces 1 (3 in total)."""
    global anywhere_map
    global edge_map
    current_size = len(atom_list)
    placed = 0
    alpha = random.randint(0,359)
    while (placed <= 359):
        alpha += 1
        oxygen_atom = Atom(current_size + 1, 'OL', 'H1A', str(added_functional_groups + 1), atom.x, atom.y, ct * 1.49 + atom.z)
        atom_list.append(oxygen_atom)
        if ((len(identify_bonds(oxygen_atom, atom_list)) == 1) and (identify_bonds(oxygen_atom, atom_list)[0][0].atom_number == atom.atom_number)):
            h = math.sin(math.radians(19)) * 0.98
            h_sp = math.cos(math.radians(19)) * 0.98
            hydrogen_atom = Atom(current_size + 2, 'HK', 'H1A', str(added_functional_groups + 1), float("{0:.3f}".format(oxygen_atom.x - math.cos(math.radians(alpha)) * h_sp)), float("{0:.3f}".format(oxygen_atom.y - math.sin(math.radians(alpha)) * h_sp)), float("{0:.3f}".format(oxygen_atom.z + ct * h)))
            atom_list.append(hydrogen_atom)
            if ((len(identify_bonds(hydrogen_atom, atom_list)) == 1) and (identify_bonds(hydrogen_atom, atom_list)[0][0].atom_number == oxygen_atom.atom_number)):
                placed = 888
                if atom in edge_map: edge_map.remove(atom)
                if atom in anywhere_map: anywhere_map.remove(atom)             
                CY = Atom(atom.atom_number, 'CY', 'H1A', oxygen_atom.residue_number, atom.x, atom.y, atom.z)
                atom_list.append(CY)
                atom_list.remove(atom)         
                del atom                                      
            else:
                placed += 5
                del hydrogen_atom
                del atom_list[current_size + 1]
                del oxygen_atom
                del atom_list[current_size + 0]
        else:
            placed += 5
            del oxygen_atom
            del atom_list[current_size + 0]
    return atom_list

def get_hydroxyl_map(atom_list):
    """Creates list of hydroxyl groups. Used when removing functional groups (eg. in rGO-PEG-NH2). Does not include CY atoms."""
    hydroxyl_map = [[atom_list[x], atom_list[x+1]] for x in range(len(atom_list)-1) if ((atom_list[x].residue_name == atom_list[x+1].residue_name == "H1A") and (atom_list[x].residue_number == atom_list[x+1].residue_number ) and (atom_list[x].atom_name != "CY" != atom_list[x+1].atom_name != atom_list[x].atom_name))]
    return hydroxyl_map

def get_epoxy_map(atom_list):
    """Creates list of epoxy groups. Used when removing functional groups (eg. in rGO-PEG-NH2). Does not include the CY atoms."""
    epoxy_map = [[atom_list[x]] for x in range(len(atom_list)) if ((atom_list[x].residue_name == "E1A") and (atom_list[x].atom_name != "CY"))]
    return epoxy_map

def get_carboxyl_map(atom_list):
    """Creates list of carboxyl groups. Used when removing functional groups (eg. in rGO-PEG-NH2). Does not include CY atoms."""
    carboxyl_map = [[atom_list[x], atom_list[x+1], atom_list[x+2], atom_list[x+3]] for x in range(len(atom_list)-3) if ((atom_list[x].residue_name == atom_list[x+1].residue_name == atom_list[x+2].residue_name == atom_list[x+3].residue_name == "C1A") and (atom_list[x].residue_number == atom_list[x+1].residue_number == atom_list[x+2].residue_number == atom_list[x+3].residue_number) and (atom_list[x].atom_name != "CY" != atom_list[x+1].atom_name != atom_list[x+2].atom_name != "CY" != atom_list[x+3].atom_name ))]
    return carboxyl_map

def remove_functional_groups(pct_C1A, pct_E1A, pct_H1A, atom_list):
    """Used to remove *pct_C1A* of present carboyl groups etc."""
    carboxyl_map = get_carboxyl_map(atom_list)
    epoxy_map = get_epoxy_map(atom_list)
    hydroxyl_map = get_hydroxyl_map(atom_list)
    remove_C1A = round(len(carboxyl_map) * pct_C1A)
    remove_E1A = round(len(epoxy_map) * pct_E1A)
    remove_H1A = round(len(hydroxyl_map) * pct_H1A)
    while (remove_C1A > 0):
        remove_C1A -= 1
        remove_group = random.choice(carboxyl_map)
        carboxyl_map.remove(remove_group)
        for element in remove_group:
            atom_list.remove(element)
            del element
    while (remove_E1A > 0):
        remove_E1A -= 1
        remove_group = random.choice(epoxy_map)
        epoxy_map.remove(remove_group)
        for element in remove_group:
            atom_list.remove(element)
            del element
    while (remove_H1A > 0):
        remove_H1A -= 1
        remove_group = random.choice(hydroxyl_map)
        hydroxyl_map.remove(remove_group)
        for element in remove_group:
            atom_list.remove(element)
            del element
    return atom_list

def find_highest_resnum(atom_list):
    """Used to find the highest / last residue number."""
    resnum = 0
    for element in atom_list:
        if (int(element.residue_number) > int(resnum)):
            resnum = int(element.residue_number)
    return resnum

def find_conn_CXCY(atom, atom_list):
    """Once a functional group is removed (such as epoxy or hydroxyl) and it is going to be replaced by a -NH-PEG-NH2 chain, one needs to know
    the exact atom to which to connect the new functional group and whether it was connected above or below. This function serves that."""
    le_list = []
    for element in identify_bonds(atom, atom_list):
        if ((element[0].atom_name == "CX") or (element[0].atom_name == "CY")):
            if (atom.z - element[0].z > 0):
                le_list.append([element[0], 1])
            else:
                le_list.append([element[0], -1])
    return le_list

def fix_sphere_m (center_x, center_y, center_z, radius, centers, radii, len_points):
    """Using the given XYZ coordinates as the center of the sphere, and the typical bond distance of the two atoms of interest
    we plot random points on a sphere. In a similar manner, we then imagine the spheres of nearby atoms (centers and radii) and 
    test the distance between each random point and nearby atoms. If too small, the points get removed.
    Len_points refers to how many random points one should plot/place."""
    
    g_x = []
    g_y = []
    g_z = []
    points = [hydrogen_coord_gen(center_x, center_y, center_z, radius) for i in range(0, len_points)]  
    x = [points[i][0] for i in range(0, len(points))]       
    y = [points[i][1] for i in range(0, len(points))]
    z = [points[i][2] for i in range(0, len(points))]

    for i in range(0, len(points)):
        check = 0
        j = 0
        while (j <= (len(centers) - 1) and (check == 0)): 
            if (calculate_3D_distance_2_centers(x[i], y[i], z[i], centers[j][0], centers[j][1], centers[j][2]) < radii[j]):
                check += 1
            j += 1
        if (check == 0):
            g_x.append(x[i])
            g_y.append(y[i])
            g_z.append(z[i])

    return g_x, g_y, g_z

def fix_sphere_h (center_x, center_y, center_z, radius, centers, radii, len_points, list_of_a):
    """Similar to above, yet takes into account through list_of_a the neighbouring skeleton atoms of the PEG-NH2 chains so that they are not missed by accident.
    This was necessary as an optimization of the *radii* or bond length as these had to be tweaked and played with slightly."""
    g_x = []
    g_y = []
    g_z = []
    points = [hydrogen_coord_gen(center_x, center_y, center_z, radius) for i in range(0, len_points)]  
    x = [points[i][0] for i in range(0, len(points))]       
    y = [points[i][1] for i in range(0, len(points))]
    z = [points[i][2] for i in range(0, len(points))]
    for i in range(0, len(points)):
        check = 0
        check_b = 0
        j = 0
        while (j <= (len(centers) - 1) and (check == 0)): 
            if (calculate_3D_distance_2_centers(x[i], y[i], z[i], centers[j][0], centers[j][1], centers[j][2]) < radii[j]):
                check += 1
            j += 1
        h = 0
        while ((check_b == 0) and (h <= len(list_of_a) -1)):
            if (calculate_3D_distance_2_centers(x[i], y[i], z[i], list_of_a[h].x, list_of_a[h].y, list_of_a[h].z) <= 1.50):               
                check_b += 1
            h += 1
        if ((check == 0) and (check_b == 0)):
            g_x.append(x[i])
            g_y.append(y[i])
            g_z.append(z[i])
    return g_x, g_y, g_z

def detect_neighbours_hydrogen(central_atom, atom_list):
    "Returns list of nighbour atoms up to 2,5A"""
    nearby_atoms_crude = [atom for atom in atom_list if ((abs(central_atom.x - atom.x) <= 2.5) and (abs(central_atom.y - atom.y) <= 2.5) and (abs(central_atom.z - atom.z) <= 2.5))]
    nearby_atoms = [atom for atom in nearby_atoms_crude if (0 < calculate_3D_distance_2_atoms(central_atom, atom) <= 2.5)]
    return nearby_atoms

def quick_sphere(center_x, center_y, center_z, radius, no_points):
    points = []
    points = [hydrogen_coord_gen(center_x, center_y, center_z, radius) for i in range(0, no_points)]  
    x = [points[i][0] for i in range(0, len(points))]       
    y = [points[i][1] for i in range(0, len(points))]
    z = [points[i][2] for i in range(0, len(points))]
    return x, y, z
    
def hydrogen_coord_gen(atomx, atomy, atomz, r):
    """Generates random XYZ points on a sphere centered at atomx, atomy, atomz with radius r"""
    theta = np.random.uniform(0, 2 * np.pi)
    cosalpha = np.random.uniform(-1, 1)
    alpha = np.arccos(cosalpha)
    x = atomx + r * np.cos(theta) * np.sin(alpha)
    y = atomy + r * np.sin(theta) * np.sin(alpha)
    z = atomz + r * np.cos(alpha)
    return [x, y, z]
       
def plot_points_on_sphere(points_x, points_y, points_z, center_x, center_y, center_z, radius):
    """Plotting function. Used to create images in the manuscript."""
    from mayavi import mlab
    mlab.figure(1, bgcolor=(1,1,1), fgcolor=(0,0,0), size=(800,800))
    return mlab.points3d(points_x, points_y, points_z, scale_factor=0.05, color=(0.25, 0.75, 0.77))

def repick(g1_x, g1_y, g1_z):
    """Repick point from the random points plotted on the sphere. Usually due to not being able to find place for a second hydrogen etc."""
    i = random.randint(0, len(g1_x) - 1)
    x = g1_x[i]
    y = g1_y[i]
    z = g1_z[i]
    return x, y, z
    
def find_CX_neighbours(list_of_atoms, atom_list):
    """As the PEG-NH2 chains can have a more or less random conformation, some of them may bend along the graphene layer.
    This function returns a list of the CX or CY graphene atoms nearby a PEG-NH2 chain."""
    my_list = []
    atom_numbers = []
    for atom in list_of_atoms:
        for element in identify_bonds(atom, atom_list):
            if (((element[0].atom_name == "CX") or (element[0].atom_name == "CY")) and (element[0].atom_number not in atom_numbers)):
                my_list.append(element[0])
                atom_numbers.append(element[0].atom_number)
    return my_list
        
def compose_listofr(atom_name, listofn):
    """Typical distances between certain atom pairs (N, O, C, H).
    Was used to tweak the distance with neighbouring spheres (see fix_sphere_m or fix_sphere_h) and improve the outcome."""
    c = 1.06
    c2 = 1.4
    listofr = []
    for x in range(len(listofn)):
        if (atom_name[0] == "N"):
            if (listofn[x].atom_name[0] == "H"):
                listofr.append(1.010*c)
            if (listofn[x].atom_name[0] == "O"):
                listofr.append(1.060*c)
            if (listofn[x].atom_name[0] == "C"):
                listofr.append(1.475*c)
            if (listofn[x].atom_name[0] == "N"):
                listofr.append(1.450*c)
        if (atom_name[0] == "O"):
            if (listofn[x].atom_name[0] == "H"):
                listofr.append(0.970*c)
            if (listofn[x].atom_name[0] == "O"):
                listofr.append(1.490*c)
            if (listofn[x].atom_name[0] == "C"):
                listofr.append(1.160*c)
            if (listofn[x].atom_name[0] == "N"):
                listofr.append(1.060*c)
        if (atom_name[0] == "C"):
            if (listofn[x].atom_name[0] == "H"):
                listofr.append(1.090*c)
            if (listofn[x].atom_name[0] == "O"):
                listofr.append(1.160*c)
            if (listofn[x].atom_name[0] == "C"):
                listofr.append(1.540*c)
            if (listofn[x].atom_name[0] == "N"):
                listofr.append(1.475*c)
        if (atom_name[0] == "H"):
            if (listofn[x].atom_name[0] == "H"):
                listofr.append(0.740*c2)
            if (listofn[x].atom_name[0] == "O"):
                listofr.append(0.970*c2)
            if (listofn[x].atom_name[0] == "C"):
                listofr.append(1.090*c2)
            if (listofn[x].atom_name[0] == "N"):
                listofr.append(1.010*c2)
    return listofr

def add_NH_PEG_NH2(GO_file, pct_C1A, pct_E1A, pct_H1A, filename1): 
    """The function used to place NH-PEG-NH2 chains with the formula: -NH-(C2H4O)2-NH2.
    Needs a GO PDB file to start from. pct_C1A, pct_E1A and pct_H1A represent how many functional groups will be removed from the GO (between 0 (none) and 1 (all)).
    After removal, all atoms are re-numbered and then the chains are placed. Each chain has 25 atoms to place so it can take a while to place one."""
    global remember_me
    global bond_list
    bond_list = bond_list_1
    atom_list = read_in_GO(GO_file)
    carboxyl_map = get_carboxyl_map(atom_list)
    epoxy_map = get_epoxy_map(atom_list)
    hydroxyl_map = get_hydroxyl_map(atom_list)
    remove_C1A = round(len(carboxyl_map) * pct_C1A)
    remove_E1A = round(len(epoxy_map) * pct_E1A)
    remove_H1A = round(len(hydroxyl_map) * pct_H1A)
    list_of_CXCY_epoxy = []
    list_of_CXCY_hydroxyl = []
    while (remove_C1A > 0):
        remove_C1A -= 1
        remove_group = random.choice(carboxyl_map)
        carboxyl_map.remove(remove_group)
        for element in remove_group:
            atom_list.remove(element)
            del element
    while (remove_E1A > 0):
        remove_E1A -= 1
        remove_group = random.choice(epoxy_map)
        epoxy_map.remove(remove_group)
        for element in remove_group:
            CXCY_epoxy = find_conn_CXCY(element, atom_list)
            if (CXCY_epoxy != []):
                list_of_CXCY_epoxy.append(CXCY_epoxy)
            atom_list.remove(element)
            del element

    while (remove_H1A > 0):
        remove_H1A -= 1
        remove_group = random.choice(hydroxyl_map)
        hydroxyl_map.remove(remove_group)
        for element in remove_group:
            CXCY_hydroxyl = find_conn_CXCY(element, atom_list)
            if (CXCY_hydroxyl != []):
                list_of_CXCY_hydroxyl.append(CXCY_hydroxyl)
            atom_list.remove(element)
            del element

    atno = 1
    resno = 1
    new_list = []
    for atom in atom_list:
        if ((atom.atom_name == "CX") or (atom.atom_name == "CY")):
            New_CX = Atom(atno, "CX", "GGG", resno, atom.x, atom.y, atom.z)
            new_list.append(New_CX)
            atno += 1 
            resno += 1

    for atom in atom_list:
        if (atom.atom_name == "C4"):
            check = False
            for atom_OJ in atom_list:
                if ((atom_OJ.atom_name == "OJ") and (atom_OJ.residue_name == "C1A") and (atom_OJ.residue_number == atom.residue_number)):
                    for atom_OK in atom_list:
                        if ((atom_OK.atom_name == "OK") and (atom_OK.residue_name == "C1A") and (atom_OK.residue_number == atom.residue_number)):
                            for atom_HK in atom_list:
                                if ((atom_HK.atom_name == "HK") and (atom_HK.residue_name == "C1A") and (atom_HK.residue_number == atom.residue_number)):
                                    New_C4 = Atom(atno + 0, "C4", "C1A", resno, atom.x, atom.y, atom.z)
                                    New_OJ = Atom(atno + 1, "OJ", "C1A", resno, atom_OJ.x, atom_OJ.y, atom_OJ.z)
                                    New_OK = Atom(atno + 2, "OK", "C1A", resno, atom_OK.x, atom_OK.y, atom_OK.z)
                                    New_HK = Atom(atno + 3, "HK", "C1A", resno, atom_HK.x, atom_HK.y, atom_HK.z)
                                    atno += 4
                                    resno += 1
                                    new_list.append(New_C4); new_list.append(New_OJ); new_list.append(New_OK); new_list.append(New_HK);
                                    check = True
                                    break
                        if (check == True):
                            break
                if (check == True):
                    break
             
                            
        elif (atom.atom_name == "OE"):  
            New_OE = Atom(atno + 0, "OE", "E1A", resno, atom.x,     atom.y,     atom.z)
            atno += 1
            resno += 1
            new_list.append(New_OE);
        elif (atom.atom_name == "OL"):
            check = False
            for atom_HK in atom_list:
                if ((atom_HK.atom_name == "HK") and (atom_HK.residue_name == "H1A") and (atom_HK.residue_number == atom.residue_number)):
                    New_OL = Atom(atno + 0, "OL", "H1A", resno, atom.x,    atom.y,    atom.z)
                    New_HK = Atom(atno + 1, "HK", "H1A", resno, atom_HK.x, atom_HK.y, atom_HK.z)
                    atno += 2
                    resno += 1
                    new_list.append(New_OL); new_list.append(New_HK);
                    check = True
                    break

    atom_list = new_list.copy()
    atoms = new_list.copy()

    bond_list = bond_list_2 + bond_list_3
    suma = 0
    quick_list = list_of_CXCY_epoxy + list_of_CXCY_hydroxyl
    
    new_quick_list = []
    for element in quick_list:
        for atom in atom_list:
            if ((atom.x == element[0][0].x) and (atom.y == element[0][0].y) and (atom.z == element[0][0].z)):
                new_quick_list.append([[atom, element[0][1]]])         
    progress = 1
    for element in new_quick_list:
        print("Progress: ", progress/(len(new_quick_list))*100, "%")
        progress +=1
        goto = 100
        points_max = 1500
        suma += 1                
        ct = element[0][1]
        resnum = find_highest_resnum(atom_list)
        current_size = len(atom_list)
        
        attempt_N1 = goto
        listofa = find_CX_neighbours([element[0][0]], atoms)
        listofn = [[atom.x, atom.y, atom.z] for atom in listofa]
        listofr = compose_listofr("N2", listofa)
        g1_x, g1_y, g1_z = fix_sphere_m(element[0][0].x, element[0][0].y, element[0][0].z, 1.475, listofn, listofr, points_max)
        while (0 <= attempt_N1 <= goto):
#            print("N1 :   ", attempt_N1)
            i = random.randint(0, len(g1_x) - 1)
            x1 = g1_x[i]
            y1 = g1_y[i]
            z1 = g1_z[i]
            while (abs(ct * z1) < element[0][0].z):
                x1, y1, z1 = repick(g1_x, g1_y, g1_z)
            N1 =  Atom(current_size + 1,   'N2', 'P1A', str(resnum + 1), float("{0:.3f}".format(x1)),         float("{0:.3f}".format(y1)),         float("{0:.3f}".format(z1)))
            atom_list.append(N1)
            if ((len(identify_bonds(N1, atom_list)) == 1) and (identify_bonds(N1,atom_list)[0][0].atom_number == element[0][0].atom_number)):
                attempt_N1 = 888
            else:
                attempt_N1 -= 1
                atom_list.remove(N1)
                del N1
                
        attempt_C1 = goto
        attempt15 = goto
        if (attempt_N1 != -1):
            listofa = find_CX_neighbours([element[0][0], N1], atoms)
            listofa.append(element[0][0])
            listofn = [[atom.x, atom.y, atom.z] for atom in listofa]
            listofr = compose_listofr("C6", listofa)
            g1_x, g1_y, g1_z = fix_sphere_m(N1.x, N1.y, N1.z, 1.475, listofn, listofr, points_max)
            while ((0 <= attempt_C1 <= goto) and (attempt_N1 != -1)):
                i = random.randint(0, len(g1_x) - 1)
                x1 = g1_x[i]
                y1 = g1_y[i]
                z1 = g1_z[i]
                while (abs(ct * z1) < N1.z):
                    x1, y1, z1 = repick(g1_x, g1_y, g1_z)
                C1 =  Atom(current_size + 2,   'C6', 'P1A', str(resnum + 1), float("{0:.3f}".format(x1)),                    float("{0:.3f}".format(y1)),                    float("{0:.3f}".format(z1)))
                atom_list.append(C1)
                if ((len(identify_bonds(C1, atom_list)) == 1) and (identify_bonds(C1, atom_list)[0][0].atom_number == N1.atom_number)):
                    attempt_C1 = 888
                    listofa = detect_neighbours_hydrogen(N1, atom_list)
                    listofa.append(element[0][0])
                    listofa.append(C1)
                    listofn = [[atom.x, atom.y, atom.z] for atom in listofa]
                    listofr = compose_listofr("H15", listofa)
                    g1_x, g1_y, g1_z = fix_sphere_h(N1.x, N1.y, N1.z, 1.01, listofn, listofr, points_max, [element[0][0], C1])
                    while (0 <= attempt15 <= goto):
                        i = random.randint(0, len(g1_x) -1)
                        x1 = g1_x[i]
                        y1 = g1_y[i]
                        z1 = g1_z[i]
                        H15 = Atom(current_size + 11,   'H15', 'P1A', str(resnum + 1), float("{0:.3f}".format(x1)),             float("{0:.3f}".format(y1)),            float("{0:.3f}".format(z1)))
                        atom_list.append(H15)
                        if ((len(identify_bonds(H15, atom_list)) == 1) and (identify_bonds(H15, atom_list)[0][0].atom_number == N1.atom_number)):
                            attempt15 = 888
                        else:
                            attempt15 -= 1
                            atom_list.remove(H15)
                            del H15  
                else:
                    attempt_C1 -= 1
                    atom_list.remove(C1)
                    del C1      

        attempt_C2 = goto
        attempt1413 = goto
        if ((attempt_C1 > -1) and (attempt_N1 > -1) and (attempt15 > -1)):
            listofa = find_CX_neighbours([element[0][0], N1, C1, H15], atoms)
            listofa.append(element[0][0])
            listofa.append(N1)
            listofa.append(H15)
            listofn = [[atom.x, atom.y, atom.z] for atom in listofa]
            listofr = compose_listofr("C5", listofa)
            g1_x, g1_y, g1_z = fix_sphere_m(C1.x, C1.y, C1.z, 1.540, listofn, listofr, points_max)
            while (0 <= attempt_C2 <= goto):
                i = random.randint(0, len(g1_x) - 1)
                x1 = g1_x[i]
                y1 = g1_y[i]
                z1 = g1_z[i]
                while (abs(ct *z1) < C1.z):
                    x1, y1, z1 = repick(g1_x, g1_y, g1_z)
                C2 =  Atom(current_size + 3,   'C5', 'P1A', str(resnum + 1), float("{0:.3f}".format(x1)),                    float("{0:.3f}".format(y1)),                    float("{0:.3f}".format(z1)))
                atom_list.append(C2)
                if ((len(identify_bonds(C2, atom_list)) == 1) and (identify_bonds(C2, atom_list)[0][0].atom_number == C1.atom_number)):
                    attempt_C2 = 888
                    listofa = detect_neighbours_hydrogen(C1, atom_list)
                    listofa.append(element[0][0])
                    listofa.append(N1)
                    listofa.append(C2)
                    listofa.append(H15)
                    listofn = [[atom.x, atom.y, atom.z] for atom in listofa]
                    listofr = compose_listofr("H14", listofa)
                    g1_x, g1_y, g1_z = fix_sphere_h(C1.x, C1.y, C1.z, 1.11, listofn, listofr, points_max, [C2, N1])
                    while (0 <= attempt1413 <= goto):
                        i = random.randint(0, len(g1_x)-1)
                        x1 = g1_x[i]
                        y1 = g1_y[i]
                        z1 = g1_z[i]   
                        H14 = Atom(current_size + 12,   'H14', 'P1A', str(resnum + 1), float("{0:.3f}".format(x1)),             float("{0:.3f}".format(y1)),            float("{0:.3f}".format(z1)))                            
                        atom_list.append(H14)
                        if ((len(identify_bonds(H14, atom_list)) == 1) and (identify_bonds(H14, atom_list)[0][0].atom_number == C1.atom_number)):
                            listofa.append(H14)
                            listofn = [[atom.x, atom.y, atom.z] for atom in listofa]
                            listofr = compose_listofr("H13", listofa)
                            g2_x, g2_y, g2_z = fix_sphere_h(C1.x, C1.y, C1.z, 1.11, listofn, listofr, points_max, [C2, N1])
                            if (len(g2_x) != 0):
                                attempt2 = goto

                                while(0 <= attempt2 <= goto):
                                    i = random.randint(0, len(g2_x)-1)
                                    x2 = g2_x[i]
                                    y2 = g2_y[i]
                                    z2 = g2_z[i]    
                               
                                    H13 = Atom(current_size + 13,   'H13', 'P1A', str(resnum + 1), float("{0:.3f}".format(x2)),             float("{0:.3f}".format(y2)),            float("{0:.3f}".format(z2)))
                                    atom_list.append(H13)
                                    if ((len(identify_bonds(H13, atom_list)) == 1) and (identify_bonds(H13, atom_list)[0][0].atom_number == C1.atom_number)):
                                        attempt2 = 888
                                        attempt1413 = 888
                                    else:
                                        attempt2 -= 1
                                        attempt1413 -= 1
                                        atom_list.remove(H13)        
                                        del H13   
                            else:
                                attempt1413 -= 1
                                atom_list.remove(H14)
                                del H14
                        else:
                            attempt1413 -= 1
                            atom_list.remove(H14)
                            del H14
                else:
                    attempt_C2 -= 1
                    atom_list.remove(C2)
                    del C2
        
        attempt_O1 = goto
        attempt1211 = goto
        if ((attempt_C2 > -1) and (attempt_C1 > -1) and (attempt_N1 > -1) and (attempt1413 > -1) and (attempt15 > -1)):
            listofa = find_CX_neighbours([element[0][0], N1, C1, H15, C2, H15, H14, H13], atoms)
            listofa.append(element[0][0])
            listofa.append(N1)
            listofa.append(C1)
            listofa.append(H15)
            listofa.append(H14)
            listofa.append(H13)
            listofn = [[atom.x, atom.y, atom.z] for atom in listofa]
            listofr = compose_listofr("O2", listofa)
            g1_x, g1_y, g1_z = fix_sphere_m(C2.x, C2.y, C2.z, 1.430, listofn, listofr, points_max)
            while (0 <= attempt_O1 <= goto):
                i = random.randint(0, len(g1_x) - 1)
                x1 = g1_x[i]
                y1 = g1_y[i]
                z1 = g1_z[i]
                while (abs(ct*z1) < C2.z):
                    x1, y1, z1 = repick(g1_x, g1_y, g1_z)
                O1 =  Atom(current_size + 4,   'O2', 'P1A', str(resnum + 1), float("{0:.3f}".format(x1)),                    float("{0:.3f}".format(y1)),              float("{0:.3f}".format(z1)))
                atom_list.append(O1)
                if ((len(identify_bonds(O1, atom_list)) == 1) and (identify_bonds(O1, atom_list)[0][0].atom_number == C2.atom_number)):
                    attempt_O1 = 888
                    listofa = detect_neighbours_hydrogen(C2, atom_list)
                    listofa.append(element[0][0])
                    listofa.append(N1)
                    listofa.append(C1)
                    listofa.append(O1)
                    listofa.append(H15)
                    listofa.append(H14)
                    listofa.append(H13)
                    listofn = [[atom.x, atom.y, atom.z] for atom in listofa]
                    listofr = compose_listofr("H12", listofa)
                    g1_x, g1_y, g1_z = fix_sphere_h(C2.x, C2.y, C2.z, 1.11, listofn, listofr, points_max, [C1, O1])
                    while (0 <= attempt1211 <= goto):
                        i = random.randint(0, len(g1_x)-1)
                        x1 = g1_x[i]
                        y1 = g1_y[i]
                        z1 = g1_z[i]   
                        H12 = Atom(current_size + 14,   'H12', 'P1A', str(resnum + 1), float("{0:.3f}".format(x1)),            float("{0:.3f}".format(y1)),            float("{0:.3f}".format(z1)))   
                        atom_list.append(H12)
                        if ((len(identify_bonds(H12, atom_list)) == 1) and (identify_bonds(H12, atom_list)[0][0].atom_number == C2.atom_number)):                            
                            listofa.append(H12)
                            listofn = [[atom.x, atom.y, atom.z] for atom in listofa]
                            listofr = compose_listofr("H11", listofa)
                            g2_x, g2_y, g2_z = fix_sphere_h(C2.x, C2.y, C2.z, 1.11, listofn, listofr, points_max, [C1, O1])
                            if (len(g2_x) != 0):
                                i = random.randint(0, len(g2_x)-1)
                                x2 = g2_x[i]
                                y2 = g2_y[i]
                                z2 = g2_z[i]
                                H11 = Atom(current_size + 15,   'H11', 'P1A', str(resnum + 1), float("{0:.3f}".format(x2)),            float("{0:.3f}".format(y2)),            float("{0:.3f}".format(z2)))            
                                atom_list.append(H11)
                                if ((len(identify_bonds(H12, atom_list)) == 1) and (identify_bonds(H12, atom_list)[0][0].atom_number == C2.atom_number) and (len(identify_bonds(H11, atom_list))==1) and (identify_bonds(H11, atom_list)[0][0].atom_number == C2.atom_number)):
                                    attempt1211 = 888
                                else:
                                    attempt1211 -= 1
                                    atom_list.remove(H12)
                                    atom_list.remove(H11)
                                    del H12
                                    del H11
                            else:
                                attempt1211 -= 1
                                atom_list.remove(H12)
                                del H12
                        else:
                            attempt1211 -= 1
                            atom_list.remove(H12)
                            del H12
                else:
                    attempt_O1 -= 1
                    atom_list.remove(O1)
                    del O1                        
        
        attempt_C3 = goto
        if ((attempt_O1 > -1) and (attempt_C2 > -1) and (attempt_C1 > -1) and (attempt_N1 > -1) and (attempt1211 > -1) and (attempt1413 > -1) and (attempt15 > -1)):
            listofa = find_CX_neighbours([element[0][0], N1, C1, C2, O1, H15, H14, H13, H12, H11], atoms)
            listofa.append(element[0][0])
            listofa.append(N1)
            listofa.append(C1)
            listofa.append(C2)
            listofa.append(H15)
            listofa.append(H14)
            listofa.append(H13)
            listofa.append(H12)
            listofa.append(H11)
            listofn = [[atom.x, atom.y, atom.z] for atom in listofa]
            listofr = compose_listofr("C4", listofa)
            g1_x, g1_y, g1_z = fix_sphere_m(O1.x, O1.y, O1.z, 1.430, listofn, listofr, points_max)
            while (0 <= attempt_C3 <= goto):
                i = random.randint(0, len(g1_x) - 1)
                x1 = g1_x[i]
                y1 = g1_y[i]
                z1 = g1_z[i]
                while (abs(ct*z1) < O1.z):
                    x1, y1, z1 = repick(g1_x, g1_y, g1_z)
                C3 =  Atom(current_size + 4,   'C4', 'P1A', str(resnum + 1), float("{0:.3f}".format(x1)),                    float("{0:.3f}".format(y1)),              float("{0:.3f}".format(z1)))
                atom_list.append(C3)
                if ((len(identify_bonds(C3, atom_list)) == 1) and (identify_bonds(C3, atom_list)[0][0].atom_number == O1.atom_number)):
                    attempt_C3 = 888
                else:
                    attempt_C3 -= 1
                    atom_list.remove(C3)
                    del C3
        
        attempt_C4 = goto
        attempt1009 = goto 
        if ((attempt_C3 > -1) and (attempt_O1 > -1) and (attempt_C2 > -1) and (attempt_C1 > -1) and (attempt_N1 > -1) and (attempt1211 > -1) and (attempt1413 > -1) and (attempt15 > -1)):
            listofa = find_CX_neighbours([element[0][0], N1, C1, C2, O1, C3, H15, H14, H13, H12, H11], atoms)
            listofa.append(element[0][0])
            listofa.append(N1)
            listofa.append(C1)
            listofa.append(C2)
            listofa.append(O1)
            listofa.append(H15)
            listofa.append(H14)
            listofa.append(H13)
            listofa.append(H12)
            listofa.append(H11)
            listofn = [[atom.x, atom.y, atom.z] for atom in listofa]
            listofr = compose_listofr("C3", listofa)
            g1_x, g1_y, g1_z = fix_sphere_m(C3.x, C3.y, C3.z, 1.540, listofn, listofr, points_max)
            while (0 <= attempt_C4 <= goto):
                i = random.randint(0, len(g1_x) - 1)
                x1 = g1_x[i]
                y1 = g1_y[i]
                z1 = g1_z[i]
                while (abs(ct*z1) < C3.z):
                    x1, y1, z1 = repick(g1_x, g1_y, g1_z)
                C4 =  Atom(current_size + 4,   'C3', 'P1A', str(resnum + 1), float("{0:.3f}".format(x1)),                    float("{0:.3f}".format(y1)),              float("{0:.3f}".format(z1)))
                atom_list.append(C4)
                if ((len(identify_bonds(C4, atom_list)) == 1) and (identify_bonds(C4, atom_list)[0][0].atom_number == C3.atom_number)):
                    attempt_C4 = 888
                    listofa = detect_neighbours_hydrogen(C3, atom_list)
                    listofa.append(element[0][0])
                    listofa.append(N1)
                    listofa.append(O1)
                    listofa.append(C1)
                    listofa.append(C2)
                    listofa.append(C4)
                    listofa.append(H15)
                    listofa.append(H14)
                    listofa.append(H13)
                    listofa.append(H12)
                    listofa.append(H11)
                    listofn = [[atom.x, atom.y, atom.z] for atom in listofa]
                    listofr = compose_listofr("H10", listofa)
                    g1_x, g1_y, g1_z = fix_sphere_h(C3.x, C3.y, C3.z, 1.11, listofn, listofr, points_max, [O1, C4])
                    while (0 <= attempt1009 <= goto):
                        i = random.randint(0, len(g1_x)-1)
                        x1 = g1_x[i]
                        y1 = g1_y[i]
                        z1 = g1_z[i]
                        H10 = Atom(current_size + 16,   'H10', 'P1A', str(resnum + 1), float("{0:.3f}".format(x1)),            float("{0:.3f}".format(y1)),            float("{0:.3f}".format(z1)))  
                        atom_list.append(H10)
                        if ((len(identify_bonds(H10, atom_list)) == 1) and (identify_bonds(H10, atom_list)[0][0].atom_number == C3.atom_number)):                            
                            listofa.append(H10)
                            listofn = [[atom.x, atom.y, atom.z] for atom in listofa]
                            listofr = compose_listofr("H09", listofa)
                            g2_x, g2_y, g2_z = fix_sphere_h(C3.x, C3.y, C3.z, 1.11, listofn, listofr, points_max, [O1, C4])
                            if (len(g2_x) != 0):
                                i = random.randint(0, len(g2_x)-1)
                                x2 = g2_x[i]
                                y2 = g2_y[i]
                                z2 = g2_z[i]
                                H9  = Atom(current_size + 17,    'H9', 'P1A', str(resnum + 1), float("{0:.3f}".format(x2)),            float("{0:.3f}".format(y2)),            float("{0:.3f}".format(z2)))
                                atom_list.append(H9)
                                if ((len(identify_bonds(H10, atom_list)) == 1) and (identify_bonds(H10, atom_list)[0][0].atom_number == C3.atom_number) and (len(identify_bonds(H9, atom_list)) == 1) and (identify_bonds(H9, atom_list)[0][0].atom_number == C3.atom_number)):
                                    attempt1009 = 888
                                else:
                                    attempt1009 -= 1      
                                    atom_list.remove(H10)
                                    atom_list.remove(H9)
                                    del H10
                                    del H9         
                            else:
                                attempt1009 -= 1
                                atom_list.remove(H10)
                                del H10
                        else:
                            attempt1009 -= 1
                            atom_list.remove(H10)
                            del H10
                else:
                    attempt_C4 -= 1
                    atom_list.remove(C4)
                    del C4

        attempt0807 = goto
        attempt_O2 = goto
        if ((attempt_C4 > -1) and (attempt_C3 > -1) and (attempt_O1 > -1) and (attempt_C2 > -1) and (attempt_C1 > -1) and (attempt_N1 > -1) and (attempt1009 > -1) and (attempt1211 > -1) and (attempt1413 > -1) and (attempt15 > -1)):
            listofa = find_CX_neighbours([element[0][0], N1, C1, C2, O1, C3, C4, H15, H14, H13, H12, H11, H10, H9], atoms)
            listofa.append(element[0][0])
            listofa.append(N1)
            listofa.append(C1)
            listofa.append(C2)
            listofa.append(O1)
            listofa.append(C3)
            listofa.append(H15)
            listofa.append(H14)
            listofa.append(H13)
            listofa.append(H12)
            listofa.append(H11)
            listofa.append(H10)
            listofa.append(H9)
            listofn = [[atom.x, atom.y, atom.z] for atom in listofa]
            listofr = compose_listofr("O1", listofa)
            g1_x, g1_y, g1_z = fix_sphere_m(C4.x, C4.y, C4.z, 1.430, listofn, listofr, points_max)
            while (0 <= attempt_O2 <= goto):
                i = random.randint(0, len(g1_x) - 1)
                x1 = g1_x[i]
                y1 = g1_y[i]
                z1 = g1_z[i]
                while (abs(ct*z1) < C4.z):
                    x1, y1, z1 = repick(g1_x, g1_y, g1_z)
                O2 =  Atom(current_size + 4,   'O1', 'P1A', str(resnum + 1), float("{0:.3f}".format(x1)),                    float("{0:.3f}".format(y1)),              float("{0:.3f}".format(z1)))
                atom_list.append(O2)
                if ((len(identify_bonds(O2, atom_list)) == 1) and (identify_bonds(O2, atom_list)[0][0].atom_number == C4.atom_number)):
                    attempt_O2 = 888
                    listofa = detect_neighbours_hydrogen(C4, atom_list)
                    listofa.append(element[0][0])
                    listofa.append(N1)
                    listofa.append(O1)
                    listofa.append(C1)
                    listofa.append(C2)
                    listofa.append(C3)
                    listofa.append(O2)
                    listofa.append(H15)
                    listofa.append(H14)
                    listofa.append(H13)
                    listofa.append(H12)
                    listofa.append(H11)
                    listofa.append(H10)
                    listofa.append(H9)
                    listofn = [[atom.x, atom.y, atom.z] for atom in listofa]
                    listofr = compose_listofr("H8", listofa)
                    g1_x, g1_y, g1_z = fix_sphere_h(C4.x, C4.y, C4.z, 1.11, listofn, listofr, points_max, [C3, O2])
                    while (0 <= attempt0807 <= goto):
                        i = random.randint(0, len(g1_x)-1)
                        x1 = g1_x[i]
                        y1 = g1_y[i]
                        z1 = g1_z[i]
                        H8  = Atom(current_size + 18,    'H8', 'P1A', str(resnum + 1), float("{0:.3f}".format(x1)),            float("{0:.3f}".format(y1)),            float("{0:.3f}".format(z1)))
                        atom_list.append(H8)
                        if ((len(identify_bonds(H8, atom_list)) == 1) and (identify_bonds(H8, atom_list)[0][0].atom_number == C4.atom_number)):                            
                            listofa.append(H8)
                            listofn = [[atom.x, atom.y, atom.z] for atom in listofa]
                            listofr = compose_listofr("H7", listofa)
                            g2_x, g2_y, g2_z = fix_sphere_h(C4.x, C4.y, C4.z, 1.11, listofn, listofr, points_max, [C3, O2])
                            if (len(g2_x) != 0):
                                i = random.randint(0, len(g2_x)-1)
                                x2 = g2_x[i]
                                y2 = g2_y[i]
                                z2 = g2_z[i]
                                H7  = Atom(current_size + 19,    'H7', 'P1A', str(resnum + 1), float("{0:.3f}".format(x2)),            float("{0:.3f}".format(y2)),            float("{0:.3f}".format(z2))) 
                                atom_list.append(H7)
                                if ((len(identify_bonds(H8, atom_list)) == 1) and (identify_bonds(H8, atom_list)[0][0].atom_number == C4.atom_number) and (len(identify_bonds(H7, atom_list)) == 1) and (identify_bonds(H7, atom_list)[0][0].atom_number == C4.atom_number)):
                                    attempt0807 = 888
                                else:
                                    attempt0807 -= 1      
                                    atom_list.remove(H8)
                                    atom_list.remove(H7)
                                    del H8
                                    del H7
                            else:
                                attempt0807 -= 1
                                atom_list.remove(H8)
                                del H8  
                        else:
                            attempt0807 -= 1
                            atom_list.remove(H8)
                            del H8
                else:
                    attempt_O2 -= 1
                    atom_list.remove(O2)
                    del O2               
        
        attempt_C5 = goto
        if ((attempt_O2 > -1) and (attempt_C4 > -1) and (attempt_C3 > -1) and (attempt_O1 > -1) and (attempt_C2 > -1) and (attempt_C1 > -1) and (attempt_N1 > -1) and (attempt0807 > -1) and (attempt1009 > -1) and (attempt1211 > -1) and (attempt1413 > -1) and (attempt15 > -1)):
            listofa = find_CX_neighbours([element[0][0], N1, C1, C2, O1, C3, C4, O2, H15, H14, H13, H12, H11, H10, H9, H8, H7], atoms)
            listofa.append(element[0][0])
            listofa.append(N1)
            listofa.append(C1)
            listofa.append(C2)
            listofa.append(O1)
            listofa.append(C3)
            listofa.append(C4)          
            listofa.append(H15)
            listofa.append(H14)
            listofa.append(H13)
            listofa.append(H12)
            listofa.append(H11)
            listofa.append(H10)
            listofa.append(H9)
            listofa.append(H8)
            listofa.append(H7)     
            listofn = [[atom.x, atom.y, atom.z] for atom in listofa]
            listofr = compose_listofr("C2", listofa)
            g1_x, g1_y, g1_z = fix_sphere_m(O2.x, O2.y, O2.z, 1.430, listofn, listofr, points_max)
            while (0 <= attempt_C5 <= goto):
                i = random.randint(0, len(g1_x) - 1)
                x1 = g1_x[i]
                y1 = g1_y[i]
                z1 = g1_z[i]
                while (abs(ct*z1) < O2.z):
                    x1, y1, z1 = repick(g1_x, g1_y, g1_z)
                C5 =  Atom(current_size + 4,   'C2', 'P1A', str(resnum + 1), float("{0:.3f}".format(x1)),                    float("{0:.3f}".format(y1)),              float("{0:.3f}".format(z1)))
                atom_list.append(C5)
                if ((len(identify_bonds(C5, atom_list)) == 1) and (identify_bonds(C5, atom_list)[0][0].atom_number == O2.atom_number)):
                    attempt_C5 = 888
                else:
                    attempt_C5 -= 1
                    atom_list.remove(C5)
                    del C5
                
        attempt_C6 = goto
        attempt0605 = goto
        if ((attempt_C5 > -1) and (attempt_O2 > -1) and (attempt_C4 > -1) and (attempt_C3 > -1) and (attempt_O1 > -1) and (attempt_C2 > -1) and (attempt_C1 > -1) and (attempt_N1 > -1) and (attempt0807 > -1) and (attempt1009 > -1) and (attempt1211 > -1) and (attempt1413 > -1) and (attempt15 > -1)):
            listofa = find_CX_neighbours([element[0][0], N1, C1, C2, O1, C3, C4, O2, C5, H15, H14, H13, H12, H11, H10, H9, H8, H7], atoms)
            listofa.append(element[0][0])
            listofa.append(N1)
            listofa.append(C1)
            listofa.append(C2)
            listofa.append(O1)
            listofa.append(C3)
            listofa.append(C4)
            listofa.append(O2)       
            listofa.append(H15)
            listofa.append(H14)
            listofa.append(H13)
            listofa.append(H12)
            listofa.append(H11)
            listofa.append(H10)
            listofa.append(H9)
            listofa.append(H8)
            listofa.append(H7)
            listofn = [[atom.x, atom.y, atom.z] for atom in listofa]
            listofr = compose_listofr("C1", listofa)
            g1_x, g1_y, g1_z = fix_sphere_m(C5.x, C5.y, C5.z, 1.540, listofn, listofr, points_max)
            while (0 <= attempt_C6 <= goto):
                i = random.randint(0, len(g1_x) - 1)
                x1 = g1_x[i]
                y1 = g1_y[i]
                z1 = g1_z[i]
                while (abs(ct*z1) < C5.z):
                    x1, y1, z1 = repick(g1_x, g1_y, g1_z)
                C6 =  Atom(current_size + 4,   'C1', 'P1A', str(resnum + 1), float("{0:.3f}".format(x1)),                    float("{0:.3f}".format(y1)),              float("{0:.3f}".format(z1)))
                atom_list.append(C6)
                if ((len(identify_bonds(C6, atom_list)) == 1) and (identify_bonds(C6, atom_list)[0][0].atom_number == C5.atom_number)):
                    attempt_C6 = 888
                    listofa = detect_neighbours_hydrogen(C5, atom_list)
                    listofa.append(element[0][0])
                    listofa.append(N1)
                    listofa.append(O1)
                    listofa.append(C1)
                    listofa.append(C2)
                    listofa.append(C3)
                    listofa.append(C4)
                    listofa.append(O2)
                    listofa.append(C6)
                    listofa.append(H15)
                    listofa.append(H14)
                    listofa.append(H13)
                    listofa.append(H12)
                    listofa.append(H11)
                    listofa.append(H10)
                    listofa.append(H9)
                    listofa.append(H8)
                    listofa.append(H7)
                    listofn = [[atom.x, atom.y, atom.z] for atom in listofa]
                    listofr = compose_listofr("H6", listofa)
                    g1_x, g1_y, g1_z = fix_sphere_h(C5.x, C5.y, C5.z, 1.11, listofn, listofr, points_max, [O2, C6])
                    while (0 <= attempt0605 <= goto):
                        i = random.randint(0, len(g1_x)-1)
                        x1 = g1_x[i]
                        y1 = g1_y[i]
                        z1 = g1_z[i]
                        H6  = Atom(current_size + 20,    'H6', 'P1A', str(resnum + 1), float("{0:.3f}".format(x1)),            float("{0:.3f}".format(y1)),            float("{0:.3f}".format(z1)))
                        atom_list.append(H6)
                        if ((len(identify_bonds(H6, atom_list)) == 1) and (identify_bonds(H6, atom_list)[0][0].atom_number == C5.atom_number)):
                            listofa = detect_neighbours_hydrogen(C5, atom_list)
                            listofa.append(H6)
                            listofn = [[atom.x, atom.y, atom.z] for atom in listofa]
                            listofr = compose_listofr("H5", listofa)
                            g2_x, g2_y, g2_z = fix_sphere_h(C5.x, C5.y, C5.z, 1.11, listofn, listofr, points_max, [O2, C6])
                            if (len(g2_x) != 0):
                                i = random.randint(0, len(g2_x)-1)
                                x2 = g2_x[i]
                                y2 = g2_y[i]
                                z2 = g2_z[i]
                                H5  = Atom(current_size + 21,    'H5', 'P1A', str(resnum + 1), float("{0:.3f}".format(x2)),             float("{0:.3f}".format(y2)),            float("{0:.3f}".format(z2)))
                                atom_list.append(H5)
                                if ((len(identify_bonds(H6, atom_list)) == 1) and (identify_bonds(H6, atom_list)[0][0].atom_number == C5.atom_number) and (len(identify_bonds(H5, atom_list)) == 1) and (identify_bonds(H5, atom_list)[0][0].atom_number == C5.atom_number)):
                                    attempt0605 = 888
                                else:
                                    attempt0605 -= 1      
                                    atom_list.remove(H6)
                                    atom_list.remove(H5)
                                    del H6
                                    del H5
                            else:
                                attempt0605 -= 1
                                atom_list.remove(H6)
                                del H6            
                        else:
                            attempt0605 -= 1
                            atom_list.remove(H6)
                            del H6
                else:
                    attempt_C6 -= 1
                    atom_list.remove(C6)
                    del C6
   
        
        attempt_N2 = goto
        attempt0403 = goto
        if ((attempt_C6 > -1) and (attempt_C5 > -1) and (attempt_O2 > -1) and (attempt_C4 > -1) and (attempt_C3 > -1) and (attempt_O1 > -1) and (attempt_C2 > -1) and (attempt_C1 > -1) and (attempt_N1 > -1) and (attempt0605 > -1) and (attempt0807 > -1) and (attempt1009 > -1) and (attempt1211 > -1) and (attempt1413 > -1) and (attempt15 > -1)):
            listofa = find_CX_neighbours([element[0][0], N1, C1, C2, O1, C3, C4, O2, C5, C6, H15, H14, H13, H12, H11, H10, H9, H8, H7, H6, H5], atoms)
            listofa.append(element[0][0])
            listofa.append(N1)
            listofa.append(C1)
            listofa.append(C2)
            listofa.append(O1)
            listofa.append(C3)
            listofa.append(C4)
            listofa.append(O2)
            listofa.append(C5)
            listofa.append(H15)
            listofa.append(H14)
            listofa.append(H13)
            listofa.append(H12)
            listofa.append(H11)
            listofa.append(H10)
            listofa.append(H9)
            listofa.append(H8)
            listofa.append(H7)
            listofa.append(H6)
            listofa.append(H5)
            listofn = [[atom.x, atom.y, atom.z] for atom in listofa]
            listofr = compose_listofr("N1", listofa)
            g1_x, g1_y, g1_z = fix_sphere_m(C6.x, C6.y, C6.z, 1.475, listofn, listofr, points_max)
            while (0 <= attempt_N2 <= goto):
                i = random.randint(0, len(g1_x) - 1)
                x1 = g1_x[i]
                y1 = g1_y[i]
                z1 = g1_z[i]
                while (abs(ct*z1) < C6.z):
                    x1, y1, z1 = repick(g1_x, g1_y, g1_z)
                N2 =  Atom(current_size + 4,   'N1', 'P1A', str(resnum + 1), float("{0:.3f}".format(x1)),                    float("{0:.3f}".format(y1)),              float("{0:.3f}".format(z1)))
                atom_list.append(N2)
                if ((len(identify_bonds(N2, atom_list)) == 1) and (identify_bonds(N2, atom_list)[0][0].atom_number == C6.atom_number)):
                    attempt_N2 = 888
                    listofa = detect_neighbours_hydrogen(C6, atom_list)
                    listofa.append(element[0][0])
                    listofa.append(N1)
                    listofa.append(O1)
                    listofa.append(C1)
                    listofa.append(C2)
                    listofa.append(C3)
                    listofa.append(C4)
                    listofa.append(O2)
                    listofa.append(C5)
                    listofa.append(H15)
                    listofa.append(H14)
                    listofa.append(H13)
                    listofa.append(H12)
                    listofa.append(H11)
                    listofa.append(H10)
                    listofa.append(H9)
                    listofa.append(H8)
                    listofa.append(H7)
                    listofa.append(H6)
                    listofa.append(H5)
                    listofa.append(N2)
                    listofn = [[atom.x, atom.y, atom.z] for atom in listofa]
                    listofr = compose_listofr("H4", listofa)
                    g1_x, g1_y, g1_z = fix_sphere_h(C6.x, C6.y, C6.z, 1.11, listofn, listofr, points_max, [C5, N2])
                    while (0 <= attempt0403 <=  goto):
                        i = random.randint(0, len(g1_x)-1)
                        x1 = g1_x[i]
                        y1 = g1_y[i]
                        z1 = g1_z[i]
                        H4  = Atom(current_size + 22,    'H4', 'P1A', str(resnum + 1), float("{0:.3f}".format(x1)),             float("{0:.3f}".format(y1)),            float("{0:.3f}".format(z1)))
                        atom_list.append(H4)
                        if ((len(identify_bonds(H4, atom_list)) == 1) and (identify_bonds(H4, atom_list)[0][0].atom_number == C6.atom_number)):
                            listofa.append(H4)                            
                            listofn = [[atom.x, atom.y, atom.z] for atom in listofa]
                            listofr = compose_listofr("H3", listofa)
                            g2_x, g2_y, g2_z = fix_sphere_h(C6.x, C6.y, C6.z, 1.11, listofn, listofr, points_max, [C5, N2])
                            if (len(g2_x) != 0):
                                i = random.randint(0, len(g2_x)-1)
                                x2 = g2_x[i]
                                y2 = g2_y[i]
                                z2 = g2_z[i]
                                H3  = Atom(current_size + 23,    'H3', 'P1A', str(resnum + 1), float("{0:.3f}".format(x2)),             float("{0:.3f}".format(y2)),            float("{0:.3f}".format(z2)))
                                atom_list.append(H3)
                                if ((len(identify_bonds(H4, atom_list)) == 1) and (identify_bonds(H4, atom_list)[0][0].atom_number == C6.atom_number) and (len(identify_bonds(H3, atom_list)) == 1) and (identify_bonds(H3, atom_list)[0][0].atom_number == C6.atom_number)):
                                    attempt0403 = 888
                                else:
                                    attempt0403 -= 1      
                                    atom_list.remove(H4)
                                    atom_list.remove(H3)
                                    del H4
                                    del H3
                            else:
                                attempt0403 -= 1
                                atom_list.remove(H4)
                                del H4            
                        else:
                            attempt0403 -= 1
                            atom_list.remove(H4)
                            del H4
                else:
                    attempt_N2 -= 1
                    atom_list.remove(N2)
                    del N2
                    
        attempt0201 = goto
        if ((attempt_N2 > -1) and (attempt_C6 > -1) and (attempt_C5 > -1) and (attempt_O2 > -1) and (attempt_C4 > -1) and (attempt_C3 > -1) and (attempt_O1 > -1) and (attempt_C2 > -1) and (attempt_C1 > -1) and (attempt_N1 > -1) and (attempt0403 > -1)  and (attempt0605 > -1) and (attempt0807 > -1) and (attempt1009 > -1) and (attempt1211 > -1) and (attempt1413 > -1) and (attempt15 > -1)):
            listofa = detect_neighbours_hydrogen(N2, atom_list)
            listofa.append(element[0][0])
            listofa.append(N1)
            listofa.append(C1)
            listofa.append(C2)
            listofa.append(O1)
            listofa.append(C3)
            listofa.append(C4)
            listofa.append(O2)
            listofa.append(C5)
            listofa.append(C6)
            listofa.append(H15)
            listofa.append(H14)
            listofa.append(H13)
            listofa.append(H12)
            listofa.append(H11)
            listofa.append(H10)
            listofa.append(H9)
            listofa.append(H8)
            listofa.append(H7)
            listofa.append(H6)
            listofa.append(H5)
            listofa.append(H4)
            listofa.append(H3)
            listofn = [[atom.x, atom.y, atom.z] for atom in listofa]
            listofr = compose_listofr("H2", listofa)
            g1_x, g1_y, g1_z = fix_sphere_h(N2.x, N2.y, N2.z, 1.03, listofn, listofr, points_max, [C6])
            while (0 <= attempt0201 <= goto):
                i = random.randint(0, len(g1_x)-1)
                x1 = g1_x[i]
                y1 = g1_y[i]
                z1 = g1_z[i]
                H2  = Atom(current_size + 24,    'H2', 'P1A', str(resnum + 1), float("{0:.3f}".format(x1)),             float("{0:.3f}".format(y1)),            float("{0:.3f}".format(z1)))                
                atom_list.append(H2)
                if ((len(identify_bonds(H2, atom_list)) == 1) and (identify_bonds(H2, atom_list)[0][0].atom_number == N2.atom_number)):
                    listofa.append(H2)
                    listofn = [[atom.x, atom.y, atom.z] for atom in listofa]
                    listofr = compose_listofr("H1", listofa)
                    g2_x, g2_y, g2_z = fix_sphere_h(N2.x, N2.y, N2.z, 1.03, listofn, listofr, points_max, [C6])
                    if (len(g2_x) != 0):
                        i = random.randint(0, len(g2_x)-1)
                        x2 = g2_x[i]
                        y2 = g2_y[i]
                        z2 = g2_z[i]
                        H1  = Atom(current_size + 25,    'H1', 'P1A', str(resnum + 1), float("{0:.3f}".format(x2)),             float("{0:.3f}".format(y2)),            float("{0:.3f}".format(z2)))
                        atom_list.append(H1)
                        if ((len(identify_bonds(H2, atom_list)) == 1) and (identify_bonds(H2, atom_list)[0][0].atom_number == N2.atom_number) and (len(identify_bonds(H1, atom_list)) == 1) and (identify_bonds(H1, atom_list)[0][0].atom_number == N2.atom_number)):
                            attempt0201 = 888
                        else:
                            attempt0201 -= 1 
                            atom_list.remove(H2)
                            atom_list.remove(H1)
                            del H2
                            del H1
                    else:
                        attempt0201 -= 1
                        atom_list.remove(H2)
                        del H2
                else:
                    attempt0201 -= 1
                    atom_list.remove(H2)
                    del H2
        
        if ((attempt_N1 > -1) and (attempt_C1 > -1) and (attempt_C2 > -1) and (attempt_O1 > -1) and (attempt_C3 > -1) and (attempt_C4 > -1) and (attempt_O2 > -1) and (attempt_C5 > -1) and (attempt_C6 > -1) and (attempt_N2 > -1) and (attempt15 > -1) and (attempt1413 > -1) and (attempt1211 > -1) and (attempt1009 > -1) and (attempt0807 > -1) and (attempt0605 > -1) and (attempt0403 > -1) and (attempt0201 > -1)):
            print("Placed.")
        else:
            print("Not placed.")
            try:
                atom_list.remove(N1)
                del N1
                atom_list.remove(C1)
                del C1
                atom_list.remove(H15)
                del H15
                atom_list.remove(C2)
                del C2
                atom_list.remove(H14)
                del H14
                atom_list.remove(H13)
                del H13
                atom_list.remove(O1)
                del O1
                atom_list.remove(H12)
                del H12
                atom_list.remove(H11)
                del H11
                atom_list.remove(C3)
                del C3
                atom_list.remove(C4)
                del C4
                atom_list.remove(H10)
                del H10
                atom_list.remove(H9)
                del H9
                atom_list.remove(C5)
                del C5
                atom_list.remove(H8)
                del H8
                atom_list.remove(H7)
                del H7
                atom_list.remove(O2)
                del O2
                atom_list.remove(H6)
                del H6
                atom_list.remove(H5)
                del H5
                atom_list.remove(C6)
                del C6
                atom_list.remove(H4)
                del H4
                atom_list.remove(H3)
                del H3
                atom_list.remove(N2)
                del N2
                atom_list.remove(H2)
                del H2
                atom_list.remove(H1)
                del H1
            except:
                print()
    writepdb(atom_list, filename1)  
    return 'done.'
    
def lw(max_no, str_obj):
    """Used to add whitespaces to fit the PDB format."""
    x = max_no - len(str_obj)
    y = 0
    string = ''
    for y in range(x):
        string = string + ' '
    return string

def writepdb(list_of_atoms, filename1):
    """PDB writer"""
    os.chdir(os.getcwd())
    with open(str(filename1), 'a') as le_file:
        for atom in list_of_atoms:
            line = "ATOM" + lw(7, str(atom.atom_number)) + str(atom.atom_number) + lw(4, str(atom.atom_name)) + str(atom.atom_name) + "  " + str(atom.residue_name) + lw(6, str(atom.residue_number)) + str(atom.residue_number) + lw(12, str(atom.x)) + str(atom.x) + lw(8, str(atom.y)) + str(atom.y) + lw(8, str(atom.z)) + str(atom.z) + "  1.00  0.00             "
            le_file.write(line + '\n')

def find_rings(atom_list):
    """Returns a list of (x1, x2, x3, x4, x5, x6) representing atoms making up a graphene ring."""    
    CX_list = [atom0 for atom0 in atom_list if ((atom0.atom_name == "CX") or (atom0.atom_name == "CY"))]
    atom_dict = {}
    for atom0 in CX_list:
        if (len(identify_bonds(atom0, atom_list)) >= 2):
            atom_dict[atom0] = {}
            for atom1 in identify_bonds(atom0, atom_list):
                if ( ((atom1[0].atom_name == "CX") or (atom1[0].atom_name == "CY")) and (len(identify_bonds(atom1[0], atom_list)) >= 2) ):
                    atom_dict[atom0][atom1[0]] = {}
                    for atom2 in identify_bonds(atom1[0], atom_list):
                        if ( ((atom2[0].atom_name == "CX") or (atom2[0].atom_name == "CY")) and (atom2[0] != atom0) and (len(identify_bonds(atom2[0], atom_list)) >= 2)):
                            atom_dict[atom0][atom1[0]][atom2[0]] = {}
                            for atom3 in identify_bonds(atom2[0], atom_list):
                                if ( ((atom3[0].atom_name == "CX") or (atom3[0].atom_name == "CY")) and (atom3[0] != atom0) and (len(identify_bonds(atom3[0], atom_list)) >= 2)):
                                    atom_dict[atom0][atom1[0]][atom2[0]][atom3[0]] = [atom3[0].atom_number]
    rings = []
    for key in atom_dict.keys():
        for key2 in atom_dict[key].keys():
            for key3 in atom_dict[key][key2].keys():
                for key4 in atom_dict[key][key2][key3].keys():
                    rings.append([key, key2, key3, key4])
    finite_rings = []
    for element in rings:
        for element2 in rings:
            if ((element[0] == element2[0]) and (element[3] == element2[3]) and (element[1] != element2[1]) and (element[1] != element2[2]) and (element[2] != element2[1]) and (element[2] != element2[2]) and (element[0] != element2[1] != element[3]) and (element[0] != element2[2] != element[3])):
                check = True
                for el in finite_rings:
                    if ((element[0] in el) and (element[1] in el) and (element[2] in el) and (element[3] in el) and (element2[0] in el) and (element2[1] in el) and (element2[2] in el) and (element2[3] in el)):
                        check = False
                if (check == True):
                    finite_rings.append([element[0], element[1], element[2], element[3], element2[1], element2[2]])
    return finite_rings

def filter_carbon_atoms(atom_list, rings):
    """Determines which atoms are available to be replaced for N_graphitic, N-Pyridinic and N_Pyrrolic atoms. """
    list_3 = []
    list_2 = []
    list_2n = []
    for atom in atom_list:
        if (check_connected(atom, identify_bonds(atom, atom_list)) == False):
            if (len(identify_bonds(atom, atom_list)) == 3):
                list_3.append(atom)
            elif (len(identify_bonds(atom, atom_list)) == 2):
                list_2.append(atom)
                for neighbour in identify_bonds(atom, atom_list):
                    if (len(identify_bonds(neighbour[0], atom_list)) == 2):
                        for ring in rings:
                            if( (atom in ring) and (neighbour[0] in ring)):
                                list_2n.append(atom)                    
    return list_3, list_2, list_2n
            
def generate_N_doping(path, N_graphitic, N_pyridinic, N_pyrrolic, filename1):
    """Generate an N-doped layer.
    N3A - residue name used for N_graphitic.
    N2A - residue name used for N_pyridinic.
    N2N - residue name used for N_pyrrolic."""
    global bond_list
    bond_list = bond_list_1 + bond_list_3
    atom_list = read_in_graphene(path)
    rings = find_rings(atom_list)
    bond_list = bond_list_1 + bond_list_3
    map_3, map_2, map_2n = filter_carbon_atoms(atom_list, rings)
    graphitic = N_graphitic 
    pyridinic = N_pyridinic
    pyrrolic = N_pyrrolic
    attempt = len(atom_list) / 10
    choices = [1, 2, 3]
    while (((N_graphitic > 0) or (N_pyridinic > 0) or (N_pyrrolic > 0)) and (attempt > 0)):
        print("Left to add: ", "N_graphitic ", N_graphitic, "N_pyridinic ", N_pyridinic, "N_pyrrolic ", N_pyrrolic)
        if (N_graphitic == 0):
            try:
                choices.remove(1)
            except:
                pass
        if (N_pyridinic == 0):
            try:
                choices.remove(2)
            except:
                pass
        if (N_pyrrolic == 0):
            try:
                choices.remove(3)
            except:
                pass
        choice = random.choice(choices)
        if (choice == 1):
            while ((N_graphitic > 0) and (len(map_3) > 0)):
                random_atom = random.choice(map_3)
                N_graphitic -= 1
                N = Atom(random_atom.atom_number, "N3", "N3A", str(graphitic - N_graphitic), float("{0:.3f}".format(random_atom.x)), float("{0:.3f}".format(random_atom.y)), float("{0:.3f}".format(random_atom.z)))
                if ((len(identify_bonds(random_atom, atom_list)) == 3) and ((identify_bonds(random_atom, atom_list)[0][0].atom_name == "CX") or (identify_bonds(random_atom, atom_list)[0][0].atom_name == "CY")) and ((identify_bonds(random_atom, atom_list)[1][0].atom_name == "CX") or identify_bonds(random_atom, atom_list)[1][0].atom_name == "CY") and ((identify_bonds(random_atom, atom_list)[2][0].atom_name == "CX") or (identify_bonds(random_atom, atom_list)[2][0].atom_name == "CY"))):
                    for ring in rings:
                        if (random_atom in ring):
                            for atom in ring:
                                try:
                                    map_3.remove(atom)
                                except:
                                    pass
                                try:
                                    map_2.remove(atom)
                                except:
                                    pass
                                try:
                                    map_2n.remove(atom)
                                except:
                                    pass
                    try:
                        atom_list.remove(random_atom)
                    except:
                        pass
                    atom_list.append(N)
                else:
                    attempt -= 1
        elif (choice == 2):
            while ((N_pyridinic > 0) and (len(map_2) > 0)): 
                random_atom = random.choice(map_2)
                N_pyridinic -= 1
                N = Atom(random_atom.atom_number, "N2", "N2A", str(pyridinic - N_pyridinic), float("{0:.3f}".format(random_atom.x)), float("{0:.3f}".format(random_atom.y)), float("{0:.3f}".format(random_atom.z)))
                if ((len(identify_bonds(random_atom, atom_list)) == 2) and ((identify_bonds(random_atom, atom_list)[0][0].atom_name == "CX") or (identify_bonds(random_atom, atom_list)[0][0].atom_name == "CY")) and ((identify_bonds(random_atom, atom_list)[1][0].atom_name == "CX") or identify_bonds(random_atom, atom_list)[1][0].atom_name == "CY") ):
                    found = False
                    for ring in rings:
                        if (random_atom in ring):
                            found = True
                            for atom in ring:
                                try:
                                    map_3.remove(atom)
                                except:
                                    pass
                                try:
                                    map_2.remove(atom)
                                except:
                                    pass
                                try:
                                    map_2n.remove(atom)
                                except:
                                    pass
                    if (found == False):
                        try:
                            map_3.remove(random_atom)
                        except:
                            pass
                        try:
                            map_2.remove(random_atom)
                        except:
                            pass
                        try:
                            map_2n.remove(random_atom)
                        except:
                            pass
                    atom_list.remove(random_atom)
                    atom_list.append(N)
                else:
                    attempt -= 1
            else: 
                attempt -= 1
        elif (choice == 3):
            while ((N_pyrrolic > 0) and (len(map_2n) > 0)):
                random_atom_1 = random.choice(map_2n)
                for neighbour in identify_bonds(random_atom_1, atom_list):
                    if (len(identify_bonds(neighbour[0], atom_list)) == 2):
                        random_atom_2 = neighbour[0]
                        break
                for ring in rings:
                    if (random_atom_1 in ring):
                        center_6 = {}
                        center_6['x'] = 0
                        center_6['y'] = 0
                        center_6['z'] = 0
                        center_4 = {}
                        center_4['x'] = 0
                        center_4['y'] = 0
                        center_4['z'] = 0
                        for atom in ring:
                            center_6['x'] += atom.x
                            center_6['y'] += atom.y
                            center_6['z'] += atom.z
                            if ((atom != random_atom_1) and (atom != random_atom_2)):
                                center_4['x'] += atom.x
                                center_4['y'] += atom.y
                                center_4['z'] += atom.z
                center_6['x'] /= 6
                center_6['y'] /= 6
                center_6['z'] /= 6
                center_4['x'] /= 4
                center_4['y'] /= 4
                center_4['z'] /= 4
                N_pyrrolic -= 1
                p = 0.6
                limit = 0.3
                if ((-limit < center_4['x'] - center_6['x'] < limit) and (-limit < center_4['y'] - center_6['y'] < limit)): 
                    N = Atom(random_atom_1.atom_number, "N1", "N2N", str(pyrrolic - N_pyrrolic), float("{0:.3f}".format(center_6['x'])), float("{0:.3f}".format(center_6['y'])), float("{0:.3f}".format(center_6['z'])))   
                elif ((-limit < center_4['x'] - center_6['x'] < limit) and (center_4['y'] - center_6['y'] < -limit)):
                    N = Atom(random_atom_1.atom_number, "N1", "N2N", str(pyrrolic - N_pyrrolic), float("{0:.3f}".format(center_6['x'])), float("{0:.3f}".format(center_6['y'] + p/2)), float("{0:.3f}".format(center_6['z'])))   
                elif ((-limit < center_4['x'] - center_6['x'] < limit) and (center_4['y'] - center_6['y'] > limit)):
                    N = Atom(random_atom_1.atom_number, "N1", "N2N", str(pyrrolic - N_pyrrolic), float("{0:.3f}".format(center_6['x'])), float("{0:.3f}".format(center_6['y'] - p/2)), float("{0:.3f}".format(center_6['z'])))                   
                elif ((center_4['x'] - center_6['x'] < -limit) and (-limit < center_4['y'] - center_6['y'] < limit)):
                    N = Atom(random_atom_1.atom_number, "N1", "N2N", str(pyrrolic - N_pyrrolic), float("{0:.3f}".format(center_6['x'] + p)), float("{0:.3f}".format(center_6['y'])), float("{0:.3f}".format(center_6['z'])))   
                elif ((center_4['x'] - center_6['x'] < -limit) and (center_4['y'] - center_6['y'] < -limit)):
                    N = Atom(random_atom_1.atom_number, "N1", "N2N", str(pyrrolic - N_pyrrolic), float("{0:.3f}".format(center_6['x'] + p)), float("{0:.3f}".format(center_6['y'] + p/2)), float("{0:.3f}".format(center_6['z'])))   
                elif ((center_4['x'] - center_6['x'] < -limit) and (center_4['y'] - center_6['y'] > limit)):
                    N = Atom(random_atom_1.atom_number, "N1", "N2N", str(pyrrolic - N_pyrrolic), float("{0:.3f}".format(center_6['x'] + p)), float("{0:.3f}".format(center_6['y'] - p/2)), float("{0:.3f}".format(center_6['z'])))               
                elif ((center_4['x'] - center_6['x'] > limit) and (-limit < center_4['y'] - center_6['y'] < limit)):
                    N = Atom(random_atom_1.atom_number, "N1", "N2N", str(pyrrolic - N_pyrrolic), float("{0:.3f}".format(center_6['x'] - p)), float("{0:.3f}".format(center_6['y'])), float("{0:.3f}".format(center_6['z'])))   
                elif ((center_4['x'] - center_6['x'] > limit) and (center_4['y'] - center_6['y'] < -limit)):
                    N = Atom(random_atom_1.atom_number, "N1", "N2N", str(pyrrolic - N_pyrrolic), float("{0:.3f}".format(center_6['x'] - p)), float("{0:.3f}".format(center_6['y'] + p/2)), float("{0:.3f}".format(center_6['z'])))   
                elif ((center_4['x'] - center_6['x'] > limit) and (center_4['y'] - center_6['y'] > limit)):
                    N = Atom(random_atom_1.atom_number, "N1", "N2N", str(pyrrolic - N_pyrrolic), float("{0:.3f}".format(center_6['x'] - p)), float("{0:.3f}".format(center_6['y'] - p/2)), float("{0:.3f}".format(center_6['z'])))   
                for ring in rings:
                    if (random_atom_1 in ring):
                        for atom in ring:
                            try:
                                map_3.remove(atom)
                            except:
                                pass
                            try:
                                map_2.remove(atom)
                            except:
                                pass
                            try:
                                map_2n.remove(atom)
                            except:
                                pass
                            for mol in identify_bonds(atom, atom_list):
                                try:
                                    map_2n.remove(mol[0])
                                except:
                                    pass
                try:
                    atom_list.remove(random_atom_1)
                    atom_list.remove(random_atom_2)
                except:
                    pass
                atom_list.append(N)
            else:
                attempt -= 1
    attempt -= 1
    writepdb(atom_list, filename1)
    print("done.")
    return 'done.'
        
def get_contour(atom_list):
    """Returns a list of atoms representing atoms on the contour/edge of the graphenic layer and of the existing holes.
    Used to prevent new holes touching existing ones."""
    initial = [atom for atom in atom_list if ((0 < len(identify_bonds(atom, atom_list)) < 3) and (check_connected(atom, identify_bonds(atom, atom_list)) == False))]
    
    extra_1 = []
    for atom in atom_list:
        neighbours = [bond[0] for bond in identify_bonds(atom, atom_list)]
        for i in neighbours:
            neighbours2 = [bond[0] for bond in identify_bonds(i, atom_list)]
            for j in neighbours2:
                if j in initial:
                    extra_1.append(atom)

    extra_2 = []
    for atom in atom_list:
        neighbours = [bond[0] for bond in identify_bonds(atom, atom_list)]
        check = 0
        for i in neighbours:
            if i in initial:
                check += 1
        if ((check == 2) and (atom not in initial)):
            extra_2.append(atom)    
    return (initial + extra_1 + extra_2)
    

def hole_cleanup(atom_list):
    """After hole creation, some small fragments (<=6 atoms) can remain entirely isolated. These are removed using this function.""" 
    joey = atom_list.copy()
    while (len(joey) != 0):
        for atom in joey:
            takein = [atom]
            source_update = takein.copy()
            check = 1
            while (check == 1):
                source = source_update.copy()
                source_update = []
                c = len(takein)
                for element in source:
                    bonds = [bond[0] for bond in identify_bonds(element, joey) if bond[0] not in takein]
                    for h in bonds:
                        takein.append(h)
                        source_update.append(h)
                if ((len(takein) == c) and (len(takein) < 6)):
                    check = 0
                    for element in takein:
                        atom_list.remove(element)
                elif (len(takein) == c):
                    check = 0
            for element in takein:
                joey.remove(element)
    return atom_list

def find_contour(hole_atoms, atom_list): 
    """Returns a list of atoms representing atoms on the contour/edge of the graphenic layer and of the existing holes.
    Used to prevent new holes touching existing ones."""
    contour_atoms = []
    extra_atoms = []
    global bond_list
    bond_list = bond_list_1
    for atom in hole_atoms:
        c = [bond[0] for bond in identify_bonds(atom, atom_list) if ((bond[0] not in hole_atoms) and (bond[0] not in contour_atoms))]
        for element in c:
            contour_atoms.append(element)
    for atom in atom_list:
        c = [bond[0] for bond in identify_bonds(atom, atom_list)]
        count = 0
        for element in c:
            if element in contour_atoms:
                count += 1
        if (count >= 2):
            extra_atoms.append(atom)
    for atom in atom_list:
        c = [bond[0] for bond in identify_bonds(atom, atom_list)]
        for element in c:
            if ((element in contour_atoms) or (element in extra_atoms)):
                for i in [bond[0] for bond in identify_bonds(element, atom_list)]:
                    if ((i in hole_atoms) and (atom not in hole_atoms) and (atom not in contour_atoms) and (atom not in extra_atoms)):
                        extra_atoms.append(atom)                        
    
    contour_atoms = contour_atoms + extra_atoms
    
    extra_atoms2 = []
    for atom in contour_atoms:
        for atom2 in contour_atoms:
            if (atom != atom2):
                c = [bond[0] for bond in identify_bonds(atom, atom_list) if ((bond in identify_bonds(atom2, atom_list)) and (bond[0] not in (contour_atoms)))]
                if (len(c) != 0):
                    extra_atoms2.append(c[0])                       
    for element in extra_atoms2:
        contour_atoms.append(element)
    return contour_atoms
                
def hole_generation(init_file, number_of_holes, atom_range, uorm, iore, cleanup, filename1):
    """Hole generation function.
    uorm - "u"nidirectional or "m"ultidirectional. 
    iore - "i"nterior (not touching edges) or "e"xterior 
    cleanup - "a" if active."""
    global atoms
    global bond_list
    bond_list = bond_list_1
    atoms = read_in_graphene(init_file)
    holes = 0
    list_of_contours = get_contour(atoms)    
    takein = {}
    while (holes < number_of_holes):
        print("Holes: ", holes, "placed out of", number_of_holes)
        if (iore == "e"):
            if (uorm == "u"):
                attempt = 0
                while (attempt < 50):
                    takein[holes] = []
                    compare = random.randint(atom_range[0], atom_range[1])
                    no_atoms = compare - 1
                    available = [atom1 for atom1 in get_map_anywhere(atoms)]
                    if (len(available) > 0):
                        random_atom = random.choice(available)
                        takein[holes].append(random_atom)
                        source = random_atom
                        while (no_atoms > 0):
                            neighbours = [bond[0] for bond in identify_bonds(source, atoms) if (bond[0] not in takein[holes])]
                            if (len(neighbours) > 0):
                                choice = random.choice(neighbours)
                                takein[holes].append(choice)
                                no_atoms -= 1
                                source = choice
                            else:
                                if (len(takein[holes]) == 1):
                                    random_atom = random.choice(available)
                                    takein[holes] = []
                                    takein[holes].append(random_atom)
                                    source = random_atom
                                    attempt += 1
                                    continue
                                else:
                                    check = 0
                                    while (check == 0):
                                        for atom in takein[holes]:
                                            neighbours_reset = [bond[0] for bond in identify_bonds(atom, atoms) if (bond[0] not in takein[holes])]
                                            if (len(neighbours_reset) != 0):
                                                source = atom
                                                check = 1
                                                break
                                        break
                                    if (check == 0):
                                        takein[holes] = []
                                        attempt += 1
                                        break
                        if (len(takein[holes]) == compare):
                            holes += 1
                            attempt = 888
                    else:
                        attempt += 1
                for key in takein.keys():
                    for atom in takein[key]:
                        if ((atom in atoms) and (atom not in list_of_contours)):
                            atoms.remove(atom)    
            elif (uorm == "m"):
                attempt = 0
                while (attempt < 50):
                    takein[holes] = []
                    compare = random.randint(atom_range[0], atom_range[1])
                    no_atoms = compare - 1
                    available = [atom1 for atom1 in get_map_anywhere(atoms)]
                    if (len(available) > 0):
                        random_atom = random.choice(available)
                        takein[holes].append(random_atom)
                        source_update = [random_atom]
                        while (no_atoms > 0):
                            source = source_update.copy()
                            source_update = []
                            check = 0
                            for source_atom in source:
                                neighbours = [bond[0] for bond in identify_bonds(source_atom, atoms) if (bond[0] not in takein[holes])]
                                if (len(neighbours) > 0):
                                    for elem in neighbours:
                                        if (no_atoms > 0):
                                            no_atoms -= 1
                                            takein[holes].append(elem)
                                            source_update.append(elem)
                                            check += 1
                            if (check == 0):
                                break
                        if (len(takein[holes]) == compare):
                            holes += 1
                            attempt = 888
                    else:
                        attempt += 1
                for key in takein.keys():
                    for atom in takein[key]:
                        if ((atom in atoms) and (atom not in list_of_contours)):
                            atoms.remove(atom)    
        elif (iore == "i"):
            if (uorm == "u"):
                attempt = 0
                while (attempt < 50):
                    takein[holes] = []
                    compare = random.randint(atom_range[0], atom_range[1])
                    no_atoms = compare - 1
                    available = [atom1 for atom1 in get_map_anywhere(atoms) if (atom1 not in list_of_contours)]
                    if (len(available) > 0):
                        random_atom = random.choice(available)
                        takein[holes].append(random_atom)
                        source = random_atom
                        while (no_atoms > 0):
                            neighbours = [bond[0] for bond in identify_bonds(source, atoms) if ((bond[0] not in list_of_contours) and (bond[0] not in takein[holes]))]
                            if (len(neighbours) > 0):
                                choice = random.choice(neighbours)
                                takein[holes].append(choice)
                                no_atoms -= 1
                                source = choice
                            else:
                                if (len(takein[holes]) == 1):
                                    random_atom = random.choice(available)
                                    takein[holes] = []
                                    takein[holes].append(random_atom)
                                    source = random_atom
                                    attempt += 1
                                    continue
                                else:
                                    check = 0
                                    while (check == 0):
                                        for atom in takein[holes]:
                                            neighbours_reset = [bond[0] for bond in identify_bonds(atom, atoms) if ((bond[0] not in list_of_contours) and (bond[0] not in takein[holes]))]
                                            if (len(neighbours_reset) != 0):
                                                source = atom
                                                check = 1
                                                break
                                        break
                                    if (check == 0):
                                        takein[holes] = []
                                        attempt += 1
                                        break                                                           
                        if (len(takein[holes]) == compare):
                            for element in find_contour(takein[holes], atoms):
                                list_of_contours.append(element)
                            holes += 1
                            attempt = 888
                    else: 
                        attempt += 1
            elif (uorm == "m"):
                attempt = 0
                while (attempt < 50):
                    takein[holes] = []
                    compare = random.randint(atom_range[0], atom_range[1])
                    no_atoms = compare - 1
                    not_in_takein = []
                    for key in takein.keys():
                        for element in takein[key]:
                            not_in_takein.append(element)
                    available = [atom1 for atom1 in get_map_anywhere(atoms) if ((atom1 not in list_of_contours) and (atom1 not in not_in_takein))]
                    if (len(available) > 0):
                        random_atom = random.choice(available)
                        timer = 50
                        while ((random_atom in list_of_contours) and (timer > 0)):
                            random_atom = random.choice(available)
                            timer -= 1
                        takein[holes].append(random_atom)
                        source_update = [random_atom]
                        while (no_atoms > 0):
                            # your 'source' is always new here
                            source = source_update.copy()
                            source_update = []
                            check = 0
                            for source_atom in source:
                                neighbours = [bond[0] for bond in identify_bonds(source_atom, atoms) if ((bond[0] not in list_of_contours) and (bond[0] not in takein[holes]))]
                                if (len(neighbours) > 0):
                                    for elem in neighbours:
                                        if ((no_atoms > 0) and (random.choice([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]) > 0)):
                                            no_atoms -= 1
                                            takein[holes].append(elem)
                                            source_update.append(elem)
                                            check += 1
                            if (check == 0):
                                break
                        if (len(takein[holes]) == compare):
                            for element in find_contour(takein[holes], atoms):
                                list_of_contours.append(element)
                            holes += 1
                            attempt = 888
                    else:
                        attempt += 1
            for key in takein.keys():
                for atom in takein[key]:
                    if ((atom in atoms) and (atom not in list_of_contours)):
                        atoms.remove(atom)
    print("Holes: ", holes, "placed out of", number_of_holes)
    if (cleanup == "a"):
        print("Hole cleanup...")
        atoms = hole_cleanup(atoms)
    
    new_list = []
    atom_list = atoms.copy()
    atno = 1
    atom_list_copy = atom_list.copy()
    for atom in atom_list_copy:
        if (atom.atom_name == "CX"):
            New_CX = Atom(atno, "CX", "GGG", atno, atom.x, atom.y, atom.z)
            atom_list.append(New_CX)
            new_list.append(New_CX)
            atom_list.remove(atom)
            del atom
            atno += 1 

    atoms = new_list.copy()
    writepdb(atoms, filename1)
    print('done.')
    return "done."

if (len(sys.argv) > 1):
    if (sys.argv[1] == "generate_PG"):
        try:
            generate_pristine_graphene(int(sys.argv[2]), int(sys.argv[3]), str(sys.argv[4]))
        except:
            print("Please type 'python GOPY.py help'")
    elif (sys.argv[1] == "generate_GO"):
        try:
            for element in sys.argv:
                print(element)
            create_GO(str(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), str(sys.argv[6]))   
        except:
            print("Please type 'python GOPY.py help'")
    elif (sys.argv[1] == "generate_rGO_PEG_NH2"):
        try:
            add_NH_PEG_NH2(str(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), str(sys.argv[6]))
        except:
            print("Please type 'python GOPY.py help'")            
    elif (sys.argv[1] == "generate_hole"):
#        try:
        hole_generation(str(sys.argv[2]), int(sys.argv[3]), [int(sys.argv[4]), int(sys.argv[5])], str(sys.argv[6]), str(sys.argv[7]), str(sys.argv[8]), str(sys.argv[9])) 
#        except:
#            print("Please type 'python GOPY.py help'")
    elif (sys.argv[1] == "generate_N_doped"):
        try:
            generate_N_doping(str(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), str(sys.argv[6]))
        except:
            print("Please type 'python GOPY.py help'")
    elif (sys.argv[1] == "help"):
        print("""Help is here. You can use GOPY in the following manner to generate graphene-based 2D PDB models:
'python GOPY.py generate_PG X Y file_to_save'                                  - Used to generate a pristine graphene layer. The X and Y dimensions are required in Å. The file_to_save
                                                                               represents the name under which the new PDB file should be saved.
                                                                               E.g. 'python GOPY.py generate_PG 30 30' will generate a 3 nm x 3 nm PG layer.
'python GOPY.py generate_GO path_to_file X Y Z file_to_save'                   - Used to generate a graphene oxide layer. The 'path_to_file' points to an
                                                                               existing pristine graphene PDB file. You may first generate a PG layer using GOPY.
                                                                               X corresponds to the desired number of carboxyl functional groups,
                                                                               Y corresponds to the number of epoxy groups and Z corresponds to the number of hydroxyl groups.
                                                                               The file_to_save represents the name under which the new PDB file should be saved.
                                                                               E.g. 'python GOPY.py /path/to/PG.pdb 30 60 60' will generate a GO layer, attempting to place
                                                                               30 carboxyl groups, 60 epoxy groups and 60 hydroxyl groups.
'python GOPY.py generate_rGO_PEG_NH2 path_to_file X Y Z file_to_save'          - Used to generate a rGO-PEG-NH2 layer. The 'path_to_file' points to an
                                                                               existing graphene oxide (GO) PDB file. X Y and Z represent the percentages of carboxyl,
                                                                               epoxy, respectively hydroxyl groups to be removed. All removed hydroxyl and epoxy groups
                                                                               that were removed will be replaced by the PEG-NH2 chains. The file_to_save
                                                                               represents the name under which the new PDB file should be saved.
                                                                               E.g. 'python GOPY.py /path/to/GO.pdb 0.6 1 1 PEG-NH2.pdb' removes 60% of carboxyl groups
                                                                               and 100% of epoxy and hydroxyl groups.
'python GOPY.py generate_hole path_to_file N R1 R2 ARG1 ARG2 C file_to_save'   - Used to generate holes in a PG layer. N represents the number of holes to be created. 
                                                                               R1, R2 represent the range expressed as a list such as "[x, y]", ARG1 should be either
                                                                               "u" (uni-directional) or "m" (multi-directional), C should be "a" if cleanup should
                                                                               be performed and file_to_save represents the name under which to save the new file.
                                                                               E.g. 'python GOPY.py /path/to/PG.pdb 10 11 20 m i a hole.pdb' reads in the PG.pdb file and
                                                                               attempts to create 10 holes with a size between 11 and 20 atoms in a multi-directional manner,
                                                                               not allowing the holes to touch edges. Cleanup is performed at the end and the file is saved 
                                                                               as hole.pdb.
'python GOPY.py generate_N_doped path_to_file 10 9 8 file_to_save'              - Used to generate an N-doped graphene layer. X Y and Z represent the number of N-graphitic,
                                                                               N-pyridinic and N-pyrrolic atoms respectively.
                                                                               E.g. 'python GOPY.py /path/to/PG.pdb 10 10 10 Ndoped.pdb' generates an N-doped molecule of the 
                                                                               PG layer given as input with 10 N-graphitic atoms, 9 N-pyridinic atoms and 8 N-pyrrolic atoms.            
            """)
    else:
        print(sys.argv[1], " is not recognized. Please type 'python GOPY.py help' for instructions.")
else:
    print("No arguments provided!")


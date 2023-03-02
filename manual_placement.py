# from my_modules.project import*
from ase.io import read
from ase.visualize import view
from ase.geometry.analysis import Analysis
from ase.atoms import Atom, Atoms
import numpy as np
from scipy.optimize import leastsq

#### CONSTANTS FOUND FROM INITIAL SCREENING ####
TOP_CO_BL = 1.1660185920095851
TOP_CNI_BL = 1.735016295455957
BRIDGE_CO_BL = 1.18634264175913
BRIDGE_CNI_BL = 1.8634753904685644
HOLLOW_CO_BL = 1.2008948988726278
HOLLOW_CNI_BL = 1.946431063991434


def get_neighbours(idx, bonds):
    neighbours = []
    for bond in bonds:
        if bond[0] == idx:
            neighbours.append(bond[1])
        if bond[1] == idx:
            neighbours.append(bond[0])
    return np.array(neighbours)

def get_norm_vec(vec):
    return vec/np.linalg.norm(vec)

def get_best_fit_plane(positions):
    XYZ = positions
    # Inital guess of the plane
    p0 = [1, 1, 1 ,1]

    def f_min(X,p):
        plane_xyz = p[0:3]
        distance = (plane_xyz*X).sum(axis=1) + p[3]
        return distance / np.linalg.norm(plane_xyz)

    def residuals(params, signal, X):
        return f_min(X, params)

    sol = leastsq(residuals, p0, args=(None, XYZ))[0]
    return sol


def get_plane(idx, positions, bonds):
    neighbours = get_neighbours(idx, bonds)
    sol = get_best_fit_plane(positions[neighbours,:])

    side = positions[idx,:]@sol[0:3] + sol[3]
    perp_dir = sol[0:3]

    if side < 0:
        perp_dir = - perp_dir
    norm_perp_dir = get_norm_vec(perp_dir)

    return norm_perp_dir


def get_bonds(template : Atoms):
    """Returns the Ni-Ni bonds"""
    anal = Analysis(template)
    bonds = anal.get_bonds("Ni","Ni")[0]
    return bonds

def get_CNibonds(template : Atoms):
    """Returns the C-Ni bonds"""
    anal = Analysis(template)
    bonds = anal.get_bonds("C","Ni")[0]
    return bonds

def get_ONibonds(template : Atoms):
    """Returns the O-Ni bonds"""
    anal = Analysis(template)
    bonds = anal.get_bonds("O","Ni")[0]
    return bonds


def get_top_atom(template : Atoms, idx, symbol, distance):
    positions = template.positions
    bonds = get_bonds(template)
    norm_perp_dir = get_plane(idx,positions,bonds )
    add_position = positions[idx,:] + norm_perp_dir*distance
    return Atom(symbol, add_position)

def add_top_atom(template : Atoms, idx, symbol, distance):
    add_atom = get_top_atom(template, idx, symbol, distance)
    template.append(add_atom)

def get_top_CO(template : Atoms, idx, distance = TOP_CNI_BL, CO_distance = TOP_CO_BL, C_first = True):
    add_O = None
    add_C = None
    if C_first:
        C_distance = distance
        O_distance = distance + CO_distance 
    else :
        C_distance = distance+ CO_distance 
        O_distance = distance 
    add_O = get_top_atom(template, idx, "O", O_distance)
    add_C = get_top_atom(template, idx, "C", C_distance)
    return Atoms([add_O, add_C])

def add_top_CO(template : Atoms, idx, distance = TOP_CNI_BL, CO_distance = TOP_CO_BL, C_first = True):
    add_atoms = get_top_CO(template, idx,distance, CO_distance, C_first )
    template.append(add_atoms[0])
    template.append(add_atoms[1])

def calculate_center(position1, position2):
    return (position1 + position2)/2

def get_center(idx1, idx2, positions):
    position1 = positions[idx1,:]
    position2 = positions[idx2,:]
    return calculate_center(position1, position2)

def get_bridge_vec(bonds, idx1, idx2, positions):
    perp_dir_1 = get_plane(idx1, positions, bonds)[0:3] # do not need constant
    perp_dir_2 = get_plane(idx2, positions, bonds)[0:3] # do not need constant
    best_dir = perp_dir_1+perp_dir_2
    best_dir = get_norm_vec(best_dir)


    position1 = positions[idx1,:]
    position2 = positions[idx2,:]
    distance_vec = position1 - position2
    best_perp_dir = best_dir - (distance_vec@best_dir)*distance_vec/(distance_vec@distance_vec)
    best_perp_dir = get_norm_vec(best_perp_dir)
    
    return best_perp_dir

def get_bridge_length(bonds, idx1, idx2, positions, distance):
    perp_dir_1 = get_plane(idx1, positions, bonds)[0:3] # do not need constant
    perp_dir_2 = get_plane(idx2, positions, bonds)[0:3] # do not need constant
    best_dir = perp_dir_1+perp_dir_2
    best_dir = get_norm_vec(best_dir)

    position1 = positions[idx1,:]
    position2 = positions[idx2,:]
    distance_vec = position1 - position2
    length = np.linalg.norm(distance_vec)
    distance_from_bond = np.sqrt(distance**2 - length**2 / 4)
    return distance_from_bond
    

def get_bridge_atom(template : Atoms, idx1, idx2, symbol, distance):
    positions = template.positions
    anal = Analysis(template)
    bonds = anal.get_bonds("Ni","Ni")[0]
    neighbours1 = get_neighbours(idx1, bonds)
    if not idx2 in neighbours1:
        raise Exception(f"{idx1} and {idx2} are not neighbours.")
    
    center = get_center(idx1, idx2, positions)
    best_perp_dir = get_bridge_vec(bonds, idx1, idx2, positions)
    distance_from_bond = get_bridge_length(bonds, idx1, idx2, positions, distance)

    position = center + distance_from_bond*best_perp_dir
    add_atom = Atom(symbol, position)
    return add_atom

def add_bridge_atom(template : Atoms, idx1, idx2, symbol, distance):
    add_atom = get_bridge_atom(template, idx1, idx2, symbol, distance)
    template.append(add_atom)

def get_bridge_CO(template : Atoms, idx1, idx2, distance = BRIDGE_CNI_BL, CO_distance = BRIDGE_CO_BL, C_first = True):
    positions = template.positions
    bonds = get_bonds(template)
    center = get_center(idx1, idx2, positions)
    best_perp_dir = get_bridge_vec(bonds, idx1, idx2, positions)
    distance_from_bond = get_bridge_length(bonds, idx1, idx2, positions, distance)
    second_add_position = center + (distance_from_bond + CO_distance)*best_perp_dir
    if C_first:
        add_C = get_bridge_atom(template, idx1, idx2, "C", distance)
        add_O = Atom("O", second_add_position)
    else :
        add_O = get_bridge_atom(template, idx1, idx2, "O", distance)
        add_C = Atom("C", second_add_position)
    
    return Atoms([add_O, add_C])


def add_bridge_CO(template : Atoms, idx1, idx2, distance = BRIDGE_CNI_BL, CO_distance = BRIDGE_CO_BL, C_first = True):
    add_atoms = get_bridge_CO(template, idx1, idx2, distance, CO_distance, C_first)
    template.append(add_atoms[0])
    template.append(add_atoms[1])


def get_hollow_atom(template : Atoms, idx1 : int, idx2 : int, idx3 : int, symbol : str, distance : float):
    positions = template.positions
    position1 = positions[idx1,:]
    position2 = positions[idx2,:]
    position3 = positions[idx3,:]
    
    vec_12 = -position1+position2
    vec_13 = -position1+position3

    perp_vec = np.cross(vec_12, vec_13)

    center = (position1 + position2 + position3)/3

    # This is a very rudementary fix
    if np.linalg.norm(center + perp_vec - positions[26,:]) < np.linalg.norm(center - perp_vec - positions[26,:]):
        perp_vec = - perp_vec

    perp_norm = get_norm_vec(perp_vec)
    vec_1tocenter = center - position1
    length = np.linalg.norm(vec_1tocenter)
    perp_distance = np.sqrt(distance**2 - length**2)
    position = center + (perp_distance)*perp_norm
    add_atom = Atom(symbol, position)
    return add_atom


def add_hollow_atom(template : Atoms, idx1, idx2, idx3, symbol, distance):
    add_atom = get_hollow_atom(template, idx1, idx2, idx3, symbol, distance)
    template.append(add_atom)

def get_hollow_CO(template : Atoms, idx1, idx2, idx3, distance = HOLLOW_CNI_BL, CO_distance = HOLLOW_CO_BL, C_first = True):
    #### Code copyed from get_hollow atom
    positions = template.positions
    position1 = positions[idx1,:]
    position2 = positions[idx2,:]
    position3 = positions[idx3,:]
    
    vec_12 = -position1+position2
    vec_13 = -position1+position3

    perp_vec = np.cross(vec_12, vec_13)
    
    center = (position1 + position2 + position3)/3

    # This is a very rudementary fix
    if np.linalg.norm(center + perp_vec - positions[26,:]) < np.linalg.norm(center - perp_vec - positions[26,:]):
        perp_vec = - perp_vec

    perp_norm = get_norm_vec(perp_vec)
    vec_1tocenter = center - position1
    length = np.linalg.norm(vec_1tocenter)
    perp_distance = np.sqrt(distance**2 - length**2)
    #### end copy
    second_add_position = center + (perp_distance + CO_distance)*perp_norm
    
    if C_first:
        add_C = get_hollow_atom(template, idx1, idx2, idx3, "C", distance)
        add_O = Atom("O", second_add_position)
    else :
        add_O = get_hollow_atom(template, idx1, idx2, idx3, "O", distance)
        add_C = Atom("C", second_add_position)
    return Atoms([add_O, add_C])

def add_hollow_CO(template : Atoms, idx1, idx2, idx3, distance = HOLLOW_CNI_BL, CO_distance = HOLLOW_CO_BL, C_first = True):
    add_atoms = get_hollow_CO(template, idx1 ,idx2, idx3, distance, CO_distance, C_first)
    template.append(add_atoms[0])
    template.append(add_atoms[1])


def classify_Csite(structure):
    "Returns number of C Ni bonds"
    bonds = get_CNibonds(structure)
    return len(bonds)

def classify_Osite(structure):
    "Returns number of C Ni bonds"
    bonds = get_ONibonds(structure)
    return len(bonds)


CO_HOLLOW_TAG = "CO_HOLOOW"
CO_BRIDGE_TAG = "CO_BRIDGE"
CO_TOP_TAG = "CO_TOP"
OC_HOLLOW_TAG = "OC_HOLOOW"
OC_BRIDGE_TAG = "OC_BRIDGE"
OC_TOP_TAG = "OC_TOP"
UNKNOWN_TAG = "UNKNOWN"
def get_tag_CO(structure):
    nCbonds = classify_Csite(structure)
    nObonds = classify_Osite(structure)
    if nCbonds == 3:
        return CO_HOLLOW_TAG
    elif nCbonds == 2:
        return CO_BRIDGE_TAG
    elif nCbonds == 1:
        return CO_TOP_TAG
    elif nObonds == 3:
        return OC_HOLLOW_TAG
    elif nObonds == 2:
        return OC_BRIDGE_TAG
    elif nObonds == 1:
        return OC_TOP_TAG
    else:
        return UNKNOWN_TAG


if __name__ == "__main__":
    #### Show all hollow structures ####
    # hollow_structures = []
    from itertools import combinations
    # for i,j,k in list(combinations(np.arange(26),3)):
    #     try:
    #         hollow_structures.append(get_hollow_CO("POSCAR",i,j,k))
    #     except Exception as e:
    #         pass
    # view(hollow_structures, block=True)


    
    # CO_add = get_top_CO("POSCAR", 0)
    # view(CO_add, block=True)
    template = read("POSCAR") 
    # add_hollow_atom(template, idx1, idx2, idx3, "C",2)
    for idx1,idx2,idx3 in list(combinations(np.arange(26),3)):
        try:
            add_hollow_atom(template, idx1, idx2, idx3, "C",2)
        except Exception as e:
            pass
        
    view(template, block=True)

    # positions = template.positions

    # core_atoms = [26,27,28,29]

    # add_top_atom(template, idx1, "C",4)

    # num_bridge = 0
    # for i in range(26):
    #     for j in range(i):
    #         try:
    #             add_bridge_CO(template, i, j,2, 1.43)
    #             num_bridge += 1
    #         except Exception as e:
    #             pass
    # print(num_bridge)
    # add_bridge_atom(template, idx1, idx2, "C", 2)
    # add_bridge_CO(template, idx1, idx2,2, 1.43)
    
    # add_top_atom(template, 30, "O", 1.43)
    # for i in range(26):
    #     add_top_atom(template, i, "C", 2)
    #     add_top_atom(template, i, "O", 3.43)




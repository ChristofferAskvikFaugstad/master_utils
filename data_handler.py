from utils.vasp_xml import *
import os
import xml.etree.ElementTree as ET
NI30BEST = "C:\\Users\\chris\\masteroppgave\\data\\dft-data\\single\\ni30\\song\\vasprun.xml"
DFTSTART = "C:\\Users\\chris\\masteroppgave\\data\\dft-data"
NEBSTART = "C:\\Users\\chris\\masteroppgave\\data\\dft-data\\neb"

# list of folders to best structures found through search
NI30COBEST = "C:\\Users\\chris\\masteroppgave\\data\\dft-data\\single\\ni30_COs\\CO12-13-22"


class EnergyForceCalculator:
    def __init__(self,energy, forces) -> None:
        self.energy = energy
        self.forces = forces
    def get_forces(self, Atoms):
        return self.forces
    def get_potential_energy(self, atoms, force_consistent = False):
        return self.energy
        
    


def read_outcar(outcar_file_path):
    with open(outcar_file_path, 'r') as f:
        lines = f.readlines()
        structures = []
        for i,line in enumerate(lines):
            if "POSCAR:" in line:
                elements = line.split()[1:]
                
            if "  direct lattice vectors " in line:
                cell_vectors = np.array([lines[i+1+j].split()[0:3] for j in range(3)], dtype=np.float64)

            if "ions per type =" in line:
                number_of_ions = [int(num_str) for num_str in line.split()[4:]]
                tot_number_atoms = sum(number_of_ions)
                elements_array = []
                for element, number in zip(elements, number_of_ions):
                    elements_array += [element]*number

            if line == " POSITION                                       TOTAL-FORCE (eV/Angst)\n":
                positions_forces = np.zeros((tot_number_atoms,6))
                for j in range(tot_number_atoms):
                    positions_forces[j,:] = np.array(lines[i+j+2].split(),dtype=np.float64)
                energy = float(lines[i+j+13].split()[-2])
                structures.append(Atoms(elements_array, positions=positions_forces[:,:3],cell=cell_vectors, pbc = True, calculator=EnergyForceCalculator(energy,positions_forces[:,3:])))

    return structures

def get_Ni30_template():
    return read_vasp_xml_final(NI30BEST)

def get_Ni30_COs_template():
    return read_vasp_xml_final("utilsCO5-25.xml")

def get_structure_path(path):
    return read_vasp_xml_final(path)


def get_structure_relpath(rel_path):
    return read_vasp_xml_final(os.path.join(DFTSTART, rel_path, "vasprun.xml"))

def get_all_structures_path(path):
    return read_vasp_xml_all(path)

def get_all_structures_relpath(relpath):
    return get_all_structures_path(os.path.join(DFTSTART, relpath, "vasprun.xml"))

def get_all_structures_sparse_path(path, sparse = 5):
    structures = []
    mytree = ET.parse(path)
    myroot = mytree.getroot()
    num_calculations = len(myroot.findall("calculation"))
    for i in range(0,num_calculations, sparse):
        structures.append(next(read_vasp_xml(path, index=i)))
    
    return structures

def get_all_structures_sparse_relpath(relpath, sparse = 5):
    return get_all_structures_sparse_path(os.path.join(DFTSTART, relpath, "vasprun.xml"), sparse = sparse)


def get_names_relfolder(folder):
    return [subfolder for subfolder in os.listdir(os.path.join(DFTSTART, folder)) if os.path.isdir(os.path.join(DFTSTART, folder, subfolder))]

def get_structures_folder(folder):
    structures = []
    paths = [os.path.join(folder, subfolder,"vasprun.xml") for subfolder in os.listdir(os.path.join(folder)) if os.path.isdir(os.path.join(folder, subfolder)) ]
    for path in paths:
        structures.append(get_structure_path(path))
    return structures

def get_structures_relfolder(folder):
    return get_structures_folder(os.path.join(DFTSTART, folder))

def get_images_relpath(rel_path):
    images = []
    paths = [os.path.join(NEBSTART, rel_path, subfolder,"vasprun.xml") for subfolder in os.listdir(os.path.join(NEBSTART,rel_path)) if os.path.isdir(os.path.join(NEBSTART,rel_path, subfolder)) ]
    for path in paths:
        images.append(get_structure_path(path))
    return images

def plot_energies_relfolder(folder):
    images = get_structures_relfolder(folder)
    names = get_names_relfolder(folder)
    energies = get_energy_images(images)
    plt.scatter(names, energies)
    plt.xticks(rotation = -90)
    plt.grid()

def get_energy_images(images):
    return np.array([image.get_potential_energy() for image in images])




def make_database_relfolder(folder, db_path,template = None, sparse = 5):
    if template == None:
        template  = get_Ni30_template()
    database = Database(db_path)
    filenames = [os.path.join(folder,name) for name in get_names_relfolder(folder)]
    for i,name in enumerate(filenames):
        for structure in get_all_structures_sparse_relpath(name, sparse= sparse):
            candidate = StandardCandidate(template=template, **structure.todict())
            candidate.calc = structure.calc
            database.store_candidate(candidate)
        print(f"Done {i} of {len(filenames)}")



def get_k_lowest_energies(k, energies, structures):
    """
    Returns the k lowest energies and their associated structures.
    
    Parameters:
        k (int): The number of lowest energies to return.
        energies (list): A list of energies.
        structures (list): A list of structures associated with the energies.
    
    Returns:
        A tuple containing two lists: the k lowest energies and their associated structures.
    """
    # Combine the energies and structures into a list of tuples
    data = list(zip(energies, structures))
    
    # Sort the list of tuples based on energy
    sorted_data = sorted(data, key=lambda x: x[0])
    
    # Extract the k lowest energies and their associated structures
    k_lowest = sorted_data[:k]
    k_energies = [d[0] for d in k_lowest]
    k_structures = [d[1] for d in k_lowest]
    
    # Return the k lowest energies and their associated structures
    return k_energies, k_structures

def extract_stop_reason(oct_file_path):
    with open(oct_file_path, 'r') as f:
        lines = f.readlines()
    for line in reversed(lines):
        if 'aborting loop' in line.lower():
            # return line.split(':')[-1].strip()
            return line
    return None

def diagnose_calculation(outcar_file_path):
    with open(outcar_file_path, 'r') as f:
        lines = f.readlines()
    
    iteration = None
    error = False
    for i,line in enumerate(lines):
        
        if "|     EEEEEEE  R     R  R     R  OOOOOOO  R     R     ###     ###     ###     |\n" == line:
            error_lines = []
            for error_line in lines[i+2:]:
                if error_line == "|       ---->  I REFUSE TO CONTINUE WITH THIS SICK JOB ... BYE!!! <----       |\n":
                    break
                else:
                    error_lines.append(error_line)
            error_msg = "".join(error_lines[:-1])
            error = True
            # return line.split(':')[-1].strip()
            
    for line in reversed(lines):
        if iteration == None:
            if "Iteration" in line:
                iteration = int(line[49:56])
    
    print(f"Used {iteration} Iterations")

    if error:
        print("The followig error occured:")
        print(error_msg)
            
def diagnose_folder(folder,i_max):
    for i in range(i_max):
        print(i)
        diagnose_calculation(f"data/dft-data/{folder}/{i}/OUTCAR")

import numpy as np

def calculate_radius_gyration(atoms  : Atoms):
    """
    Calculate the radius of gyration of an Atoms object.

    Parameters
    ----------
    atoms : ASE Atoms object
        The Atoms object to calculate the radius of gyration for.

    Returns
    -------
    float
        The radius of gyration in Å.
    """
    # Calculate the center of mass
    com = atoms.get_center_of_mass()

    # Translate the coordinates so that the center of mass is at the origin
    positions = atoms.positions - com

    rg = np.sqrt(np.sum(np.sum(positions*positions, axis = 1)*atoms.get_masses())/np.sum(atoms.get_masses()))
    return rg


from utils.manual_placement import get_CNibonds
def name_nCOstructures_in_relfolder(start_folder, nfolder, nCO : int):
    """Renames CO adsorptions structures based on the bonds C make with Ni"""
    structures = get_structures_relfolder(start_folder)
    for i in range(nfolder):
        structure = structures[i]
        CO_bonds = [[] for _ in range(nCO)]
        old_path = os.path.join(DFTSTART,start_folder, str(i))

        for bond in get_CNibonds(structure):
            CO_bonds[bond[0]-30-nCO].append(str(bond[1]))
        
        
        name = "CO" + "_CO".join(["-".join(CO_bond) for CO_bond in CO_bonds])
        new_path = os.path.join(DFTSTART,start_folder, name)
        os.popen(f"move {old_path} {new_path}")

import os
from ase.neb import NEBTools
from utils.make_POV import POVMaker
from utils.data_handler import get_images_relpath

class Reaction:
    reaction_path = "masteroppgave_tex\\Images\\Reactions"
    def __init__(self,relpath, name, title = None):
        if title == None:
             title = name
        self.title = title
        self.images = get_images_relpath(relpath)
        self.nebtool = NEBTools(self.images)
        self.name = name
        self.pov_folder = os.path.join(self.reaction_path,name)
        self.pov_maker = POVMaker(self.images, self.pov_folder)

    def view(self):
        view(self.images)
    
    def make_bandplot(self):
        self.nebtool.plot_band()
        title = plt.gca().title
        title_text = f"{self.title} "+title.get_text()
        title.set_text(title_text)
        plt.savefig(f"{self.pov_folder}/{self.name}.pdf")
        plt.show()
    
    def make_gif(self, rotation = "0x,0y,0z"):
            self.pov_maker.make_pngs( rotation=rotation)
            self.pov_maker.make_gif(name = self.name)

if __name__ == "__main__":
    # structures = read_outcar("OUTCAR")
    # view(structures, block = True)
    
    pass
    # view(get_structures_relfolder("single\\ni30_COHs"), block = True)


    # make_database_relfolder("single/ni30_COs", "CO_10.03_allrelax.db", sparse = 10)

    # pass
    # view(get_Ni30_template(), block = True)
    # view(get_images_relpath("HCO_dis/HCO5-25_CH5-24-25O12-13-25"), block = True)
    # view(get_images_relpath("HCO_dis/HCO5-25_CH5-24-25O12-13-25"), block = True)

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

from utils.vasp_xml import *
import os
import xml.etree.ElementTree as ET
NI30BEST = "C:\\Users\\chris\\masteroppgave\\data\\dft-data\\single\\ni30\\song\\vasprun.xml"
DFTSTART = "C:\\Users\\chris\\masteroppgave\\data\\dft-data"
NEBSTART = "C:\\Users\\chris\\masteroppgave\\data\\dft-data\\neb"

# list of folders to best structures found through search
NI30COBEST = "C:\\Users\\chris\\masteroppgave\\data\\dft-data\\single\\ni30_COs\\CO12-13-22"







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




def make_database_relfolder(folder, db_path,template = get_Ni30_template(), sparse = 5):
    database = Database(db_path)
    filenames = [os.path.join(folder,name) for name in get_names_relfolder(folder)]
    for i,name in enumerate(filenames):
        for structure in get_all_structures_sparse_relpath(name, sparse= sparse):
            candidate = StandardCandidate(template=template, **structure.todict())
            candidate.calc = structure.calc
            database.store_candidate(candidate)
        print(f"Done {i} of {len(filenames)}")



if __name__ == "__main__":
    # view(get_structures_relfolder("single\\ni30_COHs"), block = True)


    make_database_relfolder("single/ni30_COs", "CO_10.03_allrelax.db", sparse = 10)

    pass
    # view(get_Ni30_template(), block = True)
    # view(get_images_relpath("HCO_dis/HCO5-25_CH5-24-25O12-13-25"), block = True)
    # view(get_images_relpath("HCO_dis/HCO5-25_CH5-24-25O12-13-25"), block = True)

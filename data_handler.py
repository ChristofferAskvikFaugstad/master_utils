from utils.vasp_xml import *
from utils.vasp_outcar import *
import os
import shutil
import xml.etree.ElementTree as ET
from ase.neb import NEBTools
from utils.make_POV import POVMaker
import numpy as np


### PATHS ###
DFTSTART = os.path.join(
    "C:\\Users", "chris", "OneDrive - NTNU", "V2022", "data_masteroppgave", "dft-data"
)
# DFTSTART = "C:\\Users\\chris\\masteroppgave\\data\\dft-data"
NEBSTART = os.path.join(DFTSTART, "neb")
# NEBSTART = "C:\\Users\\chris\\masteroppgave\\data\\dft-data\\neb"

# list of folders to best structures found through search
NI30BEST = os.path.join(DFTSTART, "single\\ni30\\song\\vasprun.xml")
NI30COBEST = os.path.join(DFTSTART, "single\\ni30_COs\\CO12-13-22\\vasprun.xml")
NI302COBEST = os.path.join(
    DFTSTART, "single\\ni30_2COs\\CO5-25_CO12-13-22\\vasprun.xml"
)
NI303COBEST = os.path.join(  # 13.04
    DFTSTART, "single\\ni30_3COs\\CO11-14-19_CO12-13-22_CO1-23-24\\vasprun.xml"
)
NI304COBEST = os.path.join(
    DFTSTART, "single\\ni30_4COs\\CO11-14_CO1-23-24_CO6-13-25_CO3-17\\vasprun.xml"
)
NI305COBEST = os.path.join(
    DFTSTART,
    "single\\ni30_5COs\\CO9-16-22_CO4-14-20_CO12-13-22_CO4-15-21_CO11-14\\vasprun.xml",
)
NI306COBEST = os.path.join(  # 11.04
    DFTSTART,
    "single\\ni30_6COs\\CO10-15_CO0_CO3-6-10_CO11-14-19_CO12-13-22_CO6-13-25\\vasprun.xml",
)
NI307COBEST = os.path.join(  # 11.04
    DFTSTART,
    "single\\ni30_7COs\\CO1-23-24_CO8-17_CO5-11-23_CO12-13-22_CO11-14-19_CO7_CO10-15-24\\vasprun.xml",
)
NI308COBEST = os.path.join(  # 11.04
    DFTSTART,
    "single\\ni30_8COs\\CO5-11-23_CO0-7_CO5-12-25_CO1-4-15_CO0-8-21_CO8-17-18_CO11-14_CO1-23-24\\vasprun.xml",
)
NI309COBEST = os.path.join(  # 11.04
    DFTSTART,
    "single\\ni30_9COs\\CO5-12-25_CO6-10-24_CO9-16-22_CO1-4-14_CO14-19-20_CO3_CO2-12-16_CO1-23-24_CO10-15-18\\vasprun.xml",
)
NI3010COBEST = os.path.join(  # 11.04
    DFTSTART,
    "single\\ni30_10COs\\ni30_10COs_proper_order\\vasprun.xml",
)
NI3021COBEST = os.path.join(DFTSTART, "single\\ni30_21COs/21\\vasprun.xml")
NI3026COBEST = os.path.join(DFTSTART, "single\\ni30_26COs\\26\\vasprun.xml")
NI3031COBEST = os.path.join(DFTSTART, "single\\ni30_31COs\\31\\vasprun.xml")
NI3038COBEST = os.path.join(DFTSTART, "single\\ni30_38COs\\38\\vasprun.xml")

BEST_CO = [
    NI30COBEST,
    NI302COBEST,
    NI303COBEST,
    NI304COBEST,
    NI305COBEST,
    NI306COBEST,
    NI307COBEST,
    NI308COBEST,
    NI309COBEST,
    NI3010COBEST,
    NI3021COBEST,
    NI3026COBEST,
    NI3031COBEST,
    NI3038COBEST,
]


TEX = os.path.join("C:\\Users", "chris", "masteroppgave_tex")
TEXIMAGES = os.path.join(TEX, "Images")
TEX_REACTION = os.path.join(TEXIMAGES, "Reactions")

OUTCAR = "OUTCAR"
XML = "vasprun.xml"
NEB = "neb"


def get_teximage_path(path):
    """Returns the path joined with the the tex image folder"""
    return os.path.join(TEXIMAGES, path)


#### FUNCTIONS ####


def reader(outcar_reader: bool, final: bool):
    if outcar_reader:
        if final:
            return read_outcar_final
        else:
            return read_outcar
    else:
        if final:
            return read_vasp_xml_final
        else:
            return read_vasp_xml_all


def add_relative_path(path: str, relative: bool, relative_path=DFTSTART):
    if relative:
        return os.path.join(relative_path, path)
    else:
        return path


def add_path_ending(path, outcar_reader: bool):
    if outcar_reader:
        return os.path.join(path, "OUTCAR")
    else:
        return os.path.join(path, "vasprun.xml")


def get_Ni30_template(outcar_reader=False):
    return reader(outcar_reader, final=True)(NI30BEST)


def get_Ni30_COs_template():
    return read_vasp_xml_final("utilsCO5-25.xml")


def get_structure_path(path, outcar_reader=True, relative=True) -> Atoms:
    """Returns the final structure of the path

    Parameters
    ----------
    path : str
        The path to the folder with the calculation
    outcar_reader : bool, optional
        Determines what file should be read, by default True

    Returns
    -------
    Atoms
        Final geometry of calculation
    """
    path = add_relative_path(add_path_ending(path, outcar_reader), relative)
    return reader(outcar_reader, final=True)(path)


# def get_structure_relpath(rel_path, **kwargs):
#     """Wrapper for the get_structure_path function.
#     It ads the DFTSTART path to the rel_path supplied."""
#     return get_structure_path(os.path.join(DFTSTART, rel_path), **kwargs)


def get_all_structures_path(path, outcar_reader=True, relative=True, sparse=1):
    """Returns all the structures of the path

    Parameters
    ----------
    path : str
        The path to the folder with the calculation
    outcar_reader : bool, optional
        Determines what file should be read, by default True

    Returns
    -------
    Atoms
        Final geometry of calculation
    """
    path = add_relative_path(add_path_ending(path, outcar_reader), relative)
    return reader(outcar_reader, final=False)(path)[::sparse]


# def get_all_structures_relpath(relpath, **kwargs):
#     """Wrapper for the get_all_structures_path function.
#     It ads the DFTSTART path to the rel_path supplied."""
#     return get_all_structures_path(os.path.join(DFTSTART, rel_path), **kwargs)


# def get_all_structures_sparse_path(path, sparse=5):
#     structures = []
#     mytree = ET.parse(path)
#     myroot = mytree.getroot()
#     num_calculations = len(myroot.findall("calculation"))
#     for i in range(0, num_calculations, sparse):
#         structures.append(next(read_vasp_xml(path, index=i)))

#     return structures


# def get_all_structures_sparse_relpath(relpath, sparse=5):
#     return get_all_structures_sparse_path(
#         os.path.join(DFTSTART, relpath, "vasprun.xml"), sparse=sparse
#     )


def get_vasp_calculation_names_path(path, relative=True, outcar_reader=True):
    """Returns the subfolders in the path that contains vaspcalculations"""
    file_ending = OUTCAR if outcar_reader else XML
    path = add_relative_path(path, relative)
    return [
        subfolder
        for subfolder in os.listdir(os.path.join(path))
        if os.path.isfile(os.path.join(path, subfolder, file_ending))
    ]


# def get_vasp_calculations_names_relpath(relpath):
#     """Wrapper for the get_all_structures_path function.
#     It ads the DFTSTART path to the rel_path supplied."""
#     return get_all_structures_path(os.path.join(DFTSTART, rel_path), **kwargs)


def get_names_folder(folder):
    return [
        subfolder
        for subfolder in os.listdir(os.path.join(folder))
        if os.path.isdir(os.path.join(folder, subfolder))
    ]


def get_names_relfolder(folder):
    return [
        subfolder
        for subfolder in os.listdir(os.path.join(DFTSTART, folder))
        if os.path.isdir(os.path.join(DFTSTART, folder, subfolder))
    ]


def get_structures_folder(folder, outcar_reader=True, relative=True):
    structures = []
    paths = [
        os.path.join(folder, subfolder)
        for subfolder in get_vasp_calculation_names_path(
            folder, relative=relative, outcar_reader=outcar_reader
        )
    ]
    for path in paths:
        structures.append(
            get_structure_path(path, outcar_reader=outcar_reader, relative=relative)
        )
    return structures


def get_structures_relfolder(folder):
    return get_structures_folder(os.path.join(DFTSTART, folder))


def get_images_relpath(rel_path):
    images = []
    paths = [
        os.path.join(NEBSTART, rel_path, subfolder, "vasprun.xml")
        for subfolder in os.listdir(os.path.join(NEBSTART, rel_path))
        if os.path.isdir(os.path.join(NEBSTART, rel_path, subfolder))
    ]
    for path in paths:
        images.append(get_structure_path(path))
    return images


def plot_energies_relfolder(folder):
    images = get_structures_relfolder(folder)
    names = get_names_relfolder(folder)
    energies = get_energy_images(images)
    plt.scatter(names, energies)
    plt.xticks(rotation=-90)
    plt.grid()


def get_energy_images(images):
    """Returns a np array of energies for a list of images"""
    return np.array([image.get_potential_energy() for image in images])


def make_database_folder(
    folder, db_path, template=None, sparse=5, outcar_reader=True, relative=True
):
    if template == None:
        template = get_Ni30_template()
    database = Database(db_path)
    filenames = [os.path.join(folder, name) for name in get_names_relfolder(folder)]
    for i, name in enumerate(filenames):
        for structure in get_all_structures_path(
            name, sparse=sparse, outcar_reader=outcar_reader, relative=relative
        ):
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
    with open(oct_file_path, "r") as f:
        lines = f.readlines()
    for line in reversed(lines):
        if "aborting loop" in line.lower():
            # return line.split(':')[-1].strip()
            return line
    return None


def diagnose_calculation(outcar_file_path):
    with open(outcar_file_path, "r") as f:
        lines = f.readlines()

    iteration = None
    error = False
    for i, line in enumerate(lines):

        if (
            "|     EEEEEEE  R     R  R     R  OOOOOOO  R     R     ###     ###     ###     |\n"
            == line
        ):
            error_lines = []
            for error_line in lines[i + 2 :]:
                if (
                    error_line
                    == "|       ---->  I REFUSE TO CONTINUE WITH THIS SICK JOB ... BYE!!! <----       |\n"
                ):
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
    if "Voluntary context switches:" not in lines[-1]:
        print("The calculation did not finish")


def diagnose_folder(folder, relative=True, i_max=None):

    if i_max == None:
        paths = [
            add_relative_path(os.path.join(folder, subfolder), relative)
            for subfolder in get_vasp_calculation_names_path(folder, relative=relative)
        ]
    else:
        path = [
            add_relative_path(os.path.join(folder, subfolder), relative)
            for subfolder in range(i_max)
        ]

    for path in paths:
        print(f"Diagnosing {os.path.basename(path)}:")
        diagnose_calculation(os.path.join(path, OUTCAR))


def calculate_radius_gyration(atoms: Atoms):
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

    rg = np.sqrt(
        np.sum(np.sum(positions * positions, axis=1) * atoms.get_masses())
        / np.sum(atoms.get_masses())
    )
    return rg


from utils.manual_placement import get_CNibonds


def name_nCOstructures_in_folder(
    start_folder, nCO: int, nfolder=None, relative=True, outcar_reader=True
):
    """Renames CO adsorptions structures based on the bonds C make with Ni

    Parameters:
        start_folder (str): The folder to start in.
        nCO (int): The number of CO molecules.
        nfolder (int): The number of folders to rename following the 0, 1, 2, ... .
            If None, all folders are renamed and is unaffected by what name.

    Returns:
        None
    """
    if nfolder == None:
        subfolders = get_vasp_calculation_names_path(
            start_folder, relative=relative, outcar_reader=outcar_reader
        )
    else:
        subfolders = [str(i) for i in range(nfolder)]

    for subfolder in subfolders:
        old_path = os.path.join(start_folder, subfolder)
        structure = get_structure_path(
            os.path.join(old_path), outcar_reader=outcar_reader, relative=relative
        )
        CO_bonds = [[] for _ in range(nCO)]
        for bond in get_CNibonds(structure):
            CO_bonds[bond[0] - 30 - nCO].append(str(bond[1]))

        name = "CO" + "_CO".join(["-".join(CO_bond) for CO_bond in CO_bonds])
        new_path = os.path.join(start_folder, name)
        print(f"Renaming {subfolder} to {name}")
        shutil.move(old_path, new_path)


class Reaction:
    """A class for a reaction."""

    def __init__(self, relpath, name, title=None, relative=True):
        """Initializes a Reaction object.

        Parameters:
            relpath (str): The relative path to the folder containing the images.
            name (str): The name of the reaction.
            title (str): The title of the reaction. If None, the name is used.
        """
        if title == None:
            title = name
        self.title = title
        if relative:
            relpath = os.path.join(NEB, relpath)
        self.images = get_structures_folder(
            relpath, outcar_reader=False, relative=relative
        )
        self.nebtool = NEBTools(self.images)
        self.name = name
        self.pov_folder = os.path.join(TEX_REACTION, name)
        self.pov_maker = POVMaker(self.images, self.pov_folder)

    def view(self):
        """Views the images of the reaction."""
        view(self.images)

    def make_bandplot(self):
        """Makes a bandplot of the reaction."""
        self.nebtool.plot_band()
        title = plt.gca().title
        title_text = f"{self.title} " + title.get_text()
        title.set_text(title_text)
        save_path = os.path.join(self.pov_folder, f"{self.name}.pdf")
        plt.savefig(save_path)
        plt.show()

    def make_gif(self, rotation="0x,0y,0z"):
        """Makes a gif of the reaction.

        Parameters:
            rotation (str): The rotation of the gif.
        """
        self.pov_maker.make_pngs(rotation=rotation)
        self.pov_maker.make_gif(name=self.name)


def get_lowest_energy_subfolder_folder(folder, relative=True, outcar_reader=True):
    """Returns the name of the subfolder with the lowest energy in the folder folder."""
    subfolders = [
        os.path.join(folder, subfolder)
        for subfolder in os.listdir(folder)
        if os.path.isfile(os.path.join(folder, subfolder, "vasprun.xml"))
    ]
    energies = []
    for subfolder in subfolders:
        structure = get_structure_path(
            subfolder, relative=relative, outcar_reader=outcar_reader
        )
        energies.append(structure.get_potential_energy())
    lowest_energy = np.argmin(energies)
    return os.path.basename(subfolders[lowest_energy])


def get_lowest_energy_subfolder_relfolder(rel_folder):
    """Returns the name of the subfolder with the lowest energy in the folder rel_folder."""
    return get_lowest_energy_subfolder_folder(os.path.join(DFTSTART, rel_folder))


def copy_to_data(input_folder: str, output_folder: str):
    """Copies the folder input_folder to the folder output_folder in the DFTSTART folder

    Parameters
    ----------
    input_folder : str
        The folder to copy from
    output_folder : str
        The folder to copy to

    Returns
    -------
    None"""
    input = os.path.join("in", input_folder)
    dir = os.path.join(DFTSTART, output_folder)
    shutil.copytree(input, dir)


if __name__ == "__main__":
    from time import time

    cur_time = time()
    # structures = get_all_structures_path("in\\C2H2_newscr\\0", sparse=5)
    # structures = get_structures_folder("in\\C2H2_newscr", relative=False)
    diagnose_folder("single/ni30_3COs")
    print(f"Elapsed time : {time()-cur_time}")
    view(structures, block=True)

    pass
    # view(get_structures_relfolder("single\\ni30_COHs"), block = True)

    # make_database_relfolder("single/ni30_COs", "CO_10.03_allrelax.db", sparse = 10)

    # pass
    # view(get_Ni30_template(), block = True)
    # view(get_images_relpath("HCO_dis/HCO5-25_CH5-24-25O12-13-25"), block = True)
    # view(get_images_relpath("HCO_dis/HCO5-25_CH5-24-25O12-13-25"), block = True)

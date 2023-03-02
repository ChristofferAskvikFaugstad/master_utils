#### File with all own defines VASP utility functions for direct file handling ####
from my_modules.imports import * # importing all modlues for consistency


def get_pairs(OUTCAR_path : str):# -> list[tuple]:
    """Returns a list with tupels of all pairs as found in the OUTCAR file

    Parameters
    ----------
    OUTCAR_path : str
        The path to the OUTCAR file

    Returns
    -------
    list[tuple]
        list with all pairs of nearest neigbhours
    """
    pairs = []
    offsett = -1
    with open(OUTCAR_path, mode="r") as file:
        for line in file:
            if "nearest neighbor table" in line:
                for table_line in file:
                    if table_line != " \n":
                        split_line = [num for num in table_line.split(" ") if num != ""]
                        primary_ion = int(split_line[0])
                        secondary_ions = [int(ion_num) for ion_num in split_line[4::2]]
                        for secondary_ion in secondary_ions:
                            if secondary_ion > primary_ion:
                                pair = (primary_ion + offsett, secondary_ion + offsett)
                            else:
                                pair = (secondary_ion + offsett, primary_ion + offsett)
                            if pair not in pairs:
                                pairs.append(pair)
                    else:
                        break
        return pairs
def cohesiv_energy(final_energy, N):
    calc = Calculation.from_path("C:\\Users\\chris\\OneDrive - NTNU\\H2022\\prosjektoppgave\\ni1")
    free_ni_energy = calc.energy.to_numpy()
    print(f"Cohesive energy {(final_energy - N * free_ni_energy) / N}")


def table_coordinates(xyz_file):
    with open(xyz_file) as f:
        f.readline()
        f.readline()
        for line in f:
            print("&".join(line.split()[:4]) + "\\\\")

def everything(name, energy0, N):
    path = f"final_structures\\{name}"
    xyz_path = os.path.join(path, f"{name}.xyz")
    structure = read(xyz_path)
    calc = Calculation.from_path(path)
    energy = calc.energy.to_numpy()
    print(f"Energy above lowest {energy-energy0}")
    cohesiv_energy(energy, N)
    bond_distance(structure)
    # table_coordinates(xyz_path)
    fig = calc.dos.plot()
    fig.write_image(f"dos_{name}.pdf")
    return fig



def filter_structures(structures : List[Atoms], comparator : OFPComparator) -> List[Atoms]:
    N = len(structures)
    assert(N > 1)
    out_structures = [structures[0]]
    # index_unique = []
    for i,structure in enumerate(structures):
        any_alike = False
        for j,comp_structure in enumerate(out_structures):
            if comparator.looks_like(structure, comp_structure):
                any_alike = True
                break

        if not any_alike:
            out_structures.append(structure)
    return out_structures


def convert_pngs_to_gif(name : str  = "gif",
        save_all : bool = True,
        duration : int = 300,
        loop : int = 0,
         **kwargs) -> None:
    """In the cwd it gatters all png files and sorts on the float of the
    base name of the pngs to order

    Parameters
    ----------
    name : str, optional
        _description_, by default "gif"
    save_all : bool, optional
        _description_, by default True
    duration : int, optional
        _description_, by default 300
    loop : int, optional
        _description_, by default 0
    kwargs : 
        Sent to the PIL.Image.save function
    """
    
    from PIL import Image
    import glob

    frames =  []
    imgs = glob.glob("*.png")

    imgs = sorted(imgs, key = lambda x: float(x[:-4]) )

    for i in imgs:
        new_frame = Image.open(i)
        frames.append(new_frame)

    frames[0].save(f"{name}.gif", format='GIF',
                    append_images = frames[1:],
                    save_all = save_all,
                    duration = duration, loop = loop, **kwargs)


def bond_distance(structure, cutoff=1.1):
    N = len(structure.positions)
    positions = structure.positions
    distance = lambda x: np.linalg.norm(x, ord=2)
    n_list = build_neighbor_list(structure, cutoffs=[cutoff] * N)
    distances = [
        distance(positions[i[0]] - positions[i[1]])
        for i in n_list.get_connectivity_matrix().keys()
        if i[0] != i[1]
    ]
    avg_distance = np.average(distances)
    std_distance = np.std(distances)
    n_bonds = len(n_list.get_connectivity_matrix().keys()) - N
    print(
        f"The bond distance is {avg_distance} ({std_distance}) Å with {n_bonds} Bonds"
    )









def get_dos_xml(filename : str, spin : bool = False):# -> tuple[np.array]:
    """Finds the values for the dos from xml filename provided.


    Parameters
    ----------
    filename : str
        Name or path to xml file
    
    spin : bool
        Determines if both spin are present 

    Returns
    -------
    tuple[np.array]
        if spin False:
            energies, density of state, integrated density of state
        if spin True:
            energies, density of state up, down, integrated density of state up , down
    """

    from xml.etree.ElementTree import parse
    from numpy import zeros
    
    mytree = parse(filename)
    myroot = mytree.getroot()


    if spin:
        # Find with both spins
        calculation = myroot.find("calculation")
        dos = calculation.find("dos")
        dos_total = dos.find("total")
        dos_total_array = dos_total.find("array")
        dos_total_array_set = dos_total_array.find("set")
        dos_total_array_set_set_spin1 = dos_total_array_set[0]
        dos_total_array_set_set_spin2 = dos_total_array_set[1]


        N_points = len(dos_total_array_set_set_spin1)
        energies = zeros(N_points)

        dos_spin1 = zeros(N_points)
        dos_int_spin1 = zeros(N_points)
        for i,line in enumerate(dos_total_array_set_set_spin1):
            energies[i], dos_spin1[i], dos_int_spin1[i] = [float(num) for num in line.text.split(" ") if num != ""]

        dos_spin2 = zeros(N_points)
        dos_int_spin2 = zeros(N_points)
        for i,line in enumerate(dos_total_array_set_set_spin2):
            _, dos_spin2[i], dos_int_spin2[i] = [float(num) for num in line.text.split(" ") if num != ""]
        
        return energies, dos_spin1, dos_spin2, dos_int_spin1, dos_int_spin2
        
    else:
        # Just single spin
        calculation = myroot.find("calculation")
        dos = calculation.find("dos")
        dos_total = dos.find("total")
        dos_total_array = dos_total.find("array")
        dos_total_array_set = dos_total_array.find("set")
        dos_total_array_set_set = dos_total_array_set.find("set")
        dos_total_array_set_set


        N_points = len(dos_total_array_set_set)
        energies = zeros(N_points)
        dos = zeros(N_points)
        dos_int = zeros(N_points)
        for i,x in enumerate(dos_total_array_set_set):
            energies[i], dos[i], dos_int[i] = [float(num) for num in x.text.split(" ") if num != ""]
        
        return energies, dos, dos_int


def get_final_energies_xml(xml_filepath : str):
    """get energies from xml file

    Parameters
    ----------
    xml_filepath : str
        relative or complete paht

    Returns
    -------
    float, float, float
        free energy, energy without entropy and 0 energy
    """
    import xml.etree.ElementTree as ET
    # import matplotlib.pyplot as plt
    mytree = ET.parse(xml_filepath)
    myroot = mytree.getroot()
    dos_total = myroot.findall("calculation")[-1] #.find("structure")
    dos_total = dos_total.find("energy")#.find("crystal")#.find("energy") #.find("structure")
    free_energy = float(dos_total[0].text)
    energy_wo_entropy = float(dos_total[1].text)
    energy_0  = float(dos_total[2].text)
    
    return free_energy, energy_wo_entropy, energy_0

def get_all_energies_xml(xml_filepath : str):
    """get energies from xml file for all steps

    Parameters
    ----------
    xml_filepath : str
        relative or complete path

    Returns
    -------
    float, float, float
        free energies, energies without entropy and 0 energies
    """
    import xml.etree.ElementTree as ET
    from numpy import zeros
    # import matplotlib.pyplot as plt
    mytree = ET.parse(xml_filepath)
    myroot = mytree.getroot()
    calculation = myroot.findall("calculation")


    N_step = len(calculation)
    free_energies = zeros(N_step)
    energies_wo_entropy = zeros(N_step)
    energies_0 = zeros(N_step)
    for i, step in enumerate(calculation):
        energy = step.find("energy")
        free_energies[i] = float(energy[0].text)
        energies_wo_entropy[i] = float(energy[1].text)
        energies_0[i]  = float(energy[2].text)
    
    return free_energies, energies_wo_entropy, energies_0

def get_final_structures_xml(xml_filepath : str, elements:str = "Ni8"):
    """get structures from xml file for all steps

    Parameters
    ----------
    xml_filepath : str
        relative or complete path
    
    elements : str
        string with the elements used: Default Value ="Ni8"

    Returns
    -------
    float, float, float
        free energies, energies without entropy and 0 energies
    """
    import xml.etree.ElementTree as ET
    from numpy import asanyarray
    from ase.atoms import Atoms, Atom


    mytree = ET.parse(xml_filepath)
    myroot = mytree.getroot()
    final_step = myroot.findall("calculation")[-1]
    info = final_step.find("structure")
    basis = info.find("crystal")[0]
    positions_dat = info.find("varray")
    basis_vectors = asanyarray([asanyarray([num for num in x.text.split(" ") if num != ""], dtype = float) for x in basis])
    
    positions = []
    for x in positions_dat:
        position = asanyarray([num for num in x.text.split(" ") if num != ""], dtype = float)
        positions.append(position@basis_vectors)

    structure = Atoms(elements,positions = positions, cell = basis_vectors )

    return structure



def get_all_structures_xml(xml_filepath : str, elements:str = "Ni8"):
    """get structures from xml file for all steps

    Parameters
    ----------
    xml_filepath : str
        relative or complete path
    
    elements : str
        string with the elements used: Default Value ="Ni8"

    Returns
    -------
    float, float, float
        free energies, energies without entropy and 0 energies
    """
    import xml.etree.ElementTree as ET
    from numpy import zeros, asanyarray
    from ase.atoms import Atoms, Atom


    
    mytree = ET.parse(xml_filepath)
    myroot = mytree.getroot()
    calculation = myroot.findall("calculation")

    structures = []

    for i, step in enumerate(calculation):
        info = step.find("structure")
        basis = info.find("crystal")[0]
        positions_dat = info.find("varray")
        basis_vectors = asanyarray([asanyarray([num for num in x.text.split(" ") if num != ""], dtype = float) for x in basis])
        # print(basis_vectors)
        
        positions = []
        for x in positions_dat:
            position = asanyarray([num for num in x.text.split(" ") if num != ""], dtype = float)
            positions.append(position@basis_vectors)

        structure = Atoms(elements,positions = positions, cell = basis_vectors )

        structures.append(structure)

    return structures


def fig_w_innerplot(left: float = 0.5,
        bottom: float = 0.5,
        width :float =  0.39,
        height: float = 0.35,
        fig_width : float = plt.rcParams["figure.figsize"][0],
        fig_height : float = plt.rcParams["figure.figsize"][1]):# -> tuple[plt.Figure, plt.Axes, plt.Axes]:

    limits = [left, bottom, width, height]
    figsize = [fig_width, fig_height]
    
    fig, ax1 = plt.subplots(figsize = figsize)
    ax2 = fig.add_axes(limits)

    return fig, ax1, ax2

def fill_search_axis(
    ax1, ax2, range, free_energy,sub_range, sub_free_energy, choosen_value, xlabel, ylabel, title=""
):
    """ Fills in for the free energy search stuff, with premade template"""
    ax1.set_title(title)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    choosen_index = np.where(range == choosen_value)[0]
    ax1.scatter(
        np.delete(range, choosen_index),
        np.delete(free_energy, choosen_index),
        color=colors[0],
        marker="x",
        s=100,
        linewidth=2,
    )
    ax1.scatter(
        range[choosen_index],
        free_energy[choosen_index],
        color=colors[1],
        marker="x",
        s=100,
        linewidth=2,
    )
    ax1.grid()

    # ax2.set_title("Magnifed")
    sub_choosen_index = np.where(sub_range == choosen_value)[0]

    ax2.scatter(
        np.delete(sub_range, sub_choosen_index),
        np.delete(sub_free_energy, sub_choosen_index),
        color=colors[0],
        marker="x",
        s=100,
        linewidth=2,
    )
    ax2.scatter(
        sub_range[sub_choosen_index],
        sub_free_energy[sub_choosen_index],
        color=colors[1],
        marker="x",
        s=100,
        linewidth=2,
    )
    ax2.grid()


def plot_step_energies_xml(filepath : str,
     title: str = "",
     xlab : str = "Iteration",
     ylab : str = "Free Energy [eV]",
     zoom_index : int = 9,
     display_best : bool = True,
     save_path : str = None,
     **kwargs):
    """Plots the steps of vasp run from xml

    Parameters
    ----------
    filepath : str
        Path to xml file
    title : str, optional
        Title on main figure, by default ""
    xlab : str, optional
        label for x axis, by default "Iteration"
    ylab : str, optional
        label for y-axis, by default "Free Energy [eV]"
    display_best : bool, optional
        writes the last energy if True, by default True
    save_path : str, optional
        if provided writes the output to that file, if not 
        provided nothing is saved, by default None
    kwargs : 
        The keyword sent to fig_w_innerplot (own defined)
    """
    free_energies, _, _ = get_all_energies_xml(filepath)

    iteration = range(len(free_energies))

    last_energy = free_energies[-1]

    fig, ax1, ax2 = fig_w_innerplot(**kwargs)
    ax1.set_title(title)
    ax1.plot(iteration, free_energies, "rx")
    ax2.plot(iteration[zoom_index:], free_energies[zoom_index:], "rx")
    # ax1.grid()
    ax1.set_xlabel(xlab)
    ax1.set_ylabel(ylab)

    if display_best:
        xlab = f"{xlab} Last : {last_energy}"
    ax2.set_xlabel(xlab)

    ax2.set_ylabel(ylab)
    ax2.grid()
    if save_path != None:
        plt.savefig(save_path)
    
    plt.show()


def view_last_structure_xml(xml_filepath : str, in_notebook : bool = True, **kwargs):
    """View the last structure in the xml file

    Parameters
    ----------
    xml_filepath : str
        The source file
    in_notebook : bool, optional
        To display inline if in notebook(True) or in pop up window(False), by default True
    kwargs :
        Passed onto get_final_structures_xml
    """
    from ase.visualize import view
    structure = get_final_structures_xml(xml_filepath, **kwargs)

    if in_notebook:
        viewer = "x3d"
    else:
        viewer = "ase"

    return view(structure, viewer= viewer)
    
def view_all_structure_xml(xml_filepath : str, **kwargs):
    """View the all structure in the xml file

    Parameters
    ----------
    xml_filepath : str
        The source file
    kwargs :
        Passed onto get_all_structures_xml
    
    """
    from ase.visualize import view
    structures = get_all_structures_xml(xml_filepath, **kwargs)
        
    return view(structures)

def plot_gofee_evaluations(name, num, final_energies,fig_path):
    """plots the gofee avaluations from one run

    Parameters
    ----------
    name : str
        The link to the database
    num : int
        where to start the ploting from, remove trainingdata
    final_energies : npa
        the fully relax energies that where used as training databases
    fig_path : function 
        joins the input name with the appropriate location to store figures
    """
    properties = read_database(name)
    energies = properties["energies"][num:]
    forcesmax = properties["forcesmax"][num:]

    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    evaluations = np.arange(0, len(energies))
    scatter = ax.scatter(
        evaluations, energies, label="GOFEE evaluated candidates", c=forcesmax
    )
    ax.hlines(
        final_energies,
        0,
        200,
        label="Fully Relaxed Sturctures",
        linestyles="--",
        color=colors[1],
    )
    ax.set_ylabel("Energy [eV]")
    ax.set_xlabel("Singel Point DFT Calculations")
    ax.legend()
    cbar = plt.colorbar(scatter)
    cbar.set_label("Max Force on Atom [eV/Å] ", rotation=270, labelpad=25)
    plt.tight_layout()
    plt.savefig(fig_path("gofee_full_relax.pdf"))
    plt.grid()
    ax.set_title(
        "Ni 30 GOFEE 100 iterations 1 step with fully relaxed trajectories as training data"
    )
    plt.tight_layout()
    plt.savefig(fig_path("gofee_full_relax_t.pdf"))

    plt.show()



def read_database(file_path: str) -> dict:
    """Reads a AGOX database and return a dictionary
    with the most relevant features

    Parameters
    ----------
    file_path : str
        The filepath relative of absolute to the database

    Returns
    -------
    dict
        keys:
        -  size : The nummber of studd in it
        - energies : The energies of all candidates
        - forces : all the forces on the atoms
        - forcesmax : the largest force on one atom in size
    """
    properites_dict = {}
    
    database = Database(file_path)
    energies = database.get_all_energies()
    size = len(energies)
    properites_dict["size"] = size
    properites_dict["energies"] = energies
    all_structures = []
    forces = []
    forcesmax = np.zeros(size)
    all_structure_data = database.get_all_structures_data()
    for i in range(size):
        current_data = all_structure_data[i]
        all_structures.append(database.db_to_atoms(current_data))
        forces.append(current_data["forces"])
        forcesmax[i] = np.max(np.linalg.norm(all_structures[i].get_forces(), axis = 1))

    properites_dict["forces"] = forces
    properites_dict["forcesmax"] = forcesmax
    properites_dict["structures"] = all_structures
    return properites_dict


def parse_outfile_for_predictions(outfile: str) -> np.ndarray:
    """parse the outfile AGOX give to obtain the predicted energy of
    the best candidate of the aquistor. Does not incorperate if it is
    not convered in the evaluation. 

    Parameters
    ----------
    outfile : str
        The path to the file

    Returns
    -------
    np.ndarray
        Array with the predicted energies
    """
    predicted_energies = []
    with open(outfile) as f:
        for line in f:
            if line == "|=============================== LCBAcquisitor ===============================|\n":
                next_line = f.readline()
                predicted_energies.append(float(next_line.split(": E=")[1][:8]))
    return np.array(predicted_energies)



#########################
# Model Accuracy PLoting#
#########################
def plot_model_accuracy(
    predicted_path: str,
    db_path: str,
    num : int,
    fig_path,
    emt : bool = False,
    axis: plt.Axes = None,
    xlabel: str = "Predicted Energy [eV]",
    ylabel: str = "DFT Energy [eV]",
    cmap_label: str = "DFT Energy - Predicted Energy [eV]",
    placement: np.ndarray = None,
    title : str = "",
    save_name : str = "model_dft_comparison",
    **kwargs
):
    """Plots the dispersion between predicted and actutal energies
    for the a model

    Parameters
    ----------
    predicted_path : str
        path to outfile
    actual_energies : np.ndarray
        The actual values
    axis : plt.Axes, optional
        The axis to plot on if Not set it is taken from plt.gca(), by default None
    xlabel : str, optional
        The xlabel set, by default "Predicted Energy [eV]"
    ylabel : str, optional
        The ylabel set, by default "DFT Energy [eV]"
    placement : np.ndarray, optional
        The placment of error annotation, if not given anythin nothin is ploted , by default None
    kwargs:
        Passed onto annotate_accuracy_plot
    """
    from matplotlib.colors import TwoSlopeNorm
    fig = plt.figure(figsize=(8,8))
    if axis == None:
        axis = plt.gca()
    properties = read_database(db_path)
    actual_energies = properties["energies"][num::2]
    if not emt:
        predicted_energies = parse_outfile_for_predictions(predicted_path)
    else:
        calc = EMT()
        predicted_energies = []
        for structure in properties["structures"][num::2]:
            structure.calc = calc
            predicted_energies.append(structure.get_potential_energy())
        predicted_energies = np.array(predicted_energies)
        predicted_energies += np.min(actual_energies)-np.min(predicted_energies)
        xlabel = "EMT Energy [eV]"
        cmap_label = "DFT Energy - EMT Energy [eV]"
        save_name = "emt_dft_comparison"
    emin = np.min([actual_energies, predicted_energies])
    emax = np.max([actual_energies, predicted_energies])
    diff = actual_energies - predicted_energies
    scatter = axis.scatter(
        predicted_energies,
        actual_energies,
        c=diff,
        cmap="coolwarm",
        norm=TwoSlopeNorm(0),
    )
    axis.set_xlabel(xlabel)
    axis.set_ylabel(ylabel)

    axis.plot([emin, emax], [emin, emax], linestyle="--", color="black")
    axis.set_aspect("equal")
    axis.grid()

    cbar = plt.colorbar(scatter)
    cbar.set_label(cmap_label, rotation=270, labelpad=25)
    if placement != None:
        annotate_accuracy_plot(predicted_energies, actual_energies, placement, **kwargs)
    fig.tight_layout()
    plt.savefig(fig_path(f"{save_name}.pdf"))
    plt.grid()
    axis.set_title(
        title
    )
    fig.tight_layout()
    plt.savefig(fig_path(f"{save_name}_t.pdf"))
    plt.show()


    
def rms_error(predicted_energies, actual_energies):
    """Computes Root Mean Square error between inputs
    """
    return np.sqrt(np.mean((predicted_energies-actual_energies)*(predicted_energies-actual_energies)))

def am_error(predicted_energies, actual_energies):
    """Computes Mean Absolute error between inputs
    """
    return np.mean(np.abs(predicted_energies-actual_energies)) 

def annotate_accuracy_plot(
        predicted_energies : np.ndarray,
        actual_energies : np.ndarray,
        placement : np.ndarray,
        axis : plt.Axes = None,
        error_func = rms_error,
        label = "RMS Error :",
        unit = "eV",
        precision = 4,
        ):
    """Helper function to add annotations to model_accuracy plot

    Parameters
    ----------
    predicted_energies : np.ndarray
        predictions
    actual_energies : np.ndarray
        actual values
    placement : np.ndarray
        Where the anotations should be placed in the plot
    axis : plt.Axes, optional
        Axis to plot on, uses plt.gca() if None, by default None
    error_func : _type_, optional
        The function to compute error, by default rms_error
    label : str, optional
        The name of error, by default "RMS Error :"
    unit : str, optional
        unit of error, by default "eV"
    precision : int, optional
        Number of descimals, by default 4
    """
    
    if axis == None:
        axis = plt.gca()
    axis.annotate(f"{label} {round(error_func(predicted_energies, actual_energies),precision)} {unit}", placement)





##############################################
##############################################
# Code taken from the internet
# Cannot find site, something with pgroup
##############################################
##############################################



#!/usr/bin/env python
"""Simplified POV-RAY environment creator."""


import os
import numpy as np
from ase.data.colors import jmol_colors
from ase.data import covalent_radii

class POV:
    """Class to creat .pov files to be executed with the pov-ray raytracer.
    Initialize with

    atoms : an ase atoms object

    Keyword arguments: (all distances in Angstroms)
    ------------------
    tex : a texture to use for the atoms, either a single value or a list
        of len(atoms), default = 'vmd'
    radii : atomic radii. if a single value is given, it is interpreted as
        a multiplier for the covalent radii in ase.data. if a list of
        len(atoms) is given, it is interpreted as individual atomic radii.
        default = 1.0
    colors : a list of len(atoms) of the colors, as (r,g,b). default is
        None which will use ASE standard colors
    cameratype : type of povray camera, default='perspective'
    cameralocation : location of camera as an (x,y,z) tuple,
        default = (0., 0., 20)
    look_at : where the camera is pointed at as an (x,y,z) tuple.
        default = (0., 0., 0.)
    camera_right_up : the right and up vectors that define the image
        boundaries. The right:up ratio will also define the aspect ratio
        (width:height) of the resulting image. These two vectors should
        generally be orthogonal -- generaly on the x,y plane.
        default = [(-8.,0.,0.),(0.,6.,0.)]
    cameradirection : the initial direction vector of the camera before
        it is moved with look_at. This will also control the zoom, with
        higher values being more zoomed in. default = (0., 0., 10.)
    area_light : location and parameters of area light as [(x,y,x), color,
        width, height, Nlamps_x, Nlamps_y], default = [(20., 3., 40.),
        'White', .7, .7, 3, 3]
    background : background color, default = 'White'
    bondatoms : list of atoms to be bound together, as in
        [(index1, index2), ...], default = None
    bondradii : radii to use in drawing bonds, default = 0.1
    pixelwidth : width in pixels of the final image. Note that the height
        is set by the aspect ratio (controlled by carmera_right_up).
        default = 320
    clipplane : plane at which to clip atoms, for example "y, 0.00".
        default = None
    """

    _default_settings = {
        'tex': 'vmd',
        'radii': 1.,
        'colors': None,
        'cameratype': 'perspective',
        'cameralocation': (0., 0., 20.),
        'look_at': (0., 0., 0.),
        'camera_right_up': [(-8., 0., 0.), (0., 6., 0.)],
        'cameradirection': (0., 0., 10.),
        'area_light': [(20., 3., 40.), 'White', .7, .7, 3, 3],
        'background': 'White',
        'bondatoms': None,
        'bondradius': .1,
        'pixelwidth': 320,
        'clipplane': None,
    }

    def __init__(self, atoms, **kwargs):
        for k, v in self._default_settings.items():
            setattr(self, '_' + k, kwargs.pop(k, v))
        if len(kwargs) > 0:
            print(kwargs)
            raise TypeError('POV got one or more unexpected keywords.')
        self._atoms = atoms
        self._numbers = atoms.get_atomic_numbers()
        if self._colors is None:
            self._colors = jmol_colors[self._numbers]
        if (type(self._radii) is float) or (type(self._radii) is int):
            self._radii = covalent_radii[self._numbers] * self._radii
        if self._bondatoms is None:
            self._bondatoms = []
        self._aspectratio = (np.linalg.norm(self._camera_right_up[0]) /
                             np.linalg.norm(self._camera_right_up[1]))

    def write(self, filename, run_povray=None):
        """Writes out the .pov file for ray-tracing and also an associated
        .ini file. If filename ends in ".png" it will run povray to turn it
        into a png file. If the filename ends in ".pov" it will not. This can
        be overridden with the keyword run_povray.
        """
        if filename.endswith('.png'):
            filename = filename[:-4] + '.pov'
            if run_povray is None:
                run_povray = True
        elif filename.endswith('.pov'):
            if run_povray is None:
                run_povray = False
        else:
            raise RuntimeError('filename must end in .pov or .png')
        self._filename = filename
        filebase = filename[:-4]
        # Write the .pov file.
        f = open(filebase + '.pov', 'w')

        def w(text):
            f.write(text + '\n')

        w('#include "colors.inc"')
        w('#include "finish.inc"')
        w('')
        w('global_settings {assumed_gamma 1 max_trace_level 6}')
        w('background {color %s}' % self._background)
        w('camera {%s' % 'perspective')
        w('  location <%.2f,%.2f,%.2f>' % tuple(self._cameralocation))
        camera_right_up = self._camera_right_up
        w('  right <%.2f,%.2f,%.2f> up <%.2f,%.2f,%.2f>' %
          (camera_right_up[0][0], camera_right_up[0][1],
           camera_right_up[0][2], camera_right_up[1][0],
           camera_right_up[1][1], camera_right_up[1][2]))
        w('  direction <%.2f,%.2f,%.2f>' % tuple(self._cameradirection))
        w('  look_at <%.2f,%.2f,%.2f>}' % tuple(self._look_at))
        w('light_source {<%.2f,%.2f,%.2f> color %s' %
          tuple(list(self._area_light[0]) + [self._area_light[1]]))
        w('  area_light <%.2f,0,0>, <0,%.2f,0>, %i, %i' %
          tuple(self._area_light[2:]))
        w('  adaptive 1 jitter}')
        w('')
        w('#declare simple = finish {phong 0.7}')
        w('#declare pale = finish {ambient .5 diffuse .85 roughness .001 specular 0.200 }')
        w('#declare intermediate = finish {ambient 0.3 diffuse 0.6 specular 0.60 roughness 0.04 }')
        w('#declare vmd = finish {ambient .0 diffuse .65 phong 0.1 phong_size 40. specular 0.500 }')
        w('#declare jmol = finish {ambient .2 diffuse .6 specular 1 roughness .001 metallic}')
        w('#declare ase2 = finish {ambient 0.05 brilliance 3 diffuse 0.6 metallic specular 0.70 roughness 0.04 reflection 0.15}')
        w('#declare ase3 = finish {ambient .15 brilliance 2 diffuse .6 metallic specular 1. roughness .001 reflection .0}')
        w('#declare glass = finish {ambient .05 diffuse .3 specular 1. roughness .001}')
        w('#declare Rbond = %.3f;' % self._bondradius)
        w('')
        if self._clipplane is None:
            w('#macro atom(LOC, R, COL, FIN)')
            w('  sphere{LOC, R texture{pigment{COL} finish{FIN}}}')
            w('#end')
        else:
            w('#macro atom(LOC, R, COL, FIN)')
            w('  difference{')
            w('   sphere{LOC, R}')
            w('   plane{%s}' % self._clipplane)
            w('   texture{pigment{COL} finish{FIN}}')
            w('  }')
            w('#end')
        w('')

        for atom in self._atoms:
            w('atom(<%.2f,%.2f,%.2f>, %.2f, rgb <%.2f,%.2f,%.2f>, %s) // #%i'
              % (atom.x, atom.y, atom.z,
                 self._radii[atom.index], self._colors[atom.index][0],
                 self._colors[atom.index][1], self._colors[atom.index][2],
                 self._tex, atom.index))

        for bond in self._bondatoms:
            pos0 = self._atoms[bond[0]].position.copy()
            pos1 = self._atoms[bond[1]].position.copy()
            middle = (pos0 + pos1) / 2.
            color0 = self._colors[bond[0]]
            color1 = self._colors[bond[1]]
            w('cylinder {<%.2f,%.2f,%.2f>, <%.2f,%.2f,%.2f>, Rbond '
              'texture{pigment {rgb <%.2f,%.2f,%.2f>} finish{%s}}} '
              ' // # %i to %i' %
              (pos0[0], pos0[1], pos0[2],
               middle[0], middle[1], middle[2],
               color0[0], color0[1], color0[2], self._tex,
               bond[0], bond[1]))
            w('cylinder {<%.2f,%.2f,%.2f>, <%.2f,%.2f,%.2f>, Rbond '
              'texture{pigment {rgb <%.2f,%.2f,%.2f>} finish{%s}}} '
              ' // # %i to %i' %
              (middle[0], middle[1], middle[2],
               pos1[0], pos1[1], pos1[2],
               color1[0], color1[1], color1[2], self._tex,
               bond[0], bond[1]))
        f.close()

        # Write the .ini file.
        f = open(filebase + '.ini', 'w')
        w('Input_File_Name=%s' % os.path.split(filename)[1])
        w('Output_to_File=True')
        w('Output_File_Type=N')
        w('Output_Alpha=False')
        w('Width=%d' % self._pixelwidth)
        w('Height=%.0f' % (self._pixelwidth / self._aspectratio))
        w('Antialias=True')
        w('Antialias_Threshold=0.1')
        w('Display=False')
        w('Pause_When_Done=True')
        w('Verbose=False')

        f.close()
        if run_povray:
            self.raytrace(filename)

    def raytrace(self, filename=None):
        """Run povray on the generated file."""

        if not filename:
            filename = self._filename
        filebase = filename[:-4]
        if os.path.split(filename)[0] != '':
            pwd = os.getcwd()
            os.chdir(os.path.split(filename)[0])
        os.system('povray %s.ini' % os.path.split(filebase)[1])
        if os.path.split(filename)[0] != '':
            os.chdir(pwd)



######################################            
######################################           
# Finished with stolen part 
######################################            
######################################
from ase import Atoms

def make_image(structure : Atoms, 
    dir = None, 
    base_name = "",
    radii = 0.7,
    bounds = None,
    cameralocation=(20., 20., 20.),
    **kwargs):

    center = structure.get_center_of_mass()

    pov = POV(structure,
            pixelwidth = 2*640,
            radii = radii,
            bondatoms = bounds,
            cameralocation = cameralocation,
            look_at = center,**kwargs
            )
    if dir == None:
        file_path = f"{base_name}.pov"
    else:
        file_path = os.path.join(dir,f"{base_name}.pov")
    
    pov.write(file_path, run_povray=True)
    
    return None

def convert_pngs_to_gif(name : str  = "gif",
        save_all : bool = True,
        duration : int = 300,
        loop : int = 0,
         **kwargs) -> None:
    """In the cwd it gatters all png files and sorts on the float of the
    base name of the pngs to order

    Parameters
    ----------
    name : str, optional
        _description_, by default "gif"
    save_all : bool, optional
        _description_, by default True
    duration : int, optional
        _description_, by default 300
    loop : int, optional
        _description_, by default 0
    kwargs : 
        Sent to the PIL.Image.save function
    """

    
    from PIL import Image
    import glob

    frames =  []
    imgs = glob.glob("*.png")

    imgs = sorted(imgs, key = lambda x: float(x[:-4]) )

    for i in imgs:
        new_frame = Image.open(i)
        frames.append(new_frame)

    frames[0].save(f"{name}.gif", format='GIF',
                    append_images = frames[1:],
                    save_all = save_all,
                    duration = duration, loop = loop, **kwargs)

from my_modules.imports import*
from typing import List

def read_vasp_xml_final(xml_path : str) -> Atoms:
    """Reads the last structure from the xml path

    Parameters
    ----------
    xml_path : str
        xml path to use 

    Returns
    -------
    Atoms
        final atoms object
    """
    return next(read_vasp_xml(xml_path, index=-1))

def read_vasp_xml_all(xml_path: str) -> List[Atoms]:
    """Reads all the structures from a vasp calculation from the xml file

    Parameters
    ----------
    xml_path : str
        what path to use

    Returns
    -------
    List[Atoms]
        the list of calculation poitns
    """
    import xml.etree.ElementTree as ET

    mytree = ET.parse(xml_path)
    myroot = mytree.getroot()
    num_calculations = len(myroot.findall("calculation"))

    structures = []    
    for i in range(num_calculations):
        structures.append(next(read_vasp_xml(xml_path, index=i)))

    return structures
    
def read_vasp_xml_energies(xml_filepath : str):
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

def read_vasp_xml_energies_final(xml_filepath: str):
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

def read_vaps_xml_forces(xml_filepath : str):
    import xml.etree.ElementTree as ET
    mytree = ET.parse(xml_filepath)
    myroot = mytree.getroot()
    calculation = myroot.findall("calculation")

    forces = []
    step  = calculation[0]

    for i, step in enumerate(calculation):
        current_forces = []
        for force_xml in step.findall("varray")[0]:
            current_forces.append(np.array(force_xml.text.split(), dtype = float))
        current_forces = np.array(current_forces)
        forces.append(current_forces)
        
    return np.array(forces)

def read_vaps_xml_maxforces(xml_filepath : str):
    return np.max(np.linalg.norm(read_vaps_xml_forces(xml_filepath), axis = 2), axis = 1)


def plot_vasp_xml(xml_filepath : str):
    energies,_,_ = read_vasp_xml_energies(xml_filepath)
    max_forces = read_vaps_xml_maxforces(xml_filepath)

    fig, axis = plt.subplots(2,1,sharex=True)
    
    ax = axis[0]
    ax.plot(energies)
    ax.set_ylabel("Energy [eV]")

    ax = axis[1]
    ax.plot(max_forces)
    ax.set_ylabel("Maximum Force [eV/Å]")
    ax.set_xlabel("Number of ionic steps")
    


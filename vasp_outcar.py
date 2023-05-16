from ase.atom import Atom
from ase.atoms import Atoms
from typing import List
import numpy as np


class EnergyForceCalculator:
    """A calculator that returns the energy and forces from input, workaround to use OUTCAR files in ASE"""

    def __init__(self, energy, forces, charges, magmoms) -> None:
        self.energy = energy
        self.forces = forces
        self.charges = charges
        self.magmoms = magmoms
        self.results = {
            "energy": self.energy,
            "forces": self.forces,
            "charges": self.charges,
            "magmoms": self.magmoms,
        }

    def get_forces(self, Atoms):
        return self.forces

    def get_potential_energy(self, atoms, force_consistent=False):
        return self.energy

    def get_charges(self, atoms):
        return self.charges

    def get_magnetic_moments(self, atoms):
        return self.magmoms


def read_outcar(outcar_file_path):
    """Reads an OUTCAR file and returns a list of ASE Atoms objects with energy and forces"""
    with open(outcar_file_path, "r") as f:
        lines = f.readlines()
        structures = []
        charges = np.NaN
        magnetization = np.NaN
        for i, line in enumerate(lines):
            if "POSCAR:" in line:
                elements = line.split()[1:]

            if "  direct lattice vectors " in line:
                cell_vectors = np.array(
                    [lines[i + 1 + j].split()[0:3] for j in range(3)], dtype=np.float64
                )

            if "ions per type =" in line:
                number_of_ions = [int(num_str) for num_str in line.split()[4:]]
                tot_number_atoms = sum(number_of_ions)
                elements_array = []
                for element, number in zip(elements, number_of_ions):
                    elements_array += [element] * number

            if line == " magnetization (x)\n":
                magnetization = np.zeros(tot_number_atoms)
                for j in range(tot_number_atoms):
                    magnetization[j] = float(lines[i + j + 4].split()[-1])

            if line == " total charge\n":
                charges = np.zeros(tot_number_atoms)
                for j in range(tot_number_atoms):
                    charges[j] = float(lines[i + j + 4].split()[-1])

            if (
                line
                == " POSITION                                       TOTAL-FORCE (eV/Angst)\n"
            ):
                positions_forces = np.zeros((tot_number_atoms, 6))
                for j in range(tot_number_atoms):
                    positions_forces[j, :] = np.array(
                        lines[i + j + 2].split(), dtype=np.float64
                    )
                energy = float(lines[i + j + 13].split()[-2])
                structures.append(
                    Atoms(
                        elements_array,
                        positions=positions_forces[:, :3],
                        cell=cell_vectors,
                        pbc=True,
                        calculator=EnergyForceCalculator(
                            energy, positions_forces[:, 3:], charges, magnetization
                        ),
                    )
                )

    return structures


def read_outcar_final(outcar_file_path):
    """Reads an OUTCAR file and returns the lase ASE Atoms objects with energy and forces"""
    return read_outcar(outcar_file_path)[-1]


def read_energies_forces_neb_outar(outcar_file_path):
    with open(outcar_file_path, "r") as f:
        lines = f.readlines()
        energies = []
        forces = []
        for i, line in enumerate(lines):
            if "energy  without" in line:
                energies.append(float(line.split()[3]))
            if "max at" in line:
                forces.append(float(line.split()[4]))
    return energies, forces


if __name__ == "__main__":
    pass
    # from utils.master import *

    # read_energy_force_neb_outar(
    #     r"C:\Users\chris\OneDrive - NTNU\V2022\data_masteroppgave\dft-data\neb\10CO_coverage\CH_HtoCH2\01\OUTCAR"
    # )
    # structures = read_outcar(
    #     os.path.join(
    #         "C:\\Users",
    #         "chris",
    #         "OneDrive - NTNU",
    #         "V2022",
    #         "data_masteroppgave",
    #         "dft-data",
    #         "single/ni30_COs/CO12-13-22",
    #         "OUTCAR",
    #     )
    # )

    # from ase.visualize import view

    # view(structures, block=True)
    # # print(structures[0].get_potential_energy())
    # # print(structures[0].get_forces())
    # print(structures[0].get_charges())

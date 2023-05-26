#!/cluster/home/chrisafa/.conda/envs/agox_env/bin/python3.9
import matplotlib.pyplot as plt
import numpy as np
from utils.vasp_outcar import read_energies_forces_neb_outar
from utils.hpc_scripts.common_functions import find_numbered_folders
import os

def get_forces_image(folder):
    energies, forces = read_energies_forces_neb_outar(
        os.path.join(
            folder,
            "OUTCAR",
        )
    )
    return forces

if __name__ == "__main__":
    folders = find_numbered_folders(os.getcwd())
    n = len(folders)
    all_forces = [get_forces_image(os.path.join(os.getcwd(), folder)) for folder in folders[1:-1]]
    all_forces = np.array(all_forces)
    max_forces_iteration = np.max(all_forces, axis=0)
    i_min = np.argmin(max_forces_iteration)
    force_min = max_forces_iteration[i_min]
    print(f"Minimum force: {force_min} eV/Å on iteration {i_min} in a total of {len(max_forces_iteration)} iterations")

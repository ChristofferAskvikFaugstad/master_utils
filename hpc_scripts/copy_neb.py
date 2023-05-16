#!/usr/bin/python3
if __name__ == "__main__":
    from argparse import ArgumentParser
    from ase.io import read
    from ase.visualize import view
    import os
    from utils.vasp_outcar import (
        read_outcar_final,
        read_energies_forces_neb_outar,
        read_outcar,
    )
    from utils.data_handler import read_vasp_xml_final
    import matplotlib.pyplot as plt
    import numpy as np

    parser = ArgumentParser()
    parser.add_argument("-n")
    parser.add_argument("-i", default=-1)
    parser.add_argument("-d", default=-1)
    arguments = parser.parse_args()
    n = int(arguments.n)
    i = int(arguments.i)
    destination = arguments.d

    folders = [f"{i:02d}" for i in range(1, n)]
    image_trajectories = [
        read_outcar(os.path.join(folder, "OUTCAR")) for folder in folders
    ]
    for i, folder in enumerate(folders):
        print("Writing for folder", folder)
        image_trajectories[i].write(os.path.join(destination, folder, "POSCAR"))


#!/cluster/home/chrisafa/.conda/envs/master/bin/python3.10
if __name__ == "__main__":
    from ase.io import read
    from ase.visualize import view

    parser = ArgumentParser()
    parser.add_argument("-n")
    arguments = parser.parse_args()
    n = int(arguments.n)
    images = (
        [read(f"00/POSCAR")]
        + [read(f"{i:02d}/POSCAR") for i in range(1, n)]
        + [read(f"{n:02d}/POSCAR")]
    )
    view(images, block=True)

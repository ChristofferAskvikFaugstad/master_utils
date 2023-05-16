#!/cluster/home/chrisafa/.conda/envs/agox_env/bin/python3.9
def main(source = "."):
    from utils.hpc_scripts.common_functions import find_numbered_folders
    from utils.vasp_outcar import read_outcar
    import os
    from argparse import ArgumentParser
    
    parser = ArgumentParser()
    parser.add_argument("-i", default=-1)
    parser.add_argument("-d", default=-1)
    arguments = parser.parse_args()
    i = int(arguments.i)
    destination = arguments.d

    folders = find_numbered_folders(source)
    folders = [os.path.join(source, folder) for folder in folders]

    image_trajectories = [
        read_outcar(os.path.join(folder, "OUTCAR")) for folder in folders[1:-1]
    ]
    for i, folder in enumerate(folders[1:-1]):
        copy_from = os.path.join(folder,"OUTCAR")
        copy_to = os.path.join(destination, folder, "POSCAR")
        print(f"Copying from {copy_from} to {copy_to}")
        image_trajectories[i].write(os.path.join(destination, folder, "POSCAR"))


if __name__ == "__main__":
    main()

#!/cluster/home/chrisafa/.conda/envs/agox_env/bin/python3.9 # only valid for saga
def main(source = "."):
    from ase.io import read
    from ase.visualize import view
    from utils.hpc_scripts.common_functions import find_numbered_folders
    import os

    folders = find_numbered_folders(source)
    folders = [os.path.join(source, folder) for folder in folders]
    try:
        images = (
            [read(f"{folders[0]}/POSCAR")]
            + [read(f"{folder}/CONTCAR") for folder in folders]
            + [read(f"{folders[-1]}/POSCAR")]
        )
    except:
        images = (
            [read(f"{folders[0]}/POSCAR")]
            + [read(f"{folder}/POSCAR") for folder in folders]
            + [read(f"{folders[-1]}/POSCAR")]
        )
    view(images)


if __name__ == "__main__":
    main()
    


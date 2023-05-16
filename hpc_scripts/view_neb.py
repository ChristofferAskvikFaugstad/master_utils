def main():
    from ase.io import read
    from ase.visualize import view
    from utils.hpc_scripts.common_functions import find_numbered_folders

    folders = find_numbered_folders(".")
    try:
        images = (
            [read(f"00/POSCAR")]
            + [read(f"{folder}/CONTCAR") for folder in folders]
            + [read(f"{folders[-1]}/POSCAR")]
        )
    except:
        images = (
            [read(f"00/POSCAR")]
            + [read(f"{folder}/POSCAR") for folder in folders]
            + [read(f"{folders[-1]}/POSCAR")]
        )
    view(images)


if __name__ == "__main__":
    main()

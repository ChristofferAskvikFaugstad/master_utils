if __name__ == "__main__":
    from argparse import ArgumentParser
    from ase.io import read
    from ase.visualize import view
    import os
    from utils.vasp_outcar import read_outcar_final, read_energies_forces_neb_outar
    from utils.data_handler import read_vasp_xml_final
    import matplotlib.pyplot as plt
    import numpy as np

    parser = ArgumentParser()
    parser.add_argument("-n")
    n = int(parser.parse_args().n)

    folders = [f"{i:02d}" for i in range(n + 1)]
    image_trajectories = []
    try:
        initital = read_outcar_final(os.path.join(folders[0], "OUTCAR"))
    except:
        initial = read_vasp_xml_final(os.path.join(folders[0], "vasprun.xml"))

    try:
        final = read_outcar_final(os.path.join(folders[-1], "OUTCAR"))
    except:
        final = read_vasp_xml_final(os.path.join(folders[-1], "vasprun.xml"))

    images = [initial]
    # imagwa[read_outcar_final(os.path.join(folder, "OUTCAR")) folder in folders[1, -1]] #+ [final]

    fig, axs = plt.subplots(1, n - 1)

    for ax1, folder in zip(axs, folders):
        energies, forces = read_energies_forces_neb_outar(
            os.path.join(folders[0], "OUTCAR")
        )

        N = np.arange(len(energies))
        color = "tab:red"
        ax1.set_xlabel("Iterations")
        ax1.set_ylabel("Energy [eV]", color=color)
        ax1.plot(N, energies, color=color)
        ax1.tick_params(axis="y", labelcolor=color)

        ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

        color = "tab:blue"
        ax2.set_ylabel(
            "Max atomic force [eV/Å]", color=color
        )  # we already handled the x-label with ax1
        ax2.plot(N, forces, color=color)
        ax2.tick_params(axis="y", labelcolor=color)
    plt.show()
    view(images)

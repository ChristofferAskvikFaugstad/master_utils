import matplotlib.pyplot as plt
import numpy as np
from utils.vasp_outcar import read_energies_forces_neb_outar
import os

def find_numbered_folders(directory):
    numbered_folders = []

    for folder_name in os.listdir(directory):
        if os.path.isdir(os.path.join(directory, folder_name)):
            if folder_name.isdigit():
                numbered_folders.append(folder_name)

    return numbered_folders

def plot_folder(ax1, folder):
    energies, forces = read_energies_forces_neb_outar(
        os.path.join(
            folder,
            "OUTCAR",
        )
    )
    N = np.arange(len(energies))
    color = "tab:blue"
    ax1.set_xlabel("Iterations")
    ax1.set_ylabel("Energy [eV]", color=color)
    ax1.plot(N, energies,marker=".", color=color , linestyle = "-")
    ax1.tick_params(axis="y", labelcolor=color)

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    color = "tab:red"
    ax2.set_ylabel("Max atomic force [eV/Å]", color=color)  # we already handled the x-label with ax1
    ax2.plot(N, forces, color=color ,marker=".", linestyle = "-")
    ax2.hlines(0.05, 0, len(forces), color="tab:green")
    ax2.set_yscale("log")
    ax2.tick_params(axis="y", labelcolor=color)


def main(source = ""):
    folders = find_numbered_folders(
        source
    )
    n = len(folders)
    n_subplots = n-2
    fig, axs = plt.subplots( n_subplots, 1, figsize=(6, n_subplots*6))
    for ax1, folder in zip(axs,folders[1:-1]):
        plot_folder(ax1, os.path.join(source, folder))
    plt.tight_layout()
    plt.show()

def individual(source = ""):
    folders = find_numbered_folders(
        source
    )
    n = len(folders)
    n_subplots = n-2
    for folder in folders[1:-1]:
        ax = plt.figure(figsize=(5, 3)).add_subplot(111)
        plot_folder(ax, os.path.join(source, folder))
        ax.grid()
        plt.tight_layout()
        plt.show(block=False)
    plt.show()
    

def test():
    source = r"C:\Users\chris\OneDrive - NTNU\V2022\data_masteroppgave\dft-data\neb\C_CH2toC2H2"
    individual(source)
    



if __name__ == "__main__":
    test()
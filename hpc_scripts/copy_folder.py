#!/cluster/home/chrisafa/.conda/envs/agox_env/bin/python3.9
from argparse import ArgumentParser
import os
import shutil
import glob

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-d")
    parser.add_argument("-s", default=os.getcwd())
    destination = parser.parse_args().d

    source = parser.parse_args().s

    print(source)

    folders = glob.glob(f"{source}/*")
    for folder in folders:
        files = ["OUTCAR", "vasprun.xml"]
        for file in files:
            potential_file = os.path.join(folder, file)
            if os.path.isfile(potential_file):
                folder_basename = os.path.basename(folder)
                destination_folder = os.path.join(destination, folder_basename)
                os.makedirs(destination_folder, exist_ok=True)
                destination_file = os.path.join(destination_folder, file)
                shutil.copyfile(potential_file, destination_file)
                print(f"copied {potential_file} to {destination_file}")

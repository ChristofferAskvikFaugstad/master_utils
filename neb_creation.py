from utils.master import *
from utils.data_handler import get_vasp_calculation_names_path, get_structure_path
import glob
import shutil


class NEBCalculation:
    def __init__(
        self, num_images: int, initial: Atoms, final: Atoms, store_folder: str
    ) -> None:
        self.num_images = num_images
        self.initial = initial
        self.final = final
        self.store_folder = store_folder

    @classmethod
    def init_from_folder(
        cls, source_folder: str, store_folder: str, relative=True
    ) -> None:
        folders = get_vasp_calculation_names_path(source_folder, relative=relative)
        num_images = len(folders)

        initial = get_structure_path(os.path.join(source_folder, folders[0]))
        final = get_structure_path(os.path.join(source_folder, folders[-1]))
        neb_calculation = cls(num_images, initial, final, store_folder)
        neb_calculation.images = [
            get_structure_path(os.path.join(source_folder, folder))
            for folder in folders
        ]
        return neb_calculation

    def make_neb_folders(self):
        for i in range(self.num_images):
            folder = os.path.join(self.store_folder, f"{i:02d}")
            if not os.path.exists(folder):
                os.makedirs(folder)


from data_handler import NEBSTART


class NEBWriter:
    def __init__(
        self,
        reaction: str,
        store_folder: str,
        location: str = NEBSTART,
        time: str = "0-2:0:0",
        nodes: int = None,
        mem_per_cpu: str = "4G",
        saga: bool = True,
        **kwargs,
    ) -> None:
        self.name = reaction
        self.store_folder = store_folder
        if not os.path.exists(store_folder):
            os.makedirs(store_folder)
        self.job_folder = os.path.join(store_folder, "job")
        if not os.path.exists(self.job_folder):
            os.makedirs(self.job_folder)
        self.images_names = get_vasp_calculation_names_path(
            os.path.join(location, reaction), relative=False, outcar_reader=False
        )
        self.num_images = len(self.images_names) - 2
        self.images = [
            get_structure_path(
                os.path.join(location, reaction, image),
                outcar_reader=False,
                relative=False,
            )
            for image in self.images_names
        ]
        if nodes is None:
            nodes = self.num_images
        self.nodes = nodes
        self.location = location
        self.time = time
        self.mem_per_cpu = mem_per_cpu
        self.saga = saga
        self.kwargs = {
            "prec": '"Accurate"',
            "images": self.num_images,  # determines the number of NEB images
            "istart": 0,  # start new calculation
            "xc": '"PBE"',  # set PBE exchange correlation functional
            "ismear": 0,  # Gaussian smearing
            "sigma": 0.02,  # smearing constant
            "encut": 400,  # energy cutoff
            "nelm": 60,  # max number of SC steps default = 60
            "lorbit": 10,  # write magnetic moments
            "isym": -1,  # ingnore symmetires
            "ispin": 2,  # for spin calculations
            "ediffg": -0.05,  # force criteria
            "ibrion": 3,  # to make vasp not change
            "potim": 0,  # no vasp step
            "lcharg": False,  # toggels if CHGCAR is written, default = True
            "lwave": False,  # toggels if WAVECAR is written, default = True
            "ncore": 20 if self.saga else 32,  # number of cores to use
            "nsw": 1200,  # Number of steps
            "ichain": 0,  # indicate which TS algorithm, 0 - NEB
            "iopt": 3,  # Choose which neb optimizer
            "timestep": 0.1,  # General safe value fuond from vtst
            "maxmove": 0.1,  # General safe value fuond from vtst
            "lclimb": True,  # Climbing image NEB
            "spring": -5,  # spring constant
        }
        self.kwargs.update(kwargs)

    def write(self):
        self.write_images()
        self.write_python_script()
        self.write_slurm_script()
        self.convert_to_unix()
        self.copy_initial_and_final_data()

    def copy_initial_and_final_data(self):
        for file in ["OUTCAR", "vasprun.xml"]:
            try:
                shutil.copy(
                    os.path.join(self.location, self.name, self.images_names[0], file),
                    os.path.join(self.job_folder, self.images_names[0]),
                )
                shutil.copy(
                    os.path.join(self.location, self.name, self.images_names[-1], file),
                    os.path.join(self.job_folder, self.images_names[-1]),
                )
            except FileNotFoundError:
                pass

    def convert_to_unix(self):
        for file in glob.glob(f"{self.store_folder}/**", recursive=True):
            if os.path.isfile(os.path.join(file)):
                self.convert_file_to_unix(file)

    def write_images(self):
        for i, image in enumerate(self.images):
            if not os.path.exists(os.path.join(self.job_folder, f"{i:02d}")):
                os.makedirs(os.path.join(self.job_folder, f"{i:02d}"))
            image.write(os.path.join(self.job_folder, f"{i:02d}", "POSCAR"))

    def convert_file_to_unix(self, file_path):
        WINDOWS_LINE_ENDING = b"\r\n"
        UNIX_LINE_ENDING = b"\n"

        with open(file_path, "rb") as open_file:
            content = open_file.read()

        # Windows ➡ Unix
        content = content.replace(WINDOWS_LINE_ENDING, UNIX_LINE_ENDING)
        with open(file_path, "wb") as open_file:
            open_file.write(content)

    def write_python_script(self):
        with open(os.path.join(self.job_folder, "neb.py"), mode="w") as py:
            py.write(
                f"""from ase.io.vasp import read_vasp_xml
import os

dir = "00" 
initial_xml = "vasprun.xml"
initial = next(read_vasp_xml(os.path.join(dir, initial_xml)))

from ase.calculators.vasp import Vasp
calc = Vasp("""
            )
            for key, value in self.kwargs.items():
                py.write(f"{key}={value},\n")
            py.write(
                f""")
command = calc.make_command(calc.command)
properties = ("energy",)
from ase.calculators import calculator

system_changes = tuple(calculator.all_changes)
calc.write_input(initial, properties, system_changes)
error = calc._run(command=command, directory=calc.directory)
print(error)"""
            )

    def write_slurm_script(self):
        neb_file = os.path.join(self.store_folder, "neb.sh")
        with open(neb_file, mode="w") as slurm:
            slurm.write(
                f"""#!/bin/sh
#SBATCH --account=nn9497k
#SBATCH --job-name={self.name}
#SBATCH --time={self.time}
#SBATCH --nodes={self.nodes}
#SBATCH --output=slurm-%j.out
"""
            )
            if self.saga:
                slurm.write("#SBATCH --ntasks-per-node=40\n")
                slurm.write(f"#SBATCH --mem-per-cpu={self.mem_per_cpu}\n")
            else:
                pass
                # slurm.write("#SBATCH --ntasks-per-node=32\n")
                # slurm.write("#SBATCH --partition=bigmem\n")

            slurm.write(f"")
            slurm.write(
                f"""NAME="{self.name}"
# Load Vasp modules
module --force purge
module load StdEnv
module load VASPModules
module load VaspExtra
module load VASP/6.3.2-intel-2021b-std-wannier90-libxc-hdf5-beef-vtst-097e6cb5a78f237dc588ba9c7877f23b # for the VTST parts\n"""
            )
            if self.saga:
                slurm.write(
                    """# Loads personalized enviorment on saga\nmodule load Miniconda3/4.9.2\nexport LMOD_DISABLE_SAME_NAME_AUTOSWAP=no\nexport PS1=\$\nsource ${EBROOTMINICONDA3}/etc/profile.d/conda.sh\nconda deactivate &>/dev/null\nconda activate /cluster/home/chrisafa/.conda/envs/agox_env\n"""
                )
            else:
                slurm.write(
                    """# Loads personalized enviorment on fram\nmodule load Anaconda3/2022.05\nexport PS1=\$\nsource ${EBROOTANACONDA3}/etc/profile.d/conda.sh\nconda deactivate &>/dev/null\nconda activate /cluster/home/chrisafa/.conda/envs/master\n"""
                )
            slurm.write(
                """WORKDIR=$USERWORK/temp/$SLURM_JOB_ID
mkdir -p $WORKDIR
cp -r job/* $WORKDIR
cd ${WORKDIR}
python *.py
cat *.py
# STOREDIR=$USERWORK/temp/out/$NAME
# mkdir -p $STOREDIR
# Copy only relevant folders into dat folder
# for folder in $(ls -d $WORKDIR*);do
#     f="$(basename -- $folder)"
#     mkdir -p $USERWORK/temp/out/"$NAME"/"$f"
#     for file in OUTCAR OSZICAR vasprun.xml; do
#         cp "$folder"/"$file" $USERWORK/temp/out/"$NAME"/"$f"/"$file"
#     done;
# done
# cp $SLURM_SUBMIT_DIR/slurm-$SLURM_ARRAY_JOB_ID.out $STOREDIR/"out.out"""
            )


if __name__ == "__main__":
    # neb_writer = NEBWriter(
    #     "COtoC_O", "code/temp/COtoC_O", saga=False, nsw=500, time="2-00:00:00"
    # ).write()
    import os

    os.chdir("C:/Users/chris/masteroppgave_code/")

    names = [
        "CH_HtoCH2",  # 10CO current path seemes resonable
        "CH2_HtoCH3",  # 10CO current path seemes resonable
        "COHtoC_OH",  # 10 CO current path seemes resonable
        "HCOHtoCH_OH",  # 10 CO current path seemes resonable
        "HCOtoCH_O",  # 10 CO current path seemes resonable
    ]
    for name in names:
        neb_writer = NEBWriter(
            f"{name}_neb",
            f"code/temp/{name}10CO",
            location="in/redo_neb",
            saga=False,
            nsw=800,
            timestep=0.01,
            time="1-0:00:00",
        ).write()

    # folder = "code/temp/COtoC_O"
    # print(
    #     [
    #         file
    #         for file in glob.glob(f"{folder}/**", recursive=True)
    #         if os.path.isfile(os.path.join(file))
    #     ]
    # )
    # neb_writer.write_python_script()
    # neb_writer.write_images()
    # images = [
    #     read("code/temp/COtoC_O/job/{i:02d}/POSCAR".format(i=i)) for i in range(8)
    # ]
    # view(images, block=True)
    pass

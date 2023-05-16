from utils.master import *
from utils.data_handler import get_vasp_calculation_names_path, get_structure_path
import glob
import shutil

from utils.data_handler import NEBSTART


class CalculationWriter:
    def __init__(
        self,
        name: str,
        store_folder: str,
        time: str = "0-2:0:0",
        nodes: int = None,
        mem_per_cpu: str = "4G",
        saga: bool = True,
        **kwargs,
    ) -> None:
        self.name = name
        self.nodes = nodes
        self.time = time
        self.mem_per_cpu = mem_per_cpu
        self.saga = saga
        self.store_folder = store_folder
        if not os.path.exists(store_folder):
            os.makedirs(store_folder)
        self.job_folder = os.path.join(store_folder, "job")
        if not os.path.exists(self.job_folder):
            os.makedirs(self.job_folder)
        self.kwargs = {
            "prec": '"Accurate"',
            "istart": 0,  # start new calculation
            "xc": '"PBE"',  # set PBE exchange correlation functional
            "ismear": 0,  # Gaussian smearing
            "sigma": 0.02,  # smearing constant
            "encut": 400,  # energy cutoff
            "nelm": 60,  # max number of SC steps default = 60
            "lorbit": 10,  # write magnetic moments
            "isym": -1,  # ingnore symmetires
            "ispin": 2,  # for spin calculations
            "lcharg": False,  # toggels if CHGCAR is written, default = True
            "lwave": False,  # toggels if WAVECAR is written, default = True
            "ncore": 20 if self.saga else 32,  # number of cores to use
        }
        self.kwargs.update(kwargs)

    def get_slurm_header(self):
        header = f"""#!/bin/sh
#SBATCH --account=nn9497k
#SBATCH --job-name={self.name}
#SBATCH --time={self.time}
#SBATCH --nodes={self.nodes}
#SBATCH --output=slurm-%j.out
"""
        if self.saga:
            header += "#SBATCH --ntasks-per-node=40\n"
            header += f"#SBATCH --mem-per-cpu={self.mem_per_cpu}\n"
        else:
            pass
        return header

    def get_slurm_body(self):
        # implemented by submodules
        return None

    def convert_to_unix(self):
        for file in glob.glob(f"{self.store_folder}/**", recursive=True):
            if os.path.isfile(os.path.join(file)):
                self.convert_file_to_unix(file)

    def convert_file_to_unix(self, file_path):
        WINDOWS_LINE_ENDING = b"\r\n"
        UNIX_LINE_ENDING = b"\n"

        with open(file_path, "rb") as open_file:
            content = open_file.read()

        # Windows ➡ Unix
        content = content.replace(WINDOWS_LINE_ENDING, UNIX_LINE_ENDING)
        with open(file_path, "wb") as open_file:
            open_file.write(content)

    def write_input(self):
        # implemented by submodules
        pass

    def write_python_script(self):
        # implemented by submodules
        pass

    def write_slurm_script(self):
        with open(os.path.join(self.store_folder, "slurm.sh"), "w") as slurm:
            slurm.write(self.get_slurm_header())
            slurm.write(self.get_slurm_body())

    def write(self):
        self.write_input()
        self.write_python_script()
        self.write_slurm_script()
        self.convert_to_unix()


class NEBWriter(CalculationWriter):
    def __init__(
        self,
        reaction: str,
        store_folder: str,
        images: List[Atoms],
        time: str = "0-2:0:0",
        nodes: int = None,
        mem_per_cpu: str = "4G",
        saga: bool = True,
        **kwargs,
    ) -> None:

        self.reaction = reaction
        self.images = images
        self.num_images = len(images) - 2
        self.location = None

        if nodes is None:
            nodes = self.num_images
        spesific_kwargs = {
            "images": self.num_images,  # determines the number of NEB images
            "ediffg": -0.05,  # force criteria
            "ibrion": 3,  # to make vasp not change
            "potim": 0,  # no vasp step
            "nsw": 1200,  # Number of steps
            "ichain": 0,  # indicate which TS algorithm, 0 - NEB
            "iopt": 3,  # Choose which neb optimizer
            "timestep": 0.1,  # General safe value fuond from vtst
            "maxmove": 0.1,  # General safe value fuond from vtst
            "lclimb": True,  # Climbing image NEB
            "spring": -5,  # spring constant
        }
        large_coverage = kwargs.get("large_coverage", False)
        if large_coverage:
            spesific_kwargs.update(
                {
                    # Changes made to improve convergence for higher coverage
                    "ediff": 1e-7,  # lower ediff
                    "nedos": 50000,  # Way more DOS points
                    "algo": "'VeryFast'",  # electronic minimiztion algoirthm
                }
            )
            kwargs.pop("large_coverage")

        spesific_kwargs.update(kwargs)

        super().__init__(
            reaction,
            store_folder,
            time=time,
            nodes=nodes,
            mem_per_cpu=mem_per_cpu,
            saga=saga,
            **spesific_kwargs,
        )

    @classmethod
    def from_previous_calculation(
        cls,
        reaction: str,
        store_folder: str,
        location: str = NEBSTART,
        time: str = "0-2:0:0",
        nodes: int = None,
        mem_per_cpu: str = "4G",
        saga: bool = True,
        **kwargs,
    ) -> None:

        images_names = get_vasp_calculation_names_path(
            os.path.join(location, reaction), relative=False, outcar_reader=False
        )
        num_images = len(images_names) - 2
        images = [
            get_structure_path(
                os.path.join(location, reaction, image),
                outcar_reader=False,
                relative=False,
            )
            for image in self.images_names
        ]
        self = cls(
            reaction, store_folder, images, time, nodes, mem_per_cpu, saga, **kwargs
        )
        self.location = location
        return self

    def get_slurm_body(self):
        body = f"""NAME="{self.name}"
# Load Vasp modules
module --force purge
module load StdEnv
module load VASPModules
module load VaspExtra
module load VASP/6.3.2-intel-2021b-std-wannier90-libxc-hdf5-beef-vtst-097e6cb5a78f237dc588ba9c7877f23b # for the VTST parts\n"""
        if self.saga:
            body += """# Loads personalized enviorment on saga\nmodule load Miniconda3/4.9.2\nexport LMOD_DISABLE_SAME_NAME_AUTOSWAP=no\nexport PS1=\$\nsource ${EBROOTMINICONDA3}/etc/profile.d/conda.sh\nconda deactivate &>/dev/null\nconda activate /cluster/home/chrisafa/.conda/envs/agox_env\n"""
        else:
            body += """# Loads personalized enviorment on fram\nmodule load Anaconda3/2022.05\nexport PS1=\$\nsource ${EBROOTANACONDA3}/etc/profile.d/conda.sh\nconda deactivate &>/dev/null\nconda activate /cluster/home/chrisafa/.conda/envs/master\n"""
        body += """WORKDIR=$USERWORK/temp/$SLURM_JOB_ID
mkdir -p $WORKDIR
cp -r job/* $WORKDIR
cd ${WORKDIR}
python *.py
cat *.py
"""
        return body

    def write_input(self):
        self.write_images()
        if self.location is not None:
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

    def write_images(self):
        for i, image in enumerate(self.images):
            if not os.path.exists(os.path.join(self.job_folder, f"{i:02d}")):
                os.makedirs(os.path.join(self.job_folder, f"{i:02d}"))
            image.write(os.path.join(self.job_folder, f"{i:02d}", "POSCAR"))

    def write_python_script(self):
        with open(os.path.join(self.job_folder, "neb.py"), mode="w") as py:
            py.write(
                f"""from ase.io import read
import os

dir = "00" 
initial = "POSCAR"
initial = read(os.path.join(dir, initial))

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


class SEARCHWriter(CalculationWriter):
    def __init__(
        self,
        name: str,
        store_folder: str,
        time: str = "0-2:0:0",
        nodes: int = 4,
        mem_per_cpu: str = "4G",
        saga: bool = True,
        large_coverage: bool = False,
        **kwargs,
    ) -> None:
        spesific_kwargs = {
            "ibrion": 2,
            "ediffg": -0.02,
            "nsw": 200,
        }  # No spesific vasp kwargs yet
        if large_coverage:
            spesific_kwargs.update(
                {
                    # Changes made to improve convergence for higher coverage
                    "ediff": 1e-7,  # lower ediff
                    "nedos": 50000,  # Way more DOS points
                    "algo": "'VeryFast'",  # electronic minimiztion algoirthm
                }
            )
        spesific_kwargs.update(kwargs)

        super().__init__(
            name,
            store_folder,
            time=time,
            nodes=nodes,
            mem_per_cpu=mem_per_cpu,
            saga=saga,
            **spesific_kwargs,
        )

    def get_slurm_body(self):
        body = f"""#SBATCH --array=0-1
NAME="{self.name}"
# Load Vasp modules
module --force purge
module load StdEnv
module load VASPModules
module load VaspExtra
"""
        if self.saga:
            body += " module load VASP/6.3.2-intel-2022b-std-wannier90-libxc-hdf5-beef-0a928426e459cf2aeab3d0bf8f441c74\n"
            body += """# Loads personalized enviorment on saga\nmodule load Miniconda3/4.9.2\nexport LMOD_DISABLE_SAME_NAME_AUTOSWAP=no\nexport PS1=\$\nsource ${EBROOTMINICONDA3}/etc/profile.d/conda.sh\nconda deactivate &>/dev/null\nconda activate /cluster/home/chrisafa/.conda/envs/agox_env\n"""
        else:
            body += "module load VASP/6.3.2-intel-2021b-std-wannier90-libxc-hdf5-beef-d7238be44ec2ed23315a16cc1549a1e3 # works fine for normal\n"
            body += """# Loads personalized enviorment on fram\nmodule load Anaconda3/2022.05\nexport PS1=\$\nsource ${EBROOTANACONDA3}/etc/profile.d/conda.sh\nconda deactivate &>/dev/null\nconda activate /cluster/home/chrisafa/.conda/envs/master\n"""
        body += """WORKDIR=$USERWORK/temp/$SLURM_ARRAY_JOB_ID/$SLURM_ARRAY_TASK_ID
mkdir -p $WORKDIR
mkdir -p dat/$SLURM_ARRAY_TASK_ID
cp -r job/* $WORKDIR
cd ${WORKDIR}
python *.py -i $SLURM_ARRAY_TASK_ID
echo \n
cat *.py
"""
        return body

    def get_slurm_header(self):
        header = super().get_slurm_header()
        header = header.replace(
            "#SBATCH --output=slurm-%j.out", "#SBATCH --output=slurm-%A-%a.out"
        )
        return header

    def write_input(self):
        pass

    def write_python_script(self):
        with open(os.path.join(self.job_folder, "search.py"), mode="w") as py:
            py.write(
                """from ase.calculators.vasp import Vasp
from ase.io import read
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-i")
i = int(parser.parse_args().i)

structure = read(f"{i}.xyz")

calc = Vasp("""
            )
            for key, value in self.kwargs.items():
                py.write(f"{key}={value},\n")
            py.write(
                f""")

structure.calc = calc
# Do calculations
structure.get_potential_energy()"""
            )


class GOFEEWriter(CalculationWriter):
    def __init__(
        self,
        name: str,
        store_folder: str,
        training_database_path: str,
        gofee_kwargs: dict,
        template_path: str,
        time: str = "0-2:0:0",
        nodes: int = 2,
        mem_per_cpu: str = "4G",
        saga: bool = True,
        **kwargs,
    ) -> None:

        self.training_database_path = training_database_path
        self.template_path = template_path
        spesific_kwargs = {}  # No spesific vasp kwargs yet
        spesific_kwargs.update(kwargs)
        self.gofee_kwargs = {
            "training_database_path": f"'{os.path.basename(training_database_path)}'",
            "output_database_path": "'db.db'",
            "template_path": f"'{os.path.basename(template_path)}'",
            "iteration_start_training": 1,
            "update_interval": 10,
            "start_relax": 2,
            "sample_size": 10,
            "num_candidates": {
                0: [10, 0],  # nCOhandler  # Rattlegenerator
                10: [5, 5],
            },  # distribution between candidates,
            "kappa": 2,
            "f_max_dft": 0.01,
            "n_steps_dft": 1,
            "f_max_relax_post": 0.01,
            "n_steps_relax_post": 400,
            "N_iterations": 200,
            "seed": 0,
            "confinement_cell": "np.eye(3) * 16",
            "confinement_corner": "np.array([0, 0, 0])",
            "write_frequency": 1,  # frequency of database
        }

        self.gofee_kwargs.update(gofee_kwargs)

        super().__init__(
            name,
            store_folder,
            time=time,
            nodes=nodes,
            mem_per_cpu=mem_per_cpu,
            saga=saga,
            **spesific_kwargs,
        )

    def get_slurm_body(self):
        body = f"""NAME="{self.name}"
# Load Vasp modules
module --force purge
module load StdEnv
module load VASPModules
module load VaspExtra
module load VASP/6.3.2-intel-2021b-std-wannier90-libxc-hdf5-beef-d7238be44ec2ed23315a16cc1549a1e3 # works fine for normal\n"""
        if self.saga:
            body += """# Loads personalized enviorment on saga\nmodule load Miniconda3/4.9.2\nexport LMOD_DISABLE_SAME_NAME_AUTOSWAP=no\nexport PS1=\$\nsource ${EBROOTMINICONDA3}/etc/profile.d/conda.sh\nconda deactivate &>/dev/null\nconda activate /cluster/home/chrisafa/.conda/envs/agox_env\n"""
        else:
            body += """# Loads personalized enviorment on fram\nmodule load Anaconda3/2022.05\nexport PS1=\$\nsource ${EBROOTANACONDA3}/etc/profile.d/conda.sh\nconda deactivate &>/dev/null\nconda activate /cluster/home/chrisafa/.conda/envs/master\n"""
        body += """WORKDIR=$USERWORK/temp/$SLURM_JOB_ID
mkdir -p $WORKDIR
cp -r job/* $WORKDIR
cd ${WORKDIR}
python *.py
cat *.py
"""
        return body

    def copy_database(self):
        shutil.copy(
            self.training_database_path,
            os.path.join(
                self.job_folder, os.path.basename(self.training_database_path)
            ),
        )

    def copy_template(self):
        shutil.copyfile(
            self.template_path,
            os.path.join(self.job_folder, os.path.basename(self.template_path)),
        )

    def write_input(self):
        self.copy_database()
        self.copy_template()

    def write_python_script(self):
        with open(os.path.join(self.job_folder, "neb.py"), mode="w") as py:
            py.write(f"import numpy as np\n")
            # Write GOFEE settings
            for key, value in self.gofee_kwargs.items():
                py.write(f"{key}={value}\n")
            # Write GOFEE Imports
            py.write(
                """import matplotlib
matplotlib.use("Agg")
from agox import AGOX
from agox.databases import Database
from agox.environments import Environment
from agox.evaluators import LocalOptimizationEvaluator
from agox.generators import RattleGenerator
from agox.samplers import KMeansSampler
from agox.models import ModelGPR
from agox.acquisitors import LowerConfidenceBoundAcquisitor
from agox.postprocessors import ParallelRelaxPostprocess
from agox.collectors import ParallelCollector
from agox.collectors import StandardCollector
from agox.candidates import StandardCandidate
from ase.io import read
from ase import Atoms
# personaly made functions to interact with vaspfiles
from utils.IOdatabase import IODataBase
from utils.handlers.COHandler import nCOHandler
from utils.manual_generator import ManualGenerator, ManualObserver
# setup VASP calculator
from ase.calculators.vasp import Vasp
"""
            )

            # Write calculator
            py.write("""calc = Vasp(""")
            for key, value in self.kwargs.items():
                py.write(f"{key}={value},\n")
            py.write(f""") """)
            # Write GOFEE body
            py.write(
                """
# creat empty template
template = read(template_path)

# describe the enviorment for AGOX
environment = Environment(
    template=Atoms("", cell=confinement_cell, pbc=True),
    symbols=symbols_in_search,
    confinement_cell=confinement_cell,
    confinement_corner=confinement_corner,
)

# Database
database = Database(
    filename=output_database_path,
    order=7,
    write_frequency=write_frequency,
    initialize=True,
)

input_db = IODataBase(training_database_path)

for struc in input_db.get_all_candidates():
    candidate = StandardCandidate(template=template, **struc.todict())
    candidate.calc = struc.calc
    database.store_candidate(candidate)

# initialize default Gaussian Process Regression model
model = ModelGPR.default(environment, database)
model.iteration_start_training = iteration_start_training
model.update_interval = update_interval

# setup generator through a sampler
sampler = KMeansSampler(
    feature_calculator=model.get_feature_calculator(),
    database=database,
    sample_size=sample_size,
    order=1,
)

template = read(template_path)

nCOhandler = nCOHandler(n_CO)
manual_generator = ManualGenerator(
    nCOhandler, use_mic=False, **environment.get_confinement(), template=template
)
manual_observer = ManualObserver(nCOhandler, order=6)


rattle_generator = RattleGenerator(use_mic=False, **environment.get_confinement())
generators = [manual_generator, rattle_generator]

# controls the sampling from both generators
collector = StandardCollector(
    generators=generators,
    sampler=sampler,
    environment=environment,
    num_candidates=num_candidates,
    order=2,
)

# Acuistor function initialization
acquisitor = LowerConfidenceBoundAcquisitor(
    model_calculator=model, kappa=kappa, order=4
)


# Relaxing in surrogate potential
relaxer = ParallelRelaxPostprocess(
    model=acquisitor.get_acquisition_calculator(),
    optimizer_run_kwargs={"fmax": f_max_relax_post, "steps": n_steps_relax_post},
    constraints=environment.get_constraints(),
    order=3,
    start_relax=start_relax,
    tmp_dir="/cluster/work/users/chrisafa/r"
)


# VASP evaluation and one step
evaluator = LocalOptimizationEvaluator(
    calc,
    gets={"get_key": "prioritized_candidates"},
    optimizer_kwargs={"logfile": None},
    store_trajectory=True,
    optimizer_run_kwargs={"fmax": f_max_dft, "steps": n_steps_dft},
    order=6,
)

# collect agox modules and run algorithm
agox = AGOX(
    collector, acquisitor, relaxer, database, evaluator, manual_observer, seed=seed
)
agox.run(N_iterations=N_iterations)
"""
            )


#             py.write(
#                 f""")
# command = calc.make_command(calc.command)
# properties = ("energy",)
# from ase.calculators import calculator

# system_changes = tuple(calculator.all_changes)
# calc.write_input(initial, properties, system_changes)
# error = calc._run(command=command, directory=calc.directory)
# print(error)"""
#             )


if __name__ == "__main__":
    # neb_writer = NEBWriter(
    #     "COtoC_O", "code/temp/COtoC_O", saga=False, nsw=500, time="2-00:00:00"
    # ).write()
    import os

    os.chdir("C:/Users/chris/masteroppgave_code/")

    # names = [
    #     "CH_HtoCH2",  # 10CO current path seemes resonable
    #     "CH2_HtoCH3",  # 10CO current path seemes resonable
    #     "COHtoC_OH",  # 10 CO current path seemes resonable
    #     "HCOHtoCH_OH",  # 10 CO current path seemes resonable
    #     "HCOtoCH_O",  # 10 CO current path seemes resonable
    # ]
    # for name in names:
    # name = "CH3_HtoCH4"
    # neb_writer = NEBWriter(
    #     f"{name}_neb",
    #     f"code/temp/{name}_10COt",
    #     location="in/redo_neb",
    #     saga=False,
    #     nsw=1000,
    #     timestep=0.01,
    #     time="1-0:00:00",
    # ).write()
    # searchwriter = SEARCHWriter(
    #     "COtoC_O",
    #     "code/temp/COtoC_O_src",
    #     saga=False,
    #     nsw=1500,
    # ).write()
    # name = "20CO"
    # folder = "code/coverage/20CO/20COscr"
    # SEARCHWriter(
    #     name,
    #     folder,
    #     time="24:00:00",
    #     nodes=4,
    #     saga=False,
    #     nsw=1000,
    #     magmoms=[0.7] * 30 + [0] * 40,
    #     large_coverage=True,
    # ).write()
    GOFEEWriter(
        "10CO_gofee",
        "code/temp/10CO_gofee",
        r"code\coverage\2CO\2CO_all_14_03.db",
        {"n_CO": 10, "symbols_in_search": "'Ni30C10O10'"},
        r"code\coverage\2CO\2COgofee\job\POSCAR",
        saga=False,
        time="3-00:00:00",
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

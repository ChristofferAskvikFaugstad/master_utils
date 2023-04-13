from utils.manual_generator import *
from utils.manual_placement import *
from agox.candidates import StandardCandidate
import random
from itertools import combinations


class nCOHandler(IDHandler):
    """General class for n coadsorbed CO moleules.
    The follwing assumptions are made:
        - The handler places all sites with same probability.
        - The C always binds to Ni
        - Bridge and Hollow CO cannot share two Ni atoms
        - Top cannot share atom with bridge CO
        - Top cannot share Ni atoms with two other CO's
    """

    # Sites found using code
    avail_top = list((i,) for i in np.arange(0, 26))
    avail_bridge = [
        (0, 7),
        (0, 8),
        (0, 20),
        (0, 21),
        (1, 4),
        (1, 14),
        (1, 15),
        (1, 23),
        (1, 24),
        (2, 5),
        (2, 11),
        (2, 12),
        (2, 16),
        (2, 19),
        (3, 6),
        (3, 10),
        (3, 13),
        (3, 17),
        (3, 18),
        (4, 14),
        (4, 15),
        (4, 20),
        (4, 21),
        (5, 11),
        (5, 12),
        (5, 23),
        (6, 10),
        (6, 13),
        (6, 24),
        (6, 25),
        (7, 9),
        (7, 16),
        (7, 19),
        (7, 20),
        (8, 9),
        (8, 17),
        (8, 18),
        (8, 21),
        (9, 16),
        (9, 17),
        (9, 22),
        (10, 15),
        (10, 18),
        (10, 24),
        (11, 14),
        (11, 19),
        (11, 23),
        (12, 13),
        (12, 22),
        (12, 25),
        (13, 22),
        (13, 25),
        (14, 19),
        (14, 20),
        (14, 23),
        (15, 18),
        (15, 21),
        (15, 24),
        (16, 19),
        (16, 22),
        (17, 18),
        (17, 22),
        (18, 21),
        (19, 20),
        (20, 21),
        (23, 24),
        (23, 25),
        (24, 25),
    ]
    avail_hollow = [
        (6, 13, 25),
        (6, 24, 25),
        (8, 17, 18),
        (11, 14, 23),
        (8, 9, 17),
        (2, 16, 19),
        (10, 15, 18),
        (10, 15, 24),
        (1, 4, 14),
        (3, 6, 10),
        (1, 23, 24),
        (3, 6, 13),
        (4, 20, 21),
        (9, 16, 22),
        (4, 14, 20),
        (15, 18, 21),
        (12, 13, 22),
        (12, 13, 25),
        (11, 14, 19),
        (2, 11, 19),
        (0, 20, 21),
        (2, 5, 12),
        (4, 15, 21),
        (14, 19, 20),
        (7, 16, 19),
        (3, 17, 18),
        (6, 10, 24),
        (8, 18, 21),
        (0, 7, 20),
        (1, 14, 23),
        (5, 12, 25),
        (5, 23, 25),
        (1, 15, 24),
        (2, 5, 11),
        (7, 9, 16),
        (3, 10, 18),
        (7, 19, 20),
        (9, 17, 22),
        (5, 11, 23),
        (23, 24, 25),
        (0, 8, 21),
        (1, 4, 15),
    ]

    def __init__(self, n_CO: int, initial_candidates=[]) -> None:
        """n_CO the number of coadsorbed CO"""
        self.n_CO = n_CO
        self.avail_sites = [*self.avail_top, *self.avail_hollow, *self.avail_bridge]
        self.n_sites = len(self.avail_sites)
        self.tested_sites = []
        super().__init__(initial_candidates)

    def add_ID(self, ID):
        if not (ID in self.tested_sites):
            self.tested_sites.append(ID)

    def get_next_ID(self):
        for _ in range(5000000):  # try 5 million times to generate a new ID
            idxs = random.sample(range(0, self.n_sites), self.n_CO)
            next_ID = set([self.avail_sites[i] for i in idxs])
            if next_ID not in self.tested_sites and self.check_ID(next_ID):
                return next_ID  # return ID if not tested
        return None

    def get_additional_atoms(self, ID, template):
        ID_list = list(ID)
        assert len(ID_list) == self.n_CO

        additional_atoms = []

        for single_ID in ID_list:
            bond_type = len(single_ID)  # 1 - top, 2 - bridge, 3 - hollow
            if bond_type == 1:
                additional_atoms += get_top_CO(template, single_ID[0])
            if bond_type == 2:
                additional_atoms += get_bridge_CO(template, single_ID[0], single_ID[1])
            if bond_type == 3:
                additional_atoms += get_hollow_CO(
                    template, single_ID[0], single_ID[1], single_ID[2]
                )

        return Atoms(additional_atoms)

    def check_ID(self, ID):
        """Helper function to check that the configuration does not provide a
        very unreasonable configuration.

        Rules:
        - Bridge and Hollow CO cannot share two Ni atoms
        - Top cannot share atom with bridge CO
        - Top cannot share Ni atoms with two other CO's
        """

        # Check Bridge and hollow overlap and top and other
        for (id1, id2) in combinations(ID, 2):
            if len(id1) == 2 and len(id2) == 3:
                if id1[0] in id2 and id1[1] in id2:
                    return False
            elif len(id1) == 3 and len(id2) == 2:
                if id2[0] in id1 and id2[1] in id1:
                    return False
            elif len(id1) == 1 and len(id2) == 2:
                if id1[0] in id2:
                    return False
            elif len(id1) == 2 and len(id2) == 1:
                if id2[0] in id1:
                    return False
        # Check Top and two other overlap
        for (id1, id2, id3) in combinations(ID, 3):
            top = None
            if (len(id1) == 1) and (len(id2) > 1) and (len(id3) > 1):
                top, other1, other2 = id1, id2, id3
            if (len(id1) > 1) and (len(id2) == 1) and (len(id3) > 1):
                top, other1, other2 = id2, id1, id3
            if (len(id1) > 1) and (len(id2) > 1) and (len(id3) == 1):
                top, other1, other2 = id3, id2, id1
            if top == None:
                continue

            if (top[0] in other1) and (top[0] in other2):
                return False

        return True


def rs_test():
    import matplotlib

    matplotlib.use("Agg")

    import numpy as np
    from agox import AGOX
    from agox.databases import Database
    from agox.environments import Environment
    from agox.evaluators import LocalOptimizationEvaluator
    from agox.generators import RandomGenerator
    from agox.postprocessors import CenteringPostProcess

    from ase import Atoms

    # Calculator
    from ase.calculators.emt import EMT

    calc = EMT()

    # System & general settings:

    # Confine to one side of structure since it is symmetric
    confinement_cell = np.eye(3) * 16
    confinement_corner = np.array([0, 0, 0])

    template = read("code\\coverage\\3CO\\3COgofee\\job\\POSCAR")
    environment = Environment(
        template=Atoms("", cell=np.eye(3) * 16, pbc=True),
        # template=template,
        symbols="Ni30C2O2",
        confinement_cell=confinement_cell,
        confinement_corner=confinement_corner,
    )

    nhandler = nCOHandler(2)
    generator = ManualGenerator(
        nhandler,
        # template=template,
        **environment.get_confinement(),
        environment=environment,
        order=1
    )
    observer = ManualObserver(nhandler, order=4)

    # Database
    db_path = "db_test.db"
    database = Database(filename=db_path, order=7, write_frequency=1)

    evaluator = LocalOptimizationEvaluator(
        calc,
        gets={"get_key": "candidates"},
        optimizer_kwargs={"logfile": None},
        store_trajectory=False,
        optimizer_run_kwargs={"fmax": 0.01, "steps": 5},
        order=3,
    )

    # Run the program
    agox = AGOX(generator, database, evaluator, observer, seed=0)
    agox.run(N_iterations=20)
    view(database.get_all_candidates(), block=True)


def gofee_test():
    import numpy as np

    #### settings ####
    # General
    output_database_path = "db.db"
    template_path = "code\\coverage\\3CO\\3COgofee\\job\\POSCAR"  # read by ase.io.read
    n_CO = 4

    # gofee
    symbols_in_search = "Ni30C4O4"
    iteration_start_training = 10  # when GPR starts to trian
    update_interval = 10  # GPR update interval
    start_relax = iteration_start_training + 1
    sample_size = 15
    num_candidates = {
        0: [15, 0],  # nCOhandler  # Rattlegenerator
        10: [10, 5],
    }  # distribution between candidates
    kappa = 2
    f_max_dft = 0.01
    n_steps_dft = 1
    f_max_relax_post = 0.01
    n_steps_relax_post = 400
    N_iterations = 10
    seed = 0
    confinement_cell = np.eye(3) * 16
    confinement_corner = np.array([0, 0, 0])
    write_frequency = 1  # frequency of database

    # VASP
    nelm = 80

    import matplotlib

    matplotlib.use("Agg")

    from agox import AGOX
    from agox.databases import Database
    from agox.environments import Environment
    from agox.evaluators import LocalOptimizationEvaluator
    from agox.generators import RandomGenerator, RattleGenerator
    from agox.samplers import KMeansSampler
    from agox.models import ModelGPR
    from agox.acquisitors import LowerConfidenceBoundAcquisitor
    from agox.postprocessors import RelaxPostprocess, CenteringPostProcess
    from agox.postprocessors import ParallelRelaxPostprocess
    from agox.evaluators import SinglePointEvaluator
    from agox.collectors import StandardCollector
    from agox.collectors import ParallelCollector
    from agox.candidates import StandardCandidate
    import os

    from ase.io import read
    from ase import Atoms

    # personaly made functions to interact with vaspfiles
    from utils.IOdatabase import IODataBase
    from utils.handlers.COHandler import nCOHandler
    from utils.manual_generator import ManualGenerator, ManualObserver

    # setup VASP calculator
    from ase.calculators.emt import EMT

    calc = EMT()

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
    relaxer = RelaxPostprocess(
        model=acquisitor.get_acquisition_calculator(),
        optimizer_run_kwargs={"fmax": f_max_relax_post, "steps": n_steps_relax_post},
        constraints=environment.get_constraints(),
        order=4,
        start_relax=start_relax,
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


if __name__ == "__main__":
    gofee_test()
    # rs_test()
    # from agox.generators import RandomGenerator
    # from agox.environments import Environment
    # from ase.io import read, write
    # import numpy as np
    # from utils.manual_generator import ManualGenerator, ManualObserver

    # # from utils.handlers.COHandler import nCOHandler

    # template = read("code\\coverage\\3CO\\3COgofee\\job\\POSCAR")
    # [template.pop() for _ in range(30)]
    # # Confine to one side of structure since it is symmetric
    # confinement_cell = np.array([[16, 0, 0], [0, 16, 0], [0, 0, 16]])
    # confinement_corner = np.array([0, 0, 0])

    # environment = Environment(
    #     template=Atoms(pbc=True, cell=confinement_cell),
    #     symbols="Ni30C2O2",
    #     confinement_cell=confinement_cell,
    #     confinement_corner=confinement_corner,
    # )

    # nhandler = nCOHandler(2)
    # template = read("code\\coverage\\3CO\\3COgofee\\job\\POSCAR")
    # generator = ManualGenerator(
    #     nhandler, template=template, **environment.get_confinement()
    # )
    # observer = ManualObserver(nhandler)

    # new_candidates = []
    # while len(new_candidates) < 20:
    #     candidate = generator.get_candidates(None, environment)[0]
    #     if candidate != None:
    #         # Check that the adddition went well
    #         new_candidates.append(candidate)
    # view(new_candidates, block=True)
    pass

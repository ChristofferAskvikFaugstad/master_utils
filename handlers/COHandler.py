from utils.manual_generator import *
from utils.manual_placement import *
from agox.candidates import StandardCandidate
# The template
class COHandler(IDHandler):
    def __init__(self, initial_candidates=...) -> None:
        super().__init__(initial_candidates)
    
    def get_ID(self, candidate):
        return super().get_ID(candidate)
    
    def add_ID(self, ID):
        return super().add_ID(ID)
    
    def get_next_ID(self):
        return super().get_next_ID()
    
    def get_additional_atoms(self, ID, template):
        return super().get_additional_atoms(ID, template)




#
class HTop(IDHandler):
    maxID = 26
    
    def __init__(self, initial_candidates=[]) -> None:
        self.IDs = set()
        super().__init__(initial_candidates)

    
    
    def get_next_ID(self):
        all_IDs = set([i for i in range(self.maxID)])
        possible_IDs = list(all_IDs.difference(self.IDs))
        next_ID = possible_IDs[np.random.randint(len(possible_IDs))]
        return next_ID

    def add_ID(self, ID):
        self.IDs.add(ID)
    
    
    def get_additional_atoms(self, ID, template):
        add_atom = get_top_atom(template, ID, "H", 2)
        return add_atom


def test_in_RS():
    import matplotlib
    matplotlib.use('Agg')
    import numpy as np
    from agox import AGOX
    from agox.databases import Database
    from agox.environments import Environment
    from agox.evaluators import LocalOptimizationEvaluator
    from agox.generators import RandomGenerator
    from agox.postprocessors import CenteringPostProcess
    from utils.data_handler import get_Ni30_template
    from ase import Atoms
    ##############################################################################
    # Calculator
    ##############################################################################

    from ase.calculators.emt import EMT

    calc = EMT()

    ##############################################################################    
    # System & general settings:
    ##############################################################################
        
    template = get_Ni30_template()
    confinement_cell = np.eye(3) * 16
    confinement_corner = np.array([0, 0, 0])
    environment = Environment(template=template, symbols='H', 
        confinement_cell=confinement_cell, confinement_corner=confinement_corner)

    # Database
    db_path = 'db{}.db'.format(0) # From input argument!
    database = Database(filename=db_path, order=7, write_frequency=1, initialize= True)

    IDs = [0,1,2]

    htop = HTop()
    random_generator = ManualGenerator(htop,environment=environment,
        **environment.get_confinement(),
        sets =  {"set_key" : "candidates"},
        order = 1)

    observer = ManualObserver(htop, order = 6)


    evaluator = LocalOptimizationEvaluator(calc,
        gets =  {"get_key" : "candidates"},
        optimizer_kwargs={'logfile':None}, store_trajectory=False,
        optimizer_run_kwargs={'fmax':0.0001, 'steps':20}, order=5)



    agox = AGOX(random_generator, database, evaluator, observer, seed = 0)
    agox.run(N_iterations = 10)

    from IOdatabase import IODataBase
    iodb = IODataBase(db_path)
    view(iodb.get_all_candidates(), block = True)

def test_in_GOFEE():
    import matplotlib
    matplotlib.use('Agg')

    import numpy as np
    from agox import AGOX
    from agox.databases import Database
    from agox.environments import Environment
    from agox.evaluators import LocalOptimizationEvaluator
    from agox.generators import RandomGenerator, RattleGenerator
    from agox.samplers import KMeansSampler
    from agox.models import ModelGPR
    from agox.acquisitors import LowerConfidenceBoundAcquisitor
    from agox.postprocessors import RelaxPostprocess, CenteringPostProcess
    from agox.collectors import StandardCollector
    from data_handler import get_Ni30_template

    from ase.io.vasp import read_vasp_xml
    from ase import Atoms

    # setup VASP calculator
    from ase.calculators.emt import EMT
    calc = EMT()
    
    template = get_Ni30_template()
    confinement_cell = np.eye(3) * 16
    confinement_corner = np.array([0, 0, 0])
    environment = Environment(template=template, symbols='H', 
        confinement_cell=confinement_cell, confinement_corner=confinement_corner)

    # Database
    db_path = 'db{}.db'.format(0) # From input argument!
    database = Database(filename=db_path, order=8, write_frequency=1, initialize= True)

    IDs = [0,1,2]
    
    htop = HTop()
    # random_generator = ManualGenerator(htop,environment=environment,
    #     **environment.get_confinement(),
    #     sets =  {"set_key" : "candidates"},
    #     order = 1)

    observer = ManualObserver(htop, order = 7)

    # initialize default Gaussian Process Regression model
    model = ModelGPR.default(environment, database)

    # setup generator through a sampler
    sample_size = 20
    sampler = KMeansSampler(feature_calculator=model.get_feature_calculator(), 
        database=database, sample_size=sample_size, order=1)

    rattle_generator = RattleGenerator(use_mic = False, **environment.get_confinement())
    # random_generator = ManualGenerator(htop, use_mic = False, **environment.get_confinement())
    random_generator = RandomGenerator(use_mic = False, **environment.get_confinement())
    generators = [random_generator, rattle_generator]
    num_candidates = {0:[20, 0], 5 : [10,10]} # distribution between candidates

    # controls the sampling from both generators
    collector = StandardCollector(generators=generators, sampler=sampler,
        environment=environment, num_candidates=num_candidates, order = 2)

    # Acuistor function initialization
    acquisitor = LowerConfidenceBoundAcquisitor(model_calculator=model, 
        kappa=2, order=4)

    # Relaxing in surrogate potential
    relaxer = RelaxPostprocess(model=acquisitor.get_acquisition_calculator(), 
        constraints=environment.get_constraints(), order = 3, start_relax=0)    

    # VASP evaluation and one step
    evaluator = LocalOptimizationEvaluator(calc, 
        gets={'get_key':'prioritized_candidates'}, 
        optimizer_kwargs={'logfile':None}, store_trajectory=True,
        optimizer_run_kwargs={'fmax':0.01, 'steps':2}, order=6)


    # collect agox modules and run algorithm
    agox = AGOX(collector, acquisitor, relaxer, database, evaluator,observer, seed = 1)
    agox.run(N_iterations=100)

    from IOdatabase import IODataBase
    iodb = IODataBase(db_path)
    view(iodb.get_all_candidates(), block = True)

if __name__ == "__main__":
    from ase.visualize import view
    # from data_handler import get_Ni30_template
    # htop = HTop()
    # template = get_Ni30_template()
    # candidate = template.copy()
    # ID = htop.get_next_ID()
    # print(ID)
    # add_atom = htop.get_additional_atoms(ID, template)
    # candidate.extend(add_atom)
    # view(candidate, block = True)

    # test_in_RS()
    test_in_GOFEE()
    pass
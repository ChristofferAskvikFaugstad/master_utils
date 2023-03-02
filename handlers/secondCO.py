from utils.manual_generator import *
from utils.manual_placement import *
from agox.candidates import StandardCandidate
from utils.data_handler import get_Ni30_template

class secondCOHandler(IDHandler):
    present_CO_ID = (5,25)
    # Sites found using code
    # avail_top = list((i,) for i in np.arange(0,26))
    # avail_bridge = [(0, 7), (0, 8), (0, 20), (0, 21), (1, 4), (1, 14), (1, 15), (1, 23), (1, 24), (2, 5), (2, 11), (2, 12), (2, 16), (2, 19), (3, 6), (3, 10), (3, 13), (3, 17), (3, 18), (4, 14), (4, 15), (4, 20), (4, 21), (5, 11), (5, 12), (5, 23), (5, 25), (6, 10), (6, 13), (6, 24), (6, 25), (7, 9), (7, 16), (7, 19), (7, 20), (8, 9), (8, 17), (8, 18), (8, 21), (9, 16), (9, 17), (9, 22), (10, 15), (10, 18), (10, 24), (11, 14), (11, 19), (11, 23), (12, 13), (12, 22), (12, 25), (13, 22), (13, 25), (14, 19), (14, 20), (14, 23), (15, 18), (15, 21), (15, 24), (16, 19), (16, 22), (17, 18), (17, 22), (18, 21), (19, 20), (20, 21), (23, 24), (23, 25), (24, 25)]
    # avail_hollow = [(6, 13, 25), (6, 24, 25), (8, 17, 18), (11, 14, 23), (8, 9, 17), (2, 16, 19), (10, 15, 18), (10, 15, 24), (1, 4, 14), (3, 6, 10), (1, 23, 24), (3, 6, 13), (4, 20, 21), (9, 16, 22), (4, 14, 20), (15, 18, 21), (12, 13, 22), (12, 13, 25), (11, 14, 19), (2, 11, 19), (0, 20, 21), (2, 5, 12), (4, 15, 21), (14, 19, 20), (7, 16, 19), (3, 17, 18), (6, 10, 24), (8, 18, 21), (0, 7, 20), (1, 14, 23), (5, 12, 25), (5, 23, 25), (1, 15, 24), (2, 5, 11), (7, 9, 16), (3, 10, 18), (7, 19, 20), (9, 17, 22), (5, 11, 23), (23, 24, 25), (0, 8, 21), (1, 4, 15)]
    avail_top = [(1,),(2,)]
    avail_bridge = [(0,7),(0,8)]
    avail_hollow = [(6, 13, 25), (6, 24, 25)]
    
    tested_sites = []
    
    def __init__(self, initial_candidates=[]) -> None:
        self.avail_sites = [*self.avail_top, *self.avail_hollow, *self.avail_bridge]
        super().__init__(initial_candidates)
    
    def add_ID(self, ID):
        self.tested_sites.append(ID)
        self.avail_sites.remove(ID)
    
    def get_next_ID(self):
        return self.avail_sites[np.random.randint(0,len(self.avail_sites))]
    
    def get_additional_atoms(self, ID, template):
        bond_type = len(ID) # 1 - top, 2 - bridge, 3 - hollow
        if bond_type == 1:
            return get_top_CO(template,ID[0])
        if bond_type == 2:
            return get_bridge_CO(template,ID[0],ID[1])
        if bond_type == 3:
            return get_hollow_CO(template,ID[0],ID[1],ID[2])
        return None



def test():
    import matplotlib
    matplotlib.use('Agg')
    import numpy as np
    from agox import AGOX
    from agox.databases import Database
    from agox.environments import Environment
    from agox.evaluators import LocalOptimizationEvaluator
    from agox.generators import RandomGenerator
    from agox.postprocessors import CenteringPostProcess
    from utils.data_handler import get_Ni30_template, get_structure_relpath
    from ase import Atoms
    ##############################################################################
    # Calculator
    ##############################################################################

    from ase.calculators.emt import EMT

    calc = EMT()

    ##############################################################################    
    # System & general settings:
    ##############################################################################
        
    # template = get_Ni30_template()
    template = get_structure_relpath("single/ni30_COs/CO5-25")
    confinement_cell = np.eye(3) * 16
    confinement_corner = np.array([0, 0, 0])
    environment = Environment(template=template, symbols='CO', 
        confinement_cell=confinement_cell, confinement_corner=confinement_corner)

    # Database
    db_path = 'db{}.db'.format(0) # From input argument!
    database = Database(filename=db_path, order=7, write_frequency=1, initialize= True)


    htop = secondCOHandler()
    random_generator = ManualGenerator(htop,environment=environment,
        **environment.get_confinement(),
        sets =  {"set_key" : "candidates"},
        order = 1)

    observer = ManualObserver(htop, order = 6)


    evaluator = LocalOptimizationEvaluator(calc,
        gets =  {"get_key" : "candidates"},
        optimizer_kwargs={'logfile':None}, store_trajectory=True,
        optimizer_run_kwargs={'fmax':0.0001, 'steps':1}, order=5)
    
    from agox.evaluators.single_point import SinglePointEvaluator
    evaluator = SinglePointEvaluator(calc,
        gets =  {"get_key" : "candidates"},
        sets = {"set_key" : "evaluated_candidates"}
        )



    agox = AGOX(random_generator, database, evaluator, observer, seed = 1)
    # agox = AGOX(random_generator, database, observer, seed = 0)
    agox.run(N_iterations = 6)
    candidates = database.get_all_candidates()
    view(candidates, block = True)

if __name__ == "__main__":
    test()
    # template = get_Ni30_template()
    # add_hollow_CO(template, 6, 13, 25)
    # avail_hollow = set()
    # bonds = [bond for bond in get_bonds(template) if max(bond) < 26 and bond != (2,25)]

    # # for i, bond_i in enumerate(bonds):
    # #     for j, bond_j in enumerate(bonds):
    # #         for k, bond_k in enumerate(bonds):
    # #             if i == j or j == k or i == k:
    # #                 continue
    # #             # check if all three bonds match
    # #             bond_atoms = tuple(sorted(list(set((*bond_i,*bond_j, *bond_k)))))
    # #             if len(bond_atoms) == 3 and max(bond_atoms) < 26:
    # #                 avail_hollow.add(bond_atoms)
    # # avail_hollow = list(avail_hollow)
    # print(avail_hollow)
    # for i,j in bonds:
    #     add_bridge_atom(template, i, j, "H", 2)
    # view(template, block = True)



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
    avail_top = list((i,) for i in np.arange(0,26))
    avail_bridge = [(0, 7), (0, 8), (0, 20), (0, 21), (1, 4), (1, 14), (1, 15), (1, 23), (1, 24), (2, 5), (2, 11), (2, 12), (2, 16), (2, 19), (3, 6), (3, 10), (3, 13), (3, 17), (3, 18), (4, 14), (4, 15), (4, 20), (4, 21), (5, 11), (5, 12), (5, 23), (6, 10), (6, 13), (6, 24), (6, 25), (7, 9), (7, 16), (7, 19), (7, 20), (8, 9), (8, 17), (8, 18), (8, 21), (9, 16), (9, 17), (9, 22), (10, 15), (10, 18), (10, 24), (11, 14), (11, 19), (11, 23), (12, 13), (12, 22), (12, 25), (13, 22), (13, 25), (14, 19), (14, 20), (14, 23), (15, 18), (15, 21), (15, 24), (16, 19), (16, 22), (17, 18), (17, 22), (18, 21), (19, 20), (20, 21), (23, 24), (23, 25), (24, 25)]
    avail_hollow = [(6, 13, 25), (6, 24, 25), (8, 17, 18), (11, 14, 23), (8, 9, 17), (2, 16, 19), (10, 15, 18), (10, 15, 24), (1, 4, 14), (3, 6, 10), (1, 23, 24), (3, 6, 13), (4, 20, 21), (9, 16, 22), (4, 14, 20), (15, 18, 21), (12, 13, 22), (12, 13, 25), (11, 14, 19), (2, 11, 19), (0, 20, 21), (2, 5, 12), (4, 15, 21), (14, 19, 20), (7, 16, 19), (3, 17, 18), (6, 10, 24), (8, 18, 21), (0, 7, 20), (1, 14, 23), (5, 12, 25), (5, 23, 25), (1, 15, 24), (2, 5, 11), (7, 9, 16), (3, 10, 18), (7, 19, 20), (9, 17, 22), (5, 11, 23), (23, 24, 25), (0, 8, 21), (1, 4, 15)]
    
    
    def __init__(self, n_CO : int, initial_candidates = []) -> None:
        """ n_CO the number of coadsorbed CO"""
        self.n_CO = n_CO
        self.avail_sites = [*self.avail_top, *self.avail_hollow, *self.avail_bridge]
        self.n_sites = len(self.avail_sites)
        self.tested_sites = []
        super().__init__(initial_candidates)
    
    def add_ID(self, ID):
        if not (ID in self.tested_sites):
            self.tested_sites.append(ID)
    
    def get_next_ID(self):
        for _ in range(5000000): # try 5 million times to generate a new ID
            idxs = random.sample(range(0, self.n_sites), self.n_CO)
            next_ID = set([self.avail_sites[i] for i in idxs])
            if next_ID not in self.tested_sites and self.check_ID(next_ID):
                return next_ID # return ID if not tested
        return None
    
    def get_additional_atoms(self, ID, template):
        ID_list = list(ID)
        assert len(ID_list) == self.n_CO
        
        additional_atoms = []
        
        for single_ID in ID_list:  
            bond_type = len(single_ID) # 1 - top, 2 - bridge, 3 - hollow
            if bond_type == 1:
                additional_atoms += get_top_CO(template,single_ID[0])
            if bond_type == 2:
                additional_atoms += get_bridge_CO(template,single_ID[0],single_ID[1])
            if bond_type == 3:
                additional_atoms += get_hollow_CO(template,single_ID[0],single_ID[1],single_ID[2])
        
        return additional_atoms
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
            if (len(id1) == 2 and len(id2) == 3):
                if id1[0] in id2 and id1[1] in id2:
                    return False
            elif (len(id1) == 3 and len(id2) == 2):
                if id2[0] in id1 and id2[1] in id1:
                    return False
            elif (len(id1) == 1 and len(id2) == 2):
                if id1[0] in id2:
                    return False
            elif (len(id1) == 2 and len(id2) == 1):
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


if __name__ == "__main__":
   pass

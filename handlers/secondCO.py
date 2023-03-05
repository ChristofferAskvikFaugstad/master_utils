from utils.manual_generator import *
from utils.manual_placement import *
from agox.candidates import StandardCandidate

class SecondCOHandler(IDHandler):
    present_CO_ID = (5,25)
    # Sites found using code
    avail_top = list((i,) for i in np.arange(0,26))
    avail_bridge = [(0, 7), (0, 8), (0, 20), (0, 21), (1, 4), (1, 14), (1, 15), (1, 23), (1, 24), (2, 5), (2, 11), (2, 12), (2, 16), (2, 19), (3, 6), (3, 10), (3, 13), (3, 17), (3, 18), (4, 14), (4, 15), (4, 20), (4, 21), (5, 11), (5, 12), (5, 23), (5, 25), (6, 10), (6, 13), (6, 24), (6, 25), (7, 9), (7, 16), (7, 19), (7, 20), (8, 9), (8, 17), (8, 18), (8, 21), (9, 16), (9, 17), (9, 22), (10, 15), (10, 18), (10, 24), (11, 14), (11, 19), (11, 23), (12, 13), (12, 22), (12, 25), (13, 22), (13, 25), (14, 19), (14, 20), (14, 23), (15, 18), (15, 21), (15, 24), (16, 19), (16, 22), (17, 18), (17, 22), (18, 21), (19, 20), (20, 21), (23, 24), (23, 25), (24, 25)]
    avail_hollow = [(6, 13, 25), (6, 24, 25), (8, 17, 18), (11, 14, 23), (8, 9, 17), (2, 16, 19), (10, 15, 18), (10, 15, 24), (1, 4, 14), (3, 6, 10), (1, 23, 24), (3, 6, 13), (4, 20, 21), (9, 16, 22), (4, 14, 20), (15, 18, 21), (12, 13, 22), (12, 13, 25), (11, 14, 19), (2, 11, 19), (0, 20, 21), (2, 5, 12), (4, 15, 21), (14, 19, 20), (7, 16, 19), (3, 17, 18), (6, 10, 24), (8, 18, 21), (0, 7, 20), (1, 14, 23), (5, 12, 25), (5, 23, 25), (1, 15, 24), (2, 5, 11), (7, 9, 16), (3, 10, 18), (7, 19, 20), (9, 17, 22), (5, 11, 23), (23, 24, 25), (0, 8, 21), (1, 4, 15)]
    
    
    def __init__(self, initial_candidates=[]) -> None:
        self.avail_sites = [*self.avail_top, *self.avail_hollow, *self.avail_bridge]
        self.tested_sites = []
        super().__init__(initial_candidates)
    
    def add_ID(self, ID):
        if not (ID in self.tested_sites):
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


if __name__ == "__main__":
   pass


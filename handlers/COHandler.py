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





if __name__ == "__main__":
    pass
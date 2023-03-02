from ase.atoms import Atoms
from ase.io.pov import get_bondpairs

class COBond:
    def __init__(self, atoms : Atoms,C_index = 30,O_index = 31) -> None:
        self.C_pairs = []
        self.O_pairs = []
        self.C_index = C_index
        self.O_index = O_index

        for pairs in get_bondpairs(atoms)[:-1]:
            if pairs[0] == self.C_index:
                self.C_pairs.append(pairs[1])
            if pairs[0] == self.O_index:
                self.O_pairs.append(pairs[1])
            if pairs[1] == self.C_index:
                self.C_pairs.append(pairs[0])
            if pairs[1] == self.O_index:
                self.O_pairs.append(pairs[0])
    def __eq__(self, __o: object) -> bool:
        return (self.C_pairs == __o.C_pairs) and (self.O_pairs == __o.O_pairs)

class COSearch:
    def __init__(self, **kwargs) -> None:
        self.COBondsDict = kwargs
        self.bonds = []
    def add(self, atoms) -> bool:
        bond = COBond(atoms, **self.COBondsDict)
        accept = True
        for cur_bond in self.bonds:
            if cur_bond == bond:
                accept = False
        accept = self.filter(bond) and  accept

        if accept:
            self.bonds.append(bond)
        return accept
    def get_length(self):
        return len(self.bonds)

    def filter(self, bond : COBond) -> bool:
        if len(bond.C_pairs) + len(bond.O_pairs) < 2:
            return False
        return True
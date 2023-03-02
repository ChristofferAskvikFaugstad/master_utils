from agox.databases import Database
import numpy as np

class IODataBase(Database):
    """
    Ownmade class for database input output
    """

    def __init__(self, filename='db.db',outfile = None, initialize=False, verbose=False, write_frequency=50, call_initialize=True, store_meta_information=True, **kwargs):
        super().__init__(filename, initialize, verbose, write_frequency, call_initialize, store_meta_information, **kwargs)
        energies = self.get_all_energies()
        self.size = len(energies)
        self.all_structures = []
        self.forces = []
        self.forcesmax = np.zeros(self.size)
        all_structure_data = self.get_all_structures_data()
        for i in range(self.size):
            current_data = all_structure_data[i]
            self.all_structures.append(self.db_to_atoms(current_data))
            self.forces.append(current_data["forces"])
            self.forcesmax[i] = np.max(np.linalg.norm(self.all_structures[i].get_forces(), axis = 1))

        if outfile != None:
            self.parse_outfile_for_predictions(outfile)
        else:
            self.predictions = []


    def get_all_candidates(self):
        return self.all_structures

    def get_forces(self):
        return self.forces
    
    def get_forcesmax(self):
        return self.forcesmax

    def get_predictions(self):
        return self.predictions

    def parse_outfile_for_predictions(self, outfile: str) -> np.ndarray:
        """parse the outfile AGOX give to obtain the predicted energy of
        the best candidate of the aquistor. Does not incorperate if it is
        not convered in the evaluation. 

        Parameters
        ----------
        outfile : str
            The path to the file

        Returns
        -------
        np.ndarray
            Array with the predicted energies
        """
        predicted_energies = []
        with open(outfile) as f:
            for line in f:
                if line == "|=============================== LCBAcquisitor ===============================|\n":
                    next_line = f.readline()
                    predicted_energies.append(float(next_line.split(": E=")[1][:8]))
        self.predictions = np.array(predicted_energies)





if __name__ == "__main__":
    db = IODataBase(filename="db_ni8.db",outfile="out.txt" )
    candidates = db.get_all_candidates()




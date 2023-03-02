from agox.generators.ABC_generator import GeneratorBaseClass
from abc import ABC, abstractmethod
from agox.writer import Writer, agox_writer
from agox.observer import Observer
from agox.candidates import StandardCandidate

class IDHandler(ABC):
    IDtag = "generatorID"
    def __init__(self, initial_candidates = []) -> None:
        for candidate in initial_candidates:
            self.append(candidate)

    @abstractmethod
    def add_ID(self, ID):
        pass

    @abstractmethod
    def get_next_ID(self):
        pass
    @abstractmethod
    def get_additional_atoms(self, ID, template):
        pass
    
    def get_ID(self, candidate : StandardCandidate):
        return candidate.meta_information[self.IDtag]
    
    def append(self,candidate):
        ID = self.get_ID(candidate)
        self.add_ID(ID)

class ManualGenerator(GeneratorBaseClass):
    name = "ManualGenerator"
    generator_tries = 100
    def __init__(self, id_handler : IDHandler, **kwargs):
        self.id_handler = id_handler
        super().__init__(**kwargs)

    def get_candidates(self, sampler, environment):
        template = environment.get_template()
        candidate : StandardCandidate = template.copy()
        build_succesful = False 
        for _ in range(self.generator_tries):
            ID = self.id_handler.get_next_ID()
            add_atoms = self.id_handler.get_additional_atoms(ID, template)
            if add_atoms != None:
                candidate.extend(add_atoms)
                build_succesful = True
                break
        
        if not build_succesful:
            self.writer('Start generator failing at producing valid structure')
            return [None]
    
        candidate = self.convert_to_candidate_object(candidate, template)
        candidate.add_meta_information('description', self.name)
        candidate.add_meta_information(self.id_handler.IDtag, ID)

        return [candidate]

class ManualObserver(Observer, Writer):

    name = 'ManualObserver'

    def __init__(self, id_handler : IDHandler, order=2, gets={'get_key':'evaluated_candidates'}, 
        sets={'set_key':'evaluated_candidates'}):
        Observer.__init__(self, gets=gets, order=order, sets=sets)
        Writer.__init__(self)
        self.add_observer_method(self.basic_observer_method, order=self.order[0], 
            sets=self.sets[0], gets=self.gets[0])

        self.id_handler = id_handler

    @agox_writer
    @Observer.observer_method
    def basic_observer_method(self, state):
        self.writer(f'Iteration: {self.get_iteration_counter()}: ID Obaserver')

        candidates = state.get_from_cache(self, self.get_key)
        try:
            self.id_handler.append(candidates[-1])
        except:
            pass



        state.add_to_cache(self, self.set_key, candidates, mode='w')


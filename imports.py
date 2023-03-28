from agox.databases.database import Database
from ase.visualize import view
from ase.io.trajectory import Trajectory
from ase.neb import NEBTools
from ase.neb import NEB
from ase.io import read,write
from ase.io.vasp import *
from ase.ga.ofp_comparator import OFPComparator
from ase.visualize.plot import plot_atoms
from ase.atoms import Atoms
from agox.candidates import StandardCandidate
from ase.calculators.vasp import Vasp
from ase.calculators.emt import EMT
from ase.neighborlist import build_neighbor_list
from ase.utils.forcecurve import fit_images
from ase.build import molecule
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import (OffsetImage, AnnotationBbox)
import matplotlib.image as image
import os
from typing import List
import matplotlib.pyplot as plt
from typing import List
import dis
from agox.databases import Database
# from py4vasp import Calculation
import pandas as pd
plt.rcParams["figure.figsize"] = (12, 8)
plt.rcParams["font.size"] = 16
plt.rcParams['font.family'] = 'serif'
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
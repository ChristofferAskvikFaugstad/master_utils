from ase.io.vasp import read_vasp_xml
import os

dir = "00" 
initial_xml = "vasprun.xml"
initial = next(read_vasp_xml(os.path.join(dir, initial_xml)))

from ase.calculators.vasp import Vasp
calc = Vasp(prec="Accurate",
images=5,
istart=0,
xc="PBE",
ismear=0,
sigma=0.02,
encut=400,
nelm=60,
lorbit=10,
isym=-1,
ispin=2,
ediffg=-0.05,
ibrion=3,
potim=0,
lcharg=False,
lwave=False,
ncore=32,
nsw=1500,
ichain=0,
iopt=3,
timestep=0.05,
maxmove=0.1,
lclimb=True,
spring=-5,
)
command = calc.make_command(calc.command)
properties = ("energy",)
from ase.calculators import calculator

system_changes = tuple(calculator.all_changes)
calc.write_input(initial, properties, system_changes)
error = calc._run(command=command, directory=calc.directory)
print(error)
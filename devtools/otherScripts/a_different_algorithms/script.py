"""get paths for test data"""
from restraintmaker import test

pdb_path = test.test_system1_pdb
print(pdb_path)

from pymol import cmd
from restraintmaker.interface_Pymol import Pymol_Utitlities as u

cmd.load(pdb_path)
atom_list = u.pymol_selection_to_atom_list("all")
cmd.reinitialize()



from restraintmaker.algorithm import Filter, Optimizer
from restraintmaker.utils import Utilities as u

filtered_atoms = atom_list

## Filtering for not H
filter = Filter.ElementFilter(filtered_atoms)
filter.get_args(lambda x: ("O, N, C"))
filtered_atoms = filter.filter()
print("Filter Step1 - Elements: ", len(atom_list), "->",len(filtered_atoms))

## Filtering for rings
filter = Filter.RingFilter(filtered_atoms)
pdb_blocks = u.convert_atoms_to_pdb_molecules(filtered_atoms)
filter.get_args(lambda x: (pdb_blocks))
filtered_atoms = filter.filter()
print("Filter Step2 - Ring Atoms: ", len(atom_list), "->",len(filtered_atoms))

# Open question
number_of_restraints= 6

#Params:
methods = ["shortest"]
projections = ["pca_2d"]
cutoff_distance = 1.2


#Optimize Restraints
# Todo: fix this call
optimizer = Optimizer.TreeHeuristicOptimizer(filtered_atoms)
optimizer.get_args(lambda x: (number_of_restraints, cutoff_distance, 'shortest', 'pca_2d'))
found_restraints = optimizer.make_restraints()

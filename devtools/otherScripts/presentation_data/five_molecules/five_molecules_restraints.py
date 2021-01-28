
import glob, os, time

import pymol
from pymol import cmd
#check pdbs
molecule_dir= "."
pdbs = glob.glob(molecule_dir+"/*pdb")


#fire up pymol
pymol.finish_launching()

#load data
[cmd.load(pdb) for pdb in pdbs]
time.sleep(1)
cmd.center()
cmd.do("set_view (\
     0.999631941,    0.027128441,    0.000000000,\
    -0.027128441,    0.999631941,    0.000000000,\
     0.000000000,    0.000000000,    1.000000000,\
     0.000000000,    0.000000000,  -26.657482147,\
     0.228777885,    0.000000000,    0.001666667,\
    22.597614288,   30.717350006,  -20.000000000 )")

#Select all atoms
from restraintmaker.interface_Pymol.pymol_utilities import pymol_utitlities as u
atom_list = u.pymol_selection_to_atom_list("all")

colors = ["marine", "firebrick", "yelloworange", "hotpink", "forest"]
for ind in range(5):
    molecule2 = cmd.create("mol_"+str(ind), "resi "+str(ind))
    cmd.color(colors[ind], "mol_"+str(ind))
cmd.delete("five*")
cmd.set("grid_mode", True)

from restraintmaker.algorithm import Filter, Optimizer
from restraintmaker.utils import Utilities as u

#settings
cmd.set("sphere_scale", 0.15)

## Filtering for rings
RingFilter = Filter.RingFilter(atom_list)
pdb_blocks = u.convert_atoms_to_pdb_molecules(atom_list)
RingFilter.get_args(lambda x: (pdb_blocks))
filtered_atoms = RingFilter.filter()
print("Filter Step2 - Ring Atoms: ", len(atom_list), "->",len(filtered_atoms))
cmd.color("grey", "not id "+"+".join([str(atom.id) for atom in filtered_atoms]))

cmd.orient()
#cmd.ray(1800, 1800)
#cmd.save("molecules_selection.png")

cmd.set("grid_mode", False)
time.sleep(2)
cmd.orient()
cmd.center()
cmd.ray(1800)
cmd.save("molecules_selection_overlay.png")
exit()

#Optimizers
from restraintmaker.algorithm import Optimizer
distance_treshold = 1.5
opt = Optimizer.TreeHeuristicOptimizer(filtered_atoms)
opt.get_args(lambda x: (4, distance_treshold, 'shortest', None))
res = opt.make_restraints()

"""
Optimizer.compare_pair_optimizers(criterion=Optimizer._calculate_value_unscaled_pca_2d,
                                  atoms=filtered_atoms,
                                  opt_types=[Optimizer.TreeHeuristicOptimizer, Optimizer.TreeHeuristicOptimizer,
                                             Optimizer.TreeHeuristicOptimizer, Optimizer.TreeHeuristicOptimizer,
                                             Optimizer.BestMoleculeRingOptimizer, Optimizer.BruteForceRingOptimzer,
                                             Optimizer.BruteForceRingOptimzer],
                                  opt_args=[(4, 1.2, 'prim', None), (4, 1.2, 'cog', None),
                                            (4, 1.2, 'shortest', None), (4, 1.2, 'biased_avg', None),
                                            (4, 1.2, 'pca', 'pca_2d'), (4, 1.2, 'pca', 'None'),
                                            (4, 1.2, 'convex_hull', 'None')],
                                  out_path="tmp_res",
                                  new_dir_name="optOut")  # Name of the Ligands without path and .pdb
"""
cmd.disable()

i=200
obj="selected_res"
restraint_name = "selected_res"
sel_restraint_node = []
for restraint in res:
    atom1ID = restraint._atomA.id
    atom2ID = restraint._atomB.id

    atom1 = list(filter(lambda x: atom1ID == x.id, filtered_atoms))[0]
    atom2 = list(filter(lambda x: atom2ID == x.id, filtered_atoms))[0]

    x, y, z = ((atom1.x+atom2.x)/2, (atom1.y+atom2.y)/2, (atom1.z+atom2.z)/2)
    sres = u.Atom(elem="sres", id=str(i), name=restraint_name, x=x, y=y , z=z , chain="", resn=restraint_name, resi=str(i), alt="",b=0, label=0 )
    cmd.pseudoatom(object=obj, name=sres.name, resi=sres.id, resn=sres.resn, pos=(sres.x,sres.y,sres.z))
    sel_restraint_node.append(sres)
    i+=1

cmd.show("spheres")
cmd.hide("nonbonded")

cmd.enable("res")
cmd.color("cyan", obj)
cmd.set("sphere_scale", 0.35, obj)


cmd.disable()
cmd.enable(obj)
cmd.enable("mol_*")
atom_ids = []
resns = []
for rIND, restraint in enumerate(res):
    atom1 = restraint._atomA
    atom2 = restraint._atomB
    cmd.dist(selection1="resn "+atom1.resn+" and id "+str(atom1.id), selection2="resn "+atom2.resn+" and id "+str(atom2.id), name="selAtom")
    atom_ids.extend([atom1.id, atom2.id])

#cmd.group("sel_dist", "selAtom*")

cmd.set("sphere_scale", 0.2, obj)
cmd.set("dash_radius", 0.1)
cmd.color("lightblue", "id "+"+".join(map(str, atom_ids)))

cmd.disable("sel*")
#cmd.ray(1800, 1800)
#cmd.save("optimized_restraints_atoms_grid.png")

cmd.set("grid_mode", False)
time.sleep(2)
cmd.orient()
cmd.center()
cmd.ray(1800)
cmd.save("optimized_restraints_atoms_overlay.png")
#cmd.set("grid_mode", True)

exit()



cmd.enable("mol_*")
cmd.hide("spheres")
cmd.show("sticks")
cmd.set("dash_radius", 0.15)
#print(len(res))
#print(res)

for restraint in res:
    atom1 = restraint._atomA
    atom2 = restraint._atomB
    cmd.dist("resn "+atom1.resn+" and id "+str(atom1.id), "resn "+atom2.resn+" and id "+str(atom2.id))
    cmd.color("lightblue", "id "+str(atom1.id)+"+"+str(atom2.id)    )
cmd.group("dist", "dist*")

cmd.ray(1200)
cmd.save("optimized_restraints_atoms.png")
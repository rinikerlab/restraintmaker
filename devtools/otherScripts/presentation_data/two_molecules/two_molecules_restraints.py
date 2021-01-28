
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

#Select all atoms
from restraintmaker.interface_Pymol.pymol_utilities import pymol_utitlities as u
atom_list = u.pymol_selection_to_atom_list("all")
molecule1 = cmd.create("mol_1", "resi 0")
molecule2 = cmd.create("mol_2", "resi 1")
cmd.delete("two*")

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

from collections import defaultdict

all_atoms = [atom.id for atom in filtered_atoms]
all_atoms_mol = defaultdict(list)
for atom in filtered_atoms:
    all_atoms_mol[int(atom.resi)].append(atom.id)

mol1 = all_atoms_mol[0]
mol2 = all_atoms_mol[1]

cmd.center()
cmd.do("set_view (\
     0.810585678,   -0.306923389,    0.498748600,\
    -0.582877159,   -0.340507299,    0.737770796,\
    -0.056609299,   -0.888735831,   -0.454910129,\
     0.000000000,    0.000000000,  -25.170881271,\
     0.082305431,    0.027805567,    0.072472334,\
    19.844913483,   30.496849060,  -20.000000000 )")

cmd.color("grey", "all")
cmd.color("firebrick", "id "+"+".join(map(str, mol1)))
cmd.color("forest", "id "+"+".join(map(str, mol2)))


#cmd.ray(1200)
#cmd.save("selected_atoms.png")


#show distance_graph
cmd.hide("sticks", "all")
cmd.show("spheres", "all")

#cmd.ray(1200, async_=False)
#cmd.save("selected_atoms_all_nodes.png")
#exit()

cmd.hide("spheres", "all")
cmd.show("spheres", "id "+"+".join(map(str, mol1)))
cmd.show("spheres", "id "+"+".join(map(str, mol2)))
#cmd.ray(1200, async_=False)
#cmd.save("selected_atoms_nodes.png")


#build up edges
i=0
distance_treshold = 1.5 #A
restraints = []
for atom1 in mol1:
    for atom2 in mol2:
        n="dist"+str(i)
        d = cmd.dist(selection1="id "+str(atom1), selection2="id "+str(atom2), name=n)
        if(d>distance_treshold):
            cmd.delete(n)
        else:
            restraints.append((atom1, atom2))
            i+=1

#cmd.hide("labels")
cmd.set("dash_radius", 0.07)
cmd.set("dash_gap", 0)
cmd.set("dash_color", "orange", "dist*")
cmd.set("label_color", "black")
cmd.set("label_size", 18)

#cmd.ray(1200, async_=False)
#cmd.save("selected_atoms_graph.png")


#place new nodes
i=100
obj="res"
restraint_name = "rest"
new_restraint_node = []
for atom1ID, atom2ID in restraints:
    atom1 = list(filter(lambda x: atom1ID == x.id, filtered_atoms))[0]
    atom2 = list(filter(lambda x: atom2ID == x.id, filtered_atoms))[0]

    x, y, z = ((atom1.x+atom2.x)/2, (atom1.y+atom2.y)/2, (atom1.z+atom2.z)/2)
    tres = u.Atom(elem="res", id=str(i), name=restraint_name, x=x, y=y , z=z , chain="", resn=restraint_name, resi=str(i), alt="",b=0, label=0 )
    cmd.pseudoatom(object=obj, name=tres.name, resi=tres.id, resn=tres.resn, pos=(tres.x,tres.y,tres.z))
    new_restraint_node.append(tres)
    i+=1

cmd.set("sphere_scale", 0.15, "name "+restraint_name)
cmd.group(name="atomDistances", members="dist*")

cmd.hide("nonbonded")
cmd.show("spheres", "resn "+restraint_name)
cmd.color("deepblue", "resn "+restraint_name)
cmd.set("sphere_scale", 0.25, "resn "+restraint_name)

print("Number of possible restraints", len(new_restraint_node))
#cmd.ray(1200)
#cmd.save("restraint_nodes_withAtomdist.png")


cmd.disable("atomDistances")

#cmd.ray(1200)
#cmd.save("restraint_mols_nodes.png")

cmd.disable("mol*")
#cmd.ray(1200)
#cmd.save("restraint_nodes.png")

# Restraint Graph:
name = "r_dist"
i=0
for restraint1Ind, restraint1 in enumerate(new_restraint_node):
    for restraint2 in new_restraint_node[1+restraint1Ind:]:
        cmd.dist( selection1="resi "+str(restraint1.id)+" and resn "+restraint1.name, selection2="resi "+str(restraint2.id)+" and resn "+restraint2.name, name=name+str(i))
        i+=1

cmd.group(name="restDistances", members=name+"*")
cmd.hide("labels")
cmd.set("dash_radius", 0.01)

print("Number of possible edges", i)
#cmd.ray(1200)
#cmd.save("full_restraint_graph.png")

#Optimizers
from restraintmaker.algorithm import Optimizer

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

for rIND, restraint1 in enumerate(sel_restraint_node):
    for restraint2 in sel_restraint_node[rIND+1:]:
        cmd.dist("resi "+str(restraint2.id), "resi "+str(restraint1.id))

cmd.enable("res")
cmd.color("cyan", obj)
cmd.set("sphere_scale", 0.35, obj)
#cmd.ray(1200)
#cmd.save("optimized_restraints.png")

cmd.set("solvent_radius", 2)
cmd.set("surface_quality", 2)

cmd.set("solvent_radius", 2)
cmd.set("transparency", 0.5)
cmd.show("surface", "selected_res")

#cmd.ray(1200)
#cmd.save("convex_hull.png")
#exit()

cmd.disable()
cmd.enable(obj)
cmd.enable("mol_*")
atom_ids = []
for rIND, restraint in enumerate(res):
    atom1 = restraint._atomA
    atom2 = restraint._atomB
    cmd.dist(selection1="resn "+atom1.resn+" and id "+str(atom1.id), selection2="resn "+atom2.resn+" and id "+str(atom2.id), name="selAtom")
    atom_ids.extend([atom1.id, atom2.id])

cmd.set("sphere_scale", 0.2, obj)
cmd.set("dash_radius", 0.1)
cmd.color("lightblue", "id "+"+".join(map(str, atom_ids))+" and (mol_1 or mol_2)")

#cmd.ray(1200)
#cmd.save("optimized_restraints_atoms_overlay.png")

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
    cmd.color("lightblue", "id "+str(atom1.id)+"+"+str(atom2.id))

#cmd.ray(1200)
#cmd.save("optimized_restraints_atoms.png")
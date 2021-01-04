"""
This is a script, that starts a restraint maker session.

"""

import glob, os, time

import pymol
from pymol import cmd

molecule_dir = "ATB_molecules"

#check pdbs
pdbs = glob.glob(molecule_dir+"/*/*pdb")
pdbs = [pdb for pdb in pdbs if (not len(os.path.basename(pdb).split("_")[0]) ==0) ]
obj_names = [os.path.basename(pdb).split("_")[1] if (len(os.path.basename(pdb).split("_")[0]) ==0) else os.path.basename(pdb).split("_")[0] for pdb in pdbs]

#fire up pymol
pymol.finish_launching()

#load data
[cmd.load(pdb, obj_name) for obj_name,pdb in zip(obj_names, pdbs)]
time.sleep(1)

# set nice scene
obj_list = cmd.get_object_list()
for i, obj in enumerate(obj_list):
    cmd.alter(obj, "chain="+str(i))

cmd.save("merged.pdb")

cmd.reinitialize()

cmd.load("merged.pdb")


import restraintmaker
restraintmaker.run_plugin_gui()
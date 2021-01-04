"""
This is a script, that starts a restraint maker session.

"""

import glob

import pymol
from pymol import cmd

#possible paths
molecule_dir = "data"
pdb_path1 = molecule_dir+"/5_long_Ligands.pdb"
pdb_path2 = molecule_dir+"/5_long_Ligands.pdb"
pdb_path3 = molecule_dir+"/5_long_Ligands.pdb"
pdb_path4 = molecule_dir+"/5_long_Ligands.pdb"

pdb_singles_path1 = glob.glob(molecule_dir+"/single_ligs_bad/*.pdb")
pdb_singles_path2 = glob.glob(molecule_dir+"/single_ligs_good/*.pdb")

cmd.load(pdb_path1)


import restraintmaker
restraintmaker.run_plugin_gui()
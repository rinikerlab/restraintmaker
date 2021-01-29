"""
This is a script, that starts a restraint maker session.

"""

import glob

import pymol
from pymol import cmd

#possible paths
molecule_dir = "data"
pdb_path1 = molecule_dir+"/CHK1_5Ligands.pdb"
pdb_path2 = molecule_dir+"/BRD4_7Ligs.pdb"
pdb_path3 = molecule_dir+"/PNMT_9ligs.pdb"
pdb_path4 = molecule_dir+"/21_large_Ligands.pdb"

pdb_singles_path1 = glob.glob(molecule_dir+"/single_ligs_bad/*.pdb")
pdb_singles_path2 = glob.glob(molecule_dir+"/single_ligs_good/*.pdb")

pymol.finish_launching()

cmd.load(pdb_path1)


import restraintmaker
restraintmaker.run_plugin_gui()
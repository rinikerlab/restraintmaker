"""
This is a script, that starts a restraint maker session.

"""

import glob, os

import pymol
from pymol import cmd

#possible paths
molecule_dir = os.getcwd()
pdb_path1 = molecule_dir+"/challenging_5/challenging_5.pdb"
pdb_path2 = molecule_dir+"/easy_6/easy_6.pdb"
pdb_path3 = molecule_dir+"/flat_10/flat_10.pdb"


pymol.finish_launching()

cmd.load(pdb_path3)


import restraintmaker
restraintmaker.run_plugin_gui()
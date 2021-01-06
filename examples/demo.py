#!/usr/bin/env python
#runs apparently only on linux
import os

if __name__ == "__main__":
    from restraintmaker.restraintMaker import run_plugin_gui, _check_importing_packages

    # check Importing
    _check_importing_packages()

    # Start pymol
    import pymol
    pymol.finish_launching()

    # pymol cmd
    from pymol import cmd
    # load
    in_file_path = os.path.dirname(__file__) + "/../restraintmaker/test/test_files/ligand_system/7_veryDifferent_Ligands.pdb"
    cmd.load(in_file_path)

    # visualization
    cmd.hide()
    cmd.show("sticks")

    # view
    cmd.center()
    cmd.zoom()

    # Run Plugin
    run_plugin_gui()


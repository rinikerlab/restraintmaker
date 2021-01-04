#!/usr/bin/env python
import os
import sys

import __main__

# include python path
package_path = os.path.dirname(__main__.__file__)
sys.path.append(package_path)

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
    in_file_path = os.path.dirname(__file__) + "/test/test_files/ligand_system/7_veryDifferent_Ligands.pdb" if (
        not os.path.dirname(__file__) == "") else "/test/test_files/ligand_system/7_veryDifferent_Ligands.pdb"
    # in_file_path=os.path.dirname(__file__)+"/test/test_files/ligand_system/9_similar_Ligands.pdb" if(not os.path.dirname(__file__) == "") else "/test/test_files/ligand_system/9_similar_Ligands.pdb"

    cmd.load(in_file_path)
    # cmd.fetch("1VPR")

    # visualization
    cmd.hide()
    cmd.show("cartoon", "polymer")
    cmd.show("sticks")

    # view
    cmd.center()
    cmd.zoom()

    # Run Plugin
    run_plugin_gui()

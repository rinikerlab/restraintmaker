#!/usr/bin/env python
"""
    RestraintMaker
    Is a programm, that estimates optimal distance restraints for a given bulk of molecules, to interconnect them.

    This is the main part of restraintmaker. From here the PyMol Wizard/Plugin is triggered.

    Required Packages Needed:
        python 3
        rdkit
        pymol
        numpy
        scipy

"""

import os
import sys
import __main__
import traceback

###IF YOU WANT TO INCLUDE RESMAKER AS PLUGIN
def __init_plugin__(app=None):
    """
    Required for PyMol plugin installation

    """
    from pymol.plugins import addmenuitemqt
    addmenuitemqt('RestraintMaker', run_plugin_gui)


def run_plugin():
    """
    This is the start point for running the plugin in PyMol,
    first the package dependencies are checked and then the GUI is started
    """
    # Check imports
    try:
        _check_importing_packages()
    except ImportError as err:
        print("\nERROR: \n "
              "Importing failed! Please try to install package manually.\n")
        print("ERROR\n " + "\n ".join(err.args))
        exit(1)
    except Exception as err:
        print("\nERROR: \n "
              "A general error occured:\n " + "\n ".join(err.args))
        print("ERROR\n " + "\n ".join(err.args))
        exit(1)

    # run plugin:
    try:
        run_plugin_gui()
    except Exception as err:
        print("\nERROR: \n "
              "A general error occured:\n")
        print("ERROR\n " + "\n ".join(err.args))
        exit(1)


def run_plugin_gui():
    """
        This function fires up the plugin gui in pymol.
    """
    from pymol import cmd

    # include python path
    package_path = os.path.dirname(__main__.__file__)
    sys.path.append(package_path)

    # cmd
    import restraintmaker.interface_Pymol.RestraintMaker_Pymol as restMPym
    cmd.set_wizard(restMPym.Restraints_Wizard())


def _check_importing_packages():
    """
        This function checks if all needed packages are there.
    """

    # IMPORT PYMOL
    try:
        import pymol
    except Exception as err:
        print("WARNING: could not find PyMOL in enviroment. Try installing via anaconda3")
        if ("conda" in sys.modules):
            import conda.cli as cli
            cli.main('conda', 'install', '-y', '-c schrodinger', 'interface_Pymol')
            # start interface_Pymol
            import pymol
        else:
            raise ImportError(
                "Could not find pymol package or conda enviroment to install the package.\n " + "\n ".join(
                    err.args))

    # IMPORT RDKIT
    try:
        import rdkit
    except Exception as err:
        raise ImportError(
            "Could not find rdkit package! And also couldn't find a conda enviroment to install Tools_rdkit.\n " + "\n ".join(
                err.args))



###start pymol if main
if __name__ == "__main__":
    # Check if imports are possible
    try:
        _check_importing_packages()
    except ImportError as err:
        print("\nERROR: \n Importing failed! Please try to install package manually.\n " + "\n ".join(err.args))
        traceback.print_exc(file=sys.stderr)
        exit(1)
    except Exception as err:
        print("\nERROR: \n A general error occured  while checked.\n " + "\n ".join(err.args))
        traceback.print_exc(file=sys.stderr)
        exit(1)

    # Start pymol
    import pymol

    pymol.finish_launching()

    # run plugin:
    try:
        run_plugin_gui()
    except Exception as err:
        print("\nERROR: \n A general error occured while running.\n " + "\n ".join(err.args))
        traceback.print_exc(file=sys.stderr)
        exit(1)

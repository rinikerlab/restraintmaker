#!/usr/bin/env python
"""
.. automodule: restraintmaker
    This is the main part of restraintmaker

    Needed:
        python 3
        rdkit
        pymol
        numpy
        scipy
        (recommended anaconda3 / conda-cli)

"""

import os
import sys

import __main__

###IF YOU WANT TO INCLUDE RESMAKER AS PLUGIN
"""
def __init_plugin__(app=None):
    from interface_Pymol.plugins import addmenuitemqt
    addmenuitemqt('restraintmaker', run_plugin_gui)
"""


def run_plugin():
    """
    .. autofunction:: run_plugin
    :return:
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
    .. autofunction:: run_plugin_gui
        This function allows to run the plugin gui in pymol.
    :return:
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
    .. autofunction:: check_importing_packages
        This function checks if all needed packages are there.
    :return:
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


import traceback

###WIZARD
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

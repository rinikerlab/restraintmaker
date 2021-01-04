import traceback, sys
from restraintmaker.restraintMaker import run_plugin_gui, run_plugin, _check_importing_packages


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

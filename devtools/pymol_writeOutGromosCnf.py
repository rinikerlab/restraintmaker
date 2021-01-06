"""
funcs
"""

from PyGromos.files import coord
from generalutilities.function_libs.gromos.files.blocks import blocks
from pymol import cmd


def write_out_cnf(out_path: str, selection: str = "all") -> str:
    """with this function you can write out interface_Pymol structures to cnf :param
    out_path: :param selection: :return:

    Args:
        out_path (str):
        selection (str):
    """
    # check input
    if not out_path.endswith(".cnf"):
        raise Exception("outpath has to end with: .cnf")

    # gather data
    atom = {}
    cmd.iterate_state(1, selection=selection,
                      expression="atom.update({ ID :{ \"type\":name, \"resn\":resn, \"resi\":resi, \"coord\":[x,y,z]}})",
                      space=locals())

    # generate atom position line from dict
    pos_list = []
    for x in atom:
        pos_list.append(
            blocks.atomP(resID=atom[x]["resi"], resName=atom[x]["resn"], atomType=atom[x]["type"], atomID=x,
                         xp=atom[x]["coord"][0] / 10, yp=atom[x]["coord"][1] / 10, zp=atom[x]["coord"][0] / 10))

    # build cnf blocks
    pos = blocks.atom_pos_block(pos_list)
    title = blocks.title_block("TEST")

    # build cnf
    tmp_cnf = coord.Cnf({"POSITION": pos, "TITLE": title})
    tmp_cnf.path = out_path
    tmp_cnf.write(tmp_cnf.path)

    return out_path


if __name__ == "__main__":
    import __main__

    __main__.pymol_argv = ['interface_Pymol', "-qc"]  # , '-qc'] # Pymol: quiet and no GUI
    import pymol

    pymol.finish_launching()

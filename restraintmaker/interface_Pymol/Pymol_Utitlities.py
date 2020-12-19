"""
    contains functions that are useful to the restraintmaker_PyMOL module, but do not need to be within the Distance_restraints class
    This also includes  functions that are only usefull for debugging and testing

"""

import os
import random
import sys
import threading
import time
import typing as t

from pymol import cmd

sys.path.append(os.path.dirname(__file__) + "/..")

from restraintmaker.utils import Utilities as u
from restraintmaker.utils.Utilities import print


# ----------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------Pymol helper funcs-----------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------


def pymol_selection_to_atom_list(sele: str) -> t.List[u.Atom]:
    """
    pymol_selection_to_atom_list(...) converts a interface_Pymol selection and returns a list of restraintmaker_Logic.utilities.Atom
    :return: A list of Atoms
    :rtype: t.List[u.Atom]
    """
    atom_info: t.List[t.Dict] = []
    # Get all the information about the Atoms as a dict. I can not directly put it into an atom, because the interface_Pymol interpreter will not know what an u.Atom is.
    cmd.iterate_state(-1, selection=sele,
                      expression="atom_info.append({\"elem\":elem,\"id\":ID, \"name\":name, \"x\":x, \"y\":y, \"z\":z, \"model\": model, \"chain\":chain, \"resn\":resn, \"resi\":resi, \"alt\":alt,  \"label\":label, \"b\":b})",
                      space=locals())

    # Convert dict to Atoms
    atoms: t.List[u.Atom] = []
    for a in atom_info:
        x = u.Atom(elem=a['elem'], id=a['id'], name=a['name'], x=a['x'], y=a['y'], z=a['z'], chain=a['chain'],
                   resn=a['resn'], resi=a['resi'], alt=a['alt'], label=a['label'], b=a['b'])
        atoms.append(x)

    return (atoms)


def _atom_list_to_pymol_selection(atom_list: t.List[u.Atom]) -> str:
    """
    atom_list_to_pymol_selection converts a List of Atoms to a string that pymool can understand.

    :warning; Pymol can only handle a string with a length of ca 1000. (ca 200 atoms.) Consider using help_pymol_with_big_atoms_List
    :param atom_list: A List of Atoms
    :type atom_list: t.List[u.Atom]
    :return: A interface_Pymol selection expression, containing all those Atoms
    :rtype: str
    """
    selection_expression: str = "ID " if len(atom_list) > 0 else 'none'
    for a in atom_list:
        selection_expression += str(a.id) + ","

    if (len(selection_expression)) > 1000:
        print(
            'WARNING in _atom_list_to_pymol_selection. Pymol can not handle selection expression with more than ca 200 Atoms. Consider using help_pymol_with_big_atom_list',
            mv=3)

    return selection_expression


def colorize(atoms):
    '''
    colorize_atoms: Some tests for different ways of colouring atoms.

    When using the predefinied names (e.g. "gray") PyMol raises a Quiet Error. With this method I can test, that Pymol does accept the different
    colour numbering schemes (s for spectrum, c for complementary etc.)

    :param atoms:
    :type atoms:
    :return:
    :rtype:
    '''
    try:
        const_time_max = 10000

        print('inner1')

        def _inner_colorize_1():
            time_max = const_time_max
            time_step = time_max / (3 * len(atoms))
            while (time_max > 0):
                time_max -= time_step
                __a = atoms[random.randint(0, len(atoms) - 1)]
                __b = random.random() < 0.5
                __c = random.randint(1, 999)
                __d = ('s' if __c else 'c') + str(__c)
                cmd.color(color=__d, selection=_atom_list_to_pymol_selection([__a]))
                time.sleep(time_step / 1000)

        print('inner2')

        def _inner_colorize_2():
            time_step = const_time_max / 100
            for i in range(1, 100):
                __a = str(100 - i)
                __b = 'gray' + __a.zfill(2)
                print(__b)
                cmd.color(color=__b, selection=_atom_list_to_pymol_selection(atoms))
                time.sleep(time_step / 1000)

        print('inner3')

        def _inner_colorize_3():
            time_max = const_time_max
            time_step = (time_max / len(atoms) * 3)

            cmd.color(color='white', selection=_atom_list_to_pymol_selection(atoms))
            while (time_max > 0):
                time_max -= time_step
                __a = random.randint(0, len(atoms) - 1)
                cmd.color(color='red', selection=_atom_list_to_pymol_selection([atoms[__a]]))
                time.sleep(time_step / 1000)

        def _inner_colorize_4():
            time_step = const_time_max / 100
            cols = []
            cmd.iterate_state(state=-1, selection=_atom_list_to_pymol_selection(atoms), expression='cols.append(ID)',
                              space=locals)
            print(cols, mv=1)

        inner_funcs = [_inner_colorize_1, _inner_colorize_2, _inner_colorize_3]
        rand_i = random.randint(0, len(inner_funcs) - 1)

        my_thread = threading.Thread(target=inner_funcs[rand_i])
        # my_thread=threading.Thread(target=_inner_colorize_4())

        my_thread.start()
    except:
        pass


def create_pymol_objects_for_molecules(exclude: t.List[str] = []):
    """
    creates separate selecetions for all molecules in interface_Pymol
    :param exclude: A list of molecules that should not be considered, e.g. 'solv;, 'protein', etc
    :type exclude: t.List[str]
    :return: The list of the molecule object names
    :rtype: t.List[u.Atom]
    """
    # TODO; WRITE THE NAME OF THE MOLECULE INOT THE ATOM SOMEHOW
    # TODO: make more stable if protein is present!

    molecules: str = []

    def _append_molecule(m: str):
        if not m in molecules:
            molecules.append(m)

    """
    Deapreacitated, as chains are maybe more general!
    #try chain wise
    cmd.iterate_state(-1, selection='all', expression='_append_molecule(resi)', space=locals())

    molecule_names = []

    for m in molecules:
        if not m in exclude:
            name = 'mol_'+ m
            cmd.create(name=name, selection='resi ' + m)
            molecule_names.append(name)
    """

    cmd.iterate_state(-1, selection='all', expression='_append_molecule(chain)', space=locals())

    molecule_names = []

    i = 1
    for m in molecules:
        if not m in exclude:
            name = 'mol_' + str(i)
            cmd.create(name=name, selection='chain ' + m)
            molecule_names.append(name)
        i += 1
    return molecule_names

    # The selection containing all the Molecules must be deleted. If it is not, iterate state will duplicate each atom!
    cmd.delete(cmd.get_object_list()[0])  # Delete the object containing all molecules at once


def help_pymol_with_big_atom_list(pymol_function: t.Callable, atom_list: t.List[u.Atom], *args, **kw_args):
    '''
    help_pymol_with_big_atom_list can be wrapped around any interface_Pymol function, that takes 'selection' as an argument.

    help_pymol_with_big_atom_list  will call pymol_function several times with small enough parts of the list

    :param atom_list: A list of atoms
    :type atom_list: t.List[u.Atom]
    :param pymol_function: Any funciton from one of the interface_Pymol modules, who takes 'selection' as an argument
    :type pymol_function; t.Callable
    :param args: arguments for the interface_Pymol function
    :param kwargs; keyword argument for the interface_Pymol functiom
    :return: a list of all return values of the seperate calls #TODO; Maybe the sum?
    :rtype: t.List
    :raises: ValueError if interface_Pymol function is not passed the correct arguments
    '''

    return_list = []

    # interface_Pymol cant take long selection expressions, with more than ca 200 atoms
    # The easiest expression is 'ID x,y,z...'

    # There seems to be some variation in how long a word vcan be. The limit seems to be between 200-250 atoms => ca 1000 characters
    step = 100
    start = 0
    stop = step
    while start < len(atom_list):
        stop = min(stop, len(atom_list))
        selection_expression = _atom_list_to_pymol_selection(atom_list[start:stop])

        # CALL THE ACTUAL PYMOL FUNCTION
        pymol_function(*args, **kw_args, selection=selection_expression)

        start += step
        stop = start + step

    return return_list


def _convert_molecules_to_pdb() -> t.List[str]:
    '''
    convert_molecules_to_pdb  gives a list containing one pdb string for each molecule
    #TODO: Write/find a conversion fucntion that I can use on the u.Atoms as the Logic handler knows them
    :warning: atoms willbeassigned different IDS than they currently have in interface_Pymol!!!!!
    :return: List of pdb string for each molecule in interface_Pymol
    :rtype: t.List[str]
    '''

    pymol_objects = cmd.get_object_list()
    pdb_molecules = [cmd.get_pdbstr(obj) for obj in pymol_objects if obj[:4] == "mol_"]

    return pdb_molecules

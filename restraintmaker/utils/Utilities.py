"""
.. automodule:: Utilities
    :members:
    TODO:write DOCU!

utilities contains functions and definitions that are usefull to all modules in the restraintmaker project
"""
import builtins
import os
import typing as t
from collections import namedtuple
from enum import Enum

import __main__
if not os.path.basename(__main__.__file__).startswith("test"):  # TODO: replace by checking if display available!!!!
    from PyQt5 import QtWidgets  # make criteria to only activate if needed

# ---------------------------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------Input / Output functions-----------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------------------------------------


# Overwrite print function for debugging
verbosity_threshold = 1  #


def print(*x: str, mv: int = 0, **kargs):
    '''print is overridden for the restraintmaker project:
        The parameter mv indicates the minimal verbosity at which the statement will still be printed.

        The verobsity threshold of the program can be chagned in utilities.
        0: Print everything, 1 = Debug, 2 = Develop, 3 = Interested user, 4 = disinterested user (Error messages only)

        :param x: Arguments to be printed
        :type x: Any
        :param mv: Minimal verbosity at which this statement will still be printed.

    '''

    if mv >= verbosity_threshold:
        builtins.print(*x, **kargs)

    # Can be used find what is causing a certain output we do not want anymore
    find_print = ''
    if find_print != '':
        for output in x:
            if find_print in str(output):
                raise BaseException


def execute_at_different_verbosity(new_verbosity_threshold, func, *args, **kwargs):
    '''execute_at_different_verbosity: Wrapper function to change the program verbosity during the execution of one function.


        :warning: All Errors will be passed on to the calling function. But the verbosity_thresheold will be reset to the previous value in any case.
        :param new_verbosity_threshold: verbosity threshold to be applied duribng the execution of func
        :type new_verbosity_threshold: int
        :param func: Function which should be executed
        :type func: t.Callable
        :param args: arguments for func
        :type args: any
        :param kwargs: Key word arguments for func
        :type kwarg: Any
        :returns: return value of func
        :rtype: Any
    '''
    global verbosity_threshold
    old_verbosity_threshold = verbosity_threshold
    verbosity_threshold = new_verbosity_threshold

    try:
        return_value = func(*args, **kwargs)
        return return_value
    finally:
        verbosity_threshold = old_verbosity_threshold


# ----------------------------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------Containers / namedtuples / classes /enums, Exceptions-------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------------

# TODO: These fields are interface_Pymol based. To generalize: Just keep the ones relevant for the Program logic, and add a comment, which the GUI model can use for extra info
#  There should be a field molecule. So far I have kjust been using the resi field all over the program
#  An additional attribute comment might be good, which could be used to store additonal information by other GUI modules.
Atom = namedtuple('Atom', 'elem id name x y z chain resn resi alt b label')
# TODO: Use Pairwise distance restraint instead of connection
AtomPair = namedtuple('AtomPair', 'a1, a2, distance')  # Used in Optimizer to sotre informations more efficiently

# Not used in running code any more. Only in the old function Optimizer._connect to Molecules Kruskal
RestraintPair = namedtuple('RestraintPair', 'r1,r2,distance')


# TODO LOGIC: Combine modes and action states into one class
class ActionState:
    '''ActionState defines the different states the possible actions (filter, optimize etc.) can be in'''
    # Using the Pymol convention for numbers
    # If I use 4 for CURRENTLY_DISABLED and 1 for ALWAZS DISABLED, PYMOL WILL ONLY USE THE GREEN FONT FOR BOTH
    CURRENTLY_ENABLED = 3
    CURRENTLY_DISABLED = 1
    ALWAYS_ENABLED = 2
    ALWAYS_DISABLED = 4
    # Only usefull during debugging


class EventType(Enum):
    '''EventType defines the different events that a selection can deal with'''

    SELECT = 'select'
    MOVE = 'move'
    SIZE = 'size'
    CONFIRM = 'confirm'


class BadArgumentException(Exception):
    """
    A BadArgumentException will be raised by the Filter, Selection and Optimizer classes if their function get-args(...)
    does not produce arguments in the desired format.
    """

    def __init__(self, *args):
        super().__init__(*args)


class NoOptimalSolutionException(Exception):
    '''
    A no optimal solution excpetion will be raised by an Optimizer, if it can not find a solution fulfilling its criteria.
    '''

    def __init__(self, args):
        super().__init__(args)


def create_input_dialog(parent_window, title: str) -> t.Callable[[str], str]:
    """
    create_input_dialog returns a function which:

    Opens a pyQt5 input Dialog with the speciefed message. In this format it can be provided to the get_args functions
    :param message: Message to be dispayed in the input dialog
    :type message: str
    :return:Returns a function, that opens an input dialog, with the specified title, and is centered at the interface_Pymol window.
    :rtype: t.Callable[str,[str])
    """

    def input_dialog_with_title(message: str):
        input, ok_pressed = QtWidgets.QInputDialog.getText(parent_window, title, message)  # self.cmd.gui.get_qtwindow()
        if ok_pressed:
            return input
        else:
            return None  # Do not catch the Error here. I'd prefer to do that in get args

    return input_dialog_with_title


def create_file_open_dialog():
    """
       create_file_dialog returns a function which:
       Opens a pyQt5 file Dialog with the speciefed message. In this format it can be provided to the get_args functions
       TODO: I am doing this mainly as an execrcise in currying. Check if it is really usefull, or just foolish playing around
       :param message: Message to be dispayed in the input dialog
       :type message: str
       :return:Returns a function, that opens an input dialog, with the specified title, and is centered at the interface_Pymol window.
       :rtype: t.Callable[str,[str])
       """

    def file_dialog_with_title(message: str):
        path = QtWidgets.QFileDialog.getOpenFileName(None, 'Which File should be opened?', directory=os.getcwd())[0]
        return path

    return file_dialog_with_title


def create_file_save_dialog():
    """
       create_file_dialog returns a function which:
       Opens a pyQt5 file Dialog with the speciefed message. In this format it can be provided to the get_args functions
       TODO: I am doing this mainly as an execrcise in currying. Check if it is really usefull, or just foolish playing around
       :param message: Message to be dispayed in the input dialog
       :type message: str
       :return:Returns a function, that opens an input dialog, with the specified title, and is centered at the interface_Pymol window.
       :rtype: t.Callable[str,[str])
       """

    def file_dialog_with_title(message: str):
        path = QtWidgets.QFileDialog.getSaveFileName(None, message, directory=os.getcwd())[0]
        return path

    return file_dialog_with_title


# TODO: Might be even nicer if we do not pass a list of acceptable values but a bool function
# TODO: Pass name of the variable as arg, so we can get better Error messages ()
def check_or_convert_argument(input, desired_type, acceptable_values=None):
    '''try to convert should be used in all get_args functions to check the input from the input function and try to convert it if necessary

    :param input: input as received from the function
    :type input: t.Any
    :param desired_type:
    :type desired_type: type
    :param acceptable_values: Optional. a List of set Values the input must have.
    :type acceptable_values: t.List
    :return: converted input
    :rtype: desired_type
    :raises u.BadArgument Exception if conversion is not possible or the input is not an acceptable value
    '''

    # CONVERT
    converted_input = None
    if input.__class__ == desired_type:
        converted_input = input
    else:

        # raise Error outside of except Block to avoid carrying along hte old error.
        conversion_failed = False
        try:
            return desired_type(input)
        except (ValueError, TypeError):
            conversion_failed = True
        if conversion_failed:
            raise BadArgumentException('Failed to convert the input \'' + str(input) + '\' from ' + (
                input.__class__.__name__) + ' to ' + desired_type.__name__)

    # CHECK VALUE, IF THERE IS A LIST OF ACCEPTABLE VALUES
    if (not acceptable_values == None) and (not converted_input in acceptable_values):
        raise BadArgumentException(
            str(converted_input) + ' is not an acceptable input. Acceptable values are:' + ' '.join(
                [str(a) + ' ,' for a in acceptable_values]))

    return converted_input


# TODO get_args: I am not happpy with the create_multi_dialog solution: Now we have to 'prepare' the input function in set_x_type by checking what x(optimizer, filter...) type we want to create.
# The initial idea was that input_function is passed to get_args, which can then call it as needed, without having to know the format of its arguments
# We have to choose either:1) input function returns exactly one str. if more inputs are required, input fkt is called several times. => Does not allow to ask for all in one dialog
#                   or    :2) Abolish get_args. Each type has a list required_args. the logic handler will read that list and then get them all from the user. The Optimizer?Filter... then only has to check the type

def create_multi_dialog(title, inputs: t.List[str], options: t.Dict[str, t.List[str]] = {}, default: t.Dict = {}) -> \
t.Callable[[str], t.List[str or t.List[str]]]:
    '''
    create_multi_dialog return
    :warning: Do not call before the constructor o f the Pymol module has finished, or pyumol will hate you.
    :warning: Will not throw an error if close button is pressed. Just return an enmpty list
    :param title: Title of theDialog
    :type title: str
    :param inputs; all inputs that should be asked for
    :type  inputs: t.List[str]
    :param options: optional, specifies the options an input can have. (input as key, list of acceptable values as values)
    :type options: t.dict[str,t.List[str]]
    :param default: Specifies the options an input can have. (input as key, list of acceptable values as value0
    :type default: t.dict[str,t.Any]
    :return: A function which will open a dialog, asking for the inputs specified in inputs
    :rtype: t.Callable[[str],t.List[str or t.List[str]]]
    '''

    class MultiDialog(QtWidgets.QDialog):
        'can be used to create input functions with as many arguments as necessary on the spot'

        def __init__(self):
            super().__init__()
            self.setWindowTitle(title)

        # def __del__(self):
        #     print('HELP! I am beeing deleted',mv=1)

    def _show_multi_dialog(dummy_str):
        '''The whole function _show_multi_dialog will be returned and can then be passed to get_args'''
        answers = []
        dial = MultiDialog()
        # In python components are added to the layot, not the window
        lay = QtWidgets.QGridLayout()

        # CREATE Components and add them to the layout
        fields: t.list[QtWidgets.QLineEdit or QtWidgets.QComboBox] = []
        for i_i, input in enumerate(inputs):
            if input in options.keys():
                new_field = QtWidgets.QComboBox()
                new_field.addItems(options[input])
                if input in default.keys():
                    new_field.setCurrentText(default[input])
            else:
                new_field = QtWidgets.QLineEdit()
                if input in default.keys():
                    new_field.setText(default[input])

            fields.append(new_field)  # TODO: Acutallz I can just loop over the grid layout to get all components
            lay.addWidget(QtWidgets.QLabel(input + ':'), i_i, 1)  # Label
            lay.addWidget(new_field, i_i, 2)

        OK_Button = QtWidgets.QPushButton('OK')

        def _ok_button_pressed():
            '''iterate over all input fields and add their value to answers'''
            for field in dial.children():
                # field = lay.itemAt(i)
                if isinstance(field, QtWidgets.QLineEdit):
                    answers.append(field.text())
                elif isinstance(field, QtWidgets.QComboBox):
                    answers.append(field.currentText())
                # else: QLabel etc => Do nothing
            dial.close()
            print(answers, mv=1)

        OK_Button.clicked.connect(_ok_button_pressed)

        lay.addWidget(OK_Button, len(inputs), 2)
        dial.setLayout(lay)
        dial.exec_()  # exec interrupts program flow unitl dialog is closed
        return answers

    return _show_multi_dialog


# ----------------------------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------Functions realted to chemical concpts-------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------------


def order_atoms_by_molecule(atoms: t.List[Atom]) -> t.Dict[str, t.List[Atom]]:
    """
    :param atoms: A list of atoms
    :type atoms: t.List[Atom]
    :return: A Dictuionary: Molecules names as keys, List of atoms in the molecule as entry
    :rtype: t.Dict[t.List[Atom]]
    """
    # TODO: I am using the atom.chain attribute to identify the molecule. That is to ypmol centered. Is hould give an attribute molecule to each atom and the only use the interface_Pymol chain attribvute to set molecules.

    Molecules: t.Dict[t.List[Atom]] = {}  # Every Molecules is a list of atoms, with the Molecule attribute as key
    for a in atoms:
        if a.resi in Molecules:
            Molecules[a.resi].append(a)
        else:
            Molecules.update({a.resi: [a]})  # dict.update adds entries of one dict to another

    return [v for v in Molecules.values()]

    def convert_atoms_to_pdb_molecules(atoms: t.List[Atom]) -> t.List[str]:
        '''
        :param atoms: List of atoms
        :type atoms: t.List[Atom]
        :return: pdb strings of that molecule
        :rtype: t.List[str]
        '''

        # 1) GROUP ATOMS BT MOLECULES
        molecules: t.Dict[t.List[Atom]] = {}
        for a in atoms:
            if a.resi in molecules.keys():
                molecules[a.resi].append(a)
            else:
                molecules[a.resi] = [a]

        # 2) CONSTUCT PDB BLOCKS
        pdb_molecules: t.List[str] = []

        for m in molecules.values():
            atoms_as_lines: t.List[str] = []
            for a in m:
                atoms_as_lines.append(
                    str(a.id) + '\t' + a.name + '\t' + a.resn + '\t' + str(a.resi) + '\t' + '{:0.3f}'.format(
                        a.x) + '\t' +
                    '{:0.3f}'.format(a.y) + '\t' + '{:0.3f}'.format(
                        a.z) + '\t1.00\t0.00\t\t' + a.elem)  # TODO: What are 1.00 and o.oo for? b-value?
                # atoms_as_lines.append('{:>3}.format'(str(a.id))+'\t'+'{:3}'.format(a.name)+'\t'+a.resn+'\t'+str(a.resi)+'\t'+'{:0.3f}'.format(a.x)+'\t'+\
                #                      '{:0.3f}'.format(a.y)+'\t'+'{:0.3f}'.format(a.z)+'\t1.00\t0.00\t\t'+a.elem) #TODO: What are 1.00 and o.oo for? b-value?

            # Sort by Id: => convert str up do first space to int
            atoms_as_lines = sorted(atoms_as_lines, key=lambda x: int(x[:x.index('\t')]))
            atoms_as_lines = ['ATOM\t' + line for line in atoms_as_lines]
            molecule_as_str = '\n'.join(atoms_as_lines) + '\nEND'
            # molecule_as_str = molecule_as_str.replace('\t','    ')
            pdb_molecules.append(molecule_as_str)

            print(pdb_molecules[-1], mv=1)

        return pdb_molecules


# ----------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------Various usefull fucntions not realted to chemical concepts--------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------------


def get_all_subclasses(parent: type, include_private_sublcasses=False) -> t.List[type]:
    """
    get_all_subclasses returns ALL subclasses of a type, including grandchilder etc.
    :param parent: Parentn class
    :type obj: type
    :return: List of all subclasses
    :rtype: t.List[type]
    """

    # I refuse to use recursion in Python. I am afraid that I would melt the CPU.
    # => do it iteratively

    i: int = 0  # Index of last class that has been checked for subclasses
    all_subs = parent.__subclasses__()
    while i < len(all_subs):

        # Check if the subclass I found is already on the list, becaus OF COURSE Python supports multiple inheritance.
        # Am I not glad I did not use recursion?

        new_subs = (all_subs[i]).__subclasses__()
        for sub in new_subs:
            if not sub in all_subs:
                all_subs.append(sub)
        i += 1

    if include_private_sublcasses:
        return all_subs
    else:
        return [sub for sub in all_subs if not sub.__name__[0] == '_']


def do_nothing(dummy_string: str = ""):
    """
    do_nothing(...) does nothing.
    Can be passed to the *.get_args functions which require a function of a string as parameter. Because it is more readable than lambda x: None
    :param dummy_string:  Dummy String to satisfy format
    :type dummy_string: str
    :return: -
    :rtype: -
    """
    pass


def check_for_doubles(list) -> bool:
    '''Just used for debbuging Checks if a List contains doubles
    :param list: Any list
    :type list: t.List
    :return: Does the list contain doubles?
    :rtype: bool
    '''
    for i in range(len(list) - 1):
        for j in range(i + 1, len(list)):
            if list[i] == list[
                j]:  # WTF?1  Seems to check for ref identity, not just value identity, even though all docs say ity checks only for value!!!!!!!
                return (True)
    return False


def check_restraint_pairs_for_doubles(list):  # Also consider that a1 and a2 can be switches
    '''
    check_restraint_pairs_for_doubles checks a list of pairs for doubles. Pairs count as doubles if the order of elements is changed.
    :param list: A list of tuples
    :type list: t.List[t.Tuple]
    :return: Does the list contain doubles?
    :rtype: bool
    '''
    for i in range(len(list) - 1):
        for j in range(i + 1, len(list)):
            if (list[i].r1 == list[j].r1 and list[i].r2 == list[j].r2) or (
                    list[i].r1 == list[j].r2 and list[i].r2 == list[j].r1) or list[i].distance == list[j].distance:
                return True
    return False

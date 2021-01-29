"""
    Utilities
    utilities contains functions and definitions that are usefull to all modules in the restraintmaker project
"""
import builtins
import typing as t
from collections import namedtuple, defaultdict

"""
Verbosity
"""
verbosity_threshold = 1  #


# Overwrite print function for debugging
def print(*x: str, mv: int = 0, **kargs):
    """
    print is overridden for the restraintmaker project:
        The parameter mv indicates the minimal verbosity at which the statement will still be printed.

        The verobsity threshold of the program can be chagned in utilities.
            0 = Print everything,
            1 = Debug,
            2 = Develop,
            3 = Interested user,
            4 = disinterested user (Error messages only)

    Parameters
    ----------
    x: Any
        Arguments to be printed
    mv: int
         Minimal verbosity at which this statement will still be printed.

    Returns
    -------
    NoReturn
    """

    if mv >= verbosity_threshold:
        builtins.print(*x, **kargs)

    # Can be used find what is causing a certain output we do not want anymore
    find_print = ''
    if find_print != '':
        for output in x:
            if find_print in str(output):
                raise BaseException


def execute_at_different_verbosity(new_verbosity_threshold, func, *args, **kwargs):
    """
     Wrapper function to change the program verbosity during the execution of one function.
        @Warning All Errors will be passed on to the calling function. But the verbosity_thresheold will be reset to the previous value in any case.

    Parameters
    ----------
    new_verbosity_threshold: int
        verbosity threshold to be applied duribng the execution of func
    func: t.Callable
         Function which should be executed
    args: Any
    kwargs: Any

    Returns
    -------
    Any:
         return value of func
    """

    global verbosity_threshold
    old_verbosity_threshold = verbosity_threshold
    verbosity_threshold = new_verbosity_threshold

    try:
        return_value = func(*args, **kwargs)
        return return_value
    finally:
        verbosity_threshold = old_verbosity_threshold


"""
Exceptions
"""


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


"""
Chemical concepts related
"""
# named tuples for code
# representation of an atom in pymol
Atom = namedtuple('Atom', 'elem id name x y z chain resn resi alt b label')
# represents a pair of restraints, needed for optimizers
RestraintPair = namedtuple('RestraintPair', 'r1,r2,distance')


def order_atoms_by_molecule(atoms: t.List[Atom]) -> t.Dict[str, t.List[Atom]]:
    """

    Parameters
    ----------
    atoms : t.List[Atom]
         A list of atoms
    Returns
    -------
    t.Dict[t.List[Atom]]
        A Dictionary: Molecules names as keys, List of atoms in the molecule as entry
    """

    Molecules: t.Dict[t.List[Atom]] = {}  # Every Molecules is a list of atoms, with the Molecule attribute as key
    for a in atoms:
        if a.resi in Molecules:
            Molecules[a.resi].append(a)
        else:
            Molecules.update({a.resi: [a]})  # dict.update adds entries of one dict to another

    return [v for v in Molecules.values()]


def find_atom_by_property(atoms: t.List[Atom], property_value: any, property_name: str = "id") -> Atom:
    """
        find_atom will look for the first atom with the specified property and value

    Parameters
    ----------
    property_value : Any
        value of the property
    property_name : str, optional
        name of the property (default: id)

    Returns
    -------
    u.Atom
        returns the atom

    Raises
    ------
    ValueError
        if there is no or more than one atom with that id
    """

    if (all([hasattr(a, property_name) for a in atoms])):
        all_hits = list(filter(lambda a: a.id == property_value, atoms))
        if len(all_hits) == 0:
            raise ValueError('There is no atom with id: ' + str(id))
        elif len(all_hits) > 1:
            raise ValueError('There is more than one atom with id: ' + str(id))
        else:
            return all_hits[0]
    else:
        raise ValueError("atom does not have Property: " + property_name + " but: " + str(vars(atoms[0])))


def convert_atoms_to_pdb_molecules(atoms: t.List[Atom]) -> t.List[str]:
    """
        This function converts the atom list into pdb blocks.

    Parameters
    ----------
    atoms : t.List[Atom]
        List of atoms
    Returns
    -------
    t.List[str]
         pdb strings of that molecule
    """
    # 1) GROUP ATOMS BT MOLECULES
    molecules = defaultdict(list)
    for a in atoms:
        molecules[a.resi].append(a)

    # 2) CONSTUCT PDB BLOCKS
    #ref: https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
    pdb_format = "ATOM  {:>5d}  {:<3}{:1}{:>3}  {:1}{:>3d}{:1}   {:>7.3f}{:>7.3f}{:>7.3f}{:>5f}{:>6f}{:<3}{:>2} {:>2d}"
    dummy_occupancy= dummy_bfactor= dummy_charge = 0.0
    dummy_alt_location= dummy_chain= dummy_insertion_code= dummy_segment = ""

    pdb_molecules: t.List[str] = []
    for m_ID in sorted(molecules):
        m = molecules[m_ID]
        atoms_as_lines: t.List[str] = []
        for a in sorted(m, key= lambda x: x.id):
            atoms_as_lines.append(pdb_format.format(int(a.id), a.name, dummy_alt_location, a.resn, dummy_chain, int(a.resi), dummy_insertion_code, a.x, a.y, a.z, dummy_occupancy, dummy_bfactor, dummy_segment, a.elem, int(dummy_charge)))

        # Sort by Id: => convert str up do first space to int
        #atoms_as_lines = sorted(atoms_as_lines, key=lambda x: int(x[:x.index('\t')]))
        molecule_as_str = "TITLE "+a.resn+"\n"+'\n'.join(atoms_as_lines) + '\nEND'
        # molecule_as_str = molecule_as_str.replace('\t','    ')
        pdb_molecules.append(molecule_as_str)

        print(pdb_molecules[-1], mv=0)

    return pdb_molecules


"""
Various usefull functions
"""


def check_or_convert_argument(input, desired_type, acceptable_values=None):
    """
    try to convert should be used in all get_args functions to check the input from the input function and try to convert it if necessary

     TODO: Might be even nicer if we do not pass a list of acceptable values but a bool function
     TODO: Pass name of the variable as arg, so we can get better Error messages ()


    Parameters
    ----------
    input :  t.Any
        input as received from the function
    desired_type : type
        the needed type
    acceptable_values : t.List
        Optional. a List of set Values the input must have.

    Returns
    -------
     type
        converted input to desired_type

    Raises
    ------
    u.BadArgument_Exception
        if conversion is not possible or the input is not an acceptable value
    """

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


def get_all_subclasses(parent: type, include_private_sublcasses=False) -> t.List[type]:
    """
        get_all_subclasses returns ALL subclasses of a type, including grandchilder etc.

    Parameters
    ----------
    parent : type
        Parentn class
    include_private_sublcasses: bool, optional
        shall privates be included? (default: False)
    Returns
    -------
    t.List[type]
         List of all subclasses
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
    does nothing.
    Can be passed to the *.get_args functions which require a function of a string as parameter. Because it is more readable than lambda x: None


    Parameters
    ----------
    dummy_string: str
        Dummy String to satisfy format

    Returns
    -------
    NoReturn

    """
    pass


def check_for_doubles(list: t.List) -> bool:
    """
        Just used for debbuging Checks if a List contains doubles

    Parameters
    ----------
    list: t.List
        Any list

    Returns
    -------
    bool
        Does the list contain doubles?
    """

    for i in range(len(list) - 1):
        for j in range(i + 1, len(list)):
            if list[i] == list[
                j]:  # WTF?1  Seems to check for ref identity, not just value identity, even though all docs say ity checks only for value!!!!!!!
                return (True)
    return False


def check_restraint_pairs_for_doubles(list):  # Also consider that a1 and a2 can be switches
    """
       check_restraint_pairs_for_doubles checks a list of pairs for doubles. Pairs count as doubles if the order of elements is changed.

    Parameters
    ----------
    list : t.List[t.Tuple]
        A list of tuples

    Returns
    -------
    bool
         Does the list contain doubles?
    """
    for i in range(len(list) - 1):
        for j in range(i + 1, len(list)):
            if (list[i].r1 == list[j].r1 and list[i].r2 == list[j].r2) or (
                    list[i].r1 == list[j].r2 and list[i].r2 == list[j].r1) or list[i].distance == list[j].distance:
                return True
    return False

"""
The module Filters contains functionality of filtering selections.
"""

import os
import sys
import typing as t

sys.path.append(os.path.dirname(__file__))

from restraintmaker.utils import Utilities as u
# from Utilities import print
from restraintmaker.tools_Rdkit import \
    Rdkit_Functions  # TODO: Copy the functions and fragments we actually need to either utilites or Filter directly

from rdkit import Chem


class _Filter():
    """
    .. autoclass::    _Filter
        This is the private parent class to all filter Classes.
    """

    def __init__(self, atoms: t.List[u.Atom]):
        """
        :param atoms: List of atoms to be filtered
        :type atoms: t.List[u.Atom]
        :param args: Extra Arguments the subclasses might need
        :type args: t.List[str]
        """
        self.atoms: t.List[u.Atom] = atoms

    def get_args(self, input_function: t.Callable):
        """
        get_args(...) should be called by every subclass of _Filter, right after it is created.

         It uses an input function to find all arguments required by a certain instance of a _Filter
        :param input_function: A function that can get the input
        :type input_function: function(str)
        :param message: A string that is displayed if if input function needs to communicate with the user
        :type message: String
        :return: -
        :rtype: -
        :raises BadArgumentException: If the input function does not provide all arguments in the necessary format
        :raises NotImplementedError: If _Filter.get_args is called directly.
        """
        raise NotImplementedError("Call of function get_args(...) of abstract parent class _Filter")

    def filter(self) -> t.List[u.Atom]:
        '''
        filter must be overridden by every subclass of _Filter. Returns a filtered list of atoms
        :return: Filtered List of Atoms
        :rtype:  t.List[u.Atom]
        :raises NotImplementedErrror: If _Filter.filter is called directly
        '''

        raise NotImplementedError("Call of function filter() of abstract parent class _Filter")


class PropertyFilter(_Filter):
    # TODO: At the moment this works for exactly one value of the criterion. For some criteria it would make more sense, to give a list of acceptable values, or a continous range
    '''
    ..class PropertyFilter:  Filters by any attribute of atom
    '''

    def __init__(self, atoms: t.List[u.Atom]):
        """
        :param atoms: List of atoms that are filtered over
        :type: t.List[u.Atom]
        """
        super().__init__(atoms)

        # Attributes that will be set in get_args
        self.criterion: str = None  # Filter criterion
        self.value: t.Any = None  # Desired Value. Is introduced as string and then cast to datatype of the chosen criterion

    # TODO STRUCTURE: get_args calls input function twice. This way it is not possible to make shortcut sublcasses, that pas it standard values for the criterion directly. => Rethink. Maybe if I introduce provider functions in Pymol I can fix that

    def get_args(self, input_function: t.Callable):
        """
               get_args(...) Needs to be called, before filter() can be used. For the Property_Filter it finds out which
                property should be filter-criterion and then which value it should have.
               :param input_function: A function that will find the input needed ofr the filter.
               :type input_function: function(str)
               :param message: A string that is displayed if if input function needs to communicate with the user
               :type message: String
               :return: -
               :rtype: -
               :raises: BadArgumentException if the input function does not provide all arguments in the necessary format
        """
        input = input_function('Which atom-property should be the filter critrion?')
        self.criterion: str = u.check_or_convert_argument(input, str)
        if not self.criterion in u.Atom._fields:
            raise u.BadArgumentException(
                "The filter criterion \' " + self.criterion + " \' provided by the input function is not an attribute of an Atom.")

        input = input_function('What value should ' + self.criterion + ' have?')
        try:
            self.value = u.check_or_convert_argument(input, self.atoms[0].__getattribute__(self.criterion).__class__)
        except IndexError:  # Happens when len(self.atoms) is 0. No special handling necessary. It just means no atoms will be selected.
            pass

    def filter(self) -> t.List[u.Atom]:
        """\
        Returns a list of atoms, filtered using the provided criterion and value

        :return: Filtered atom list
        :rtype: t.List[u.Atom]
        '''
        """
        return list(filter(lambda x: x.__getattribute__(self.criterion) == self.value, self.atoms))


# TODO STRUCTURE: Should really be a sublclass of property filter, but that does not work at the moment, because property Filter calls its input function twice
# => See TODO GET_ARGS
class ElementFilter(_Filter):
    """
    ..class Element_Filter: Filter for certain Element
    """

    def __init__(self, atoms: t.List[u.Atom]):
        """
        :param atoms: List of atoms that are filtered over
        :type: t.List[u.Atom]
        """
        super().__init__(atoms)

    def get_args(self, input_function: t.Callable = lambda *args: 'C'):
        """
        get_args(...) Needs to be called, before filter() can be used.

         For the Element_Filter it finds out which element should be filtered
        :param input_function: A function that will find the input needed ofr the filter.
        :type input_function: function(str)
        :param message: A string that is displayed if if input function needs to communicate with the user
        :type message: String
        :return: -
        :rtype: -
        """
        input: str = input_function("Which Element do do you want to filter for?")
        self.in_path = u.check_or_convert_argument(input, str)
        self.elem = input

    def filter(self) -> t.List[u.Atom]:
        """
        Returns a filtered version of the class attribute atoms, filtered by the class attribute element


         :return: Filtered atom list
        :rtype: t.List[u.Atom]
        """

        # Use the built in python function filter(expression, iterable)
        return list(filter(lambda x: x.elem == self.elem, self.atoms))


class RingFilter(_Filter):
    '''...class Ring Filter chooses all Atoms that are part of a Ring
    :warning: using tools_Rdkit. Needs conversion to tools_Rdkit molecules. => Needs the pdb as input
    '''

    def __init__(self, atoms: t.List[u.Atom]):
        """
        :param atoms: List of atoms to be filtered
        :type atoms: t.List[u.Atom]
        :param args: Extra Arguments the subclasses might need
        :type args: t.List[str]
        """

        self.atoms: t.List[u.Atom] = atoms

        # Attributes set in get_args
        self.molecules_rdk: t.List[Chem.Mol]

    def get_args(self, input_function: t.Callable):
        """
        get_args(...) Needs to be called, before filter() can be used. It uses an input function to find all arguments required by a certain instance of a _Filter.

        For RingFilter that is a List of all Molecules in the pdb format (CAREFULL: The Molecules need to contain ALL atoms in the Molecules, not just the selected one.)


        :param input_function: A function that can get the input
        :type input_function: function(str)
        :param message: A string that is displayed if if input function needs to communicate with the user
        :type message: String
        :return: -
        :rtype: -
        :raises: BadArgumentException if the input function does not provide all arguments in the necessary format
        """

        input = input_function("A List of Molecules in the PDB FORMAT:")
        input = u.check_or_convert_argument(input, list)
        molecules_pdb = [u.check_or_convert_argument(i, str) for i in input]
        self.molecules_rdk = Rdkit_Functions.parse_pdb_blocks_to_rdkit(molecules_pdb)

        # TODO: What Kind of Error will I get if tools_Rdkit conversion fails?/ Catch it here

    def filter(self):
        '''
        filter should be overridden by every subclass of _Filter. Returns a filtered List of atoms

        In the case of Ring Filter: All atoms which are in rings
        :return: Filtered atom list
        :rtype: t.List[u.Atom]
        '''
        ids_of_atoms_in_rings = Rdkit_Functions.ring_atom_filter(selected=None, mols=self.molecules_rdk,
                                                                 selected_mols=None)
        return [a for a in self.atoms if a.id in ids_of_atoms_in_rings]

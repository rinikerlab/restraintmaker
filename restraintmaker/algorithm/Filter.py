"""
The module Filters contains functionality of filtering selections by criteria like in ring structure or element.

"""

import typing as t
from rdkit import Chem

from restraintmaker.tools_Rdkit import Rdkit_Functions
from restraintmaker.utils import Utilities as u


class _Filter():

    def __init__(self, atoms: t.List[u.Atom]):
        """_Filter
        This is the private parent class to all filter Classes.

        Parameters
        ----------
        atoms : t.List[u.Atom]
            List of atoms to be filtered
        args : t..List[str]
            Extra Arguments the subclasses might need

        """
        self.atoms: t.List[u.Atom] = atoms

    def get_args(self, input_function: t.Callable):
        """
            get_args(...) should be called by every subclass of _Filter, right after it is created.

            It uses an input function to find all arguments required by a certain instance of a _Filter

        Parameters
        ----------
        input_function : function(str)
            A function that can get the input
        message : str
            A string that is displayed if if input function needs to communicate with the user

        Returns
        -------
        NoReturn

        Raises
        ------
        BadArgumentException
             If the input function does not provide all arguments in the necessary format
        NotImplementedError
             If _Filter.get_args is called directly.
        """
        raise NotImplementedError("Call of function get_args(...) of abstract parent class _Filter")

    def filter(self) -> t.List[u.Atom]:
        """
            filter must be overridden by every subclass of _Filter. Returns a filtered list of atoms

        Returns
        -------
        t.List[u.Atom]
            Filtered List of Atoms
        Raises
        ------
        NotImplementedErrror
            If _Filter.filter is called directly
        """

        raise NotImplementedError("Call of function filter() of abstract parent class _Filter")


class PropertyFilter(_Filter):

    def __init__(self, atoms: t.List[u.Atom]):
        """
        Filters by any attribute of atom

        Parameters
        ----------
        atoms : t.List[u.Atom]
            List of atoms that are filtered over
        """

        super().__init__(atoms)

        # Attributes that will be set in get_args
        self.criterion: str = None  # Filter criterion
        self.value: t.Any = None  # Desired Value. Is introduced as string and then cast to datatype of the chosen criterion

    def get_args(self, input_function: t.Callable):
        """
            get_args(...) Needs to be called, before filter() can be used. For the Property_Filter it finds out which
                property should be filter-criterion and then which value it should have.

            Todo: think about multiple property selections

        Parameters
        ----------
        input_function : function(str)
            A function that will find the input needed ofr the filter.
        message : str
            A string that is displayed if if input function needs to communicate with the user

        Returns
        -------
        NoReturn

        Raises
        ------
        BadArgumentException
            if the input function does not provide all arguments in the necessary format

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
        """
                Returns a list of atoms, filtered using the provided criterion and value

        Returns
        -------
        t.List[u.Atom]
            Filtered atom list
        """

        return list(filter(lambda x: x.__getattribute__(self.criterion) == self.value, self.atoms))


class ElementFilter(_Filter):
    def __init__(self, atoms: t.List[u.Atom]):
        """
            Filter for certain Element(s)

        Parameters
        ----------
        atoms: t.List[u.Atom]
            List of atoms that are filtered over
        """
        super().__init__(atoms)

    def get_args(self, input_function: t.Callable = lambda *args: 'C'):
        """
            get_args(...) Needs to be called, before filter() can be used.
             For the Element_Filter it finds out which element should be filtered

        Parameters
        ----------
        input_function: function(str)
            A function that will find the input needed ofr the filter.
        message: str
            A string that is displayed if if input function needs to communicate with the user

        Returns
        -------
        NoReturn

        """

        input: str = input_function(
            "Which Element do do you want to filter for?\n (select multiple elements seperated by , ; e.g. 'O, N')")

        if (len(input) > 1):
            input = list(map(lambda x: x.strip(), input.split(",")))
        else:
            input = [input]

        self.in_path = u.check_or_convert_argument(input, str)
        self.elem = input

    def filter(self) -> t.List[u.Atom]:
        """
            Returns a filtered version of the class attribute atoms, filtered by the class attribute element

        Returns
        -------
        t.List[u.Atom]
            Filtered atom list
        """

        # Use the built in python function filter(expression, iterable)
        return list(filter(lambda x: x.elem in self.elem, self.atoms))


class RingFilter(_Filter):

    def __init__(self, atoms: t.List[u.Atom]):
        """
            Ring Filter chooses all Atoms that are part of a Ring
            Warnings: using tools_Rdkit. Needs conversion to tools_Rdkit molecules. => Needs the pdb as input

        Parameters
        ----------
        atoms : t.List[u.Atom]
            List of atoms to be
        args : t.List[str]
            Extra Arguments the subclasses might need
        """
        self.atoms: t.List[u.Atom] = atoms

        # Attributes set in get_args
        self.molecules_rdk: t.List[Chem.Mol]

    def get_args(self, input_function: t.Callable):
        """
            get_args(...) Needs to be called, before filter() can be used. It uses an input function to find all arguments required by a certain instance of a _Filter.

            For RingFilter that is a List of all Molecules in the pdb format (CAREFULL: The Molecules need to contain ALL atoms in the Molecules, not just the selected one.)

            TODO: What Kind of Error will I get if tools_Rdkit conversion fails?/ Catch it here

        Parameters
        ----------
        input_function: function(str)
            A function that can get the input
        message: str
            A string that is displayed if if input function needs to communicate with the user

        Returns
        -------
        NoReturn

        Raises
        ------
        BadArgumentException
            if the input function does not provide all arguments in the necessary format

        """

        input = input_function("A List of Molecules in the PDB FORMAT:")
        input = u.check_or_convert_argument(input, list)
        molecules_pdb = [u.check_or_convert_argument(i, str) for i in input]
        self.molecules_rdk = Rdkit_Functions.parse_pdb_blocks_to_rdkit(molecules_pdb)

    def filter(self):
        """
                filter should be overridden by every subclass of _Filter. Returns a filtered List of atoms
        In the case of Ring Filter: All atoms which are in rings

        Returns
        -------
        t.List[u.Atom]
            Filtered atom list
        """

        ids_of_atoms_in_rings = Rdkit_Functions.ring_atom_filter(selected=None, mols=self.molecules_rdk,
                                                                 selected_mols=None)
        print(ids_of_atoms_in_rings)
        return [a for a in self.atoms if a.id in ids_of_atoms_in_rings]

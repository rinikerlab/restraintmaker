"""
    The module Selections contains the selection logics. A selected atom will be considered for restraint placement.
"""

import math
import typing as t

from restraintmaker.tools_Rdkit import Rdkit_Functions
from restraintmaker.utils import Utilities as u
# from restraintmaker.utils.Utilities import print


class _Selection():

    def __init__(self, all_atoms: t.List[u.Atom]):
        """Selection
        This is the private parent class to all selection Classes.
        A selection contains a list of Atoms, and provides functions for their manipulation

        Parameters
        ----------
        all_atoms : t.List[u.Atom]
            all possible selectable atoms

        """
        self.atoms = []
        self.has_finished = False
        # If we keep it this way, we could think about just passsing indices, instead of the atoms themselves in _update_selected
        self.all_atoms = all_atoms

    def get_args(self, input_function: t.Callable) -> bool:
        """get_args(...)
        should be overridden by every subclass of _Selection. It will assign all necessart variables using input function.


        Parameters
        ----------
        input_function: t.Callable
            a function that will provide the argumetns for the selection in the necessary format
        Returns
        -------
        noReturn

        Raises
        -------
         BadArgumentException
            if the input function does not provide all arguments in the necessary format
         NotImplementedError
            If _Selection.get_args is called directly.
        """

        raise (NotImplementedError("Direct call of function get_args of abstract parent class _Selection"))

    # All update functions. They will be called by the GUI-program (restraintmaker_PyMOL) to inform the selection of user actions.
    # The different update functions here are formulated as generally as possible. It is up to the GUI-Program to decide exactly which is called when
    # Instead of explicit event handling the Selections can just override all _update_functions corresponding to events they want to be notifies for

    # ATOMS_SELECTED: Atom marked/highlighted/clicked in GUI program
    #   Recommended GUI-event: Click ON an Atom/GUI specific pick/select/add event
    # TODO REORG: Only allow one atom to be selected at a time?
    def _update_select(self, new_atoms: t.List[u.Atom]):
        """
                _update_select should be overridden by every subclass of _Selection that wants to be notified when new atoms have been selected/marked in the GUI

        Parameters
        ----------
        new_atoms : t.List[u.Atom]
             List of atom selection that has been added since last updated

        Returns
        -------

        """

        pass

    # SIZE_CHANGE: The user indicated that they want to change the size or increase something
    #   Recommended GUI-event: Keyboard arrow / mouse wheel
    def _update_size(self, increase: bool):
        """
             _update_size should be overridden by every subclass of _Selection who has some value, such as a size, that can be changed by the user during the selection.
            Warnings: This method merely shows that the user WANTS to change the size. The actual changes are made by the selection itself.

        Parameters
        ----------
        increase: bool
             Should the size change or not?

        Returns
        -------
        NoReturn

        """

        pass

    # MOVE: The user moved something
    #   Recommended GUI-event: Keyboard arrows / mouse moved
    def _update_move(self, new_x: float, new_y: float, new_z: float):
        """
             _update_move should be overridden by every subclass of _Selection that diplays a postion or similar, that the user can change.

        Parameters
        ----------
        new_x : float
            New x Position
        new_y : float
            New y Position
        new_z : float
            New z Position

        Returns
        -------
        NoReturn

        """

        pass

    # CONFIRM: The user indicated, that he wants to accept the current state of the selection
    #   Recommended GUI Events: Enter, Double-Click
    def _update_confirm(self):
        """
            should be overriden by everz subclass of _Selection that needs some kind of confirmation before it finsihes.

        Returns
        -------
        NoReturn

        """

    def reset(self):
        """
            reset resets the atom list to [] and has_finished to false.
        Usefull to reuse the same type of selection from the GUI program.
        This leaves the arguments particular to each Selection Type unchanged  => This way we do not have to create a new selection and call get_arguments again.
        Subclasses can override reset() if some parameters should be kept.

        Returns
        -------
        NoReturn
        """

        self.atoms = []
        self.has_finished = False


class LimitedSelection(_Selection):
    def __init__(self, all_atoms: t.List[u.Atom]):
        """
            A selection containing a fix maximal number of atoms

        Parameters
        ----------
        all_atoms : t.List[u.Atom]
            all possible selectable atoms
        """

        super().__init__(all_atoms)

    def get_args(self, input_function: t.Callable):
        """
            should be overridden by every subclass of _Selection. It will assign all necessart variables using input function.
            For Limited_Filter those are: max_size: int - Maximal number of atoms in selection.

        Parameters
        ----------
        input_function: t.Callable(str)
             a function that will provide the arguments for the selection in the necessary format: int

        Returns
        -------
        NoReturn

        Raises
        ------
        BadArgumentException
            if the input function does not provide all arguments in the necessary format
        """

        input = input_function("Maximal Number of Atoms in Limited Selection (int):")
        self.max_size = u.check_or_convert_argument(input, int)

    def _update_select(self, new_atoms: t.List[u.Atom]):
        """
            _update_select should be overridden by every subclass of _Selection that wants to be notified when new atoms have been selected/marked in the GUI
            LimitedSelection.update () only accepts one molecule in new_atoms at a time.

        Parameters
        ----------
        new_atoms : t.List
            List of atom selection that has been added since last updated

        Returns
        -------
        NoReturn

        """

        if self.has_finished:
            return

        if (len(new_atoms) > 1):
            print(
                'WARNING: LimitedSelection.update(...) has been called with more than one atom in its new_atoms List.',
                mv=4)
            print("   Only the fist Atom of the List is added: ", new_atoms[0], mv=4)

        new_atom = new_atoms[0]

        if new_atom in self.atoms:
            self.atoms.remove(new_atom)
        else:
            self.atoms.append(new_atom)
            self.has_finished = len(self.atoms) >= self.max_size


class SingleAtomSelection(LimitedSelection):
    def __init__(self, all_atoms: t.List[u.Atom]):
        """
            A selection containing only one atom. A Shortcut to create a Limited_Selection with max_atoms+1
            This is usefull for GUI programs, that use refresh, once a selection is full: Just allow to pick atoms.

        Parameters
        ----------
        all_atoms
        """
        super().__init__(all_atoms)

    def get_args(self, input_function: t.Callable):
        """
        should be overridden by every subclass of _Selection. It will assign all necessart variables using input function.
        For Pairwise_Filter there is no additional arguments. get_args(...) still needs to be called because it sets the fixed arguments directly.
        Parameters
        ----------
        input_function : t.Callable(str)
            a function that will provide the arguments for the selection in the necessary format. Here: Not necessesary. Just pass an 'empty' function or any function at all. It will not be called.

        Returns
        -------
        NoReturn

        Raises
        ------
        BadArgumentException
            if the input function does not provide all arguments in the necessary format
        """

        super().get_args(lambda *args: 1)


class PairSelection(LimitedSelection):
    def __init__(self, all_atoms: t.List[u.Atom]):
        """
            A Selection containing up to 2 Atoms. A Shortcut to create a Limited_Selection with max_atoms = 2

        Parameters
        ----------
        all_atoms : t.List[u.Atom]
            all possible selectable atoms

        """
        super().__init__(all_atoms)

    def get_args(self, input_function: t.Callable):
        """
        should be overridden by every subclass of _Selection. It will assign all necessart variables using input function.

        For Pairwise_Filter there is no additional arguments. get_args(...) still needs to be called because it sets the fixed arguments directly.

        Parameters
        ----------
        input_function : t.Callable (str)
            a function that will provide the arguments for the selection in the necessary format. Here: Not necessesary. Just pass an 'empty' function or any function at all. It will not be called.

        Returns
        -------
        NoReturn

        Raises
        ------
        BadArgumentException
            if the input function does not provide all arguments in the necessary format
        """

        super().get_args(lambda *args: 2)


class UniversalSelection(_Selection):

    def __init__(self, all_atoms: t.List[u.Atom]):
        """
            Universial Selections can be used to select all present atoms with one click. (Actually it 2 Clicks)

        Parameters
        ----------
        all_atoms : t.List[u.Atom]
            all possible selectable atoms
        """
        super().__init__(all_atoms)

    def get_args(self, input_function: t.Callable):
        """
            should be overridden by every subclass of _Selection. It will assign all necessary variables using input function.
            In Universal selection, the selection will be set directly by this function to all atoms.

        Parameters
        ----------
        input_function : t.Callable (str)
             a function that will provide the arguments for the selection in the necessary format. Here: Not necessesary. Just pass an 'empty' function or any function at all. It will not be called.

        Returns
        -------
        NoReturn

        Raises
        ------
        BadArgumentException
            if the input function does not provide all arguments in the necessary format

        """
        self.update()

    def update(self):
        """
            will be called by get_args, to initially set all atoms as selected. => The Selection will start the first time it will be checked

        Returns
        -------
        NoReturn
        """

        self.atoms = self.all_atoms
        self.has_finished = True

    def _update_select(self, new_atoms: t.List[u.Atom]):
        """
            Dummy no selection possible

        Parameters
        ----------
        new_atoms : t.List[u.Atom]
            not used here, because all atoms are selected

        Returns
        -------
        NoReturn
        """
        pass

    def _update_move(self, new_x: float, new_y: float, new_z: float):
        """
            Dummy no selection possible

        Parameters
        ----------
        new_atoms : t.List[u.Atom]
            not used here, because all atoms are selected

        Returns
        -------
        NoReturn
        """
        pass

    def _update_size(self, increase: bool):
        """
            Dummy no selection possible

        Parameters
        ----------
        new_atoms : t.List[u.Atom]
            not used here, because all atoms are selected

        Returns
        -------
        NoReturn
        """
        pass

    def _update_confirm(self):
        """
            Dummy no selection possible

        Parameters
        ----------
        new_atoms : t.List[u.Atom]
            not used here, because all atoms are selected

        Returns
        -------
        NoReturn
        """
        pass


class SphericalSelection(_Selection):
    """..class SphericalSelection:"""

    def __init__(self, all_atoms: t.List[u.Atom]):
        """
             A spherical selection selects all atoms within a given sphere. The Handling Program (restraintmaker_PyMol needs to provide a function that displays and chagnes the sphere

        Parameters
        ----------
        all_atoms : t.List[u.Atom]
            all possible selectable atoms

        """
        super().__init__(all_atoms)

        # atrributes that will be set in get_args(...)

    def get_args(self, input_function: t.Callable[[str or float], str]):
        """
            should be overridden by every subclass of _Selection. It will assign all necessart variables using input function.

        Parameters
        ----------
        input_function : t.Callable (str)
             a function that will provide the arguments for the selection in the necessary format. Here: Not necessesary. Just pass an 'empty' function or any function at all. It will not be called.

        Returns
        -------
        NoReturn

        Raises
        ------
        BadArgumentException
             if the input function does not provide all arguments in the necessary format
        """

        # attributes that will be changed by the GUI
        self.radius: float = 3  # Standard length. => Handling program needs to adjust it
        self.x = 0
        self.y = 0
        self.z = 0
        self.zoom_factor = 1.05

        pass

    def _update_size(self, increase: bool):
        """
        should be overridden by every subclass of _Selection who has some value, such as a size, that can be changed by the user during the selection.

         Carefull: This method merely shows that the user WANTS to change the size. The actual changes are made by the selection itself.
         For SphericalSelection this will change the radius of the spehere

        Parameters
        ----------
        increase : bool
            should the radius increase(True) or decrease (False)

        Returns
        -------
        NoReturn

        """

        self.radius *= self.zoom_factor if increase else 1 / self.zoom_factor
        self.select_atoms_within_sphere()

    def _update_move(self, new_x: float, new_y: float, new_z: float):
        """
         should be overridden by every subclass of _Selection that diplays a postion or similar, that the user can change.

        Parameters
        ----------
        new_x : float
            New x Position
        new_y : float
            New y Position
        new_z : float
            New z Position

        Returns
        -------

        """

        self.x = new_x
        self.y = new_y
        self.z = new_z
        self.select_atoms_within_sphere()

    def _update_confirm(self):
        """
        should be overriden by everz subclass of _Selection that needs some kind of confirmation before it finsihes.
        For Spherical_Selection, this indicates, that the user has moved to Sphere to where he wants it.

        Returns
        -------
        NoReturn

        """

        self.has_finished = True

    def atom_within_sphere(self, a: u.Atom) -> bool:
        """
            get all atoms within the sphere

        Parameters
        ----------
        a: u.Atom
            the atom, to be checked if it is inside the selection sphere.

        Returns
        -------
        bool:
            atom is in sphere?

        """
        dis = math.sqrt(math.pow(a.x - self.x, 2) + math.pow(a.y - self.y, 2) + math.pow(a.z - self.z, 2))
        return dis <= self.radius

    def select_atoms_within_sphere(self):
        """
        update selection, select all atoms inside the sphere

        Returns
        -------
        NoReturn

        """

        self.atoms = list(filter(self.atom_within_sphere, self.all_atoms))


class PaintSelection(SphericalSelection):

    def __init__(self, all_atoms: t.List[u.Atom]):
        """
            Works like a Spherical selection, but does keep Atoms Selected, once the Selection Spher has touched them

        Parameters
        ----------
        all_atoms : t.List[u.Atom]
            all possible selectable atoms
        """

        super().__init__(all_atoms)

    # get args inherited from Spherical selection
    # update functions inherited from Spherical Selection

    def select_atoms_within_sphere(self):
        """
            get all atoms within the sphere
            Override select_atoms_within_sphere, so it will keep old atoms
            This works, because when the inherited update functions are called, they will still call the overriden select_atoms_within_sphere

        Parameters
        ----------
        a: u.Atom
            the atom, to be checked if it is inside the selection sphere.

        Returns
        -------
        bool:
            atom is in sphere?

        """
        atoms_in_sphere = list(filter(self.atom_within_sphere, self.all_atoms))
        # Filter out doubles
        self.atoms += list(filter(lambda a: not a in self.atoms, atoms_in_sphere))


class MCS_Selection(_Selection):

    def __init__(self, all_atoms):
        """
             Selects all atoms, which are in the MCS of at least on pair of molecules

        Parameters
        ----------
        all_atoms : t.List[u.Atom]
            all possible selectable atoms
        """
        super().__init__(all_atoms)
        # Attributes set in get args:
        self.rdkit_molecules = None

    def get_args(self, input_function: t.Callable):
        """
            should be overridden by every subclass of _Selection. It will assign all necessart variables using input function.
            For MCS_Selection the argument is a list of the Molecules in the pdb format.

        Parameters
        ----------
        input_function : t.Callable (str)
            a function that will provide the arguments for the selection in the necessary format. Here: Not necessesary. Just pass an 'empty' function or any function at all. It will not be called.

        Returns
        -------
        NoReturn

        Raises
        ------
        BadArgumentException:
            if the input function does not provide all arguments in the necessary format

        """

        input = input_function("A List of Molecules in the PDB FORMAT:")
        input = u.check_or_convert_argument(input, list)
        molecules_pdb = [u.check_or_convert_argument(i, str) for i in input]
        self.molecules_rdk = Rdkit_Functions.parse_pdb_blocks_to_rdkit(molecules_pdb)

        self.update()

    def update(self):
        """
         general update function

        Returns
        -------
        NoReturn

        """
        mcs_ids = Rdkit_Functions.mcs_selection(self.molecules_rdk)
        self.atoms = list(filter(lambda a: a.id in mcs_ids, self.all_atoms))
        self.has_finished = True

    def _update_select(self, new_atoms: t.List[u.Atom]):
        """
            Dummy no selection possible

        Parameters
        ----------
        new_atoms : t.List[u.Atom]
            not used here, because all atoms are selected

        Returns
        -------
        NoReturn
        """
        pass

    def _update_move(self, new_x: float, new_y: float, new_z: float):
        """
            Dummy no selection possible

        Parameters
        ----------
        new_atoms : t.List[u.Atom]
            not used here, because all atoms are selected

        Returns
        -------
        NoReturn
        """
        pass

    def _update_size(self, increase: bool):
        """
            Dummy no selection possible

        Parameters
        ----------
        new_atoms : t.List[u.Atom]
            not used here, because all atoms are selected

        Returns
        -------
        NoReturn
        """
        pass

    def _update_confirm(self):
        """
            Dummy no selection possible

        Parameters
        ----------
        new_atoms : t.List[u.Atom]
            not used here, because all atoms are selected

        Returns
        -------
        NoReturn
        """
        pass

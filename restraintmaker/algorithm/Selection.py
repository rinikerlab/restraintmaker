"""
    The module Selections contains the selectoin logics.
"""

import math
import typing as t

from restraintmaker.tools_Rdkit import Rdkit_Functions
from restraintmaker.utils import Utilities as u
from restraintmaker.utils.Utilities import print


class _Selection():
    """
    ..class:    Selection
        This is the private parent class to all selection Classes.
        A selection contains a list of Atoms, and provides functions for their manipulation
    """

    def __init__(self, all_atoms: t.List[u.Atom]):
        self.atoms = []
        self.has_finished = False
        # TODO STRUC: Consider if it is realle necessary to give all_atoms as argument to _Selection. Only Spherical selection, really needs them.
        # If we keep it this way, we could think about just passsing indices, instead of the atoms themselves in _update_selected
        self.all_atoms = all_atoms

    def get_args(self, input_function: t.Callable) -> bool:
        """
        get_args(...) should be overridden by every subclass of _Selection. It will assign all necessart variables using input function.

        :param input_function: a function that will provide the argumetns for the selection in the necessary format
        :type input_function: t.Callable
        :return: -
        :rtype: -
        :raises BadArgumentException: if the input function does not provide all arguments in the necessary format
        :raises NotImplementedError: If _Selection.get_args is called directly.
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

        :param new_atoms: List of atom selection that has been added since last updated
        :type new_atoms: t.List
        :param display_function:
        :type display_function:
        :return: -
        :rtype: -
        """
        pass

    # SIZE_CHANGE: The user indicated that they want to change the size or increase something
    #   Recommended GUI-event: Keyboard arrow / mouse wheel
    def _update_size(self, increase: bool):
        """
         _update_size should be overridden by every subclass of _Selection who has some value, such as a size, that can be changed by the user during the selection.

        :warning: This method merely shows that the user WANTS to change the size. The actual changes are made by the selection itself.
        :param increase:  Should the size change or not?
        :type bool
        :return: -
        :rtype: -
        """
        pass

    # MOVE: The user moved something
    #   Recommended GUI-event: Keyboard arrows / mouse moved
    def _update_move(self, new_x: float, new_y: float, new_z: float):
        '''
         _update_move should be overridden by every subclass of _Selection that diplays a postion or similar, that the user can change.
        :param new_x: New x Position
        :type new_x: float
        :param new_y: New y Position
        :type new_y: float
        :param new_z: New x Position
        :type new_z: float
        :return: -
        :rtype: None
        '''
        pass

    # CONFIRM: The user indicated, that he wants to accept the current state of the selection
    #   Recommended GUI Events: Enter, Double-Click
    def _update_confirm(self):
        '''
        _update_confirm(...) should be overriden by everz subclass of _Selection that needs some kind of confirmation before it finsihes.
        :return: -
        :rtype: None
        '''

    def reset(self):
        """
        reset resets the atom list to [] and has_finished to false.

        Usefull to reuse the same type of selection from the GUI program.
        This leaves the arguments particular to each Selection Type unchanged  => This way we do not have to create a new selection and call get_arguments again.
        Subclasses can override reset() if some parameters should be kept.
        :return: -
        :rtype: -
        """
        self.atoms = []
        self.has_finished = False


class LimitedSelection(_Selection):
    """
    .. class Limited Selection: A selection containing a fix maximal number of atoms
    """

    def __init__(self, all_atoms: t.List[u.Atom]):
        """
        :param max_atoms: Maximal number of atoms in Selection
        :type max_atoms: int
        """
        super().__init__(all_atoms)

    def get_args(self, input_function: t.Callable):
        """
        get_args(...) should be overridden by every subclass of _Selection. It will assign all necessart variables using input function.
        For Limited_Filter those are: max_size: int - Maximal number of atoms in selection.
        :param input_function: a function that will provide the arguments for the selection in the necessary format: int
        :type input_function: t.Callable (str)
        :return: -
        :rtype: -
        :raises BadArgumentException: if the input function does not provide all arguments in the necessary format
        """

        input = input_function("Maximal Number of Atoms in Limited Selection (int):")
        self.max_size = u.check_or_convert_argument(input, int)

    # TODO: Change to onlz accepting a single atom nstead of a list
    def _update_select(self, new_atoms: t.List[u.Atom]):
        """
        _update_select should be overridden by every subclass of _Selection that wants to be notified when new atoms have been selected/marked in the GUI
        LimitedSelection.update () only accepts one molecule in new_atoms at a time.

        :param new_atoms: List of atom selection that has been added since last updated
        :type new_atoms: t.List
        :param display_function:
        :type display_function:
        :return: -
        :rtype: -
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
    """
    ..class Single Atom Selection: A selection containing only one atom. A Shortcut to create a Limited_Selection with max_atoms+1
    This is usefull for GUI programs, that use refresh, once a selection is full: Just allow to pick atoms.
    """

    def __init__(self, all_atoms: t.List[u.Atom]):
        super().__init__(all_atoms)

    def get_args(self, input_function: t.Callable):
        """
        get_args(...) should be overridden by every subclass of _Selection. It will assign all necessart variables using input function.
        For Pairwise_Filter there is no additional arguments. get_args(...) still needs to be called because it sets the fixed arguments directly.
        :param input_function: a function that will provide the arguments for the selection in the necessary format. Here: Not necessesary. Just pass an 'empty' function or any function at all. It will not be called.
        :type input_function: t.Callable (str)
        :return: -
        :rtype: -
        :raises BadArgumentException: if the input function does not provide all arguments in the necessary format

        """

        super().get_args(lambda *args: 1)

    def _update_select(self, new_atoms: t.List[u.Atom]):
        """
        _update_select should be overridden by every subclass of _Selection that wants to be notified when new atoms have been selected/marked in the GUI
        PairSelection.update(...) just calls super().update(...) (TODO CLEAN: In that case I can remove it during the final cleanup)
        :param display_function: A function of the calling program that can be used to display the current selection
        :type display_function: t.Callable
        :return: -
        :rtype: -
        """
        super()._update_select(new_atoms)


class PairSelection(LimitedSelection):
    """
    A Selection containing up to 2 Atoms. A Shortcut to create a Limited_Selection with max_atoms = 2
    """

    def __init__(self, all_atoms: t.List[u.Atom]):
        super().__init__(all_atoms)

    def get_args(self, input_function: t.Callable):
        """
        get_args(...) should be overridden by every subclass of _Selection. It will assign all necessart variables using input function.

        For Pairwise_Filter there is no additional arguments. get_args(...) still needs to be called because it sets the fixed arguments directly.
        :param input_function: a function that will provide the arguments for the selection in the necessary format. Here: Not necessesary. Just pass an 'empty' function or any function at all. It will not be called.
        :type input_function: t.Callable (str)
        :return: -
        :rtype: -
        :raises BadArgumentException: if the input function does not provide all arguments in the necessary format
        """

        super().get_args(lambda *args: 2)

    # is_full is inherited from Limited_Selection

    # TODO STRUC: Either only accept a single atom as input, or chech here, that we do not receive to many atoms
    def _update_select(self, new_atoms: t.List[u.Atom]):
        """
        _update_select should be overridden by every subclass of _Selection that wants to be notified when new atoms have been selected/marked in the GUI
       
        PairSelection.update(...) just calls super().update(...) (TODO CLEAN: In that case I can remove it during the final cleanup)
        :param display_function: A function of the calling program that can be used to display the current selection
        :type display_function: t.Callable
        :return: -
        :rtype: -
        """
        super()._update_select(new_atoms)


class UniversalSelection(_Selection):
    '''Universial Selections can be used to select all present atoms with one click. (Actually it 2 Clicks)'''

    # TODO STRUCTURE: At the moment one still has to select an atom to get to update() via do_select(), for the selection to do its work
    #   Instead: Call update directly from get_args or introduce a update_creation method, which will be called when selection is created?
    def __init__(self, all_atoms: t.List[u.Atom]):
        super().__init__(all_atoms)

    def get_args(self, input_function: t.Callable):
        """
            get_args(...) should be overridden by every subclass of _Selection. It will assign all necessart variables using input function.

            For Universal_Selection there is no additional arguments.
            :param input_function: a function that will provide the arguments for the selection in the necessary format. Here: Not necessesary. Just pass an 'empty' function or any function at all. It will not be called.
            :type input_function: t.Callable (str)
            :return: -
            :rtype: -
            :raises BadArgumentException: if the input function does not provide all arguments in the necessary format
\
        """

        pass

    def update(self):
        """
        Universal_Selection.update(...) will be called by any other update Function. => The Selection will start the first time it will be checked
        :param display_function: A function of the calling program that can be used to display the current selection
        :type display_function: t.Callable
        :return: -
        :rtype: -
        """
        self.atoms = self.all_atoms
        self.has_finished = True

    def _update_select(self, new_atoms: t.List[u.Atom]): self.update()

    def _update_move(self, new_x: float, new_y: float, new_z: float): self.update()

    def _update_size(self, increase: bool):   self.update()

    def _update_confirm(self): self.update()


class SphericalSelection(_Selection):
    """..class SphericalSelection: A spherical selection selects all atoms within a given sphere. The Handling Program (restraintmaker_PyMol needs to provide a function that displays and chagnes the sphere"""

    def __init__(self, all_atoms: t.List[u.Atom]):
        super().__init__(all_atoms)

        # atrributes that will be set in get_args(...)

        # attributes that will be changed by the GUI
        self.radius: float = 5  # Standard length. => Handling program needs to adjust it
        self.x = 0
        self.y = 0
        self.z = 0

    def get_args(self, input_function: t.Callable[[str or float], str]):
        """
            get_args(...) should be overridden by every subclass of _Selection. It will assign all necessart variables using input function.

            For SphericalSelection there is no additional arguments.
            :param input_function: a function that will provide the arguments for the selection in the necessary format. Here: Not necessesary. Just pass an 'empty' function or any function at all. It will not be called.
            :type input_function: t.Callable (str)
            :return: -
            :rtype: -
            :raises BadArgumentException: if the input function does not provide all arguments in the necessary format
\
        """
        # TODO: Allow get args to set a starting postion and size, instead of hardcoding it.
        pass

    def _update_size(self, increase: bool):
        '''
        _update_size should be overridden by every subclass of _Selection who has some value, such as a size, that can be changed by the user during the selection.

         Carefull: This method merely shows that the user WANTS to change the size. The actual changes are made by the selection itself.
         For SphericalSelection this will change the radius of the spehere
        :param increase: should the radius increase(True) or decrease (False)
        :type increase: bool
        :return: -
        :rtype: -
        '''

        # TODO STRUC: Make zoom-factyor a class varaible, so it could in theory be changed from outside
        zoom_factor = 1.05
        self.radius *= zoom_factor if increase else 1 / zoom_factor
        self.select_atoms_within_sphere()

    def _update_move(self, new_x: float, new_y: float, new_z: float):
        '''
            _update_moved(...) should be overridden by every subclass of _Selection that diplays a postion or similar, that the user can change.
           :param new_x: New x Position
           :type new_x: float
           :param new_y: New y Position
           :type new_y: float
           :param new_z: New x Position
           :type new_z: float
           :return: -
           :rtype: -
           '''
        self.x = new_x
        self.y = new_y
        self.z = new_z
        self.select_atoms_within_sphere()

    def _update_confirm(self):
        '''
        _update_confirm(...) should be overriden by everz subclass of _Selection that needs some kind of confirmation before it finsihes.
        For Spherical_Selection, this indicates, that the user has moved to Sphere to where he wants it.
        :return: -
        :rtype: None
        '''
        self.has_finished = True

    def atom_within_sphere(self, a: u.Atom) -> bool:
        dis = math.sqrt(math.pow(a.x - self.x, 2) + math.pow(a.y - self.y, 2) + math.pow(a.z - self.z, 2))
        return dis <= self.radius

    def select_atoms_within_sphere(self):
        self.atoms = list(filter(self.atom_within_sphere, self.all_atoms))


class PaintSelection(SphericalSelection):
    '''..class Paint Selection: Works like a Spherical selection, but does keep Atoms Selected, once the Selection Spher has touched them'''

    def __init__(self, all_atoms: t.List[u.Atom]):
        super().__init__(all_atoms)

    # get args inherited from Spherical selection
    # update functions inherited from Spherical Selection

    # Override select_atoms_within_sphere, so it will keep old atoms
    # This works, because when the inherited update functions are called, they will still call the overriden select_atoms_within_sphere
    def select_atoms_within_sphere(self):
        atoms_in_sphere = list(filter(self.atom_within_sphere, self.all_atoms))
        # Filter out doulbes
        self.atoms += list(filter(lambda a: not a in self.atoms, atoms_in_sphere))


class MCS_Selection(_Selection):
    '''..class MCS Selection: Selects all atoms, which are in the MCS of at least on pair of molecules'''

    def __init__(self, all_atoms):
        super().__init__(all_atoms)
        # Attributes set in get args:
        self.rdkit_molecules = None

    def get_args(self, input_function: t.Callable):
        """
        get_args(...) should be overridden by every subclass of _Selection. It will assign all necessart variables using input function.

        For MCS_Selection the argument is a list of the Molecules in the pdb format.
        :param input_function: a function that will provide the arguments for the selection in the necessary format. Here: Not necessesary. Just pass an 'empty' function or any function at all. It will not be called.
        :type input_function: t.Callable (str)
        :return: -
        :rtype: -
        :raises BadArgumentException: if the input function does not provide all arguments in the necessary format
        \
        """
        input = input_function("A List of Molecules in the PDB FORMAT:")
        input = u.check_or_convert_argument(input, list)
        molecules_pdb = [u.check_or_convert_argument(i, str) for i in input]
        self.molecules_rdk = Rdkit_Functions.parse_pdb_blocks_to_rdkit(molecules_pdb)

    def _update(self):
        mcs_ids = Rdkit_Functions.mcs_selection(self.molecules_rdk)
        self.atoms = list(filter(lambda a: a.id in mcs_ids, self.all_atoms))
        self.has_finished = True

    # TODO: CREATION EVENT OR CAL UPDATE FROM GET ARGS. See also Universal Selection
    def _update_size(self, increase: bool): self._update()

    def _update_select(self, new_atoms: t.List[u.Atom]): self._update()

    def _update_move(self, new_x: float, new_y: float, new_z: float): self._update()

    def _update_confirm(self): self._update()

"""
.. automodule:: Restraint
    :members:
Defines how selected or optimized atoms are linked together.

"""
import os
import sys

sys.path.append(os.path.dirname(__file__))

import typing as t
from restraintmaker.algorithm import Selection
from restraintmaker.utils import Utilities as u


class _Restraint():
    """
    .. autoclass:: RestraintType
        This is the private parent class to all selection Classes.
    """

    # Selections types from whose atoms this Restraint can be cosnsturcted
    accepted_selection_types = []

    def __init__(self, atoms):
        '''
        :type atoms:
        '''
        self.atoms = atoms


# TODO: Rename to distance restraint
class Pair_Restraint(_Restraint):
    """
    .. autoclass:: Pair_Restraint
        This class is defining the rules for a Gromos Distance Restraint
    """

    # 'static' variables
    atom_limit: int = 2
    accepted_selection_types: t.List[type] = [Selection.PairSelection]

    # TODO: Do not initializw with list of len 2 but 2 single atoms
    def __init__(self, atoms: t.List[u.Atom]):
        '''
        :param atoms:
        :type atoms: t.List[u.Atom]
        '''
        super().__init__(atoms)
        # TODO: Attribute distance


class CoM_Restraint(_Restraint):
    # TODO Implement
    """
    .. autoclass:: CoM_Restraint
        This class is a type Class, defining the properties of CoM restraints
    """

    accepted_selection_types = [Selection.PairSelection, Selection.LimitedSelection, Selection.SphericalSelection,
                                Selection.PaintSelection, Selection.UniversalSelection]


class Position_restraint(_Restraint):
    # TODO: Implement
    """
    .. autoclass:: Position_restraint
        This is type class, that should predefine position_Restraint
    """
    accepted_selection_types = [Selection.SingleAtomSelection]

    atom_limit: int = 1

# Todo: implement Pairwise type and forward to Slection. (Filter and opt?)

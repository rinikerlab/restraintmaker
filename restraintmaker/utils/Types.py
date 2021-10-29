"""
    Defines how Restraints are selected or optimized from atoms.
    Todo: implement Pairwise type and forward to Slection. (Filter and opt?)
    Todo: Molecule class?
"""

import typing as t
import numpy as np

from restraintmaker.algorithm import Selection
from restraintmaker.utils.Utilities import Atom
#from restraintmaker.algorithm import Optimizer


class _Restraint():
    # Selections types from whose atoms this Restraint can be cosnsturcted
    accepted_selection_types = []
    #accepted_optimizer_types = []


    def __init__(self, atoms):
        """
            This is the private parent class to all selection Classes.

        Parameters
        ----------
        atoms
        """
        self.atoms = atoms

    def _calc_distance(self, atomA, atomB) -> float:
        """
            Calculate the distance between two atoms and their 3D positions.

        Parameters
        ----------
        atomA : u.Atom
        atomB : u.Atom


        Returns
        -------
        float
            returns the distance
        """
        return np.sqrt((atomA.x - atomB.x) ** 2 + (atomA.y - atomB.y) ** 2 + (atomA.z - atomB.z) ** 2)


class Position_restraint(_Restraint):
    accepted_selection_types: t.List[Selection._Selection] = [Selection.SingleAtomSelection, Selection.AllSelection, Selection.MCS_Selection, Selection.PaintSelection, Selection.SphericalSelection]
    accepted_imports: t.List = []
    #accepted_optimizer_types: t.List[Optimizer._Optimizer] = []
    optimizable: bool = False
    atom_limit: int = 1

    def __init__(self, atomA: Atom, reference_atom: Atom = None):
        """
            This is type class, that should predefine position_Restraint

        Parameters
        ----------
        atomA :  t.List[u.Atom]

        """

        self.atomA = atomA

        if (reference_atom is None):
            self._reference_atom = atomA
        else:
            self._reference_atom = reference_atom

        self._distance_to_reference_position = self._calc_distance(self.atomA, self.reference_atom)
        atoms = [atomA]
        super().__init__(atoms)

    @property
    def atomA(self) -> Atom:
        return self.atomA

    @property
    def reference_atom(self) -> Atom:
        return self._reference_atom

    @property
    def distance_to_reference_position(self) -> float:
        return self._distance_to_reference_position


class Distance_Restraint(_Restraint):
    # 'static' variables
    atom_limit: int = 2
    accepted_selection_types: t.List[Selection._Selection] = [Selection.PairSelection, Selection.AllSelection, Selection.MCS_Selection, Selection.SphericalSelection, Selection.PaintSelection]
    optimizable: bool = True #accepted_optimizer_types: t.List = [Optimizer.TreeHeuristicOptimizer]

    def __init__(self, atomA: Atom, atomB: Atom):
        """
            This class is defining the rules for a Gromos Distance Restraint

        Parameters
        ----------
        atoms :  t.List[u.Atom]
        """

        self._atomA = atomA
        self._atomB = atomB
        atoms = [self._atomA, self._atomB]

        self._distance = self._calc_distance(atomA, atomB)

        super().__init__(atoms)

    def __str__(self):
        return str(self.__class__.__name__)+"\t"+str(self._atomA.id)+ "\t"+str(self._atomB.id)+"\t"+str(self._distance)

    @property
    def atomA(self) -> Atom:
        return self._atomA

    @property
    def atomB(self) -> Atom:
        return self._atomB

    @property
    def distance(self) -> float:
        return self._distance


class __CoM_Restraint(_Restraint):
    #accepted_selection_types = [Selection.PairSelection, Selection.LimitedSelection, Selection.SphericalSelection,
    #                            Selection.PaintSelection, Selection.UniversalSelection]

    def __init__(self, atoms: t.List[Atom]):
        """
            This class is a type Class, defining the properties of CoM restraints
            TODO: Implement

        Parameters
        ----------
        atoms :  t.List[u.Atom]
        """

        raise NotImplemented("WIP")
        # super().__init__(atoms)

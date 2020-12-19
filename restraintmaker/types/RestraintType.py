"""
.. automodule:: Filters
    :members:
    TODO:write DOCU!

Defines how selected or optimized atoms are linked together.

"""


class _RestraintType():
    """
    ..autoclass:    RestraintType
        This is the private parent class to all selection Classes.
    """

    pass


class Pair_Restraint(_RestraintType):
    """
    ..autoclass:    Pair_Restraint
        This class is defining the rules for a
    """
    atom_limit: int = 2

    settings: dict
    k: float
    l: float

    pass


class CoM_Restraint(_RestraintType):
    """
    ..autoclass:    CoM_Restraint
        This class is a type Class, defining the properties of CoM restraints
    """

    atom_limit: int = 4


class Position_restraint(_RestraintType):
    """
    ..autoclass:    Position_restraint
        This is type class, that should predefine position_Restraint
    """

    atom_limit: int = 1

# Todo: implement Pairwise type and forward to Slection. (Filter and opt?)

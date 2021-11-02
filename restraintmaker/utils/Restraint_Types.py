import typing as t

from restraintmaker.algorithm import Selection, Optimizer
from restraintmaker.io import Importer, Exporter
from restraintmaker.utils import Restraints
from restraintmaker.utils.Utilities import get_all_subclasses



class _Restraint_Type():
    accepted_selection_types_res: t.List[Selection._Selection]
    accepted_selection_types_atoms: t.List[Selection._Selection]
    accepted_optimizer_types: t.List[Optimizer._Optimizer]
    accepted_exporter_types: t.List[Optimizer._Optimizer]
    accepted_importer_types: t.List[Optimizer._Optimizer]
    restraint: Restraints._Restraint



class PositionRestraintType(_Restraint_Type):
    accepted_selection_types_res: t.List[Selection._Selection] = [Selection.SingleAtomSelection, Selection.AllSelection, Selection.MCSSelection, Selection.SphericalSelection]
    accepted_selection_types_atoms: t.List[Selection._Selection]  = []
    accepted_optimizer_types: t.List[Optimizer._Optimizer] = []
    accepted_exporter_types: t.List[Exporter._Export_Position_Restraints] = get_all_subclasses(Exporter._Export_Position_Restraints)
    accepted_importer_types: t.List[Importer._Importer_Position_Restraints] = get_all_subclasses(Importer._Importer_Position_Restraints)
    restraint:type(Restraints.PositionRestraint) = Restraints.PositionRestraint
    atom_limit: int = 1


class DistanceRestraintType(_Restraint_Type):
    # 'static' variables
    accepted_selection_types_atoms: t.List[Selection._Selection] = [Selection.SingleAtomSelection, Selection.AllSelection, Selection.MCSSelection, Selection.SphericalSelection]
    accepted_selection_types_res: t.List[Selection._Selection] = [Selection.PairSelection]
    accepted_optimizer_types: t.List[Optimizer._Optimizer] = [Optimizer.GreedyGraphOptimizer, Optimizer.BruteForceRingOptimzer]
    accepted_exporter_types: t.List[Exporter._Export_Distance_Restraints] = get_all_subclasses(Exporter._Export_Distance_Restraints)
    accepted_importer_types: t.List[Importer._Importer_Distance_Restraints] = get_all_subclasses(Importer._Importer_Distance_Restraints)
    restraint:type(Restraints.DistanceRestraint) = Restraints.DistanceRestraint
    atom_limit: int = 2
"""
    The optimizers of this package are contained in this module. The optimizers can be used too....

"""

import typing as t
import numpy as np
from collections import namedtuple

from scipy.spatial import ConvexHull
from scipy.spatial.qhull import QhullError
from scipy.special import binom

from restraintmaker.tools_Rdkit import Rdkit_Functions as rdkit_functions
from restraintmaker.utils import Utilities as u, Types
from restraintmaker.utils.Utilities import print


class _Optimizer():
    """
    .. autoclass::   _Optimizer
    """

    def __init__(self, atoms: t.List[u.Atom]):
        """
            This is the private parent class to all Optimizer Classes.

        Parameters
        ----------
        atoms :  t.List[u.Atom]
            list of atoms to be considered

        """
        self.atoms = atoms

    def get_args(self):
        """
            to be called by every Optimizer after its creation, t It uses an input function to find all arguments required by a certain instance of a _Optimize

        Returns
        -------
        NoReturn

        """

        raise NotImplementedError("Direct call of abstract function  _Optimizer.get_args(...)")

    def make_restraints(self) -> t.List[Types._Restraint]:
        """
                make_restraints() return a List of restraints that have been optimally set, depending on the criteria defined by the specific Optimizer used.

        Returns
        -------
        t.List[RestraintType.RestraintType]
            The list of optimized Restraints

        Raises
        ------
        NoOptimalSolutionException
            if no Solution fulfilling all criteria can be found.
        """

        raise NotImplementedError("Direct call of abstract function  _Optimizer.make_restraints()")


class _MoleculeRingOptimizer(_Optimizer):
    '''Abstract parent class for all Molecule Ring Optimziers. Molecule Ring Optimziers set their restraints by finding the best set of restraints between pairs of Molecules. In the end everz molecule is connected to two other molecules'''

    def __init__(self, atoms):
        super().__init__(atoms)
        self.Molecules: t.List[t.List[u.Atom]] = u.order_atoms_by_molecule(
            self.atoms)  # Every Molecules is a list of atoms, with the Molecule attribute as key

        # get_args: inherited from _Optimizer

        # to be set in get_args of children:
        self.cutoff = -1
        self.arrange_algo = ''  # criterion for how molecules should be sorted

    def make_restraints(self) -> t.List[Types._Restraint]:
        """
        make_restraints() returns a List of restraints that have been optimally set, depending on the criteria defined by the specific Optimizer used:

        All MoleculeRingOptimizers inherit make_restraint form _Molecule RingOptimizer. Only the function connect_two_Molecules differs. Connects Molecules pairwise, forming a topological ring. (Not all Molecules are connected to each other. Each is only connected to 2 neighbours).
        The Molecules are connected in alphanumerical order of their name. No Optimization.
        The molecules are connected, by applying a Pairwise distance restraint to n pairs of Atoms.
        These Atoms are chosen heuristically by finding distance restraints that are as far apart as possilbe. (Non optimal solution)
        :return: A List of Restraints
        :rtype: t.list[RestraintType.PairRestraint]
        :raises: NoOptimalSolutionException if no Solution fulfilling all criteria can be found. For MoleculeRingOptimizer: Not enough restraints shorter than cutoff

        """
        restraints: t.List[Types._Restraint] = []

        if len(self.Molecules) < 2:
            raise u.BadArgumentException(
                'A Molecule Ring Optimizer needs at least 2 Molecules to connect')
        elif len(self.Molecules) == 2:
            restraints = self.connect_two_molecules(self.Molecules[0], self.Molecules[1], self.n)
        elif self.arrange_algo == 'None':  # Str 'None', not None:

            for i_m in range(len(self.Molecules) - 1):
                print()
                print('\n Connecting molecules #' + str(i_m + 1) + ' and #' + str(i_m + 2), mv=3)
                restraints.extend(self.connect_two_molecules(self.Molecules[i_m], self.Molecules[i_m + 1], self.n))
                # Connect last molecule to first
                # TODO: Parallellize: Every moleucle pair gets one thread

                # Connect last Molecules to first

            print('\n Connecting molecules #' + str(len(self.Molecules)) + ' and #' + str(1), mv=3)
            restraints.extend(self.connect_two_molecules(self.Molecules[-1], self.Molecules[0], self.n))
        # An Optimization algorithm has been given
        else:

            # 0) Choose the algo from a predefined list
            # Do not just use the value returned by connect_two_molecules, because I want to be able to use a different criterion here.
            # TODO: Now the scaling is automatically done by taking calculate_value of all atom positions. For other criteria, it might be better to define a function scale, along with
            # calculate_value

            # TODO: Move to get-args and make calculate_value an attribute
            calculate_value = u.do_nothing
            if self.arrange_algo == 'pca_2d':
                calculate_value = _calculate_value_unscaled_pca_2d
                # print('Warning: Optimization of Molecule pairs is done using unscaled pca values',mv=4)
            elif self.arrange_algo == 'convex_hull':
                calculate_value = _calculate_value_convex_hull
            else:
                raise u.BadArgumentException(
                    self.arrange_algo + 'is not a known algorithm to arrange the molecules')

            # 1) Calclulate ALL PAIRS and their value

            # Save restraints and the value of a molecule pair in a list
            # Also create a list indicting which molecule pair is at what index of the list, to save time later on

            restraints_list = []  # List containing all restraints of each molecule pair
            value_list = []  # List containing all values of molecules pairs
            molecule_pair_list = []  # Needed later. Construct it here, as we are doing the loop anywaz

            for i_m1 in range(len(self.Molecules) - 1):
                for i_m2 in range(i_m1 + 1, len(self.Molecules)):  # Cant use enumerate because i_m2 would start from 0
                    m1, m2 = self.Molecules[i_m1], self.Molecules[i_m2]
                    print('\n Connecting molecules #' + str(i_m1 + 1) + ' and #' + str(i_m2 + 1), mv=3)
                    try:
                        pairwise_restraints = self.connect_two_molecules(m1, m2, self.n)

                        atom_positions_of_both_molecules = [(a.x, a.y, a.z) for a in m1 + m2]
                        scaling_factor = 1 / calculate_value(atom_positions_of_both_molecules)
                        restr_pos = [find_middle(r) for r in pairwise_restraints]
                        unscaled_value = calculate_value(restr_pos)
                        value = unscaled_value * scaling_factor

                    except NoOptimalSolutionException:
                        value = -1
                        print('Warning: Failed to connect molecule #' + str(i_m1 + 1) + ' and #' + str(i_m2 + 1), mv=4)

                    if (value > 0):
                        restraints_list.append(pairwise_restraints)
                        value_list.append(value)
                        molecule_pair_list.append((i_m1, i_m2))

                    # print(i_m1,i_m2,mv=0)

            # 2) Build a Ring using modified Kruskal

            i_chosen_pairs = maximal_weight_ring(edge_list=molecule_pair_list, value_list=value_list,
                                                 n_nodes=len(self.Molecules))
            chosen_restraints = []
            for i in i_chosen_pairs:
                chosen_restraints.extend(restraints_list[i])

            return chosen_restraints

        return restraints

    def connect_two_molecules(self, m1: t.List[u.Atom], m2: t.List[
        u.Atom], n) -> t.List[Types._Restraint]:
        """     connect_two_molecules places n restraint between 2 molecules, using criteria depending on the Optimizer

               :param m1: Molecule 1
               :type m1:  t.List[u.Atom]
               :param m2: Molecule 2
               :type m2: t.List[u.Atom]
               :return: Restraints between those atoms
               :rtype: t.List[RestraintType._RestraintType]:
               :raises NotImplementedError if _MoleculeRingOptimizer.connect_two_molecules is called directly.
        """
        raise NotImplementedError(
            'Direct call of abstract parent function _MoleculeRingOptimizer.connect_two_molecules')


class TreeHeuristicOptimizer(_MoleculeRingOptimizer):
    """
    ..class MoleculeRingOptimizer_0_1. Connects Molecules pairwise, forming a topological ring. (Not all Molecules are connected to each other. Each is only connected to 2 neighbours).
    First the all Molecules pairs are connected. Then a ring is formed by using the best pairs.
    These restraints between two molecules are chosen heuristically by finding distance restraints that are as far apart as possilbe. (Non optimal solution)
    """

    # #TODO: Speed is far from optimal, because some of the taks done by the utility functions are redundant:
    #     1)  We do not need to construct the whole spanning tree. We could stop, after we have found best connections
    #     2)  The info about which atom belongs to which molecule is lost when we call get_all_short_connections(...), and get_all_short_connections checks it again
    #     3)  The info about connectiviy in the tree is lost when it is returned as a list as has to be tediouly reconstructed while walking the tree

    def __init__(self, atoms):
        super().__init__(atoms)

        # Attributes inherited from _MoleculeRingOptiimzer:
        # self.Molecules
        self.algo = None
        self.n: int = 0
        self.cutoff = 0

        # Check for duplicates:
        for a1 in self.atoms[:-1]:
            for a2 in self.atoms[self.atoms.index(a1) + 1:]:
                if a1.id == a2.id:
                    print('WARNING: Optimizer has duplicates in its atom list:  i=', self.atoms.index(a1))

        if len(self.Molecules) < 2:
            raise u.BadArgumentException(
                str(self.__class__) + " must be initialized with atoms from at least 2 Molecules.")

    # TODO STRUC: Only Optimizer uses argument unpacking after the input function. Exporters, Filters etc call the input_functions several times if the y need more than one arg
    # =>Unify that!
    def get_args(self, input_function: t.Callable):
        '''
               get_args needs to be called by every Optimizer after its creation, t It uses an input function to find all arguments required by a certain instance of a _Filter
              :param input_function: A function that can get the input
              :type input_function: t.Callable[str]
              :param message: A string that is displayed if input function needs to communicate with the user
              :type message: str
              :return: -
              :rtype: -
               :raises: BadArgumentException if the input function does not provide all arguments in the necessary format
        '''
        inputs = input_function(
            'List of all args for TreeHeuristicOptimizer')  # TODO: See TODO at u.create_multi_dialog
        self.n = u.check_or_convert_argument(inputs[0], int)
        if self.n < 2:
            raise u.BadArgumentException(
                "A MoleculeRingOptimizer_0_1 needs to set at least 2 connection per Molecules pair. ")
        self.cutoff = u.check_or_convert_argument(inputs[1], float)
        self.algo = u.check_or_convert_argument(inputs[2], str, ['prim', 'shortest', 'cog', "biased_avg"])
        self.arrange_algo = u.check_or_convert_argument(inputs[3], str, ['None', 'convex_hull', 'pca_2d'])

    def connect_two_molecules(self, m1: t.List[u.Atom], m2: t.List[
        u.Atom], n) -> t.List[Types._Restraint]:
        """
        connect_two_molecules places n restraint between 2 molecules, using criteria depending on the Optimizer

        For three heurstic Optimizer that is by greedly picking restraints according to a criterion depending on the already chosen restraints

           :param m1: Molecule 1
           :type m1:  t.List[u.Atom]
           :param m2: Molecule 2
           :type m2: t.List[u.Atom]
           :return: Restraints between those atoms
           :rtype: t.List[RestraintType._RestraintType]:
           :raises NotImplementedError if _MoleculeRingOptimizer.connect_two_molecules is called directly.

        """

        return maximal_spanning_tree_greedy(get_all_short_pair_restraints(m1 + m2, self.cutoff), n, self.algo)


class BruteForceRingOptimzer(_MoleculeRingOptimizer):
    '''A Molecule Ring Optimzier that connects TwoMolecules by finding Optimizing some criterion by brute force
    '''

    # TODO: Think: Could we use a version of volume optimizer that does not do pairwise connections?

    def __init__(self, atoms):
        super().__init__(atoms)
        self.Molecules: t.List[t.List[u.Atom]] = u.order_atoms_by_molecule(
            self.atoms)  # Every Molecules is a list of atoms, with the Molecule attribute as key

        # attributes to be set in get_args
        self.n = 0  # Number of restrained Atoms per molecule
        self.algo = ''

        # Check for duplicates:
        for a1 in self.atoms[:-1]:
            for a2 in self.atoms[self.atoms.index(a1) + 1:]:
                if a1.id == a2.id:
                    print('WARNING: Optimizer has duplicates in its atom list:  i=', self.atoms.index(a1))

        if len(self.Molecules) < 2:
            raise u.BadArgumentException(
                str(self.__class__) + " must be initialized with atoms from at least 2 Molecules.")

    # TODO STRUC:Either use argumetn unpacking for all get arg functions or for none
    def get_args(self, input_function: t.Callable):
        '''
               get_args needs to be called by every Optimizer after its creation, t It uses an input function to find all arguments required by a certain instance of a _Filter
              :param input_function: A function that can get the input
              :type input_function: t.Callable[str]
              :param message: A string that is displayed if input function needs to communicate with the user
              :type message: str
              :return: -
              :rtype: -
               :raises: BadArgumentException if the input function does not provide all arguments in the necessary format
        '''
        # TODO: Consitency Change the string for self.algo to pca_2d as well
        inputs = input_function('List of all args for BruteForceOptimizer')  # TODO: See TODO at u.create_multi_dialog
        self.n = u.check_or_convert_argument(inputs[0],
                                             int)  # TODO: Add condition for acceptable values, once check_or_convert accepts bool functions as cirterion
        self.cutoff = u.check_or_convert_argument(inputs[1], float)
        self.algo = u.check_or_convert_argument(inputs[2], str, ['convex_hull', 'pca'])
        self.arrange_algo = u.check_or_convert_argument(inputs[3], str, ['None', 'convex_hull', 'pca_2d'])

        if self.algo == 'convex_hull' and self.n < 4:
            raise u.BadArgumentException(
                "A VolumeOptimizer needs to set at least 4 connections per Molecules pair. ")

    def connect_two_molecules(self, m1: t.List[u.Atom], m2: t.List[
        u.Atom], n) -> t.List[Types._Restraint]:
        """
               connect_two_molecules places n restraint between 2 molecules, using criteria depending on the Optimizer

               For three BruteForeceRingOptimizer that is by evaluating a criterion for every possible set of rsetraints and choosing the best.
                  :param m1: Molecule 1
                  :type m1:  t.List[u.Atom]
                  :param m2: Molecule 2
                  :type m2: t.List[u.Atom]
                  :return: Restraints between those atoms
                  :rtype: t.List[RestraintType._RestraintType]:
                  :raises NotImplementedError if _MoleculeRingOptimizer.connect_two_molecules is called directly.

               """

        # TODO: Each function call will cost time => For the Brute force Optimzer it might be better to have a single class for each case, instead of different functions
        # 1) set the get value function and prepare the restraint positions (e.g ifthey need to be moved, normalized etc)

        # TODO IMPORTANT, PROGRAM STRUCTURE: To generalize to different restraints and more sophisticated criteria, whihc do not only depend on atom position: Give restraints a property position, whichis set at creation.
        # Then: Do not pass restraint positions but restraints themselves to calculate_value, which could then use the involved atoms as well.

        print(self.algo, mv=1)
        if self.algo == 'convex_hull':
            calculate_value = _calculate_value_convex_hull
        elif self.algo == 'pca':
            calculate_value = _calculate_value_unscaled_pca_2d
        else:
            raise u.BadArgumentException(
                self.algo + 'is not an known criterion for optimization.')

        potential_restraints = get_all_short_pair_restraints(m1 + m2, self.cutoff)
        if len(potential_restraints) < n:
            raise NoOptimalSolutionException('There are not enough restraints to connect the two molecules')

        # Iterate (recursively) over all combinations of self.n restraints from the potential restraints
        # TODO: Generalize the algorithm to maximize anything using lambda functions
        count_precision_errors = 0  # Precision Errors happen, when the points are to close to one plane for calculating QHull
        progres = 0
        progres_step = 100 / len(potential_restraints)

        def find_max_val_recursively(max_depth: int, current_depth: int = 0, restraint_indices=[]):
            # #print(restraint_indices)
            # if current_depth < max_depth-1:
            #     print(current_depth,end='-')
            #     if current_depth == max_depth - 2:
            #         print('---')

            # print(restraint_indices,mv=1)

            if current_depth == max_depth:

                # Find positions of the restraints from their indices and return the volume of their Convex hull
                points = [restraint_positions[i] for i in restraint_indices]
                try:
                    max_val = calculate_value(points)
                except QhullError:  # TODO: Catch corresponding math errors for the other calculate value function. consider creating a CalculationError that the differen functions can throw if necessarz. I think we save ime if we only have one try block instead oftwo neste ones
                    max_val, best_ir = -1, []
                    nonlocal count_precision_errors
                    count_precision_errors += 1

                return max_val, restraint_indices

            else:
                if current_depth == 1:
                    nonlocal progres
                    print('>{:3d}'.format(int(progres)) + '%', mv=3, end=" ")
                    progres += progres_step

                max_val = 0
                best_restraint_indices = []

                start = restraint_indices[-1] + 1 if len(restraint_indices) > 0 else 0
                stop = len(potential_restraints) - (max_depth - current_depth - 1)
                for i_r in range(start, stop):
                    new_val, new_best_restraint_indices = find_max_val_recursively(max_depth=max_depth,
                                                                                   current_depth=current_depth + 1,
                                                                                   restraint_indices=restraint_indices + [
                                                                                       i_r])

                    if new_val > max_val:
                        max_val = new_val
                        best_restraint_indices = new_best_restraint_indices

                return max_val, best_restraint_indices

        restraint_positions = [find_middle(r) for r in potential_restraints]

        maximal_value, best_indices = find_max_val_recursively(max_depth=self.n)
        print("", mv=3)

        print('Value of the Optimal Solution: ' + '{:05.2f}'.format(maximal_value), mv=3)
        print('Checked  ' + str(int(binom(len(potential_restraints), self.n)) - count_precision_errors) + ' of ' + str(
            int(binom(len(potential_restraints), self.n))) + ' sets of restraints.', mv=3)
        print('Encountered ' + str(count_precision_errors) + ' precision_errors', mv=3)

        return [potential_restraints[i_r] for i_r in best_indices]


# END OF class BruteForceRingOptimizer -------------------------------------------------------------

class MetaMoleculeRingOptimizer(_MoleculeRingOptimizer):
    '''class MetaMoleculeRingOptimizer: Meta Optimizer which compares the solution of different Optimizers for each Molecule pair  and picks the best

        Warning: As implemented now, it can NOT be generalized to other Optimizer types than RingOptimizers.
        Warning: The complete solution of the BestMoleculesRingOPtimizer is NOT equal to any of the Optimizers it test. It will pick the best solution for each pair of molecules, leading to a different ring than any of the single Optimizers
        Warning: As Implemented now, it can not deal with future Molecule Ring Optimizers, which might have additional arguments
    '''

    # TODO GENERALIZE: 1) Allow to indicate which Optimizers to include, instead of hardcoding all of them
    #                 2) Allow to give the arguments for each single Optimizer, instead of using the same ones for each Optimizer(except of course criterion/update function.)
    # Would be perfect if we could pass Optimizers as arguments
    def __init__(self, atoms):
        super().__init__(atoms)

        # Arguments set in get_args
        self.sub_optimizers: t.List[_MoleculeRingOptimizer] = []
        self.criterion: t.Callable[[t.Tuple[float]], float] = None

        # Arguments inherited from _MoleculeRingOptimzier, to be set in ger_args
        self.cutoff = -1
        self.arrange_algo = ''  # criterion for how molecules should be sorted
        self.n = ''  # Number of restraints

    def get_args(self, input_function: t.Callable):
        '''
        :return:
        :rtype:
        '''

        inputs = input_function('List of all args for BestMoleculeRingOptimizer')
        self.n = u.check_or_convert_argument(inputs[0], int)
        if self.n < 2:
            raise u.BadArgumentException(
                "A MoleculeRingOptimizer_0_1 needs to set at least 2 connection per Molecules pair. ")
        self.cutoff = u.check_or_convert_argument(inputs[1], float)
        self.criterion = u.check_or_convert_argument(inputs[2], str, ['convex_hull',
                                                                      'pca'])  # would be nice to pass a function instead of just a description here
        self.arrange_algo = u.check_or_convert_argument(inputs[3], str, ['None', 'convex_hull', 'pca_2d'])

        # Create sub optimizers
        tree_algos = ["prim", "cog", "shortest",
                      "biased_avg"]  # TODO: Remove Hardcoded algorithms and allow to pass a list of all algos of interest

        for algorithm in tree_algos:
            new_optimizer = TreeHeuristicOptimizer(self.atoms)
            new_optimizer.get_args(lambda _: (self.n, self.cutoff, algorithm,
                                              None))  # TODO: Allow to use different kinds of Molecules Ring Optimizers, does not have to be TreeOptimizer
            self.sub_optimizers.append(new_optimizer)

    def connect_two_molecules(self, m1: t.List[u.Atom], m2: t.List[
        u.Atom], n):
        # TODO:Move this block to get_args of _MoleculeRingOptimzer and make calculate value an atribute of MoelculeRing Optimzer.
        # At the moment we have to put this block into make restraints and connect two molecules
        if self.arrange_algo == 'pca_2d':
            calculate_value = _calculate_value_unscaled_pca_2d
            # print('Warning: Optimization of Molecule pairs is done using unscaled pca values',mv=4)
        elif self.arrange_algo == 'convex_hull':
            calculate_value = _calculate_value_convex_hull
        else:
            raise u.BadArgumentException(
                self.arrange_algo + 'is not a known algorithm to arrange the molecules')

        best_restraint_set = None
        value_of_best_set = -1
        best_algorithm = ""
        for sub_optimizer in self.sub_optimizers:
            new_restraint_set = sub_optimizer.connect_two_molecules(m1, m2, n)
            new_value = calculate_value([find_middle(r) for r in new_restraint_set])

            if new_value > value_of_best_set:
                value_of_best_set = new_value
                best_restraint_set = new_restraint_set
                best_algorithm = sub_optimizer.algo

        print("Best algorithm: ", best_algorithm, mv=3)
        return best_restraint_set


# END OF class BestMoleculeRingOptimizer -------------------------------------------------------------


# --------------------Static Utilitiy functions ->

# TODO: change connections to Pairwise Restraints
def _get_all_short_connections(atoms: t.List[u.Atom], max_dis: float) -> t.List[
    Types.Distance_Restraint]:
    """
    Old. Only here to keep get_maximal_tree_kruskal_running. Delete when not needed anymore
    get_all_short_connections(...) finds all pairs of atoms within a certain cutoff distance of each other. Atoms withing the same molecules can are not connected
    :param atoms: List of atoms
    :type atoms: t.List[u.Atom]
    :param max_dis: maximal distance of two connected atoms in nm
    :type max_dis: float
    :return: A list of all connections, shorter than max_dis
    :rtype: t.List[u.connection]
    """

    max_dis_sq: float = max_dis * max_dis

    # Check for all pairs of atoms if they are within the specified cutoff of each other. Try to minimize necessary calculations using shortcut function of and and or
    # TODO: Could be quicker if we know ahead which are the prefferred directions of the molecules
    # TODO: Instead of first checking if distance in all three spacial directions is , cutoff it might be quicker to first check x distance, then x^2 + y^2 distance and then to add z^2
    # => Decrease cost in case of success a bit, increas cost in case of failure after first square. Decrease chance to onlzy fail after the last calculation
    connections: t.List[Types.Distance_Restraint] = []
    for i_a1 in range(0, len(atoms) - 1):
        a1 = atoms[
            i_a1]  # Quicker to loop over index of a1 and then acess it once for all inner loops, than findix index of a1 once for every onner loop
        for a2 in atoms[i_a1 + 1:]:
            if a1.chain != a2.chain and abs(a1.x - a2.x) <= max_dis and abs(a1.y - a2.y) <= max_dis and (
                    a1.z - a2.z) <= max_dis:
                dis_sq = (a1.x - a2.x)**2 + (a1.y - a2.y)**2 + (a1.z - a2.z)**2
                if dis_sq <= max_dis_sq:
                    connections.append(Types.Distance_Restraint(atomA=a1, atomB=a2))

    return connections


def get_all_short_pair_restraints(atoms: t.List[u.Atom], max_dis: float) -> t.List[
    Types.Distance_Restraint]:
    """
    get_all_short_pair_restraints: finds all pairs of atoms within a certain cutoff distance of each other. Atoms withing the same molecules can are not connected
    :param __atoms: List of atoms
    :type __atoms: t.List[u.Atom]
    :param max_dis: maximal distance of two connected atoms in nm
    :type max_dis: float
    :return: A list of all connections, shorter than max_dis
    :rtype: t.List[RestraintType.Pair_Restraint]:
    """
    # TODO: SPEED: first sort bymolecules and then only compare the molecule list, to avaiod all forbidden combinations of atoms within the samemolecules

    max_dis_sq: float = max_dis * max_dis

    # Check for all pairs of atoms if they are within the specified cutoff of each other. Try to minimize necessary calculations using shortcut function of and and or
    # TODO: Instead of first checking if distance in all three spacial directions is , cutoff it might be quicker to first check x distance, then x^2 + y^2 distance and then to add z^2
    # => Decrease cost in case of success a bit, increas cost in case of failure after first square. Decrease chance to onlzy fail after the last calculation
    __atom_pairs: t.List[Types.Distance_Restraint] = []
    for i_a1 in range(0, len(atoms) - 1):
        a1 = atoms[
            i_a1]  # Quicker to loop over index of a1 and then acess it once for all inner loops, than findix index of a1 once for every onner loop
        for a2 in atoms[(i_a1 + 1):]:
            if a1.resi != a2.resi and abs(a1.x - a2.x) <= max_dis and abs(a1.y - a2.y) <= max_dis and (
                    a1.z - a2.z) <= max_dis:
                dis_sq = (a1.x - a2.x)**2 + (a1.y - a2.y)**2 + (a1.z - a2.z)**2
                if dis_sq <= max_dis_sq:
                    __atom_pairs.append(Types.Distance_Restraint(a1, a2))

    # Convert the atom pairs to restraints
    return __atom_pairs


# TODO: Delete connect_two_molecules_kruskal, maximum_spanning_tree_kruskal
def _not_debugged_maximum_spanning_tree_kruskal(cons: t.List[
    Types.Distance_Restraint]) -> t.List[u.RestraintPair]:
    """
    Finds a maximal spanning tree in which the nodes represent distance restraints and edges the distance between their middles.
    The 'tree' will be provided in the form of a list. Connectivity info is not explicitly given. Due to the priority queue implementation of Kruskal's Algorithm the list is guaranteed to be sorted by descending distance.

    :WARNING: But it is not guaranteed that consecutive edges are connected!!!
    :WARNING: Involves calculating distances between all Restraint pairs. => O(n^2). Only use on preselected lists of connections
    :warning: old and not supported. Only kept here to keep connect two Molecules_kruskal_running
    :param cons: List of connections (Pairwise Distance Restraints)
    :type cons: t.List[u.connections]
    :return: The tree in the form of Pairs of connections and their distances
    :rtype: .list[u.RestraintPair]
    """

    # Step 1) Build List of all distances between the restraints: WARNING: Expensive (O(n!)). => Provide a preselected list of distance restraints

    # TODO: Could save time by implementing it as a heap, while building the Priority Quue, instead of sorting the list afterwards: Sorted probably O(n log n), heap O(n)

    def _check_circle(_edge_list: t.List[u.RestraintPair],
                      _new_edge: Types.Distance_Restraint) -> bool:
        """
        KEEP AS INNER FUNCTION OF maximum_spanning_tree_Kruskal: The function saves time, by NOT checking the whole Graph for circles, but just the part connected to the new edge.
        Technically it will return true if the new edge closes a circle or if the new edge is added to a subgraph that already contains a circle
        """
        # TODO: Clean this function. 1) Generalize: A node is just an index of anything, and an edge a tuple of two inidces.
        #                           2) Clean: Do not start at both edges of the new edge. Just start at one node of the new edge and see if you can reach the other.

        all_edges = [_new_edge] + _edge_list[
                                  :]  # CAREFULL: 1) Add _new_edge before, not after _edgeList. It is important, that the algorithm starts on the new edge, to avoid searching through a disconnected subgraph
        #          2) Do not work on _edge_list directly!
        # CARFULL: By starting on all_edges[0], ie new edge

        nodes: t.List = [all_edges[0].r1]

        ind = -1  # is the index of the first, node whose neighbours have not been added to the list yet
        old_len = 0  # length of nodes[] in last round
        while old_len < len(nodes):  # 'New nodes have been discovered in last iteration'
            ind += 1
            old_len = len(nodes)
            # STEP 1) Add all neighbours of the node to the list, and remove the edge that was passed to that neighbour from the list of edges
            for e in all_edges:
                if e.r1 == nodes[ind]:
                    nodes.append(e.r2)
                    all_edges.remove(e)
                elif e.r2 == nodes[ind]:
                    nodes.append(e.r1)
                    all_edges.remove(e)

            # STEP 2) Check if any of the points added are a point that was already visited,
            for old_node in nodes[:ind]:
                for new_node in nodes[ind + 1:]:
                    if old_node == new_node:
                        return True
                        # TODO:There is a lot of redundant testing here: When ind moves by one we will still check many of the combinations we already checked

        return False

    # ---------------------------------------------------------------------------------------------------------------------------------

    queue: t.List[u.RestraintPair] = []
    for i_r1 in range(len(cons) - 1):
        r1 = cons[i_r1]
        for r2 in cons[i_r1 + 1:]:
            queue.append(u.RestraintPair(r1, r2, np.sqrt(
                ((r1.atomA.x + r1.atomB.x - r2.atomA.x - r2.atomB.x) / 2)**2 + (
                    (r1.atomA.y + r1.atomA.y - r2.atomA.y - r2.atomB.y) / 2)**2 + (
                    (r1.atomA.z + r1.atomB.z - r2.atomA.z - r2.atomB.z) / 2)**2)))

    list.sort(queue, key=lambda x: x.distance, reverse=True)  # reversed = true: Descending order
    # Step 2) Build the list by moving through priority Queue, accepting objects that connect 2 nodes, at least one of which is not in the tree. Remove edges, that would connect two nodes already in the tree
    # TODO: Use a proper heap priority queue instead of a list +> Reduce time complexity
    # TODO: As we are working with a fully connected graph prims algorithms would  be better. (Always choose an edge connecting one node in the tree with one not in the tree)
    # That would require firstly: Implementationof queue as binary heap, update of edge weights after every addation of
    # => Problem: Would need a priority queue(implemented as heap) over vertices. Priority of each vertice is updated when a new edge is added to the tree. In our case that is relatively expensive. And as our n are not soooo high that might actually be a stronger effect than decreased complexity
    # The difference in complexity is relatively small: Kruskal e log e, e log v Prim, well implemented

    # Some Internet research: Prim & Kruskal are proven to give the optimal solution => We can be sure to get the same result both ways, EXCPET if there are different optimal solutions (i.e. ther are edges with same weight.)

    edges: t.list[u.RestraintPair] = []

    for new_edge in queue:
        # TODO: save some time by turning boolof check_circle,so I do not need not operator, and by seven len cons -1
        if not _check_circle(edges, new_edge):
            edges.append(new_edge)
            if (len(edges) == len(
                    cons) - 1):  # A spanning tree over a fully connected Graph with v vertices has v-1 edges. (Gartehaagproblem)
                break

    else:  # Moved through queue without finding enough edges. Impossible for a fully connected input graph
        raise u.BadArgumentException(
            "maximum_spanning_tree_kruskal(...) did not find enough edges for a spanning tree even though that should be mathematically impossible.  Input might be corrupted. Does it contain an Atom restrained to itself?")

    return edges


def maximal_spanning_tree_greedy(potential_restraints: t.List[Types.Distance_Restraint], n: int, priority_algo: str) -> \
        t.List[Types.Distance_Restraint]:
    """
    Finds a maximal spanning tree in which the nodes represent distance restraints and edges the distance between their middles.
    The 'tree' will be provided in the form of a list. Connectivity info is not explicitly given. Due to the working of Prims Algorithm it is guaranteed, that every vortex in the list is connected to one vortex before it in the list
    WARNING: But it is not guaranteed that the edges are ordered by size
    WARNING: Involves calculating distances between all Restraint pairs. => O(n!). Only use on preselected lists of connections
    :param cons: List of connections (Pairwise Distance Restraints)
    :type cons: t.List[u.connections]
    :return: The tree in the form of Pairs of connections and their distances
    :param n: how many restraints should be set?
    :type n: int
    :param priority_algo: Which algorithm should be used to update priorities of not yet chosen nodes?
    :type priority_algo: str
    :rtype: .list[u.RestraintPair]
    """

    # TODO: Implement different options for distance calculation within this function.
    # I think passing a function get_priority as argument would become really dirty really, quickly. But I can define some options as inner functions
    # BUT: Do I slow it down thisway. What is the cost of a function call?

    # Options i want: 1) Prim: Distance to farthest node in tree not yet in tree
    #                2) closest: Distance to closeest node in maximal_spanning_tree_prim)
    #                3) Average: Average distance to all nodes already in tree

    # Have to define priority_node as class for now, even though it is totally overkill. I need it to be mutable, so I can update priorities. But now I actually have to use opeator overloading to make it heapq compatible =>
    # TODO: Find a better solution

    class priority_node:
        '''..class priority_node: Used to assign priorities to nodes of a tree. A node can be any object. For speed reasons the node does not save the object itself, but only a index
            indicating the position of the 'node' in a list.
            The < (lt) operator is overloaded, to only compare by priority

            '''

        def __init__(self, priority: float, index: int):
            self.priority = priority
            self.index = index  # index in list of potential restraints
            # self.restraint = restraint

        # This is super dirty, because in python I cant even overload the , operator JUST for two p_nodes. This will now be used if i compare a priority_node to anything.
        # (But not if I compare anything to a priority_node )
        # => TODO: Kill it with fire
        def __lt__(self, other):
            if other.__class__ != priority_node:
                raise TypeError(
                    'The < (lt) operator is overloaded for priority_node. It should only be called tocompare it to another priority_node')

            return self.priority < other.priority

    # Give each node its initial priority: The dstance to the FARTHEST node
    # Read carefully: I sacrificed readability to mimize having to call certain attributes to often. Use _ as . operator

    # TODO: Move those inner functions out, so they can be reused in othre places
    # Define the different Algorithms

    # TODO: The current setup is a bit confusing:
    # _update_priority_(v) accesses the newly chosen node from within the function. Would be cleaner to give list of chosen restraints as arg as well

    def _update_priority_prim(v: priority_node):
        """The priority of the node is the MAXIMAL distance it has to any node in the tree """

        if m_sq[v.index][chosen_v.index] > v.priority:
            v.priority = m_sq[chosen_v.index][v.index]

    def _update_priority_shortest(v: priority_node):
        """The priority of the_node is the SHORTEST distance to any node in the tree"""

        # Shorter means -m > priority
        if m_sq[v.index][
            chosen_v.index] < v.priority or v.priority == 0:  # TODO Find a more general way to initialize, so I do not have to check for 0
            v.priority = m_sq[v.index][chosen_v.index]

    def _update_priority_cog(v: priority_node):
        """Choose the node, farthest from the average position of all chosen restraints """

        # TODO: Would be Quicker if we did not have to calculate center from scratch, but could keep it saved outside of fkt and update it
        # In return we make no use of the distance matrix => Maybe make a different function for this

        x_cog = 0
        y_cog = 0
        z_cog = 0

        for r in chosen_nodes:
            x_cog += pos_x[r.index]
            y_cog += pos_y[r.index]
            z_cog += pos_z[r.index]

        x_cog /= len(chosen_nodes)
        y_cog /= len(chosen_nodes)
        z_cog /= len(chosen_nodes)

        x_node = pos_x[v.index]
        y_node = pos_y[v.index]
        z_node = pos_z[v.index]

        # TODO: Speed up by doing intermediate checks (if x_dis>dis_sq do not calulate further)

        # v.prioirity changes everytime
        v.priority = np.sqrt((x_cog - x_node)**2 + (y_cog - y_node)**2 + (z_cog - z_node)**2)

    # TODO: Implement the avg method: Average of distance to all chosen restraints, not distance to cog (average_position)

    def _update_priority_biased_avg(v: priority_node):
        """Choose the node, farthest from the average position of all chosen restraints """

        # TODO: Ithink we could get a similar result quicker if we first simply summed the squares ofthe distances and then  took an cubic or bigger root?

        bias_exponent = 1 / 2  # Must be between 0 and 1
        corrected_exponent = bias_exponent * 1 / 2  # Corrected Exponent can be used on the matrix containing the squared distances directlz
        sum_of_biased_distances = 0

        for old_node in chosen_nodes:
            sum_of_biased_distances += (m_sq[old_node.index][v.index], corrected_exponent)**2

        v.priority = sum_of_biased_distances

    # TODO: Implement the avg method: Average of distance to all chosen restraints, not distance to cog (average_position)

    update_priority: t.Callable = None
    if priority_algo == 'prim':
        update_priority = _update_priority_prim
    elif priority_algo == 'shortest':
        update_priority = _update_priority_shortest
    elif priority_algo == 'cog':
        update_priority = _update_priority_cog
    elif priority_algo == 'biased_avg':
        update_priority = _update_priority_biased_avg
    else:
        raise u.BadArgumentException(
            priority_algo + ' is not an acceptable algorithm for maximal_spanning_tree_heuristic(...). Accepted values are: \'prim\' and \'shortest\'')

    # STEP 0) Check if there are enough restraints
    if len(potential_restraints) < n:
        raise NoOptimalSolutionException("There are not enough possible restraints to connect the two molecules")

    # STEP 1) Find positions of all restraints
    pos_x: t.List = [(r.atoms[0].x + r.atoms[1].x) / 2 for r in potential_restraints]
    pos_y: t.List = [(r.atoms[0].y + r.atoms[1].y) / 2 for r in potential_restraints]
    pos_z: t.List = [(r.atoms[0].z + r.atoms[1].z) / 2 for r in potential_restraints]

    # STEP 2) Find pairwise distances of all restraints

    m_sq: t.List[t.List[float]] = []  # Save squared distances. No need to waste time taking the sqrt
    # TODO SPEED: saving coords of r1 instead of calling them each time
    # TODO: Is there really no better 2d data structure than that? Maybe I implement it myself as array using multiplication to get 2d

    # Create the empty matrix
    empty_line = []
    for i in range(len(potential_restraints)):
        empty_line.append(None)
    for i in range(len(potential_restraints)):
        m_sq.append(
            empty_line[:])  # Slicing necessary: attach a copy of new_line, do NOT attach the samelist several times

    # Fill matrix with distance_sq values
    for i_r1 in range(len(potential_restraints) - 1):
        line = []
        for i_r2 in range(i_r1 + 1, len(potential_restraints)):
            dis_sq = ((pos_x[i_r1] - pos_x[i_r2])**2 + (pos_y[i_r1] - pos_y[i_r2])**2 + (
                pos_z[i_r1] - pos_z[i_r2])**2)
            m_sq[i_r1][i_r2] = dis_sq
            m_sq[i_r2][i_r1] = dis_sq
        m_sq.append(line)  # append, not extend.
    #
    # Print Matrix to check format
    # for i in range(len(m_sq)):
    #     print('')
    #     for j in range(len(m_sq[i])):
    #         print(str(m_sq[i][j]),end= ' ')

    # STEP 3) Set up priority queue and Choose the first node: One of the nodes, which is part of the longest edge

    first_node = priority_node(priority=0, index=-1)
    for i_line in range(len(potential_restraints) - 1):
        for i_col in range(i_line + 1, len(potential_restraints)):
            if m_sq[i_line][i_col] > first_node.priority:
                first_node = priority_node(priority=m_sq[i_line][i_col], index=i_line)

    priority_q = [priority_node(0, i_r) for i_r in range(len(potential_restraints))]
    chosen_nodes: t.list[priority_node] = [first_node]
    priority_q.remove(
        priority_q[first_node.index])  # Only works here because priority_q is still a copy of potential restraints
    # ANCHOR 20Mar M: Check if correct ind is removed

    # STEP 4) Update priority of all  nodes, according to distance to FARTHEST (for pure Prim algo) node that is already in tree

    while len(chosen_nodes) < n:
        # Update priorities
        chosen_v = chosen_nodes[-1]  # Last node that was selected
        for v in priority_q:
            update_priority(v)

        new_node = max(priority_q)
        priority_q.remove(new_node)
        chosen_nodes.append(new_node)

    return [potential_restraints[x.index] for x in chosen_nodes]


class NoOptimalSolutionException(Exception):
    '''
    A NoOptimalSolutionException is thrown, when an Optimizer reaches the end of its procedure without having found a solution fulfilling its criteria.
    '''

    def __init__(self, args):
        super().__init__(args)


priority_node = namedtuple('priority_node', 'priority node')


# --------------------------------------------------------------------Small algorithm functions----------------
# CALCULATE_VALUE: Calculate value of a set of points
def _calculate_value_convex_hull(restraint_pos):
    return ConvexHull(restraint_pos).volume


def _calculate_value_unscaled_pca_2d(restraint_pos):
    '''Quantity to optimize is the 2d-volume of the PCA of the restraint positions'''

    # TODO: Test!!!!!!
    # TODO: Think about other possibilities: 3d-volume / covariance / pca of all involved atoms instead of restraint
    eigen_vals = rdkit_functions._calc_pca_without_scaling(restraint_pos, dims=2, verbose=False)
    return eigen_vals[0] * eigen_vals[1]  # TODO GEneralize for other dimensions


def calculate_value_pca_relative_to_pca_of_both_molecules(restraint_pos):
    eigen_vals = rdkit_functions._calc_pca_without_scaling(restraint_pos, dims=2, verbose=False)


def _calculate_value_scaled_pca_2d(restraint_pos):
    # TODO: Scaling before pca is not really reasonable, because we do prefer a big relative variance in the biggest dimensions
    # TODO: Function returns a lot of stuff we do not need.
    raise NotImplementedError(
        'Scaled PCA is not implemented. WARNING: The tools_Rdkit function calc_pca only CENTERS the atoms, it does not scale!!!!')
    eigen_vals = rdkit_functions._calc_pca(restraint_pos, dims=2, verbose=False)[2]
    return eigen_vals[0] * eigen_vals[1]  # TODO GEneralize for other dimensions


# ------------utilites------------------------------------------
# TODO: make get postion a function of restraint
def find_middle(r: Types.Distance_Restraint):
    return (r.atoms[0].x + r.atoms[1].x) / 2, (r.atoms[0].y + r.atoms[1].y) / 2, (
            r.atoms[0].z + r.atoms[1].z) / 2


def maximal_weight_ring(edge_list: t.List[t.Tuple[int, int]], value_list: t.List[float], n_nodes: int):
    '''
    maximal_weight_ring chooses from a list of edges the ones that form the biggest ring including all nodes.
    Nodes are considered as ints (indices) and edges as tuples of nodes
    :warning: if the input is NOT a fully connected graph there might be no solution
    :param edge_list: list of all edges as tuples of ints
    :type edge_list: t.List[t.Tuple[int,int]]
    :param value_list: A list of all edge weights
    :type value_list: t.List[float]
    :param n_nodes: Number of nodes
    :type n_nodes: int
    :returns: Indices of the chosen edges
    :rtype t.List[ind]
    '''

    # 0)
    if len(edge_list) != len(value_list):
        raise u.BadArgumentException(
            'The lengths of the list inidcating edges weights and edges must be equal')
    if not len(edge_list) == n_nodes * (n_nodes - 1) / 2:
        print('WARNING: Method maximal_weight_ring has not been tested for not fully connected graphs', mv=4)

    def _edge_is_acceptable(__i_1, __i_2):
        '''warning: Call BEFORE  new edge is added'''
        # 1) Does not cause a fork (i.e create a node with 3 edges)
        if len(neighbours_of_node[__i_1]) == 2 or len(neighbours_of_node[__i_2]) == 2:
            # print(__i_1,__i_2,'causes fork',mv=1)
            return False
        # 2 Does not close a Ring. Much easier here than in general Kruskal: Because there are no forks: Just follow the chain from one node of the new edge and see if you can reach the other node

        # Start at node1 and see if i can reach node 2

        if len(neighbours_of_node[__i_1]) > 0:
            i_current_node = neighbours_of_node[__i_1][0]
            coming_from = 0 if neighbours_of_node[i_current_node][0] == __i_1 else 1
        else:
            i_current_node, coming_from = -1, -1

        # Follow the chain until you reach the other node of the new edge, or a loose end
        # print(' ',mv=1)
        while (i_current_node >= 0):
            if i_current_node == __i_2:  # Edge closes a ring
                # print(__i_1, __i_2, 'closes ring', mv=1)
                return False

            if (len(neighbours_of_node[i_current_node]) == 2):  # Chain goes on
                i_next_node = neighbours_of_node[i_current_node][0 if coming_from == 1 else 1]
                coming_from = 0 if neighbours_of_node[i_next_node][0] == i_current_node else 1
            else:
                i_next_node, coming_from = -1, -1

            # print(i_current_node, '-.>', i_next_node, mv=1)
            i_current_node = i_next_node

        return True

    # 1) List of indices used to acess edges sorted by values
    sorted_by_value = sorted(list(range(len(edge_list))), key=lambda i: value_list[i], reverse=True)

    # 2) Iterate over edges in order of sorted list
    neighbours_of_node = [[] for dummy_var in range(
        n_nodes)]  # Carefull: [[]]*n_nodes will fill the list with the SAME empty list n_nodes times
    i_chosen_edges = []
    # TODO: Instead of iterating of ind, loop while total connections not big enough
    for ind in sorted_by_value:  # CARFEULL: ind is reused in while loop to find last edge
        i_m1 = edge_list[ind][0]
        i_m2 = edge_list[ind][1]
        #        print(neighbours_of_node,mv=1)

        # Condition is stricter but easier to calculate than Kruskal
        # ) 1) There can be no  node with more than two edges (because in the end we want to have a ring
        # ) 2) We can not close the ring (because the last edge we add, will close it)
        #   => Because no node can have more than 2 edges that means there must always be at least one node with one edge. (Actually there must be two,but it is impossible that there is just one as long as no rings are allowed.

        if _edge_is_acceptable(i_m1, i_m2):
            neighbours_of_node[i_m1].append(i_m2)
            neighbours_of_node[i_m2].append(i_m1)
            i_chosen_edges.append(ind)
            print('ADDING', i_m1 + 1, i_m2 + 1, '(value =  ' + '{:05.2f}'.format(value_list[ind]) + ')', mv=1)
            if len(i_chosen_edges) == n_nodes - 1:
                last_position = sorted_by_value.index(ind)
                break
        else:
            print('DISCARDING', i_m1 + 1, i_m2 + 1, mv=1)

    else:  # For loop quit without break\
        raise ValueError(
            'Kruskal-like algorithm faild to connect all molecules. This can only happen if there is a serious bug in the Kruskal algorithm, or the input is not a fully connected graph.')

    # TODO: Clean: Start a new for loop, which starts where the one above ended. Do not check in loop if ind is bigger than possible, but just check with a for else statement that an edge was found
    # ADD the last edge: Does close ring, but is not allowed to fork
    found_last_edge = False
    print('----------', mv=1)
    for ind in sorted_by_value[last_position + 1:]:
        i_m1 = edge_list[ind][0]
        i_m2 = edge_list[ind][1]
        if not (len(neighbours_of_node[i_m1]) == 2 or len(neighbours_of_node[i_m2]) == 2):
            print('ADDING', i_m1 + 1, i_m2 + 1, '(value =  ' + '{:05.2f}'.format(value_list[ind]) + ')', mv=1)
            i_chosen_edges.append(ind)
            break
        else:
            print('DISCARDING', i_m1 + 1, i_m2 + 1, mv=1)
    else:  # No break => Did not find last edge
        raise ValueError(
            'Kruskal-like algorithm faild to connect all nodes. This can only happen if there is a serious bug in the Kruskal algorithm, or the input has a bad format')

    return i_chosen_edges


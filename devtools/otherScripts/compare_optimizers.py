
# ----------------------------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------Functions to test & compare Optimizers----------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------------



import typing as t
from scipy.spatial.qhull import QhullError

import os
import time
from restraintmaker.io import Exporter
from restraintmaker.utils import Utilities as u, Types
from restraintmaker.utils.Utilities import print
from restraintmaker.algorithm import Optimizer as o

def compare_pair_optimizers(criterion: t.Callable[[t.Tuple[float]], float], atoms: t.List[
    u.Atom],
                            opt_types: t.List[t.Type[o._MoleculeRingOptimizer]], opt_args: t.List[t.Tuple], out_path,
                            new_dir_name):
    '''
    compare_pair_optimizers runs different optimizers over the same set of atoms, saves the restraints for each PAIR of Molecules (not just the ones chosen for the ring), assigns them a value and saves all values acheived with different optimizers
    :param criterion: A function which assigns a value to a set of restraints, depending on their positions
    :type criterion: criterion: t.Callable[[t.Tuple[t.float]],float]
    :param atoms: List of atoms
    :type atoms: t.List[u.atom]
    :param opt_types: types of Optimzers which have to be compared. Must be Molecules Ring Optimziers
    :type opt_types: t.List[t.Type[_MoleculeRingOptimizer]]
    :param opt_args: args for the different Optimizers. Needs to have same length a opt, types
    :type opt_args: t.List[t.Tuple]
    :param out_path: path where a folder containing all the output should be saved.
    :type out_path: str or None
    :raises ValueError if arguments have a bad format
    :raises OSError if the directory can not be created
    :return:
    :rtype:
    '''

    if not len(opt_types) == len(opt_args):
        raise ValueError('There needs to be one set of arguments for each optimizer')
    out_dir = out_path + '/' + new_dir_name
    os.mkdir(out_dir)
    # Let the OSError pass if the dir can not be made

    # Arrange as dict of list:  outer dict: key = numbers of both molecules, value = list of values created by different optimzers
    values_of_mol_pairs_by_all_optimizers = {}
    # Create dict entries  and Folders for each Molecule pair
    n_mols = len(u.order_atoms_by_molecule(atoms))
    for m1 in range(1, n_mols):
        for m2 in range(m1 + 1, n_mols + 1):
            os.mkdir(out_dir + '/restraints_' + str(m1) + '-' + str(m2))
            values_of_mol_pairs_by_all_optimizers.update({str(m1) + '_' + str(m2): [None] * len(opt_types)})

    for ind, (cur_opt_type, cur_opt_args) in enumerate(zip(opt_types, opt_args)):

        my_Optimizer = cur_opt_type(atoms)
        my_Optimizer.get_args(lambda _: cur_opt_args)

        exporter_called = 0

        # Optimize all Molecule pairs
        for i_m1 in range(len(my_Optimizer.Molecules) - 1):
            for i_m2 in range(i_m1 + 1,
                              len(my_Optimizer.Molecules)):  # Cant use enumerate because i_m2 would start from 0
                print(
                    '_______________________________________________________________________________________________________________________________________',
                    mv=3)
                print('  testing ' + cur_opt_type.__name__ + '( ' + ', '.join(
                    str(a) for a in cur_opt_args) + ' ) on Molecules mol' + str(i_m1 + 1) + ' and mol' + str(i_m2 + 1),
                      mv=3)
                m1, m2 = my_Optimizer.Molecules[i_m1], my_Optimizer.Molecules[i_m2]
                filename = my_Optimizer.__class__.__name__ + '_' + '_'.join(
                    [str(a) for a in cur_opt_args]) + '.disres'
                try:
                    pairwise_restraints = u.execute_at_different_verbosity(new_verbosity_threshold=5,
                                                                           func=my_Optimizer.connect_two_molecules,
                                                                           m1=m1, m2=m2, n=4)

                    # Export Restraints
                    my_Exporter = Exporter.Gromos_Distance_Restraint_Exporter(pairwise_restraints)
                    my_Exporter.get_args(
                        lambda _: out_dir + '/restraints_' + str(i_m1 + 1) + '-' + str(i_m2 + 1) + '/' + filename)
                    my_Exporter.export_restraints(pairwise_restraints)
                    exporter_called += 1
                    print('Exporting', str(i_m1 + 1), str(i_m2 + 1), 'expoerter-called:', str(exporter_called), mv=2)

                    # Calculate Value

                    values_of_mol_pairs_by_all_optimizers[str(i_m1 + 1) + '_' + str(i_m2 + 1)][ind] = criterion(
                        [o.find_middle(r) for r in pairwise_restraints])
                except u.NoOptimalSolutionException:
                    my_file = open(out_dir + '/restraints_' + str(i_m1 + 1) + '-' + str(i_m2 + 1) + '/' + filename, 'w')
                    my_file.write('#NO Restraints could be generated. (Not enough atoms within cutoff distance')
                    my_file.close()
                    values_of_mol_pairs_by_all_optimizers[str(i_m1 + 1) + '_' + str(i_m2 + 1)][ind] = -1
                except QhullError:
                    # Restraints were set and saved. But the value could not be evaluated
                    values_of_mol_pairs_by_all_optimizers[str(i_m1 + 1) + '_' + str(i_m2 + 1)][ind] = -1
                    print("WARNING: Optimal solution could not be assigned a value (Precision Error)")

    # Write overview File containig all Values

    my_file = open(out_dir + '/values_overview.txt', 'w')

    # Some general infp
    my_file.write('#' + time.strftime('%d/%m/%Y %H:%M:%S') + '\n')
    my_file.write('#Criterion for comparison of Optimizers: ' + criterion.__name__ + '\n')
    # Write all Optimzers that are compared
    for ind, (o_type, o_args) in enumerate(zip(opt_types, opt_args)):
        my_file.write(
            '# o' + str(ind + 1) + ': ' + o_type.__name__ + ' ( ' + ', '.join(str(a) for a in o_args) + ' ) \n')

    my_file.write('\n')

    # Top line pf table: Names of the Optimizers
    my_file.write('mols \t' + '\t'.join('o' + str(ind + 1) for ind in range(len(opt_types))) + '\n')

    # Write table
    for mol_pair, values in zip(values_of_mol_pairs_by_all_optimizers.keys(),
                                values_of_mol_pairs_by_all_optimizers.values()):
        line = mol_pair.replace('_', ' & ') + '\t'
        line += '\t'.join('{:5.1f}'.format(v) for v in values)
        my_file.write(line + '\n')

    my_file.close()

import os
import pickle
import unittest

import restraintmaker.algorithm.Optimizer as Optimizer

# Input = 7_veryDifferentLigands. Unfiltered
test_file_dir = os.path.dirname(__file__) + "/test_files"
test_files_optimizers = test_file_dir + "/optimizers"

in_atom_list_path = test_files_optimizers + "/in_atom_list.p"
out_result_cog_MAD_optimizer_path = test_files_optimizers + "/out_cog_MAD_optimizer.p"
out_result_longestShortest_MAD_optimizer_path = test_files_optimizers + "/out_longestShortest_MAD_optimizer.p"
out_result_prim_optimizer_path = test_files_optimizers + "/out_prim_optimizer.p"
out_result_bruteForce_optimizer_path = test_files_optimizers + "/out_bruteForce_optimizer.p"

input = pickle.load(open(in_atom_list_path, "rb"), fix_imports=True)
file_handle_prim = open(out_result_prim_optimizer_path, "rb")
file_handle_cog_MAD = open(out_result_cog_MAD_optimizer_path, "rb")
file_handle_longestShortest_MAD = open(out_result_longestShortest_MAD_optimizer_path, "rb")
file_handle_bruteForce = open(out_result_bruteForce_optimizer_path, "rb")



class test_Optimizer(unittest.TestCase):

    def test_prim_optimizer(self):
        my_Optimizer = Optimizer.TreeHeuristicOptimizer(input)
        my_Optimizer.get_args(lambda x: (3, 1.2, 'prim', 'None'))  # nRes, distance, algo, ordering of chain
        found_restraints = my_Optimizer.make_restraints()

        #with open(out_result_prim_optimizer_path, "wb") as solution_file:
        #    pickle.dump(found_restraints, solution_file, fix_imports=True)

        expected_prim_soultion = pickle.load(file_handle_prim, fix_imports=True)
        self.check_restraint_results(found_restraints=found_restraints, expected_restraints=expected_prim_soultion)

    def test_shortest_optimizer(self):
        my_Optimizer = Optimizer.TreeHeuristicOptimizer(input)
        my_Optimizer.get_args(lambda x: (3, 1.2, 'shortest', 'None'))
        found_restraints = my_Optimizer.make_restraints()

        #with open(out_result_longestShortest_MAD_optimizer_path, "wb") as solution_file:
        #    pickle.dump(found_restraints, solution_file, fix_imports=True)

        expected_longestShortest_MAD_solution = pickle.load(file_handle_longestShortest_MAD, fix_imports=True)
        self.check_restraint_results(found_restraints=found_restraints,
                                     expected_restraints=expected_longestShortest_MAD_solution)

    def test_cog_optimizer(self):
        my_Optimizer = Optimizer.TreeHeuristicOptimizer(input)
        my_Optimizer.get_args(lambda x: (3, 1.2, 'cog', "None"))
        found_restraints = my_Optimizer.make_restraints()

        #with open(out_result_cog_MAD_optimizer_path, "wb") as solution_file:
        #   pickle.dump(found_restraints, solution_file, fix_imports=True)

        expected_cog_MAD_solution = pickle.load(file_handle_cog_MAD, fix_imports=True)
        self.check_restraint_results(found_restraints=found_restraints, expected_restraints=expected_cog_MAD_solution)

    def test_bruteForce_optimizer(self):
        my_Optimizer = Optimizer.BruteForceRingOptimzer(input)
        my_Optimizer.get_args(lambda x: (3, 1.2, 'pca', "None"))
        found_restraints = my_Optimizer.make_restraints()

        #with open(out_result_bruteForce_optimizer_path, "wb") as solution_file:
        #    pickle.dump(found_restraints, solution_file, fix_imports=True)

        expected_bruteForce_solution = pickle.load(file_handle_bruteForce, fix_imports=True)
        self.check_restraint_results(found_restraints=found_restraints,
                                     expected_restraints=expected_bruteForce_solution)


    def check_restraint_results(self, found_restraints, expected_restraints):
        # check if same ammount of resis was found
        if (len(found_restraints) != len(expected_restraints)):
            raise Exception(" The expected restraint ammount was " + str(
                len(expected_restraints)) + " but " + __name__ + " implementation found " + str(len(found_restraints)))

        # check if all atoms res ar the same:
        found = []
        not_found = []
        for found_res in found_restraints:
            setted_value = False

            found_atoms_tuple = found_res.atoms
            for expected in expected_restraints:
                expected_atom_tuple = expected.atoms
                if (found_atoms_tuple[0] == expected_atom_tuple[0]):
                    if (found_atoms_tuple[1] == expected_atom_tuple[1]):
                        found.append(True)
                        setted_value = True

                elif (found_atoms_tuple[0] == expected_atom_tuple[1]):
                    if (found_atoms_tuple[1] == expected_atom_tuple[0]):
                        found.append(True)
                        setted_value = True
            if (not setted_value):
                found.append(False)
                not_found.append(found_atoms_tuple)
        if (all(found)):
            pass
        else:
            raise Exception("Test failed, could not find all restraints! " + str(
                len(not_found)) + " were missing.\n Input: \n" + str(not_found) + "\n\nExpected:\n" + str(
                expected_restraints))

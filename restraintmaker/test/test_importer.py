"""
..module: Tests
    This file is testing, if the Exporter class.
     It checks if:
        * output is correct file format.
        * class is constructable
"""

import os
import pickle
import sys
import unittest

sys.path.append(os.path.dirname(__file__) + "../")

import restraintmaker.io.Importer as Importer

test_file_dir = os.path.dirname(__file__) + "/test_files"
test_files_exporter = test_file_dir + "/IO"

in_atom_list_path = test_files_exporter + "/in_atom_list.p"
in_gromos_disres = test_files_exporter + "/in_solution_disres.dat"
in_solution_restraints_path = test_files_exporter + "/in_restraints.p"

all_atoms = pickle.load(open(in_atom_list_path, "rb"), fix_imports=True)
solution_restraint = pickle.load(open(in_solution_restraints_path, "rb"), fix_imports=True)


class test_cnf(unittest.TestCase):
    """
    ..class: This class is testing the Importer Class
    """

    def test_Importer_Gromos_construct(self):
        gromos_importer = Importer.Gromos_Pair_Restraint_Importer(all_atoms)

    def test_Importer_Gromos_import_file_notFound(self):
        pass

    def test_Importer_Gromos_import_disresDat(self):
        # Atom Ifds of all Pair restraints
        gromos_importer = Importer.Gromos_Pair_Restraint_Importer(all_atoms)
        gromos_importer.get_args(lambda x: (in_gromos_disres))
        disres = gromos_importer.import_restraints()
        self.check_restraint_results(disres, solution_restraint)

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


if __name__ == '__main__':
    unittest.main()

"""
    Test Importer
    This file is testing, if the Exporter class.
     It checks if:
        * output is correct file format.
        * class is constructable
"""

import os
import pickle
import unittest
from pymol import cmd

import restraintmaker.io.Importer as Importer
from restraintmaker.interface_Pymol.pymol_utilities import pymol_utitlities


class test_Importer(unittest.TestCase):
    """
        This class is testing the Importer Class

    """

    test_file_dir = os.path.dirname(__file__) + "/test_files"

    test_files_io = test_file_dir + "/IO"
    in_disres1 = test_files_io+"/in_5ligs_disres.dat"
    in_disres2 = test_files_io + "/in_7ligs_disres.dat"

    test_files_structures = test_file_dir +"/test_systems"
    in_pdb1 = test_files_structures+"/systemA/CHK1_5Ligs.pdb"
    in_pdb2 = test_files_structures + "/systemB/BRD4_7Ligs.pdb"

    def setUp(self) -> None:
        cmd.load(self.in_pdb1)
        self.all_atoms1 = pymol_utitlities.pymol_selection_to_atom_list("all")
        cmd.reinitialize()

        cmd.load(self.in_pdb2)
        self.all_atoms2 = pymol_utitlities.pymol_selection_to_atom_list("all")
        cmd.reinitialize()


    def test_Importer_Gromos_construct(self):
        gromos_importer = Importer.Gromos_Distance_Restraint_Importer(self.all_atoms1)
        print(vars(gromos_importer))

    def test_Importer_Gromos_import_file_notFound(self):
        try:
            gromos_importer = Importer.Gromos_Distance_Restraint_Importer(self.all_atoms1)
            gromos_importer.get_args(lambda x: self.in_disres1+"df")
            print(vars(gromos_importer))
            raise Exception("There should be an Error!")
        except IOError as err:
            print("JUHU")


    def test_Importer_Gromos_import_disresDat(self):
        # Atom Ifds of all Pair restraints
        print("SubTest1")
        gromos_importer = Importer.Gromos_Distance_Restraint_Importer(self.all_atoms1)
        gromos_importer.get_args(lambda x: (self.in_disres1))
        disres1 = gromos_importer.import_restraints()
        print("\n".join(map(str, disres1)))
        print()

        print("SubTest2")
        gromos_importer = Importer.Gromos_Distance_Restraint_Importer(self.all_atoms2)
        gromos_importer.get_args(lambda x: (self.in_disres2))
        disres2 = gromos_importer.import_restraints()

        print("\n".join(map(str, disres2)))
        print()


    def _check_restraint_results(self, found_restraints, expected_restraints):
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

"""
    Tests
    This file is testing, if the Exporter class.
     It checks if:
        * output is correct file format.
        * class is constructable
"""
import os, pickle
import unittest

from restraintmaker.io import Exporter




class test_Exporter(unittest.TestCase):
    """
    ..autoclass:: This class is testing the Exporters
    """

    test_file_dir = os.path.dirname(__file__) + "/test_files"

    test_files_IO = test_file_dir + "/IO"
    in_disresObj_path1 = test_files_IO + "/in_5ligs_disres.p"
    in_disresObj_path2 = test_files_IO + "/in_7ligs_disres.p"
    out_gromos_disres_path1 = test_files_IO + "/out_5ligs_disres.dat"
    out_gromos_disres_path2 = test_files_IO + "/out_7ligs_disres.dat"
    in_solution_disres_path1 = test_files_IO + "/in_5ligs_disres.dat"
    in_solution_disres_path2 = test_files_IO + "/in_7ligs_disres.dat"

    def setUp(self) -> None:
        file_handle1 = open(self.in_disresObj_path1, "rb")
        file_handle2 = open(self.in_disresObj_path2, "rb")
        self.disres1 = pickle.load(file_handle1, fix_imports=True)
        self.disres2 = pickle.load(file_handle2, fix_imports=True)

        file_handle1.close()
        file_handle2.close()

    def test_Exporter_Gromos_construct(self):
        gromos_Exporter = Exporter.Gromos_Distance_Restraint_Exporter(self.disres1)

    def test_Exporter_Gromos_getargs(self):
        exporter = Exporter.Gromos_Distance_Restraint_Exporter(self.disres1)
        exporter.get_args(lambda x: self.out_gromos_disres_path1)

    def test_Exporter_Gromos_export_disresDat(self):
        exporter = Exporter.Gromos_Distance_Restraint_Exporter(self.disres1)
        exporter.get_args(lambda x: (self.out_gromos_disres_path1))
        out = exporter.export_restraints()

        with open(self.in_solution_disres_path1, "r") as solution:
            expected_filelines = solution.readlines()
            expected_str = "".join(expected_filelines)
        self.assertEqual(str(out), expected_str, msg="file text is not expected.")

        exporter = Exporter.Gromos_Distance_Restraint_Exporter(self.disres2)
        exporter.get_args(lambda x: (self.out_gromos_disres_path2))
        out = exporter.export_restraints()

        with open(self.in_solution_disres_path2, "r") as solution:
            expected_filelines = solution.readlines()
            expected_str = "".join(expected_filelines)
        self.assertEqual(str(out), expected_str, msg="file text is not expected.")


if __name__ == '__main__':
    unittest.main()

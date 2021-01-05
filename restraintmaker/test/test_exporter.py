"""
    Tests
    This file is testing, if the Exporter class.
     It checks if:
        * output is correct file format.
        * class is constructable
"""
import os, sys
import pickle, unittest

from restraintmaker.io import Exporter

test_file_dir = os.path.dirname(__file__) + "/test_files"
test_files_exporter = test_file_dir + "/IO"

in_data_disres_path = test_files_exporter + "/in_restraints.p"
in_solution_disres_path = test_files_exporter + "/in_solution_disres.dat"
out_gromos_disres_path = test_files_exporter + "/out_disres.dat"

distance_restraints = pickle.load(open(in_data_disres_path, "rb"), fix_imports=True)

class test_Exporter(unittest.TestCase):
    """
    ..autoclass:: This class is testing the Exporters
    """

    def test_Exporter_Gromos_construct(self):
        gromos_Exporter = Exporter.GromosPairRestraintExporter(distance_restraints)

    def test_Exporter_Gromos_getargs(self):
        exporter = Exporter.GromosPairRestraintExporter(distance_restraints)
        exporter.get_args(lambda x: (out_gromos_disres_path))

    def test_Exporter_Gromos_export_disresDat(self):
        exporter = Exporter.GromosPairRestraintExporter(distance_restraints)
        exporter.get_args(lambda x: (out_gromos_disres_path))
        out = exporter.export_restraints()

        expected_filelines = open(in_solution_disres_path, "r").readlines()
        expected_str = "".join(expected_filelines)

        self.assertEqual(str(out), expected_str, msg="file text is not expected.")


if __name__ == '__main__':
    unittest.main()

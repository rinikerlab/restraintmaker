"""
    Importer
    This module contains the functions, to import Files.

"""

import os, json
import typing as t
import warnings

from restraintmaker.io import Files
from restraintmaker.utils import Utilities as u, Restraints
from restraintmaker.utils.Utilities import print


class _Importer():


    def __init__(self, all_atoms: t.List[u.Atom]):
        """
            Importer Base class

        Parameters
        ----------
        all_atoms :  t.List[u.Atom]
            List of all Atoms: Needed because the importer will only read the ids of the atoms. It then/
            has to look up this id to find the rest of the informations
        """
        self.all_atoms = all_atoms

    def get_args(self, input_function: t.Callable):
        """
          get_args(...) should be overridden by every subclass of Exporter. It will assign all necessary varaibles using input_function

         For Gromos_Importer: in_path

        Parameters
        ----------
        input_function : t.Callable[[str],t.Any]
            a function that will provide the arguments for the selection in the necessary format.

        Returns
        -------
        NoReturn

        Raises
        ------
        BadArgument
            Excpetion if the input function does not deliver the arguments as needed

        NotImplementedError
            if _Importer.get_args is called directly
        """
        raise NotImplementedError('Direct call of method get_args of abstract parent class _Importer')

    def import_restraints(self, verbose: bool = False) -> t.List[dict]:
        """
            This function reads in a distance restraint file. Here in parent
            class, it serves as an interface for the child Classes.

        Parameters
        ----------
        in_path : str
            path to the disres file
        verbose : bool
             shall I be loud and noisy?

        Returns
        -------
        t.List[dict]
            returns a dict containing the atom ids for each disres. (has to be
            translated in interface_Pymol

        """
        raise NotImplementedError("This function is not implemented yet.")

class _Importer_Distance_Restraints(_Importer):
    file_ending:str
    def get_args(self, input_function: t.Callable):
        """
        should be overridden by every subclass of Exporter. It will assign all necessary varaibles using input_function

         For Gromos_Importer: in_path

        Parameters
        ----------
        input_function : t.Callable[[str], t.Any]
             a function that will provide the arguments for the selection in the necessary format.

        Returns
        -------
        NoReturn

        Raises
        ------
        BadArgument
            Excpetion if the input function does not deliver the arguments as needed
        NotImplementedError
            if _Importer.get_args is called directly
        """

        input = input_function('Which file should the restrants be read from? ')
        self.in_path = u.check_or_convert_argument(input, str)
        if(not self.in_path.endswith(self.file_ending)):
            warnings.warn("I could not find the correct file ending. However I will try importing! File ending should be: "+str(self.file_ending))
        if self.in_path == '' or self.in_path == 'None':
            raise u.BadArgumentException("Empty filename. (Unless you actually wanted to call your file 'None'. \n"
                                         "In which case you have to blame Python's promiscuous type conversion.  And yourself, for not using file extensions.)")

        if (not os.path.isfile(self.in_path)):
            raise IOError("Could not find File in import path: " + str(self.in_path))
    pass

class import_Gromos_Distance_Restraints(_Importer_Distance_Restraints):
    file_ending:str = "disres"
    def import_restraints(self, verbose: bool = False) -> t.List[dict]:
        """
        This function reads in a gromos distance restraint file. TODO: Read
            in additional Settings (r0, w0, etc...)

        TODO: WE could try to read thefile in get_args to get the Error there already.
         IDEA: READ THE WHOLE TEXT IN GET ARGS AND ONLY CONVERT IT HERE. OS ALLERRORS HAPPEN IN GET ARGS

        Parameters
        ----------
        verbose

        Returns
        -------
        t.List[dict]
            returns a dict containing the atom ids for each disres. (has to be
            translated in interface_Pymol


        """

        restraint_objects = []  # Define before try_block, so we can return the empty list in case of errors
        try:
            # readFiles
            disres_file = Files.Gromos_files.disres(self.in_path)
            if verbose: print(disres_file)
            if verbose: print("READ: " + "".join(disres_file["TITLE"]))

            # Analogous new version
            for restraint in disres_file.distance_res_spec_block.RESTRAINTS:
                if verbose: print(restraint)
                atom1 = u.find_atom_by_property(self.all_atoms, restraint.atom1i)
                atom2 = u.find_atom_by_property(self.all_atoms, restraint.atom2i)
                new_restraint = Restraints.DistanceRestraint(atomA=atom1, atomB=atom2)
                restraint_objects.append(new_restraint)

            # PRINT RESULTS
            for i_r, r in enumerate(restraint_objects):
                print('PairRestraint', str(i_r), ': Atom1_ID:', str(r.atoms[0].id), 'Atom2_ID:', str(r.atoms[1].id))

        except FileNotFoundError:
            print("Error: Could not find the file: \'" + self.in_path + "\'", mv=4)
        except KeyError:
            print("Error: Could not read the file: " + self.in_path,
                  mv=4)
        except ImportError:
            print("BAD ERROR: FAILED TO IMPORT THE MODULES NECESSARY FOR IMPORTING AND EXPORTING GROMOS FILES!", mv=4)

        if(len(restraint_objects) == 0):
            warnings.warn("Could not find any Restraint in file: "+str(self.in_path))
        return restraint_objects

class import_Gromacs_Distance_Restraints(_Importer_Distance_Restraints):
    file_ending:str = "itp"

    def import_restraints(self, verbose: bool = False) -> t.List[dict]:
        """
        This function reads in a gromos distance restraint file. TODO: Read
            in additional Settings (r0, w0, etc...)

        TODO: WE could try to read thefile in get_args to get the Error there already.
         IDEA: READ THE WHOLE TEXT IN GET ARGS AND ONLY CONVERT IT HERE. OS ALLERRORS HAPPEN IN GET ARGS

        Parameters
        ----------
        verbose

        Returns
        -------
        t.List[dict]
            returns a dict containing the atom ids for each disres. (has to be
            translated in interface_Pymol


        """

        restraint_objects = []  # Define before try_block, so we can return the empty list in case of errors
        try:
            # readFiles
            atoms_list = Files.Gromacs_files.Gromacs_distanceRestraint_file.load(in_path=self.in_path)

            # Analogous new version
            for restraint in atoms_list:
                if True: print(restraint, mv=5)
                atom1 = u.find_atom_by_property(self.all_atoms, restraint[0])
                atom2 = u.find_atom_by_property(self.all_atoms, restraint[1])
                new_restraint = Restraints.DistanceRestraint(atomA=atom1, atomB=atom2)
                restraint_objects.append(new_restraint)

            # PRINT RESULTS
            for i_r, r in enumerate(restraint_objects):
                print('PairRestraint', str(i_r), ': Atom1_ID:', str(r.atoms[0].id), 'Atom2_ID:', str(r.atoms[1].id))

        except FileNotFoundError:
            print("Error: Could not find the file: \'" + self.in_path + "\'", mv=4)
        except KeyError:
            print("Error: Could not read the file: " + self.in_path,
                  mv=4)
        except ImportError:
            print("BAD ERROR: FAILED TO IMPORT THE MODULES NECESSARY FOR IMPORTING AND EXPORTING GROMOS FILES!", mv=4)

        if(len(restraint_objects) == 0):
            warnings.warn("Could not find any Restraint in file: "+str(self.in_path))
        return restraint_objects


class import_JSON_Distance_Restraints(_Importer_Distance_Restraints):
    file_ending:str = "json"

    def import_restraints(self, verbose: bool = False) -> t.List[dict]:
        """
        This function reads in a gromos distance restraint file. TODO: Read
            in additional Settings (r0, w0, etc...)

        TODO: WE could try to read thefile in get_args to get the Error there already.
         IDEA: READ THE WHOLE TEXT IN GET ARGS AND ONLY CONVERT IT HERE. OS ALLERRORS HAPPEN IN GET ARGS

        Parameters
        ----------
        verbose

        Returns
        -------
        t.List[dict]
            returns a dict containing the atom ids for each disres. (has to be
            translated in interface_Pymol


        """

        restraint_objects = []  # Define before try_block, so we can return the empty list in case of errors
        try:
            tmp_file_content = json.load(open(self.in_path, "r"))

            # Analogous new version
            for i,restraint in tmp_file_content.items():
                if verbose: print(restraint)
                atom1 = u.find_atom_by_property(self.all_atoms, restraint['a1']['id'])
                atom2 = u.find_atom_by_property(self.all_atoms, restraint['a2']['id'])
                new_restraint = Restraints.DistanceRestraint(atomA=atom1, atomB=atom2)
                restraint_objects.append(new_restraint)

            # PRINT RESULTS
            for i_r, r in enumerate(restraint_objects):
                print('PairRestraint', str(i_r), ': Atom1_ID:', str(r.atoms[0].id), 'Atom2_ID:', str(r.atoms[1].id))

        except FileNotFoundError:
            print("Error: Could not find the file: \'" + self.in_path + "\'", mv=4)
        except KeyError:
            print("Error: Could not read the file: " + self.in_path,
                  mv=4)
        except ImportError:
            print("BAD ERROR: FAILED TO IMPORT THE MODULES NECESSARY FOR IMPORTING AND EXPORTING JSON FILES!", mv=4)

        if(len(restraint_objects) == 0):
            warnings.warn("Could not find any Restraint in file: "+str(self.in_path))
        return restraint_objects


class _Importer_Position_Restraints(_Importer):
    def get_args(self, input_function: t.Callable):
        """
        should be overridden by every subclass of Exporter. It will assign all necessary varaibles using input_function

         For Gromos_Importer: in_path

        Parameters
        ----------
        input_function : t.Callable[[str], t.Any]
             a function that will provide the arguments for the selection in the necessary format.

        Returns
        -------
        NoReturn

        Raises
        ------
        BadArgument
            Excpetion if the input function does not deliver the arguments as needed
        NotImplementedError
            if _Importer.get_args is called directly
        """

        input = input_function('Which file should the restrants be read from? ')
        self.in_path = u.check_or_convert_argument(input, str)

        if(self.in_path.endswith(self.file_ending)):
            warnings.warn("I could not find the correct file ending. However I will try importing! File ending should be: "+str(self.file_ending))

        if self.in_path == '' or self.in_path == 'None':
            raise u.BadArgumentException("Empty filename. (Unless you actually wanted to call your file 'None'. \n"
                                         "In which case you have to blame Python's promiscuous type conversion.  And yourself, for not using file extensions.)")

        if (not os.path.isfile(self.in_path)):
            raise IOError("Could not find File in import path: " + str(self.in_path))


class import_Gromos_Position_Restraints(_Importer_Position_Restraints):
    file_ending:str = "por"

    def import_restraints(self, verbose: bool = True) -> t.List[dict]:
        """
        This function reads in a gromos distance restraint file. TODO: Read
            in additional Settings (r0, w0, etc...)

        TODO: WE could try to read thefile in get_args to get the Error there already.
         IDEA: READ THE WHOLE TEXT IN GET ARGS AND ONLY CONVERT IT HERE. OS ALLERRORS HAPPEN IN GET ARGS

        Parameters
        ----------
        verbose

        Returns
        -------
        t.List[dict]
            returns a dict containing the atom ids for each disres. (has to be
            translated in interface_Pymol


        """

        restraint_objects = []  # Define before try_block, so we can return the empty list in case of errors
        try:
            # readFiles
            from pygromos.files.coord.posres import Position_Restraints

            posres_file = Position_Restraints(self.in_path)
            if verbose: print(posres_file,mv=5)

            # Analogous new version
            for restraint in posres_file.POSRESSPEC:
                if verbose: print(restraint)
                atom1 = u.find_atom_by_property(self.all_atoms, restraint.atomID)
                new_restraint = Restraints.PositionRestraint(atomA=atom1)
                restraint_objects.append(new_restraint)

            # PRINT RESULTS
            for i_r, r in enumerate(restraint_objects):
                print('PairRestraint', str(i_r), ': Atom1_ID:', str(r.atomA.id))

        except FileNotFoundError:
            print("Error: Could not find the file: \'" + self.in_path + "\'", mv=4)
        except KeyError:
            print("Error: Could not read the file: " + self.in_path,
                  mv=4)
        except ImportError:
            print("BAD ERROR: FAILED TO IMPORT THE MODULES NECESSARY FOR IMPORTING AND EXPORTING GROMOS FILES!", mv=4)

        if (len(restraint_objects) == 0):
            warnings.warn("Could not find any Restraint in file: " + str(self.in_path))

        return restraint_objects


class import_Gromacs_Position_Restraints(_Importer_Position_Restraints):
    file_ending:str = "itp"

    def import_restraints(self, verbose: bool = False) -> t.List[dict]:
        """
        This function reads in a gromos distance restraint file. TODO: Read
            in additional Settings (r0, w0, etc...)

        TODO: WE could try to read thefile in get_args to get the Error there already.
         IDEA: READ THE WHOLE TEXT IN GET ARGS AND ONLY CONVERT IT HERE. OS ALLERRORS HAPPEN IN GET ARGS

        Parameters
        ----------
        verbose

        Returns
        -------
        t.List[dict]
            returns a dict containing the atom ids for each disres. (has to be
            translated in interface_Pymol


        """

        restraint_objects = []  # Define before try_block, so we can return the empty list in case of errors
        try:
            # readFiles
            atoms_list = Files.Gromacs_files.Gromacs_positionRestraint_file.load(in_path=self.in_path)

            # Analogous new version
            for restraint in atoms_list:
                if verbose: print(restraint, mv=5)
                atom1 = u.find_atom_by_property(self.all_atoms, restraint)
                new_restraint = Restraints.PositionRestraint(atomA=atom1)
                restraint_objects.append(new_restraint)

            # PRINT RESULTS
            for i_r, r in enumerate(restraint_objects):
                print('PairRestraint', str(i_r), ': Atom1_ID:', str(r.atomA.id))

        except FileNotFoundError:
            print("Error: Could not find the file: \'" + self.in_path + "\'", mv=4)
        except KeyError:
            print("Error: Could not read the file: " + self.in_path,
                  mv=4)
        except ImportError:
            print("BAD ERROR: FAILED TO IMPORT THE MODULES NECESSARY FOR IMPORTING AND EXPORTING GROMACS FILES!", mv=4)

        if(len(restraint_objects) == 0):
            warnings.warn("Could not find any Restraint in file: "+str(self.in_path))

        return restraint_objects


class import_JSON_Position_Restraints(_Importer_Position_Restraints):
    file_ending:str = "json"

    def import_restraints(self, verbose: bool = False) -> t.List[dict]:
        """
        This function reads in a gromos distance restraint file. TODO: Read
            in additional Settings (r0, w0, etc...)

        TODO: WE could try to read thefile in get_args to get the Error there already.
         IDEA: READ THE WHOLE TEXT IN GET ARGS AND ONLY CONVERT IT HERE. OS ALLERRORS HAPPEN IN GET ARGS

        Parameters
        ----------
        verbose

        Returns
        -------
        t.List[dict]
            returns a dict containing the atom ids for each disres. (has to be
            translated in interface_Pymol


        """

        restraint_objects = []  # Define before try_block, so we can return the empty list in case of errors
        try:
            tmp_file_content = json.load(open(self.in_path, "r"))

            # Analogous new version
            for i,restraint in tmp_file_content.items():
                if verbose: print(restraint)
                atom1 = u.find_atom_by_property(self.all_atoms, restraint['a1']['id'])
                atomRef = u.find_atom_by_property(self.all_atoms, restraint['aR']['id'])
                new_restraint = Restraints.PositionRestraint(atomA=atom1, reference_atom=atomRef)
                restraint_objects.append(new_restraint)

            # PRINT RESULTS
            for i_r, r in enumerate(restraint_objects):
                print('PositionRestraint', str(i_r), ': Atom1_ID:', str(r.atoms[0].id), 'AtomR_ID:', str(r.reference_atom.id))

        except FileNotFoundError:
            print("Error: Could not find the file: \'" + self.in_path + "\'", mv=4)
        except KeyError:
            print("Error: Could not read the file: " + self.in_path,
                  mv=4)
        except ImportError:
            print("BAD ERROR: FAILED TO IMPORT THE MODULES NECESSARY FOR IMPORTING AND EXPORTING JSON FILES!", mv=4)

        if(len(restraint_objects) == 0):
            warnings.warn("Could not find any Restraint in file: "+str(self.in_path))
        return restraint_objects


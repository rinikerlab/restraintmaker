"""
.. automodule:: Importer
    :members:
    This module contains the functions, to import Files.

"""

import os
import sys
import typing as t

import restraintmaker.utils.Types
import restraintmaker.utils.Utilities
import restraintmaker.utils.program_states

sys.path.append(os.path.dirname(__file__))

from restraintmaker.utils import Utilities as u, Types
from restraintmaker.utils.Utilities import print


class _Importer():
    """
    ..class: _Importer
        private parent class for all Importers.
    """

    def __init__(self, all_atoms: t.List[restraintmaker.utils.Utilities.Atom]):
        '''

        :param all_atoms: List of all Atoms: Needed because the importer will only read the ids of the atoms. It then/
            has to look up this id to find the rest of the information
        :type all_atoms:  t.List[u.Atom]
        '''
        self.all_atoms = all_atoms

    def get_args(self, input_function: t.Callable):
        '''
         get_args(...) should be overridden by every subclass of Exporter. It will assign all necessary varaibles using input_function

         For Gromos_Importer: in_path

         :param input_function: a function that will provide the arguments for the selection in the necessary format.
         :type input_function: t.Callable[[str],t.Any]
         :return: -
         :rtype: None
         :raises BadArgument Excpetion if the input function does not deliver the arguments as needed
         :raises NotImplementedError if _Importer.get_args is called directly
        '''

        raise NotImplementedError('Direct call of method get_args of abstract parent class _Importer')

    def import_restraints(in_path: str, verbose: bool = False) -> t.List[dict]:
        """
        ..function: import_restraints
            This function reads in a distance restraint file. Here in parent
            class, it serves as an interface for the child Classes.

        Args:
            in_path (str): path to the disres file
            verbose (bool): shall I be loud and noisy?

        Returns:
            returns a dict containing the atom ids for each disres. (has to be
            translated in interface_Pymol

            t.List[dict]:
        """

        raise NotImplementedError("This function is not implemented yet.")


# TODO: Need to adjust this if we rename Pair Restraint to Distance Restraint
# TODO: Do we really need a separte importert & exporter for each restraint type/ Couldnt we write different restraint types into the same file? => GROMOS_IMPORTER would be enough
class Gromos_Pair_Restraint_Importer(_Importer):
    """
    ..class: Gromos_Importer
        This class contains functions, for reading in Gromos Files.
    """

    def __init__(self, all_atoms: t.List[restraintmaker.utils.Utilities.Atom]):
        super().__init__(all_atoms)

        # attributes to be set in get_args
        self.in_path = ''

    def get_args(self, input_function: t.Callable):
        '''
         get_args(...) should be overridden by every subclass of Exporter. It will assign all necessary varaibles using input_function

         For Gromos_Importer: in_path

         :param input_function: a function that will provide the arguments for the selection in the necessary format.
         :type input_function: t.Callable[[str],t.Any]
         :return: -
         :rtype: None
         :raises BadArgument Excpetion if the input function does not deliver the arguments as needed
        :raises NotImplementedError if _Importer.get_args is called directly

        '''

        input = input_function('Which file should the restrants be read from? ')
        self.in_path = u.check_or_convert_argument(input, str)
        if self.in_path == '' or self.in_path == 'None':
            raise restraintmaker.utils.Utilities.BadArgumentException("Empty filename. (Unless you actually wanted to call your file 'None'. \n"
                                         "In which case you have to blame Python's promiscuous type conversion.  And yourself, for not using file extensions.)")
        if (not os.path.isfile(self.in_path)):
            raise IOError("Could not find File in import path: " + str(self.in_path))

    def import_restraints(self, verbose: bool = False) -> t.List[dict]:
        """
        ..function: import_restraints
            This function reads in a gromos distance restraint file. TODO: Read
            in additional Settings (r0, w0, etc...)

        Args:
            in_path (str): path to the disres file
            verbose (bool): shall I be loud and noisy?

        Returns:
            returns a dict containing the atom ids for each disres. (has to be
            translated in interface_Pymol

            t.List[dict]:

        :raises: TODO
        """

        restraint_objects = []  # Define before try_block, so we can return the empty list in case of errors
        try:
            # TODO CLEAN IMPORT: WE should do this import in the beginning. It is extremley annoying, if you  go through the whole program an the get an error because this import fails. Sam e with exporter
            import Files
            # readFiles

            disres_file = Files.Gromos_files.disres(self.in_path)
            if verbose: print(disres_file)
            if verbose: print("READ: " + "".join(disres_file["TITLE"]))

            # Analogous new version
            for restraint in disres_file.distance_res_spec_block.RESTRAINTS:
                if verbose: print(restraint)
                atom1 = self.find_atom_by_id(restraint.atom1i)
                atom2 = self.find_atom_by_id(restraint.atom2i)
                new_restraint = Types.Distance_Restraint(atoms=[atom1, atom2])
                restraint_objects.append(new_restraint)

            # PRINT RESULTS
            for i_r, r in enumerate(restraint_objects):
                print('PairRestraint', str(i_r), ': Atom1_ID:', str(r.atoms[0].id), 'Atom2_ID:', str(r.atoms[1].id))


        except FileNotFoundError:
            print("Error: Could not find the file: \'" + self.in_path + "\'", mv=4)
        except KeyError:
            print("Error: Could not read the file: " + self.in_path,
                  mv=4)  # TODO: WE could try to read thefile in get_args to get the Error there already. IDEA: READ THE WHOLE TEXT IN GET ARGS AND ONLY CONVERT IT HERE. OS ALLERRORS HAPPEN IN GET ARGS
        except ImportError:
            print("BAD ERROR: FAILED TO IMPORT THE MODULES NECESSARY FOR IMPORTING AND EXPORTING GROMOS FILES!", mv=4)

        return restraint_objects

    # TODO:Move to utilites, generalize for all properties
    def find_atom_by_id(self, id: int) -> restraintmaker.utils.Utilities.Atom:
        '''
        find_atom will look for the atom with the specified id

        :param id: id of the atom
        :type id: int
        :raises: ValueError if there is no or more than one atom with that id
        '''

        all_hits = list(filter(lambda a: a.id == id, self.all_atoms))
        if len(all_hits) == 0:
            raise ValueError('There is no atom with id: ' + str(id))
        elif len(all_hits) > 1:
            raise ValueError('There is more than one atom with id: ' + str(id))
        else:
            return all_hits[0]

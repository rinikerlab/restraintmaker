"""
    The module Exporter is providing the required code for exporting files from PyMol
"""
import typing as t
from math import sqrt

from restraintmaker.io import Files
from restraintmaker.types import Restraint
from restraintmaker.utils import Utilities as u
from restraintmaker.utils.Utilities import print


class _Exporter():
    """
  ..class: _exporter This class is the private parent
    class, that is giving the interface for submethods.
    """

    def __init__(self, restraints: t.List[Restraint._Restraint]):
        '''
        :param restraints: Restraints to be saved
        :type restraints: Restraint
        '''
        self.restraints = restraints

    def get_args(self, input_function: t.Callable[[str], t.Any]):
        '''
        get_args(...) should be overridden by every subclass of Exporter. It will assign all necessary attributes using input_function

        :param input_function: a function that will provide the arguments for the selection in the necessary format.
        :type input_function: t.Callable[[str],t.Any]
        :return: -
        :rtype: None
        :raises: u.BadArgumentException

        '''
        raise NotImplementedError("Direct call of method get_args of abstract parent class _Exporter.")

    def export_restraints(settings: t.Dict[str, t.Any], verbose: bool = False):
        '''
        export_restraints must be overridden by every subclass of Exporter. Writes the restraints into a file.

        :param verbose: print progress (True) or not (False)
        :type verbose: bool
        :return: -
        :rtype: None
        :raises FileNotFoundError
        '''
        raise NotImplementedError("Direct call of method export_disres of abstract parent class _Exporter.")


# TODO: If I rename Pair restraint to distrance restraint this must be adjusted here
class GromosPairRestraintExporter(_Exporter):

    # TODO: Think about: Do I really need a seperate Gromos_exporter for every kind of DisRes, or can I use one exporter that
    # accepts all gromos compatible restraints
    # TODO:COMMENT - the restraint file definitions are sadly very different!
    def __init__(self, restraints: t.List[Restraint.Pair_Restraint]):
        '''
        :param restraints: Restraints to be saved
        :type restraints: Restraint
        '''
        for r in restraints:
            if not isinstance(r, Restraint.Pair_Restraint):
                raise TypeError('Gromos_Pair_Restriant_Exporter only accepts Pair restraints as input')
        super().__init__(restraints)

        # attributes to be specified in get_args:
        self.out_path = None

    def get_args(self, input_function: t.Callable[[str], t.Any]):
        '''
        get_args(...) should be overridden by every subclass of Exporter. It will assign all necessary varaibles using input_function

        For Gromos_Exporter: out_path

        :param input_function: a function that will provide the arguments for the selection in the necessary format.
        :type input_function: t.Callable[[str],t.Any]
        :return: -
        :rtype: None
        :raises: u.BadArgumentException
        '''

        # Error checking can only be done when we try to acess the file
        self.out_path = u.check_or_convert_argument(input_function('Name of the output File:'), str)
        if self.out_path == '' or self.out_path == 'None':
            raise u.BadArgumentException("Empty filename. (Unless you actually wanted to call your file None. \n"
                                         "In which case you have to blame Python's promiscuous type conversion.  And yourself, for not using file extensions.)")

    # TODO CLEAN: At the moment there is no way to change the default argument, because it is essential for the program logic that.

    def export_restraints(self, verbose: bool = False) -> str:
        # TODO: Go through all the functions it calls to see where we would get a file error
        '''
        export_restraints must be overridden by every subclass of Exporter. Writes the restraints into a file.

        For Gromos Exporter it will be in a gromos compatible format.

        :param verbose: print progress (True) or not (False)
        :type verbose: bool
        :return: -
        :rtype: None
        :raises FileNotFoundError
        '''

        ##disres_file settings  #todo: realise these parameters as settings dict, which the user can provide as input
        KDISH = 0.1
        KDISC = 0.153
        fullharm = 1
        deviationdist = None
        general_dist = None

        # build_up clean disres file;
        # clean_restrains = data

        def build_pair_distance_restraints(self, restrain_atoms: t.List[dict], fullharm: int = 1, deviationdist=None,
                                           general_dist=None) -> list:

            #   Can do it in two steps: 1)type = fullharm, deviationdist, generaldist, 2) value
            # check input
            if (fullharm == 0):
                if (deviationdist != None or general_dist != None):
                    raise IOError("Please use one option, fullharm, ddist or generaldist.")

            elif (deviationdist != None):

                if (fullharm == 0 or general_dist != None):
                    raise IOError("Please use one option, fullharm, ddist or generaldist.")

            elif (general_dist != None):
                if (deviationdist != None or fullharm == 0):
                    raise IOError("Please use one option, fullharm, ddist or generaldist.")

            if (deviationdist == None):
                deviationdist = 0.0

            restraint_dict = []
            # TODO use zipto iterate over a1 and a2 directly

            for r in self.restraints:
                a1 = r.atoms[0]
                a2 = r.atoms[1]

                # TODO:Make distance an attribute of RestraintPair and set it once at creation. BUT THEN WE HAVE TO MAKE SUERE ATOMS ARE NOT CHANGE DAFTER THAT
                distance_A = sqrt((float(a1.x) - float(a2.x)) ** 2 + (float(a1.y) - float(a2.y)) ** 2 + (
                            float(a1.z) - float(a2.z)) ** 2)
                distance_nm = round(distance_A, 2) / 10
                comment = "##\t" + a1.resn + "/" + a1.name + " " + str(
                    a1.id) + " - " + a2.resn + "/" + a2.name + " " + str(a2.id) + "\n"
                distance = r.atoms
                new_entry = Files.Gromos_blocks.atom_pair_distanceRes(i1=a1.id, j1=0, k1=0, l1=0, type1=0, i2=a2.id,
                                                                      j2=0, k2=0, l2=0, type2=0, r0=distance_nm, w0=1.0,
                                                                      rah=fullharm,
                                                                      comment=comment)
                # print(new_entry)
                restraint_dict.append(new_entry)

            return restraint_dict

        disres_out = build_pair_distance_restraints(self, restrain_atoms=self.restraints, fullharm=fullharm,
                                                    deviationdist=deviationdist)

        ##WRITE out disres.dat
        if verbose: print("generate out_dict")
        disres_out_dict = {"KDISH": KDISH, "KDISC": KDISC,
                           "RESTRAINTHEADER": "i  j  k  l  type    i  j  k  l  type    r0    w0    rah".split(),
                           # header
                           "RESTRAINTS": disres_out}
        if verbose: print("generate top_disres obj")
        disres_file = Files.Gromos_files.disres()
        if verbose: print("top_disres obj add:")
        disres_file.add_block(blocktitle="TITLE", content="generated disres file for BRD4 with PYMOL wizard\n",
                              verbose=True)
        # disres.write('/home/rhinerc/Desktop/testingRestraints_afterTitelBlock.disres')
        disres_file.add_block(blocktitle="DISTANCERESSPEC", content=disres_out_dict, verbose=verbose)
        # disres_file.write('/home/rhinerc/Desktop/testingRestraints_afterDISRESBLOCK.disres')
        # print(disres_out_dict,mv=1)
        disres_file.write(self.out_path)
        if verbose: print("wrote to: " + self.out_path)

        return str(disres_file)

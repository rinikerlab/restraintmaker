"""
    The module Exporter is providing the required code for exporting files from PyMol
    If you want to add a new exporter, simply just derive it from _Exporter and the code will automatically add it to the GUI.
"""
import typing as t
import json
from math import sqrt
import numpy as np

from restraintmaker.io import Files

from restraintmaker.utils import Utilities as u, Restraints
from restraintmaker.utils.Utilities import print


class _Exporter():
    """
  ..class: _exporter
    """

    def __init__(self, restraints: t.List[Restraints._Restraint]):
        """
        This class is the private parent
        class, that is giving the interface for submethods.

        Parameters
        ----------
        restraints : Restraints
             Restraints to be saved
        """
        self.restraints = restraints

    def get_args(self, input_function: t.Callable[[str], t.Any]):
        """
             should be overridden by every subclass of Exporter. It will assign all necessary attributes using input_function

        Parameters
        ----------
        input_function:  t.Callable[[str],t.Any]
            a function that will provide the arguments for the selection in the necessary format.

        Returns
        -------
        NoReturn

        Raises
        ------
        u.BadArgumentException
        """
        raise NotImplementedError("Direct call of method get_args of abstract parent class _Exporter.")

    def export_restraints(settings: t.Dict[str, t.Any], verbose: bool = False):
        """
                export_restraints must be overridden by every subclass of Exporter. Writes the restraints into a file.

        Parameters
        ----------
        verbose : bool
            print progress (True) or not (False)

        Returns
        -------
        NoReturn

        """

        raise NotImplementedError("Direct call of method export_disres of abstract parent class _Exporter.")


class _Export_Distance_Restraints(_Exporter):
    def __init__(self, restraints: t.List[Restraints.DistanceRestraint]):
        """
            This is a exporting class for Gromos Distance Restraints

        Parameters
        ----------
        restraints : Restraints
            Restraints to be saved
        """
        for r in restraints:
            if not isinstance(r, Restraints.DistanceRestraint):
                raise TypeError('Distance_Pair_Restraint_Exporter only accepts Pair restraints as input')
        super().__init__(restraints)

        # attributes to be specified in get_args:
        self.out_path = None

    def get_args(self, input_function: t.Callable[[str], t.Any]):
        """should be overridden by every subclass of Exporter. It will assign all necessary varaibles using input_function

        For Gromos_Exporter: out_path

        Parameters
        ----------
        input_function:  t.Callable[[str],t.Any]
            a function that will provide the arguments for the selection in the necessary format.

        Returns
        -------
        NoReturn

        Raises
        ------
         u.BadArgumentException
        """

        # Error checking can only be done when we try to acess the file
        self.out_path = u.check_or_convert_argument(input_function('Name of the output File:'), str)
        if self.out_path == '' or self.out_path == 'None':
            raise u.BadArgumentException(
                "Empty filename. (Unless you actually wanted to call your file None. \n"
                "In which case you have to blame Python's promiscuous type conversion.  And yourself, for not using file extensions.)")
        self.out_path = ".".join(self.out_path.split(".")[:-1]) if(len(".".join(self.out_path.split(".")[:-1]))>0) else self.out_path

    pass


class export_Gromos_Distance_Restraints(_Export_Distance_Restraints):
    def export_restraints(self, verbose: bool = True) -> str:
        """
            export_restraints must be overridden by every subclass of Exporter. Writes the restraints into a file.
            For Gromos Exporter it will be in a gromos compatible format.

            todo: realise these parameters as settings dict, which the user can provide as input

        Parameters
        ----------
        verbose : bool
            print progress (True) or not (False)

        Returns
        -------

        """

        ##disres_file settings
        KDISH = 0.1
        KDISC = 0.153
        fullharm = 1
        deviationdist = None
        general_dist = None

        # build_up clean disres file;
        def build_pair_distance_restraints(self, restrain_atoms: t.List, fullharm: int = 1, deviationdist=None,
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
                                                                      j2=0, k2=0, l2=0, type2=0, r0=np.floor(distance_nm*100)/100, w0=1.0,
                                                                      rah=fullharm,
                                                                      comment=comment)
                # print(new_entry)
                restraint_dict.append(new_entry)

            return restraint_dict

        disres_out = build_pair_distance_restraints(self, restrain_atoms=self.restraints, fullharm=fullharm,
                                                    deviationdist=deviationdist)

        ##WRITE out disres.dat
        print("generate out_dict", mv=0)
        disres_out_dict = {"KDISH": KDISH, "KDISC": KDISC,
                           "RESTRAINTHEADER": "i  j  k  l  type    i  j  k  l  type    r0    w0    rah".split(),
                           # header
                           "RESTRAINTS": disres_out}
        print("generate top_disres obj", mv=0)
        disres_file = Files.Gromos_files.disres()
        print("top_disres obj add:", mv=0)
        disres_file.add_block(blocktitle="TITLE", content="generated disres file with restraintmaker\n",
                              verbose=True)
        disres_file.add_block(blocktitle="DISTANCERESSPEC", content=disres_out_dict, verbose=verbose)
        disres_file.write(self.out_path+".disres")
        print("wrote to: " + self.out_path+".disres", mv=4)

        return str(disres_file)


class export_Gromacs_Distance_Restraints(_Export_Distance_Restraints):
    def export_restraints(self, verbose: bool = True) -> str:
        """
            export_restraints must be overridden by every subclass of Exporter. Writes the restraints into a file.
            For Gromacs Exporter it will be in a gromacs compatible format.


        Parameters
        ----------
        verbose : bool
            print progress (True) or not (False)

        Returns
        -------

        """
        gmx_disres_file = Files.Gromacs_files.Gromacs_distanceRestraint_file()
        gmx_disres_file.build_disres(self.restraints)
        gmx_disres_file.write(out_path=self.out_path+".itp")
        return str(gmx_disres_file)


class export_JSON_Distance_Restraints(_Export_Distance_Restraints):
    def export_restraints(self, verbose: bool = True) -> str:
        """
            export_restraints must be overridden by every subclass of Exporter. Writes the restraints into a file.
            For Gromacs Exporter it will be in a gromacs compatible format.


        Parameters
        ----------
        verbose : bool
            print progress (True) or not (False)

        Returns
        -------

        """
        restrains_d = {}
        for i, res in enumerate(self.restraints):
            a1 = res.atoms[0]
            a2 = res.atoms[1]
            d = res._distance/10

            restrains_d.update({i: {"a1":a1._asdict(),
                                    "a2":a2._asdict(),
                                    "d":d}})
        json.dump(restrains_d, open(self.out_path+".json", "w"),  indent=4)
        return str(restrains_d)


class _Export_Position_Restraints(_Exporter):
    def __init__(self, restraints: t.List[Restraints.PositionRestraint]):
        """
            This is a exporting class for Gromos Distance Restraints

        Parameters
        ----------
        restraints : Restraints
            Restraints to be saved
        """
        for r in restraints:
            if not isinstance(r, Restraints.PositionRestraint):
                raise TypeError('Position_Restriant_Exporter only accepts Position restraints as input')
        super().__init__(restraints)

        # attributes to be specified in get_args:
        self.out_path = None

    def get_args(self, input_function: t.Callable[[str], t.Any]):
        """should be overridden by every subclass of Exporter. It will assign all necessary varaibles using input_function

        For Gromos_Exporter: out_path

        Parameters
        ----------
        input_function:  t.Callable[[str],t.Any]
            a function that will provide the arguments for the selection in the necessary format.

        Returns
        -------
        NoReturn

        Raises
        ------
         u.BadArgumentException
        """

        # Error checking can only be done when we try to acess the file
        self.out_path = u.check_or_convert_argument(input_function('Name of the output File:'), str)
        if self.out_path == '' or self.out_path == 'None':
            raise u.BadArgumentException(
                "Empty filename. (Unless you actually wanted to call your file None. \n"
                "In which case you have to blame Python's promiscuous type conversion.  And yourself, for not using file extensions.)")
        self.out_path = ".".join(self.out_path.split(".")[:-1]) if(len(".".join(self.out_path.split(".")[:-1]))>0) else self.out_path



    pass


class export_Gromos_Position_RESTRAINTS(_Export_Position_Restraints):

    def export_restraints(self, verbose: bool = True) -> str:
        """
            export_restraints must be overridden by every subclass of Exporter. Writes the restraints into a file.
            For Gromos Exporter it will be in a gromos compatible format.

            todo: realise these parameters as settings dict, which the user can provide as input

        Parameters
        ----------
        verbose : bool
            print progress (True) or not (False)

        Returns
        -------

        """

        from pygromos.files.coord.posres import Position_Restraints
        from pygromos.files.coord.refpos import Reference_Position
        from pygromos.files.blocks.coord_blocks import TITLE, POSRESSPEC, REFPOSITION, atomP

        title_block = TITLE(content="Position Restraints Generated with RestraintMaker.")
        restraints = [atomP(r.atomA.resi, r.atomA.resn, r.atomA.name,atomID=r.atomA.id,
                            xp=r.atomA.x, yp=r.atomA.y, zp=r.atomA.z) for r in self.restraints]
        posres_block = POSRESSPEC(None)
        posres_block.content=restraints

        path_prefix = ".".join(self.out_path.split(".")[:-1]) if(len(".".join(self.out_path.split(".")[:-1]))>0) else self.out_path

        prF = Position_Restraints(None)
        prF.add_block(block=title_block)
        prF.add_block(block=posres_block)
        prF.write(path_prefix+".por")


        title_block = TITLE(content="Position Restraints Generated with RestraintMaker.")
        restraints = [atomP(r.reference_atom.resi, r.reference_atom.resn, r.reference_atom.name,atomID=r.reference_atom.id,
                            xp=r.reference_atom.x, yp=r.reference_atom.y, zp=r.reference_atom.z) for r in self.restraints]
        refpos_block = REFPOSITION(None)
        refpos_block.content=restraints
        rpF = Reference_Position(None)
        rpF.add_block(block=title_block)
        rpF.add_block(block=refpos_block)
        rpF.write(path_prefix+".rpr")


        return str(prF)


class export_Gromacs_Position_RESTRAINTS(_Export_Position_Restraints):
    def export_restraints(self, verbose: bool = True) -> str:
        """
            export_restraints must be overridden by every subclass of Exporter. Writes the restraints into a file.
            For Gromacs Exporter it will be in a gromacs compatible format.


        Parameters
        ----------
        verbose : bool
            print progress (True) or not (False)

        Returns
        -------

        """
        gmx_disres_file = Files.Gromacs_files.Gromacs_positionRestraint_file()
        gmx_disres_file.build_posres(self.restraints)
        gmx_disres_file.write(out_path=self.out_path+".itp")
        return str(gmx_disres_file)


class export_JSON_Position_RESTRAINTS(_Export_Position_Restraints):
    def export_restraints(self, verbose: bool = True) -> str:
        """
            export_restraints must be overridden by every subclass of Exporter. Writes the restraints into a file.
            For Gromacs Exporter it will be in a gromacs compatible format.


        Parameters
        ----------
        verbose : bool
            print progress (True) or not (False)

        Returns
        -------

        """
        restrains_d = {}
        for i, res in enumerate(self.restraints):
            a1 = res.atomA
            aR = res.reference_atom
            d = res.distance_to_reference_position / 10

            restrains_d.update({i: {"a1": a1._asdict(),
                                    "aR": aR._asdict(),
                                    "dAR": d}})
        json.dump(restrains_d, open(self.out_path+".json", "w"), indent=4)
        return str(restrains_d)
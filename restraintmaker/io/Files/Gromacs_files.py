"""
    Gromacs_files
    TODO:write DOCU! and code :)
"""
import numpy as np
from typing import List, Tuple
from collections import namedtuple


gmx_disres = namedtuple("distanceRestraint", ["i", "j","type", "index", "typeb", "low", "up1", "up2", "fac", "comment"])
gmx_posres = namedtuple("positionRestraint", ["i", "r", "funct", "fc", "comment"])

class Gromacs_distanceRestraint_file():
    header:List[str] = ["i", "j", "type", "index", "type'", "low", "up1", "up2", "fac"]
    block_name = "distance_restraints"
    disres = []
    std_low:float = 0.0
    std_up1:float = 0.0
    std_up2:float = 0.15
    std_index:int = 1


    def __init__(self, disres_funct:int=1, disres_factor:int=1, standard_up1:bool=False, standard_up2:bool=False):
        self.disres_funct = disres_funct
        self.disres_functb = disres_funct

        self.disres_factor = disres_factor

        self.standard_up1=standard_up1
        self.standard_up2=standard_up2
        pass

    def __str__(self):
        msg = ["[ "+self.block_name+" ]"]
        msg.append(";"+" ".join(self.header))
        format_str = "{:>5} {:>5}      {:>1}      {:>1}      {:>1}      {:>2.1f}     {:>2.1f}     {:>2.1f}      {:>2.1f}; {}"
        for res in self.disres:
            msg.append(format_str.format(res.i, res.j, res.type, res.index, res.typeb, res.low, res.up1, res.up2, res.fac, res.comment))

        return "\n".join(msg)

    def write(self, out_path:str)->str:
        out_file = open(out_path, "w")
        out_file.write(str(self))
        out_file.close()
        return out_path

    @staticmethod
    def load(in_path:str)->List[Tuple[int, int]]:
        """
        This function is a very simple loading function, that gets the atom IDs of selected restratins in the in_file.
        """
        in_file = open(in_path, "r")
        lines = in_file.readlines()

        block_start = False
        atoms_list = []
        for line in lines:
            if(line.strip().startswith("[") and Gromacs_distanceRestraint_file.block_name in line and not block_start):
                block_start = True
            elif(line.strip().startswith("[") and block_start):
                block_start = False
                break
            elif(block_start):
                if(line.strip().startswith(";")):
                    continue
                else:
                    res_list = line.split()
                    atoms = (int(res_list[0]), int(res_list[1]))
                    atoms_list.append(atoms)
        return atoms_list

    def build_disres(self, selected_restraints: List):

        restraints = []
        for i, r in enumerate(selected_restraints):
            a1 = r.atoms[0]
            a2 = r.atoms[1]

            distance_A = np.sqrt((float(a1.x) - float(a2.x)) ** 2 + (float(a1.y) - float(a2.y)) ** 2 + (
                    float(a1.z) - float(a2.z)) ** 2)
            distance_nm = round(distance_A, 2) / 10

            r0 = self.std_low

            if(self.standard_up1):
                r1 = self.std_up1
            else:
                r1 = distance_nm

            if(self.std_up2):
                r2 = self.std_up2
            else:
                r2 = distance_nm+self.std_up2



            res = gmx_disres(i=a1.id, j=a2.id, type=self.disres_funct, index=self.std_index, typeb=self.disres_functb,
                             low=self.std_low, up1=self.std_up1, up2=self.std_up2,
                             fac=self.disres_factor, comment=a1.resn+"/"+a1.name+" - "+a2.resn+"/"+a2.name)
            restraints.append(res)

        self.disres = restraints


class Gromacs_positionRestraint_file():
    header:List[str] = ["ai", "funct", "fc"]
    block_name = "position_restraints"
    posres = []
    std_fc:float = tuple([1000,1000,1000])
    std_funct:float = 0.0

    def __init__(self,  std_funct:int=1, std_fc:float=tuple([1000,1000,1000])):
        self.std_funct = std_funct
        self.std_fc=std_fc
        pass

    def __str__(self):
        msg = ["[ "+self.block_name+" ]"]
        msg.append(";"+" ".join(self.header))

        format_str = "{:>4} {:>5} {:>8} {:>8} {:>8} ;{}"
        for res in self.posres:
            msg.append(format_str.format(res.i, res.funct, res.fc[0], res.fc[1], res.fc[2], res.comment))

        return "\n".join(msg)

    def write(self, out_path:str)->str:
        out_file = open(out_path, "w")
        out_file.write(str(self))
        out_file.close()
        return out_path

    @staticmethod
    def load(in_path:str)->List[Tuple[int, int]]:
        """
        This function is a very simple loading function, that gets the atom IDs of selected restratins in the in_file.
        """
        in_file = open(in_path, "r")
        lines = in_file.readlines()

        block_start = False
        atoms_list = []
        for line in lines:
            if(line.strip().startswith("[") and Gromacs_positionRestraint_file.block_name in line and not block_start):
                block_start = True
            elif(line.strip().startswith("[") and block_start):
                block_start = False
                break
            elif(block_start):
                if(line.strip().startswith(";")):
                    continue
                else:
                    res_list = line.split()
                    atom = int(res_list[0])
                    atoms_list.append(atom)
        return atoms_list



    def build_posres(self, selected_restraints: List):

        restraints = []
        for i, r in enumerate(selected_restraints):
            a1 = r.atomA
            aR = r.reference_atom

            distance_A = np.sqrt((float(a1.x) - float(aR.x)) ** 2 + (float(a1.y) - float(aR.y)) ** 2 + (
                    float(a1.z) - float(aR.z)) ** 2)
            distance_nm = round(distance_A, 2) / 10


            res = gmx_posres(i=a1.id, r=aR.id, funct=self.std_funct, fc=self.std_fc,
                             comment=a1.resn+"/"+a1.name+" to "+aR.resn+"/"+aR.name)
            restraints.append(res)

        self.posres = restraints

"""
.. automodule:: Gromos_blocks
    :members:
    TODO:write DOCU!
"""

import typing as t
from enum import Enum


# enum
##TOP
class distant_Restraint_Type(Enum):
    half_harmonic_repulsive = -1
    full_harmonic_distance_restraint = 0
    half_harmonic_attractive = 1


class geometric_code(Enum):
    real_atom = 0
    virtual_H_atom_aliphaticC = 1
    virtual_H_atom_aromaticC = 2
    virtual_H_atoms_goc_aliph = 3
    virtual_H_atom_aliphaticC_strange = 4
    pseudo_H_atom_goc_H_atoms_CH3 = 5
    pseudo_H_atoms_goc_of_two_CH3 = 6
    pseudo_H_atoms_goc_of_three_CH3 = 7
    virtual_atoms_cog = -1
    virtual_atoms_com = -2


# FIELDS
class _generic_field():
    def __str__(self):
        return self.to_string()

    def to_string(self):
        raise NotImplementedError("to string method needs to be implemented!")


##COORD
class atomP(_generic_field):
    def __init__(self, resID, resName, atomType, atomID, xp, yp, zp):
        """
        Args:
            resID:
            resName:
            atomType:
            atomID:
            xp:
            yp:
            zp:
        """
        self.resID = resID
        self.resName = resName
        self.atomType = atomType
        self.atomID = atomID
        self.xp = xp
        self.yp = yp
        self.zp = zp

    def to_string(self) -> str:
        return "{: >5} {: >3}  {: <5}{: >7}{: 15.9f}{:15.9f}{:15.9f}\n".format(self.resID, self.resName, self.atomType,
                                                                               self.atomID, self.xp,
                                                                               self.yp, self.zp)


class atomV(_generic_field):
    def __init__(self, resID, resName, atomType, atomID, xv, yv, zv):
        """
        Args:
            resID:
            resName:
            atomType:
            atomID:
            xv:
            yv:
            zv:
        """
        self.resID = resID
        self.resName = resName
        self.atomType = atomType
        self.atomID = atomID
        self.xv = xv
        self.yv = yv
        self.zv = zv

    def to_string(self) -> str:
        return "{: >5} {: >3}  {: <5}{: >7}{: 15.9f}{:15.9f}{:15.9f}\n".format(self.resID, self.resName, self.atomType,
                                                                               self.atomID, self.xv,
                                                                               self.yv, self.zv)


class lattice_shift(_generic_field):
    def __init__(self, atomID, x, y, z):
        """
        Args:
            atomID:
            x:
            y:
            z:
        """
        self.atomID = atomID
        self.x = x
        self.y = y
        self.z = z

    def to_string(self) -> str:
        return "{:>10}{:>10}{:>10}\n".format(self.x, self.y, self.z)


##TOP
class atom_pair_distanceRes(_generic_field):
    def __init__(self, i1: int, j1: int, k1: int, l1: int, type1: (geometric_code or int),
                 i2: int, j2: int, k2: int, l2: int, type2: (geometric_code or int),
                 r0: float, w0: float, rah: (distant_Restraint_Type or int), comment: str = ""):
        """
            ..function: Constructor of Gromos atom_pair_distanceRes

        Args:
            i1 (int): id of atom i of first molecule
            j1 (int): id of atom j of first molecule
            k1 (int): id of atom k of first molecule
            l1 (int): id of atom l of first molecule
            type1: geometric restraintype of first molecule
            i2 (int): id of atom i of second molecule
            j2 (int): id of atom j of second molecule
            k2 (int): id of atom k of second molecule
            l2 (int): id of atom l of second molecule
            type2:
            r0 (float): radius_0 of restraint
            w0 (float): weighting of restraint
            rah: restraint_type
            comment (str):
        """
        try:
            self.atom1i = int(i1)
            self.atom1j = int(j1)
            self.atom1k = int(k1)
            self.atom1l = int(l1)

            if (type(type1) is geometric_code):
                self.atom1ic = type1
            elif (type(type1) is int or (type(type1) is str and str(type1).isdigit())):
                self.atom1ic = geometric_code(int(type1))
            else:
                raise ValueError("geometric index atom1ic in atom_pair_distanceRes unknown\n" + str(type1))

            self.atom2i = int(i2)
            self.atom2j = int(j2)
            self.atom2k = int(k2)
            self.atom2l = int(l2)

            if (type(type2) is geometric_code):
                self.atom2ic = type2
            elif (type(type2) is int or (type(type2) is str and str(type2).isdigit())):
                self.atom2ic = geometric_code(int(type2))
            else:
                raise ValueError("geometric index atom2ic in atom_pair_distanceRes unknown\n" + str(type2))

            self.radius_0 = float(r0)
            self.weight = float(w0)
            if (type(rah) is int or (type(rah) is str and str(rah).isdigit())):
                self.disResType = distant_Restraint_Type(int(rah))
            elif (type(rah) is distant_Restraint_Type):
                self.disResType = rah
            else:
                raise ValueError("DisresType in atom_pair_distanceRes unknown\n" + str(rah))
            self.comment = comment
        except:
            raise IOError("COULD NOT convert a parameter for distancerestraint field into correct form!")

    def to_string(self):
        if (len(self.comment) > 0 and not self.comment.endswith("\n")):
            self.comment += "\n"
        return self.comment + "{:>5} {:>5} {:>5} {:>5} {:>5}    {:>5} {:>5} {:>5} {:>5} {:>5}  {:10.5f} {:10.5f} {:>3}\n".format(
            self.atom1i, self.atom1j, self.atom1k, self.atom1l, self.atom1ic.value, self.atom2i, self.atom2j,
            self.atom2k,
            self.atom2l, self.atom2ic.value, self.radius_0, self.weight, self.disResType.value)


# BLOCKS
##genericblock:
class _generic_gromos_block:
    def __init__(self, name: str, used: bool):  # content:str,
        """
        Args:
            name (str):
            used (bool):
        """
        self.used = used
        self.name = name
        self.line_seperator = " \n"
        self.field_seperator = " \t "
        # self.content = content

    def block_to_string(self) -> str:
        result = self.name + self.line_seperator
        result += self.field_seperator.join(self.content)
        result += "END" + self.line_seperator
        return result

    def get_name(self):
        return self.name

    def __str__(self):
        return self.block_to_string()


##MISC
class timestep_block(_generic_gromos_block):
    def __init__(self, t: float, step: int, subcontent=False):
        """
        Args:
            t (float):
            step (int):
            subcontent:
        """
        super().__init__(used=True, name="TIMESTEP")
        self.t = t
        self.step = step
        self.subcontent = subcontent

    def block_to_string(self) -> str:
        result = self.name + "\n"
        result += " \t{} \t{} \n".format(self.step, self.t)
        result += "END\n"
        return result


class title_block(_generic_gromos_block):
    order = [[["content"]]]

    def __init__(self, content):
        """
        Args:
            content:
        """
        super().__init__(used=True, name="TITLE")
        self.content = content

    def block_to_string(self) -> str:
        result = ""
        result += str(self.name) + "\n"
        result += "".join(self.content)
        if ("\t >>> Generated with python lib function_libs utilities. (riniker group)\n" not in self.content):
            result += "\t >>> Generated with python lib function_libs utilities. (riniker group)\n"
            result += "\t >>> line_seperator: " + repr(self.line_seperator) + "\t field_seperator: " + repr(
                self.field_seperator) + "\n"
        result += "END\n"
        return result


##COORD Blocks
class atom_pos_block(_generic_gromos_block):
    def __init__(self, content: t.List[atomP]):
        """
        Args:
            content:
        """
        if (content):
            super().__init__(used=True, name="POSITION")
            self.content = content
        else:
            raise ValueError("atom_pos block Constructor is missing correct argument combination!")

    def block_to_string(self) -> str:
        result = self.name + "\n"
        for x in self.content:
            result += x.to_string()
        result += "END\n"

        return result


class atom_vel_block(_generic_gromos_block):
    def __init__(self, content: t.List[atomV]):
        """
        Args:
            content:
        """
        super().__init__(used=True, name="VELOCITY")
        self.content = content

    def block_to_string(self) -> str:
        result = self.name + "\n"
        for x in self.content:
            result += x.to_string()
        result += "END\n"

        return result


class atom_ref_pos_block(_generic_gromos_block):

    def __init__(self, content: list):
        """
        Args:
            content (list):
        """
        if (content):
            super().__init__(used=True, name="REFPOSITION")
            self.content = content
        else:
            raise ValueError("atom_pos block Constructor is missing correct argument combination!")

    def block_to_string(self) -> str:
        result = self.name + "\n"
        for x in self.content:
            result += x.to_string()
        result += "END\n"

        return result


class lattice_shifts_block(_generic_gromos_block):
    def __init__(self, content: t.List[lattice_shift]):
        """
        Args:
            content:
        """
        super().__init__(used=True, name="LATTICESHIFTS")
        self.content = content

    def block_to_string(self) -> str:
        result = self.name + "\n"
        for x in self.content:
            result += x.to_string()

        result += "END\n"

        return result


class genbox_block(_generic_gromos_block):
    def __init__(self, pbc: Enum, length: list, angles: list, euler: list, origin: list):
        """
        Args:
            pbc (Enum):
            length (list):
            angles (list):
            euler (list):
            origin (list):
        """
        super().__init__(used=True, name="GENBOX")
        self.pbc = pbc
        self.length = length
        self.angles = angles
        self.euler = euler
        self.origin = origin

    def block_to_string(self) -> str:
        result = self.name + "\n"
        result += "{:>5}\n".format(str(self.pbc.value))
        result += "{:>15.9f}{:>15.9f}{:>15.9f}\n".format(self.length[0], self.length[1], self.length[2])
        result += "{:>15.9f}{:>15.9f}{:>15.9f}\n".format(self.angles[0], self.angles[1], self.angles[2])
        result += "{:>15.9f}{:>15.9f}{:>15.9f}\n".format(self.euler[0], self.euler[1], self.euler[2])
        result += "{:>15.9f}{:>15.9f}{:>15.9f}\n".format(self.origin[0], self.origin[1], self.origin[2])

        result += "END\n"

        return result


##TOP
class distance_res_spec_block(_generic_gromos_block):
    def __init__(self, KDISH: int, KDISC: int, RESTRAINTHEADER: list, RESTRAINTS: list):
        """
        Args:
            KDISH (int):
            KDISC (int):
            RESTRAINTHEADER (list):
            RESTRAINTS (list):
        """
        super().__init__(used=True, name="DISTANCERESSPEC")
        self.KDISH = KDISH
        self.KDISC = KDISC
        self.RESTRAINTHEADER = RESTRAINTHEADER
        self.RESTRAINTS = RESTRAINTS

    def block_to_string(self) -> str:
        result = self.name + self.line_seperator
        result += "# KDISH" + self.field_seperator + "KDISC" + self.line_seperator
        result += self.field_seperator + str(self.KDISH) + self.field_seperator + str(self.KDISC) + self.line_seperator
        result += "#{:>4} {:>5} {:>5} {:>5} {:>5}    {:>5} {:>5} {:>5} {:>5} {:>5}  {:>10} {:>10} {:>3}\n".format(
            *self.RESTRAINTHEADER)
        for x in self.RESTRAINTS:
            result += x.to_string()
        result += "END\n"
        return result

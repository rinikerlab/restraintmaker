"""
.. automodule:: Gromos Files
    :members:
    Description:
    in this lib, gromos topo file mainpulating functions are gathered
    Author: Benjamin Schroeder
    TODO:write DOCU!
"""

# imports
from typing import *

from restraintmaker.io.Files import Gromos_blocks as blocks
from restraintmaker.utils.Utilities import print


# file Classes
class general_gromos_file():
    path: str
    content: Dict[str, Any]
    required_blocks: List[str] = ["TITLE"]
    blocksset: List[dict] = []

    def __init__(self, in_path: (str or dict) = None):

        """
        Args:
            in_path:
        """
        self.blocksset = []
        self.block_names = {"TITLE": "title_block"}
        if type(in_path) is str:
            self.path = in_path
        elif (type(in_path) == dict):
            self.content = in_path
            self.path = None
        elif (in_path == None):
            print("WARNING!: generated empty gromosFileobj!")
        else:
            raise IOError("Gromos File Input path is not str, dict or None!")

    def __str__(self):
        text = ""
        for block in self.blocksset:
            text += self.__getattribute__(self.block_names.get(block)).block_to_string()
        return text

    def add_block(self, blocktitle: str, content: dict, verbose: bool = False):
        """
        Args:
            blocktitle (str):
            content (dict):
            verbose (bool):
        """
        if blocktitle in self.block_names:
            if blocktitle == "TITLE":
                self.title_block = blocks.title_block(content)
                self.blocksset.append(blocktitle)
            else:
                try:
                    content = {k.split("(")[0]: v for k, v in content.items()}
                    content = {k.split(":")[0]: v for k, v in content.items()}
                    content = {k.split(" ")[0]: v for k, v in content.items()}

                    self.__setattr__(self.block_names.get(blocktitle),
                                     blocks.__getattribute__(self.block_names.get(blocktitle))(**content))
                    if verbose:
                        print("++++++++++++++++++++++++++++++")
                        print("New Block: Adding " + blocktitle + " block")
                        print(content)
                    self.blocksset.append(blocktitle)
                except:
                    print("Error while adding new value - can not resolve value names in \'" + blocktitle + "\' block!")
                    print("Content is " + str(content))
                    print("Block knows " + str(
                        (blocks.__getattribute__(self.block_names.get(blocktitle)).__init__.__code__.co_varnames)[1:]))
                    exit(1)
            if verbose:
                print("Block " + blocktitle + " added to gromos File object.")

        else:

            print(" Warning while adding new block \'" + blocktitle + "\' - can not resolve \'" + blocktitle + "\'\n" +
                  "  -> Block will be skipped!")

        return 0

    def _read_gromos_block(self, lines: list):
        """
        Args:
            lines (list):
        """
        first_key = True
        blocks = {}
        key = ""
        block_lines = []
        for line in lines:
            if (first_key and not line.startswith("#")):
                key = line.strip()
                first_key = False

            elif (line.strip().startswith("END")):
                blocks.update({key: block_lines})
                first_key = True
                block_lines = []
            else:
                block_lines.append(line)
        return blocks

    def write(self, path):
        """
        Args:
            path:
        """
        file = open(path, "w")
        file.write(self.__str__())
        file.flush()
        file.close()
        print("New file generated: " + path)
        return path


class disres(general_gromos_file):
    required_blocks: List[str] = ["TITLE", "DISTANCERESPEC"]

    def __init__(self, in_path: (str or dict) = None):
        """
        Args:
            in_path:
        """
        super().__init__(in_path=in_path)
        self.block_names.update({"DISTANCERESSPEC": "distance_res_spec_block"})

        if ("path" in dir(self) and self.path != None):
            self.read_disres_file(in_path)

    def parse_disres_file(self, path: str) -> dict:
        # specific parser for subblocks of disres Files
        """
        Args:
            path (str):
        """

        def _read_disres_subblock(blocks):
            result_data = {}
            for block in blocks:
                if (block == "TITLE"):
                    result_data.update({block: blocks[block]})

                elif (block == "DISTANCERESSPEC"):
                    subcontent = {}
                    restrains = []
                    current_block = blocks[block]

                    # readout KDISH or KDISC
                    keys = current_block[0].replace("#", "").strip().split()
                    values = current_block[1].split()
                    subcontent.update({key: values[keys.index(key)] for key in keys})

                    # read list header:
                    line_header = current_block[2].replace("#", "").split()
                    ##unify keys:
                    key_dict = {"i": 1, "j": 1, "k": 1, "l": 1, "type": 1}
                    renamed_header = []
                    for x in line_header:
                        if (x in key_dict):
                            renamed_header.append(x + str(key_dict[x]))
                            key_dict[x] += 1
                        else:
                            renamed_header.append(x)
                    line_header = renamed_header

                    # read restraints
                    for line in current_block[3:]:
                        if (not line.startswith("#") and len(line.split()) == len(line_header)):
                            values = line.split()
                            restrains.append({key: values[line_header.index(key)] for key in line_header})
                            # print(restrains)
                        elif (line.startswith("#")):
                            continue
                        else:
                            print("WARNING! could not Read in :" + line)
                            continue
                    subcontent.update({"RESTRAINTHEADER": line_header})
                    subcontent.update({"RESTRAINTS": restrains})
                    result_data.update({block: subcontent})
                else:
                    raise IOError("DISRES parser does not know block: " + str(block) + "\n with content: " + "\n".join(
                        blocks[block]))
            return result_data

        with open(path, "r") as infile:
            lines = infile.readlines()
            # parse the coarse gromos_block structure
            blocks = self._read_gromos_block(lines)
            # parse data of the blocks
            data = _read_disres_subblock(blocks)
        return data

    def read_disres_file(self, path: str):
        # parse file into dicts
        """
        Args:
            path (str):
        """
        data = self.parse_disres_file(path)

        # convert distance_res lines to objects
        data["DISTANCERESSPEC"]["RESTRAINTS"] = list(
            map(lambda x: blocks.atom_pair_distanceRes(**x), data["DISTANCERESSPEC"]["RESTRAINTS"]))

        # add blocks as attribute to objects
        for key, sub_content in data.items():
            print(sub_content)
            self.add_block(blocktitle=key, content=sub_content)

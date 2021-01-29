#%%
import os, sys, glob
from collections import defaultdict
from pygromos.gromos.gromosPP import GromosPP
import restraintmaker

#CHANGE HERE
gromos_bin_path = "/home/bschroed/code/gromosPlsPls/gromos++/installed/bin"
restraintmaker_path = os.path.abspath(os.path.dirname(restraintmaker.__file__)+"/..")

#%%
gromPP = GromosPP(gromos_bin_path)
atb_dirs = restraintmaker_path+"/devtools/otherScripts/b_ATB_solvationFreeEnergies/ATB_molecules"
mstate_dir = restraintmaker_path+"/devtools/otherScripts/b_ATB_solvationFreeEnergies/sets/multistate"

sys.path.append(atb_dirs+"/..")
import utils_test_set_ATB as util

gromos_out_dir = restraintmaker_path+"/devtools/otherScripts/b_ATB_solvationFreeEnergies/sets/gromos_files"
ifp_path = gromos_out_dir+"/54A7.ifp"
mtb_path = gromos_out_dir+"/omni_unitedatom_atb_solv.mtb"
resn_lib_path = gromos_out_dir+"/resn_lib.lib"

#translation lib
from pygromos.files.otherfiles import residue_library

long_short = {}
mol_names = util.translate
for ID in mol_names:
    long_short.update({mol_names[ID]["short"]:[mol_names[ID]["orig_name"]]})
resn_lib = residue_library.residue_library()
resn_lib.RESIDUENAMELIB.pdb_top.update(long_short)
resn_lib.write(resn_lib_path)

#
for mol_dir in os.listdir(atb_dirs):
    mtb_path = glob.glob(atb_dirs+"/"+mol_dir+"/*.mtb")[0]
    mol_name = os.path.basename(mtb_path).split("_")[0]
    top_path = atb_dirs+"/"+mol_dir+"/"+mol_name+".top"
    gromPP.make_top(out_top_path=top_path, in_sequence=mol_name,
                    in_building_block_lib_path=mtb_path, in_parameter_lib_path=ifp_path, use_arg_file=False)




"""
# Single Tops
out_dir = root_dir+"/all_thirteen"
sequence = ["G27", "M09", "801", "6J2", "S00", "TVV", "E1V", "6KE", "M03",
            "M21", "F31", "T06","T07", "O71", "G07", "G20", "TP8"]
out_top_path = out_dir+"/out_all_thirteen.top"

gromPP.make_top(out_top_path=out_top_path, in_sequence=" ".join(sequence),
                in_building_block_lib_path=mtb_path, in_parameter_lib_path=ifp_path)
"""
#pdbs:
#gromPP.pdb2gromos(in_lib_path=)



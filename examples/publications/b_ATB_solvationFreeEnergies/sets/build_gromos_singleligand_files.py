#%%
import os, sys, glob
from collections import defaultdict
from pygromos.gromos.gromosPP import GromosPP
from pygromos.gromos.gromosXX import GromosXX
from pygromos.files import imd

import restraintmaker

#CHANGE HERE
gromos_bin_path = "/home/bschroed/Documents/code/gromosPP/installed/bin"
restraintmaker_path = os.path.abspath(os.path.dirname(restraintmaker.__file__)+"/..")

control_dict = {
    "gen_resn_lib": False,
    "gen_single_tops": False, #Buggy!
}
#%%
gromPP = GromosPP(gromos_bin_path)
atb_dirs = restraintmaker_path+"/devtools/otherScripts/b_ATB_solvationFreeEnergies/ATB_molecules"

sets_dir = restraintmaker_path+"/devtools/otherScripts/b_ATB_solvationFreeEnergies/sets"
mstate_dir = sets_dir+"/multistate"
pairwise_dir = sets_dir+"/pairwise"

if(not os.path.exists(mstate_dir)):
    os.mkdir(mstate_dir)

if (not os.path.exists(pairwise_dir)):
    os.mkdir(pairwise_dir)

sys.path.append(atb_dirs+"/..")
import utils_test_set_ATB as util


#RESNlib
resn_lib_path = sets_dir+"/resn_lib.lib"
if(control_dict['gen_resn_lib']):
    #translation lib
    from pygromos.files.otherfiles import residue_library

    long_short = {}
    mol_names = util.translate
    for ID in mol_names:
        long_short.update({mol_names[ID]["short"]:[mol_names[ID]["orig_name"]]})
    resn_lib = residue_library.residue_library()
    resn_lib.RESIDUENAMELIB.pdb_top.update(long_short)
    resn_lib.write(resn_lib_path)

#Generate TOPOS - THIS STEP IS BUGGY don't use, rather manual?
if (control_dict['gen_single_tops']):
    for mol_dir in os.listdir(atb_dirs):
        print(mol_dir)
        mtb_path = glob.glob(atb_dirs+"/"+mol_dir+"/*.mtb")[0]
        ifp_path = glob.glob(atb_dirs+"/"+mol_dir+"/*.ifp")[0]

        mol_name = "_"+os.path.basename(mtb_path).split("_")[1] if(not mtb_path.startswith("_")) else os.path.basename(mtb_path).split("_")[0]
        top_path = atb_dirs+"/"+mol_dir+"/"+mol_name+".top"

        gromPP.make_top(out_top_path=top_path, in_sequence=mol_name, in_solvent="H2O",
                        in_building_block_lib_path=mtb_path, in_parameter_lib_path=ifp_path, use_argfile=False)


#Systems
##get all_single file_tops:
all_tops = glob.glob(atb_dirs+"/*/*top")
state_all_tops={os.path.basename(value).split(".")[0]: value for value in all_tops}
all_mstate_sys = glob.glob(mstate_dir+"/*")

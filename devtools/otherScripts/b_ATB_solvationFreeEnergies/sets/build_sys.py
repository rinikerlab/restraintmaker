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
    "gen_multi_state": False,
    "gen_pairwise": True,
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

##multi_state different_sets
import utils_test_set_ATB as utils

if(control_dict['gen_multi_state']):
    for name,state_set in utils.multistate_ligand_sets.items():
        print(name, state_set)
        system_tops = [state_all_tops[x] for x in state_set]

        #build eds system
        out_dir = list(filter(lambda x: name+"_" in x, all_mstate_sys))[0]
        out_prefix_path = out_dir+"/"+os.path.basename(out_dir)
        out_top_path, out_ptp_path = gromPP.prep_eds(in_top_paths=system_tops, number_of_eds_states=len(system_tops), out_file_path=out_prefix_path)

        #generate cnf
        in_pdb = glob.glob(out_dir+"/*.pdb")[0]
        out_cnf_path = gromPP.pdb2gromos(in_pdb_path=in_pdb, in_top_path=out_top_path, in_lib_path=resn_lib_path,
                                         out_cnf_path=out_prefix_path+".cnf")
        print(out_cnf_path)

        #build posres/refpos
        from pygromos.files.coord import cnf
        cnf_file = cnf.Cnf(out_cnf_path)
        cnf_file.write_refpos(out_prefix_path+".rfp")
        cnf_file.write_possrespec(out_prefix_path+".por", residues=list(cnf_file.residues.keys()))

## PairWise Systems
if(control_dict['gen_pairwise']):
    all_states = utils.multistate_ligand_sets["all"]
    all_combos = os.listdir(pairwise_dir)
    all_combos = [x if(len(x)==2) else (x[0], "_"+x[2]) for x in list(map(lambda x: x.split("_"), all_combos))]
    print(all_combos)
    gromXX = GromosXX("/home/bschroed/Documents/code/gromosXX/installed/bin")
    for state_a, state_b in all_combos:
            system_tops = [state_all_tops[x] for x in [state_a, state_b]]

            name = state_a+"_"+state_b
            out_dir = pairwise_dir+"/"+name
            out_prefix_path = out_dir+"/"+name
            print(name, [state_a, state_b])

            if(not os.path.exists(out_dir)):
                os.mkdir(out_dir)

            #build dualTop system
            out_top_path, out_ptp_path = gromPP.prep_eds(in_top_paths=system_tops, number_of_eds_states=len(system_tops), out_file_path=out_prefix_path)

            #generate cnf
            in_pdb = glob.glob(out_dir+"/*.pdb")[0]
            out_cnf_path = gromPP.pdb2gromos(in_pdb_path=in_pdb, in_top_path=out_top_path, in_lib_path=resn_lib_path,
                                             out_cnf_path=out_prefix_path+".cnf")

            #build posres/refpos
            from pygromos.files.coord import cnf
            cnf_file = cnf.Cnf(out_cnf_path)
            refpos_path,_ = cnf_file.write_refpos(out_prefix_path+".rfp")
            posresspec_path,_ = cnf_file.write_possrespec(out_prefix_path+".por", residues=list(cnf_file.residues.keys()))

            print(out_cnf_path)

            #vac - emin
            from pygromos.data.imd_templates import template_emin
            from pygromos.files.blocks import imd_blocks

            in_imd = imd.Imd(template_emin)
            in_imd.STEP.NSTLIM = 10
            in_imd.SYSTEM.NSM = 0
            in_imd.FORCE.NRE =[12,28]
            in_imd_path = in_imd.write(out_path=out_dir+"/vac_emin.imd")

            gromXX.md_run(in_coord_path=out_cnf_path, in_topo_path=out_top_path, in_imd_path=in_imd_path,out_prefix=out_dir+"/emin",
                          in_pert_topo_path=out_ptp_path, in_refpos_path=refpos_path, in_posresspec_path=posresspec_path)
                          #in_disres_path=out_prefix_path+".disres")

            break

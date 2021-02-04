"""
This is a script, that starts a restraint maker session.

"""

import glob, os, time

import pymol
#pymol.finish_launching()
from pymol import cmd

#Alignment
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS

#project
import utils_test_set_ATB as utils


#Distance restraints
from restraintmaker.interface_Pymol.pymol_utilities import pymol_utitlities as up
from restraintmaker.algorithm import Filter, Optimizer
from restraintmaker.io import Exporter
from restraintmaker.utils import Utilities as u

#check pdbs
distance_treshold = 1.5
nrestraints = 4

molecule_dir= os.path.abspath("../../ATB_molecules")

all_pdbs = glob.glob(molecule_dir + "/*/*pdb")
state_all_pdbs = {os.path.basename(value).split(".")[0].replace("_unitedatom_optimised_geometry", ""): value for value in all_pdbs}

import restraintmaker
restraintmaker_path = os.path.abspath(os.path.dirname(restraintmaker.__file__)+"/..")
sets_dir = restraintmaker_path+"/devtools/otherScripts/b_ATB_solvationFreeEnergies/sets"
mstate_dir = sets_dir+"/multistate"

print(state_all_pdbs)
print(state_all_pdbs.keys())

#FUNCTION
## Aligning was cumbersome!
def align_mols_mcs_all(system_pdbs, align_to=0):

    ##Load mols
    mols = [Chem.MolFromPDBFile(pdb, removeHs=False) for pdb in system_pdbs]
    num_states = len(system_pdbs)
    ##Align with mcs
    ref = mols[align_to]
    print("ref:\t", align_to, os.path.basename(system_pdbs[align_to]))

    for mol2ID,mv in enumerate(mols):
        if(mol2ID==align_to):
            continue
        print("move:\t", mol2ID, os.path.basename(system_pdbs[mol2ID]))

        if(mol2ID == 4):
            mcs = rdFMCS.FindMCS([ref, mv], ringMatchesRingOnly=True)
        elif(mol2ID == 10):
            ref = mols[3]
            mcs = rdFMCS.FindMCS([ref, mv], ringMatchesRingOnly=True,ringCompare=rdFMCS.RingCompare.PermissiveRingFusion, atomCompare=rdFMCS.AtomCompare.CompareAny)
        elif(mol2ID == 12):
            ref = mols[10]
            mcs = rdFMCS.FindMCS([ref, mv], ringMatchesRingOnly=True)
        else:
            ref = mols[align_to]
            mcs = rdFMCS.FindMCS([ref, mv], ringMatchesRingOnly=True, atomCompare=rdFMCS.AtomCompare.CompareAnyHeavyAtom)

        smart = mcs.smartsString #"[#6&R]1-&@[#6&R](-&!@[#1&!R])-&@[#6&R]-&@[#6&R]-&@[#6&R]-&@[#6&R]" #mcs.smartsString
        patt = Chem.MolFromSmarts(smart)  # smartsString
        print("patternMol: ", mcs.smartsString)
        refMatch = ref.GetSubstructMatch(patt)
        print("refMatch:\t", refMatch)
        mvMatch = mv.GetSubstructMatch(patt)
        print("mvMatch:\t", mvMatch)

        try:
            AllChem.AlignMol(mv, ref, atomMap=list(zip(mvMatch, refMatch)),reflect=True)
        except Exception as err:
            raise err

    ##write out
    out_dir = "align"
    if (not os.path.exists(out_dir)):
        os.mkdir(out_dir)
    path_prefix = out_dir + "/aligned_"
    aligned_pdb_paths = []
    for mol1ID, in_pdb in enumerate(system_pdbs):
        base_name = os.path.basename(in_pdb)
        tmp_out = path_prefix + base_name
        Chem.MolToPDBFile(mols[mol1ID], tmp_out)
        aligned_pdb_paths.append(tmp_out)

    return aligned_pdb_paths

def align_mols_mcs_challenging(system_pdbs, align_to=3):

    ##Load mols
    mols = [Chem.MolFromPDBFile(pdb, removeHs=False) for pdb in system_pdbs]
    num_states = len(system_pdbs)
    ##Align with mcs
    print("ref:\t", align_to, os.path.basename(system_pdbs[align_to]))

    for mol2ID,mv in enumerate(mols):
        if(mol2ID==align_to):
            continue
        print("move:\t", mol2ID, os.path.basename(system_pdbs[mol2ID]))
        if(mol2ID == 4):
            print("AHA")
            ref= mols[2]
            mcs = rdFMCS.FindMCS([ref,mv], atomCompare=rdFMCS.AtomCompare.CompareAnyHeavyAtom)
        else:
            ref = mols[align_to]
            mcs = rdFMCS.FindMCS([ref, mv], ringMatchesRingOnly=True)

        smart = mcs.smartsString #"[#6&R]1-&@[#6&R](-&!@[#1&!R])-&@[#6&R]-&@[#6&R]-&@[#6&R]-&@[#6&R]" #mcs.smartsString
        patt = Chem.MolFromSmarts(smart)  # smartsString
        print("patternMol: ", mcs.smartsString)
        refMatch = ref.GetSubstructMatch(patt)
        print("refMatch:\t", refMatch)
        mvMatch = mv.GetSubstructMatch(patt)
        print("mvMatch:\t", mvMatch)

        try:
            AllChem.AlignMol(mv, ref, atomMap=list(zip(mvMatch, refMatch)),reflect=True)
        except Exception as err:
            raise err

        mv = mols[0]
        ref = mols[3]
        mcs = rdFMCS.FindMCS([ref, mv], ringMatchesRingOnly=True)

        smart = mcs.smartsString #"[#6&R]1-&@[#6&R](-&!@[#1&!R])-&@[#6&R]-&@[#6&R]-&@[#6&R]-&@[#6&R]" #mcs.smartsString
        patt = Chem.MolFromSmarts(smart)  # smartsString
        print("patternMol: ", mcs.smartsString)
        refMatch = ref.GetSubstructMatch(patt)
        print("refMatch:\t", refMatch)
        mvMatch = mv.GetSubstructMatch(patt)
        print("mvMatch:\t", mvMatch)

        try:
            AllChem.AlignMol(mv, ref, atomMap=list(zip(mvMatch, refMatch)),reflect=True)
        except Exception as err:
            raise err

    ##write out
    out_dir = "align"
    if (not os.path.exists(out_dir)):
        os.mkdir(out_dir)
    path_prefix = out_dir + "/aligned_"
    aligned_pdb_paths = []
    for mol1ID, in_pdb in enumerate(system_pdbs):
        base_name = os.path.basename(in_pdb)
        tmp_out = path_prefix + base_name
        Chem.MolToPDBFile(mols[mol1ID], tmp_out)
        aligned_pdb_paths.append(tmp_out)

    return aligned_pdb_paths

def align_mols_mcs(system_pdbs, align_to=0):

    ##Load mols
    mols = [Chem.MolFromPDBFile(pdb, removeHs=False) for pdb in system_pdbs]
    num_states = len(system_pdbs)
    ##Align with mcs
    ref = mols[align_to]
    print("ref:\t", align_to, os.path.basename(system_pdbs[align_to]))

    for mol2ID,mv in enumerate(mols):
        if(mol2ID==align_to):
            continue
        print("move:\t", mol2ID, os.path.basename(system_pdbs[mol2ID]))

        mcs = rdFMCS.FindMCS([ref, mv], ringMatchesRingOnly=True)

        smart = mcs.smartsString
        patt = Chem.MolFromSmarts(smart)  # smartsString
        print("patternMol: ", mcs.smartsString)
        refMatch = ref.GetSubstructMatch(patt)
        print("refMatch:\t", refMatch)
        mvMatch = mv.GetSubstructMatch(patt)
        print("mvMatch:\t", mvMatch)

        try:
            AllChem.AlignMol(mv, ref, atomMap=list(zip(mvMatch, refMatch)))
        except Exception as err:
            raise err

    ##write out
    out_dir = "align"
    if (not os.path.exists(out_dir)):
        os.mkdir(out_dir)
    path_prefix = out_dir + "/aligned_"
    aligned_pdb_paths = []
    for mol1ID, in_pdb in enumerate(system_pdbs):
        base_name = os.path.basename(in_pdb)
        tmp_out = path_prefix + base_name
        Chem.MolToPDBFile(mols[mol1ID], tmp_out)
        aligned_pdb_paths.append(tmp_out)

    return aligned_pdb_paths


############################################3
#DO
for set_name_prefix in  utils.multistate_ligand_sets:

    print(set_name_prefix)
    #if(set_name_prefix in ["all", "challenging", "singles", "flat"]):
    #    continue

    state_set = utils.multistate_ligand_sets[set_name_prefix]
    system_pdbs = [state_all_pdbs[x] for x in state_set]
    num_states = len(system_pdbs)
    set_name = set_name_prefix+"_"+str(num_states)
    print("set: ", set_name,"\n\t", system_pdbs)

    #outdir:
    out_dir_path = mstate_dir + "/" + set_name
    if(not os.path.exists(out_dir_path)):
        os.mkdir(out_dir_path)

    # Algin
    if(set_name_prefix in ["all"] ):
        aligned_pdb_paths_list = align_mols_mcs_all(system_pdbs, align_to=5)
    elif(set_name_prefix in ["challenging"]):
        aligned_pdb_paths_list =align_mols_mcs_challenging(system_pdbs, align_to=2)
    else:
        aligned_pdb_paths_list = align_mols_mcs(system_pdbs, align_to=0)


    # fire up pymol
    pymol.finish_launching()

    # load data
    [cmd.load(pdb, obj_name) for obj_name, pdb in zip(map(lambda x: x.replace("_", "T"), state_set), aligned_pdb_paths_list)]
    time.sleep(1)

    # set nice scene
    obj_list = cmd.get_object_list()
    for i, obj in enumerate(obj_list):
        cmd.alter(obj, "chain=" + str(i))
        cmd.alter(obj, "resi=" + str(i))

        if (obj == "TO71"):
            cmd.alter(obj, "resn=\"O71\"")

        cmd.alter(obj, "resn=resn.replace(\"_\", \"T\")")

    #exit()
    out_pdb = out_dir_path+"/"+set_name+".pdb"
    cmd.save(out_pdb)
    cmd.center()
    cmd.ray(1200)
    cmd.png(out_dir_path+"/"+set_name+"_overlay.png")
    cmd.set("grid_mode", 1)
    cmd.ray(1200)
    cmd.png(out_dir_path+"/"+set_name+"_grid.png")

    cmd.reinitialize()
    cmd.load(out_pdb)

    ## GET ATOMS
    atom_list = up.pymol_selection_to_atom_list("all")

        ## Filtering for rings
    try:
        RingFilter = Filter.RingFilter(atom_list)
        #pdb_blocks = u.convert_atoms_to_pdb_molecules(atom_list)
        pdb_blocks = [cmd.get_pdbstr(pdb_id) for pdb_id in cmd.get_object_list()]
        RingFilter.get_args(lambda x: (pdb_blocks))
        filtered_atoms = RingFilter.filter()
    except Exception as err:
        print("failed ", err.args)

        ## Optimizers
    try:
        opt = Optimizer.TreeHeuristicOptimizer(filtered_atoms)
        opt.get_args(lambda x: (nrestraints, distance_treshold, 'shortest', None))
        res = opt.make_restraints()
        ## Export
        exporter = Exporter.Gromos_Distance_Restraint_Exporter(restraints=res)
        exporter.get_args(lambda x: out_dir_path + "/" + set_name + ".disres")
        exporter.export_restraints()
    except:
        print("manual mode:")
        import restraintmaker
        restraintmaker.run_plugin_gui()

    break
exit()
"""
This is a script, that starts a restraint maker session.

"""

import glob, os, time

import pymol
pymol.finish_launching()
from pymol import cmd

#RDKIT - Mol Alignment
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS

#Distance restraints
from restraintmaker.interface_Pymol.pymol_utilities import pymol_utitlities as up
from restraintmaker.algorithm import Filter, Optimizer
from restraintmaker.io import Exporter
from restraintmaker.utils import Utilities as u

#check pdbs
molecule_dir= "../../ATB_molecules"
pdbs = glob.glob(molecule_dir+"/*/*pdb")

aligned_mol="aligned.pdb"
orig_pdbs = [pdb for pdb in pdbs ]
print(orig_pdbs)

obj_names = [os.path.basename(pdb).split("_")[1] if (len(os.path.basename(pdb).split("_")[0]) ==0) else os.path.basename(pdb).split("_")[0] for pdb in pdbs]
all_ligs = obj_names

from itertools import combinations
all_combos = list(combinations(all_ligs,2))

M030_TI_pairwise = list(filter(lambda x: "M030" in x, all_combos))

REEDS_sub_combos = ['G277', "_O6T", "M097", "6KET", "F313",]
all_add_reeds = list(combinations(REEDS_sub_combos,2))

#selected_pairs = M030_TI_pairwise
selected_pairs = all_add_reeds
#selected_pairs = [('M030', '_O6T'), ('M030', '_O71'), ('M030', '_O70'), ('M030', "_P8I")]
#selected_pairs = [("M030", "G078")]
#selected_pairs = [("M097", "G277")]
#selected_pairs = [("M030", "G209"), ("8018", "M030")]
#selected_pairs = (("6KET", "_O6T"), ("M097", "_O6T"),("G277", "_O6T"), ("F313", "G277"), ("F313", "_O6T"))
selected_pairs = (('M030', '_O6T')) # list(filter( lambda x: not "G209" in x, M030_TI_pairwise))


print(M030_TI_pairwise)



def vis( res, path_prefix, c=["firebrick", "forest", 'purple', 'salmon', "gold", "red"]):


    for ind, r in enumerate(res):
        for a in r.atoms:
            cmd.show("spheres", "id " + str(a.id))
            cmd.set("sphere_color", c[ind], "id " + str(a.id))
        cmd.distance("d" + str(ind), "id " + str(r.atoms[0].id), "id " + str(r.atoms[1].id))
        cmd.hide("labels")
        cmd.set("grid_slot", -2, "d" + str(ind))

    cmd.ray(1200)
    time.sleep(2)
    cmd.set("grid_mode", 0)
    cmd.png(path_prefix + "_restraints.png")
    #cmd.set("grid_mode", 1)
    cmd.move("zoom", -10)


    cmd.ray(1200,)
    time.sleep(4)
    cmd.png(path_prefix + "_grid_restraints.png")
    cmd.set("grid_mode", 0)
    cmd.hide("spheres")
    cmd.delete("d*")



for indA, molA in enumerate(orig_pdbs):
    for indB, molB in enumerate(orig_pdbs):
        molA_name_partI = os.path.basename(molA).split(".")[0]
        molB_name_partI = os.path.basename(molB).split(".")[0]

        molA_name = "_"+molA_name_partI.split("_")[1] if(molA_name_partI.startswith("_")) else molA_name_partI.split("_")[0]
        molB_name = "_"+molB_name_partI.split("_")[1] if(molB_name_partI.startswith("_")) else molB_name_partI.split("_")[0]

        out_prefix = molA_name+"_"+molB_name
        out_dir = "pictures/"+out_prefix  #os.getcwd()+"/"+out_prefix

        print(molA_name)
        print("\t", molB_name)

        #if((molA_name, molB_name) in selected_pairs or (molB_name, molA_name) in selected_pairs ):
        #    cmd.reinitialize()
        #    time.sleep(1)
        #    print(out_prefix)
        #else:
        #    continue

        if(not os.path.exists(out_dir)):
            os.mkdir(out_dir)

        print(out_dir)

        # ALIGN THE TWO MOLS
        ##Load mols
        mols = [Chem.MolFromPDBFile(pdb, removeHs=False) for pdb in [molA, molB]]

        ##Align with mcs
        ref = mols[0]
        mv = mols[1]
        #mcs = rdFMCS.FindMCS([ref, mv], completeRingsOnly=True, matchValences=True, ringMatchesRingOnly=True) # G078
        mcs = rdFMCS.FindMCS([ref, mv],  completeRingsOnly=True, ringMatchesRingOnly=True)

        smartsString = mcs.smartsString
        #from rdkit.Chem import MCS
        #smartsString = MCS.FindMCS(mols, atomCompare="any").smarts

        patt = Chem.MolFromSmarts(smartsString)  # smartsString
        refMatch = ref.GetSubstructMatch(patt)
        mvMatch = mv.GetSubstructMatch(patt)

        try:
            AllChem.AlignMol(mv, ref, atomMap=list(zip(mvMatch, refMatch)))
        except Exception as err:
            print(err.args)
            pass

        ##write out
        out_pdb_path = out_dir + "/" + out_prefix + ".pdb"
        out_text = ""
        for mol in mols:
            out_text+=Chem.MolToPDBBlock(mol)


        file_out = open(out_pdb_path, "w")
        file_out.write(out_text)
        file_out.close()


        ###############################
        #BUILD DISRES

        distance_treshold=1.0
        nrestraints = 4

            ##load data
        cmd.load(out_pdb_path)
        time.sleep(1)


        # set nice scene
        obj_list = cmd.get_object_list()
        cmd.set("pdb_retain_ids", 0)
        offset = 0
        for i, obj in enumerate(obj_list):
            cmd.alter(obj, "chain="+str(i)) #fix chain
            cmd.alter(obj, "resi="+str(i))  #fix res numbers
            cmd.alter(obj, "ID=ID+"+str(offset))    #generate unique IDS
            offset = len(cmd.get_model(obj).atom)

            if(obj =="O71"):
                cmd.alter(obj, "resn=\"O71\"") #a very special case

            cmd.alter(obj, "resn=resn.replace(\"_\", \"T\")")   #replace underscores
        cmd.sync()
        cmd.save(out_pdb_path)
        print(out_pdb_path)

        ## GET ATOMS
        atom_list = up.pymol_selection_to_atom_list("all")

            ## Filtering for rings
        try:
            RingFilter = Filter.RingFilter(atom_list)
            pdb_blocks = [cmd.get_pdbstr(obj) for obj in cmd.get_object_list()] #u.convert_atoms_to_pdb_molecules(atom_list)
            RingFilter.get_args(lambda x: (pdb_blocks))
            filtered_atoms = RingFilter.filter()

            print(filtered_atoms)
        except Exception as err:
            print("failed ", err.args)
            continue

            ## Optimizers
        opt = Optimizer.TreeHeuristicOptimizer(filtered_atoms)
        opt.get_args(lambda x: (nrestraints, distance_treshold, 'shortest', None))
        res = opt.make_restraints()

            ## Export
        exporter = Exporter.Gromos_Distance_Restraint_Exporter(restraints=res)
        exporter.get_args(lambda x: out_dir+"/"+out_prefix+".disres")
        exporter.export_restraints()


        ## additional:
        #import restraintmaker
        #restraintmaker.run_plugin_gui()
        ##
        cmd.bg_color("white")
        cmd.set("stick_radius", 0.1)
        cmd.set("sphere_scale", 0.2)
        cmd.color("vanadium", "elem C")
        cmd.color("copper", "ID " + "+".join([str(a.id) for a in filtered_atoms]))

        cmd.ray(1200,800)
        cmd.png(out_dir+"/"+out_prefix+"_overlay.png")
        cmd.set("grid_mode", "1")
        cmd.ray(1200,800)
        cmd.png(out_dir+"/"+out_prefix+"_grid.png")
        cmd.set("grid_mode", "0")
        vis(res,out_dir+"/"+out_prefix)


print("fini")
exit()

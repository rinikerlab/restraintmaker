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
obj_names = [os.path.basename(pdb).split("_")[1] if (len(os.path.basename(pdb).split("_")[0]) ==0) else os.path.basename(pdb).split("_")[0] for pdb in pdbs]

all_ligs = ["_O6T", "G277", "M097", "6KET", "F313"]
from itertools import combinations
all_combos = list(combinations(all_ligs,2))
print(all_combos)

print(obj_names)
all_combos = [("8018", "6J29")]
#
#("M097", "F313")]#, ("G277", "6KET"), ("_O6T", "F313"), ("_O6T", "6KET"), ("6KET", "F313")]

for indA, molA in enumerate(orig_pdbs):

    for indB, molB in enumerate(orig_pdbs):
        molA_name_partI = os.path.basename(molA).split(".")[0]
        molB_name_partI = os.path.basename(molB).split(".")[0]

        molA_name = "_"+molA_name_partI.split("_")[1] if(molA_name_partI.startswith("_")) else molA_name_partI.split("_")[0]
        molB_name = "_"+molB_name_partI.split("_")[1] if(molB_name_partI.startswith("_")) else molB_name_partI.split("_")[0]

        out_prefix = molA_name+"_"+molB_name
        out_dir = os.getcwd()+"/"+out_prefix

        print(molA_name, molB_name)

        if((molA_name, molB_name) in all_combos):
            cmd.reinitialize()
            print(out_prefix)
        else:
            continue

        if(not os.path.exists(out_dir)):
            os.mkdir(out_dir)

        print(out_dir)

        # ALIGN THE TWO MOLS
        ##Load mols
        mols = [Chem.MolFromPDBFile(pdb, removeHs=False) for pdb in [molA, molB]]

        ##Align with mcs
        ref = mols[0]
        mv = mols[1]
        mcs = rdFMCS.FindMCS([ref, mv])
        patt = Chem.MolFromSmarts(mcs.smartsString)  # smartsString
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

        print(out_text)
        print(out_pdb_path)

        file_out = open(out_pdb_path, "w")
        file_out.write(out_text)
        file_out.close()


        ###############################
        #BUILD DISRES

        distance_treshold=1.5
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

        cmd.ray(1200,800)
        cmd.png(out_dir+"/"+out_prefix+"_overlay.png")
        cmd.set("grid_mode", "1")
        cmd.ray(1200,800)
        cmd.png(out_dir+"/"+out_prefix+"_grid.png")
        cmd.set("grid_mode", "0")

exit()

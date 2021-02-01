"""
This is a script, that starts a restraint maker session.

"""

import glob, os, time

import pymol
#pymol.finish_launching()
from pymol import cmd

#RDKIT - Mol Alignment
from rdkit import Chem
from rdkit.Chem import AllChem
# from rdkit.Chem import rdFMCS
from rdkit.Chem import MCS as rdFMCS

#Distance restraints
from restraintmaker.interface_Pymol.pymol_utilities import pymol_utitlities as up
from restraintmaker.algorithm import Filter, Optimizer
from restraintmaker.io import Exporter
from restraintmaker.utils import Utilities as u

#check pdbs
molecule_dir= "../../ATB_molecules"
pdbs = glob.glob(molecule_dir+"/*/*pdb")
aligned_mol="aligned.pdb"
orig_pdbs = [pdb for pdb in pdbs  ]
obj_names = [os.path.basename(pdb).split("_")[1] if (len(os.path.basename(pdb).split("_")[0]) ==0) else os.path.basename(pdb).split("_")[0] for pdb in pdbs]


print(obj_names)
for indA, molA in enumerate(orig_pdbs):
    for indB, molB in enumerate(orig_pdbs):
        molA_name_partI = os.path.basename(molA).split(".")[0]
        molB_name_partI = os.path.basename(molB).split(".")[0]

        molA_name = "_"+molA_name_partI.split("_")[1] if(molA_name_partI.startswith("_")) else molA_name_partI.split("_")[0]
        molB_name = "_"+molB_name_partI.split("_")[1] if(molB_name_partI.startswith("_")) else molB_name_partI.split("_")[0]

        out_prefix = molA_name+"_"+molB_name
        out_dir = os.getcwd()+"/"+out_prefix

        if(#molB_name == "M030" or not molA_name == "M030" ):
            not (#(molA_name == "8018" and molB_name == "G078") or
                 #(molA_name == "8018" and molB_name == "6J29") or
                (molA_name == "6J29" and molB_name == "G078"))):
            continue
        else:
            print(out_prefix)


        if(not os.path.exists(out_dir)):
            os.mkdir(out_dir)

        #print(out_dir)

        # ALIGN THE TWO MOLS
        ##Load mols
        mols = [Chem.MolFromPDBFile(pdb) for pdb in [molA, molB]]

        ##Align with mcs
        ref = mols[0]
        mv = mols[1]
        mcs = rdFMCS.FindMCS([ref, mv], ringMatchesRingOnly=True)
        patt = Chem.MolFromSmarts(mcs.smarts)  # smartsString
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

        #print(out_text)

        file_out = open(out_pdb_path, "w")
        file_out.write(out_text)
        file_out.close()


        ###############################
        #BUILD DISRES

        distance_treshold=1.5

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

        ## additional:
        import restraintmaker

        restraintmaker.run_plugin_gui()
        exit()

        ## GET ATOMS
        atom_list = up.pymol_selection_to_atom_list("all")

            ## Filtering for rings
        try:
            RingFilter = Filter.RingFilter(atom_list)
            pdb_blocks = u.convert_atoms_to_pdb_molecules(atom_list)
            RingFilter.get_args(lambda x: (pdb_blocks))
            filtered_atoms = RingFilter.filter()
        except Exception as err:
            print("failed ", err.args)
            continue

            ## Optimizers
        opt = Optimizer.TreeHeuristicOptimizer(filtered_atoms)
        opt.get_args(lambda x: (4, distance_treshold, 'shortest', None))
        res = opt.make_restraints()

            ## Export
        exporter = Exporter.Gromos_Distance_Restraint_Exporter(restraints=res)
        exporter.get_args(lambda x: out_dir+"/"+out_prefix+".disres")
        exporter.export_restraints()

        ##
        cmd.sync()
        cmd.reinitialize()
        time.sleep(1)

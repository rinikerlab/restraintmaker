"""
This is a script, that starts a restraint maker session.

"""

import glob, os, time

import pymol
pymol.finish_launching()

from pymol import cmd
#check pdbs
molecule_dir= "../../ATB_molecules"
pdbs = glob.glob(molecule_dir+"/*/*pdb")
aligned_mol="aligned.pdb"
orig_pdbs = [pdb for pdb in pdbs  ]
obj_names = [os.path.basename(pdb).split("_")[1] if (len(os.path.basename(pdb).split("_")[0]) ==0) else os.path.basename(pdb).split("_")[0] for pdb in pdbs]

print(obj_names)

#orig_pdbs = [orig_pdbs[0], orig_pdbs[1]]+ [orig_pdbs[i] for i in range(4,15+1)] #The simple thirteen
#obj_names = [obj_names[0]]+ [obj_names[i] for i in range(3,15+1)]

#orig_pdbs = [orig_pdbs[0], orig_pdbs[2], orig_pdbs[3], orig_pdbs[11], orig_pdbs[12]] #The challenging five

#orig_pdbs = [orig_pdbs[i] for i in range(4, 11)] + [orig_pdbs[0], orig_pdbs[1]] #The nine singles

#orig_pdbs = [orig_pdbs[i] for i in [0,1,3,4,5,7,8,9,10,11]] #The challenging five

#orig_pdbs = [orig_pdbs[i] for i in [0,1,4,7,8]] #The easy five

#Alignment
from rdkit import Chem
from rdkit.Chem import AllChem
#from rdkit.Chem import rdFMCS
from rdkit.Chem import MCS as rdFMCS

##Load mols
num_orig_pdbs = len(orig_pdbs)
mols = [Chem.MolFromPDBFile(pdb) for pdb in orig_pdbs]
print(len(mols), print(len(pdbs)))
print(mols)
##Align with mcs
for mol1ID in range(0,num_orig_pdbs):
    ref = mols[mol1ID]
    for mol2ID in range(mol1ID+1, num_orig_pdbs):
        mv = mols[mol2ID]

        mcs = rdFMCS.FindMCS([ref, mv], ringMatchesRingOnly=True)
        patt = Chem.MolFromSmarts(mcs.smarts)  # smartsString
        refMatch = ref.GetSubstructMatch(patt)
        mvMatch = mv.GetSubstructMatch(patt)
        #print(mvMatch, refMatch)
        try:
            AllChem.AlignMol(mv, ref, atomMap=list(zip(mvMatch, refMatch)))
        except Exception as err:
            print(err.args)
            pass
    break

##write out
out_dir = "align"
if(not os.path.exists(out_dir)):
    os.mkdir(out_dir)
path_prefix = out_dir+"/aligned_mol"

pdbs = []
for mol1ID in range(num_orig_pdbs):
    tmp_out = path_prefix+"_"+str(mol1ID)+".pdb"
    Chem.MolToPDBFile(mols[mol1ID], tmp_out)
    pdbs.append(tmp_out)

#fire up pymol
pymol.finish_launching()

#load data
[cmd.load(pdb, obj_name) for obj_name,pdb in zip(obj_names, pdbs)]
time.sleep(1)

# set nice scene
obj_list = cmd.get_object_list()
for i, obj in enumerate(obj_list):
    cmd.alter(obj, "chain="+str(i))
    cmd.alter(obj, "resi="+str(i))

    if(obj =="O71"):
        cmd.alter(obj, "resn=\"O71\"")

    cmd.alter(obj, "resn=resn.replace(\"_\", \"T\")")

cmd.save("merged.pdb")

#For Restraint generation
cmd.reinitialize()
cmd.load("merged.pdb")

import restraintmaker
restraintmaker.run_plugin_gui()
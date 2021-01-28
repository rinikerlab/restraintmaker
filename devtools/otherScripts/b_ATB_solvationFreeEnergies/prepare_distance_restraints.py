"""
This is a script, that starts a restraint maker session.

"""

import glob, os, time

import pymol
pymol.finish_launching()

from pymol import cmd
#check pdbs
molecule_dir= "ATB_molecules"
pdbs = glob.glob(molecule_dir+"/*/*pdb")
orig_pdbs = [pdb for pdb in pdbs if (not len(os.path.basename(pdb).split("_")[0]) ==0) ]
obj_names = [os.path.basename(pdb).split("_")[1] if (len(os.path.basename(pdb).split("_")[0]) ==0) else os.path.basename(pdb).split("_")[0] for pdb in pdbs]

#Alignment

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS

def alignLigs(ligands, lig1, lig2):
  print(str(lig1) + " " + str(lig2))
  print(mcs)
  #print(Chem.MolToSmiles(patt))
  ref = ligands[lig1]
  refMatch = ref.GetSubstructMatch(patt)
  #print(refMatch)
  mv = ligands[lig2].GetSubstructMatch(patt)
  #print(mv)
  AllChem.AlignMol(ligands[lig2], ref, atomMap=list(zip(mv, refMatch)))

num_orig_pdbs = len(orig_pdbs)
mols =  [Chem.MolFromPDBFile(pdb) for pdb in orig_pdbs]

mcs = rdFMCS.FindMCS([mols[0], mols[1]], ringMatchesRingOnly=False,)
patt = Chem.MolFromSmarts(mcs.smartsString)
print(mcs.smartsString)
print(orig_pdbs)


from rdkit.Chem import Subshape
from rdkit.Chem.Subshape import SubshapeAligner
from rdkit.Chem.Subshape import SubshapeBuilder


for mol1ID in range(num_orig_pdbs):
    ref = mols[mol1ID]
    refMatch = ref.GetSubstructMatch(patt)
    for mol2ID in range(mol1ID+1, num_orig_pdbs):
        mv = mols[mol2ID]
        mvMatch = mv.GetSubstructMatch(patt)
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
    cmd.alter(obj, "resn=\"mol" + str(i)+"\"")

exit()


cmd.save("merged.pdb")


#molecule_dir = "ATB_molecules"

cmd.reinitialize()

cmd.load("merged.pdb")
exit()
import restraintmaker
restraintmaker.run_plugin_gui()
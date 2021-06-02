"""
This is a script, that starts a restraint maker session.

"""
import datetime

import tqdm
from itertools import permutations
import glob, os, time

import numpy as np
import pandas as pd
import pymol
from pymol import cmd
from scipy.spatial import ConvexHull

#Distance restraints
from restraintmaker.interface_Pymol.pymol_utilities import pymol_utitlities as up
from restraintmaker.algorithm import Filter, Optimizer


#check pdbs
pairwise_files = "/home/bschroed/Documents/code/restraintmaker/devtools/otherScripts/b_ATB_solvationFreeEnergies/sets/pairwise"

pairs = glob.glob(pairwise_files+"/*")
pairs = list(filter(lambda x: "M030" in x, pairs))

pymol.finish_launching()
volumes=  {}
for ind, pair_path in enumerate(pairs):
    cmd.reinitialize()
    cmd.bg_colour("white")
    cmd.set("label_color", "black")
    if(os.path.isfile(pair_path)):
        continue

    pair_str = os.path.basename(pair_path)

    molA_name = pair_str.split("_")[0]
    molB_name = pair_str.split("_")[1] if(len(pair_str.split("_")) == 2) else  "_"+ pair_str.split("_")[2]
    print(molA_name, molB_name)

    out_pdb_path = pair_path+"/"+pair_str+".pdb"

    ###############################
    #BUILD DISRES

    distance_treshold=1.5
    nrestraints = 4

    ##load data
    cmd.set("retain_order", 1)
    cmd.load(out_pdb_path)
    cmd.show("lines", "all")
    cmd.hide("sticks")
    cmd.set("sphere_scale", 0.2)
    time.sleep(2)


    # set nice scene
    obj_list = cmd.get_object_list()
    cmd.create("res1", "resn "+molA_name[:3])
    cmd.create("res2", "resn "+molB_name[:3].replace("_", "T"))

    cmd.delete(obj_list[0])
    obj_list = cmd.get_object_list()
    cmd.set("pdb_retain_ids", 0)

    first = True
    for i, obj in enumerate(obj_list):
        cmd.alter(obj, "chain="+str(i+1)) #fix chain
        cmd.alter(obj, "resi="+str(i+1))  #fix res numbers
        if(first):
            #cmd.alter(obj, "ID=ID")    #generate unique IDS
            first = False
        else:
            pass
           # cmd.alter(obj, "ID=ID+"+str(1))    #generate unique IDS

        cmd.alter(obj, "resn=resn.replace(\"_\", \"T\")")   #replace underscores

    cmd.label("all", "ID")
    #cmd.set("grid_mode")
    cmd.sync()

    ## GET ATOMS
    atom_list = up.pymol_selection_to_atom_list("all")

    ## Filtering for rings
    try:
        RingFilter = Filter.RingFilter(atom_list)
        pdb_blocks = [cmd.get_pdbstr(obj) for obj in cmd.get_object_list()]  # u.convert_atoms_to_pdb_molecules(atom_list)
        RingFilter.get_args(lambda x: (pdb_blocks))
        filtered_atoms = RingFilter.filter()
    except Exception as err:
        print(len(pdb_blocks))
        raise err
        print("failed ", err.args)
        continue

    # Optimizers

    ## Greedy
    startt = datetime.datetime.now()
    opt = Optimizer.TreeHeuristicOptimizer(filtered_atoms)
    opt.get_args(lambda x: (nrestraints, distance_treshold, 'shortest', None))
    res = opt.make_restraints()
    endt = datetime.datetime.now()
    delta_t_greed = (endt-startt).seconds

    """
    ### Visualization
    c = ["firebrick", "forest", 'purple', 'salmon']

    for ind,r in enumerate(res):
        for a in r.atoms:
            cmd.show("spheres", "id "+str(a.id))
            cmd.set("sphere_color", c[ind], "id " + str(a.id))
        cmd.distance("d"+str(ind), "id "+str(r.atoms[0].id), "id "+str(r.atoms[1].id))
        cmd.set("grid_slot", -2, "d"+str(ind))
    """
    #cmd.ray(1200)
    #time.sleep(4)
    #cmd.save(pair_str+"_greedy_restraints.png")


    #cmd.set("grid_mode", 1)
    #cmd.move("zoom", -20)

    #cmd.ray(1200)
    #time.sleep(4)
    #cmd.save(pair_str+"_grid_greedy_restraints.png")
    #cmd.set("grid_mode", 0)

    ###measure convex Hull:
    atoms = []
    for r in res:
        a1 = r.atoms[0]
        a2 = r.atoms[1]
        cog_atom = np.array([a1.x+a2.x, a1.y+a2.y, a1.z+a2.z])/2
        atoms.append(cog_atom)

    d_greed = np.round(np.sum([np.sum([np.sqrt(np.sum(np.square([a1[0]-a2[0], a1[1]-a2[1], a1[2]-a2[2]]))) for a2 in atoms[ind:]]) for ind, a1 in enumerate(atoms)]),5)
    print(d_greed)
    print([[np.sqrt(np.sum(np.square([a1[0]-a2[0], a1[1]-a2[1], a1[2]-a2[2]]))) for a2 in atoms[ind:]] for ind, a1 in enumerate(atoms)])

    v_greed = ConvexHull(atoms).volume
    print(v_greed, Optimizer._calculate_value_convex_hull(atoms))
    ## Brute Force
    opt = Optimizer.BruteForceRingOptimzer(filtered_atoms)
    opt.get_args(lambda x: (4, 1.5, 'convex_hull', 'None'))
    res = opt.make_restraints()

    c = ["wheat", "palegreen", 'lightblue', 'lightorange']
    for ind,r in enumerate(res):
        for a in r.atoms:
            cmd.show("spheres", "id "+str(a.id))
            cmd.set("sphere_color",c[ind], "id " + str(a.id))
        cmd.distance("d"+str(ind), "id "+str(r.atoms[0].id), "id "+str(r.atoms[1].id))
        cmd.set("grid_slot", -2, "d"+str(ind))


    #cmd.ray(1200)
    #time.sleep(4)
    #cmd.save(pair_str+"_bruteF_restraints.png")

    #cmd.set("grid_mode", 1)
    #cmd.move("zoom", -20)

    #cmd.ray(1200)
    #time.sleep(4)
    #cmd.save(pair_str+"_grid_bruteF_restraints.png")
    #cmd.set("grid_mode", 0)

    ###measure convex Hull:
    atoms = []
    for r in res:
        a1 = r.atoms[0]
        a2 = r.atoms[1]
        cog_atom = np.array([a1.x+a2.x, a1.y+a2.y, a1.z+a2.z])/2
        atoms.append(cog_atom)
    v_bruteF = ConvexHull(atoms).volume

    atoms = []
    all_possible_restraints = []
    all_res = {}
    t = []
    ind = 0
    for i,a1 in enumerate(filtered_atoms):
        for a2 in filtered_atoms[i:]:
            d =  np.sqrt(np.sum(np.square(np.array([a1.x-a2.x, a1.y-a2.y, a1.z-a2.z]))))
            if(a1 != a2 and d < 1.5 and a1.resn != a2.resn):
                cog_atom = np.array([a1.x + a2.x, a1.y + a2.y, a1.z + a2.z]) / 2

                all_possible_restraints.append(cog_atom)
                atoms_ind = [a1.id, a2.id]
                all_res.update({ind: {"atoms":(a1, a2), "atomsind":atoms_ind, "cog": cog_atom}})
                ind +=1

    v_mol = ConvexHull(all_possible_restraints).volume

    #Real brute force
    startt = datetime.datetime.now()
    print("possible restraints ", len(all_res))
    nRestraints = 4

    molecule_ids = []
    cmd.iterate("res1", "molecule_ids.append(ID)", space={"molecule_ids":molecule_ids})
    all_res_ind = list(all_res.keys())
    ind_combinations = permutations(all_res_ind, nRestraints)
    cleaned_combinations = {}
    bruteF_distance_d = {}
    selected = []
    iterator = tqdm.tqdm(ind_combinations, leave=False)
    for res_comb in iterator:
        unique_atoms =set(sorted(np.unique(np.concatenate(list(map(lambda x: all_res[x]['atomsind'], res_comb))))))
        if(len(unique_atoms) == nRestraints *2 and not unique_atoms in selected):
            selected.append(unique_atoms)
            cogs = np.array(list(map(lambda x: all_res[x]['cog'], res_comb)))
            d = np.round(np.sum([np.sum([np.sqrt(np.sum(np.square([a1[0]-a2[0], a1[1]-a2[1], a1[2]-a2[2]]))) for a2 in cogs[ind:]]) for ind, a1 in enumerate(cogs)]),5)
            cleaned_combinations.update({res_comb: round(ConvexHull(cogs).volume,4)})
            bruteF_distance_d.update({res_comb: d})

    max_key = max(cleaned_combinations, key=lambda x: cleaned_combinations[x])
    max_val = cleaned_combinations[max_key]
    max_atoms = list(map(lambda x: all_res[x]['atomsind'], max_key))
    vol_dist = list(cleaned_combinations.values())
    print("all 4 res_combinations: ", len(cleaned_combinations))
    print("maxV: ", max_key, max_val, max_atoms)
    greed_rank = round(len([x for x in cleaned_combinations.values() if(x > v_greed)])/len(cleaned_combinations),2)*100
    print("Rank: ", greed_rank)
    endt = datetime.datetime.now()
    deltat_bruteF = (endt - startt).seconds


    print(molecule_ids)
    c = ["wheat", "palegreen", 'lightblue', 'lightorange']
    cmd.hide("spheres")
    for i,r in enumerate(max_atoms):
        for a in r:
            cmd.show("spheres", "id "+str(a))
            cmd.set("sphere_color",c[i], "id " + str(a))
        cmd.distance("d"+str(ind), "id "+str(r[0]), "id "+str(r[1]))
        cmd.set("grid_slot", -2, "d"+str(ind))

    cmd.ray(1200)
    time.sleep(4)
    cmd.save(pair_str+"_mbruteF_restraints.png")

    cmd.set("grid_mode", 1)
    cmd.move("zoom", -20)

    cmd.ray(1200)
    time.sleep(4)
    cmd.save(pair_str+"_grid_mbruteF_restraints.png")
    cmd.set("grid_mode", 0)
    #time.sleep(3)

    for cog_atom in all_possible_restraints:
        cmd.pseudoatom(object=obj, name="Pseudo", resi=12, resn="Pseudo", pos=tuple(cog_atom))

    cmd.show("nb_spheres", "resn Pseudo")
    cmd.color("lightblue", "resn Pseudo")

    print("bruteF: ", v_bruteF,
          "greed", v_greed,
          "vMol", v_mol,
          "rank", greed_rank)

    if(v_bruteF > v_mol):
        print("bRUTEF was bigger than tot space?")
        exit(0)
    if( v_greed > v_bruteF):
        print("greed was bigger than bruteF?")
        exit(0)


    volumes.update({pair_str: {"brutF": np.round(v_bruteF,4), "greed": np.round(v_greed,4), "mol": np.round(v_mol,4), 'rank':greed_rank, "myBruteF": np.round(max_val,4), "myBruteF_t":deltat_bruteF, "greed_t":delta_t_greed, "myBruteF_dist":vol_dist,
                               "bruteF_d":list(bruteF_distance_d.values()), "greed_d": d_greed }})
    pd.DataFrame(volumes).T.to_csv("./tmp_volume_ana.tsv", sep="\t")
    #exit()

print("fini")
exit()

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
from pygromos.utils import bash
#Distance restraints
from restraintmaker.interface_Pymol.pymol_utilities import pymol_utitlities as up
from restraintmaker.algorithm import Filter, Optimizer

#check pdbs
pairwise_files =os.getcwd()+"/data"  #/home/bschroed/Documents/code/restraintmaker/devtools/otherScripts/c_Protein_FreeEnergyCalculation/data"

pairs = glob.glob(pairwise_files+"/*")

pymol.finish_launching()
data=  {}
distance_treshold = 1.5

def analysis(res,file_name, approach, data, duration)->dict:
    ###measure convex Hull:
    atoms = []
    for r in res:
        a1 = r.atoms[0]
        a2 = r.atoms[1]
        cog_atom = np.array([a1.x + a2.x, a1.y + a2.y, a1.z + a2.z]) / 2
        atoms.append(cog_atom)

    d_greed = np.round(np.sum(
        [np.sum([np.sqrt(np.sum(np.square([a1[0] - a2[0], a1[1] - a2[1], a1[2] - a2[2]]))) for a2 in atoms[ind:]]) for
         ind, a1 in enumerate(atoms)]), 5)
    if(len(res)>3):
        v_greed = ConvexHull(atoms).volume
    else:
        v_greed = np.NAN

    data[file_name].update({ approach+"_volume": v_greed, approach+"_distance": d_greed, approach+"_t": duration})

    return data

def vis( res, path_prefix,nres=4, c=["firebrick", "forest", 'purple', 'salmon', "gold", "red", "marine", "deepolive"]):
    cmd.set("dash_radius", 0.1)
    for ind, r in enumerate(res):
        for a in r.atoms:
            cmd.show("spheres", "id " + str(a.id))
            cmd.set("sphere_color", c[ind%nres], "id " + str(a.id))
        cmd.distance("d" + str(ind), "id " + str(r.atoms[0].id), "id " + str(r.atoms[1].id))
        cmd.hide("labels")
        cmd.set("grid_slot", -2, "d" + str(ind))

    cmd.zoom()
    cmd.ray(1200)
    time.sleep(2)
    cmd.set("grid_mode", 0)
    cmd.png(path_prefix + "_restraints.png")
    cmd.set("grid_mode", 1)
    cmd.move("zoom", -10)

    cmd.ray(1200,)
    time.sleep(4)
    cmd.png(path_prefix + "_grid_restraints.png")
    cmd.set("grid_mode", 0)
    cmd.hide("spheres")
    cmd.delete("d*")

from restraintmaker.io.Exporter import Gromos_Distance_Restraint_Exporter


def write_out(res, outpath):
    exp = Gromos_Distance_Restraint_Exporter(res)
    exp.get_args(lambda x: outpath)
    exp.export_restraints()

for nrestraints in [4]: #[4, 6, 8]:
    for ind, pair_path in enumerate(pairs[1:]):
        cmd.reinitialize()
        cmd.bg_colour("white")
        cmd.set('ray_opaque_background', 'off')
        cmd.set("label_color", "black")
        cmd.set('stick_radius', '0.1')
        print(pair_path)
        if(not os.path.isfile(pair_path)):
            continue

        file_name = os.path.basename(pair_path).replace(".pdb", "")
        approach_path=bash.make_folder(os.getcwd()+"/disres_n" + str(nrestraints))



        data.update({file_name:{}})
        ###############################
        #BUILD DISRES
        ##load data
        cmd.set("retain_order", 1)
        cmd.load(pair_path)
        #cmd.show("lines", "all")
        cmd.show("sticks")
        cmd.set("sphere_scale", 0.2)
        time.sleep(2)

        resns = []
        cmd.iterate("all", "v.append(resn)", space={"v":resns})
        resns = list(set(resns))

        # set nice scene
        obj_list = cmd.get_object_list()
        [cmd.create("res"+str(i), "resn "+r) for i, r in enumerate(resns)]

        cmd.delete(obj_list[0])
        obj_list = cmd.get_object_list()
        cmd.set("pdb_retain_ids", 0)

        """
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
        """
        #cmd.label("all", "ID")
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

        print(filtered_atoms)
        cmd.show('spheres', "ID "+"+".join([str(a.id) for a in filtered_atoms]))
        exit()

        # Optimizers
        ## Greedy
        ### COG
        startt = datetime.datetime.now()
        opt = Optimizer.TreeHeuristicOptimizer(filtered_atoms)
        opt.get_args(lambda x: (nrestraints, distance_treshold, 'cog', 'convex_hull'))
        res = opt.make_restraints()
        endt = datetime.datetime.now()
        duration = (endt - startt).seconds

        ####Vis:
        vis(res=res, path_prefix=approach_path+"/"+file_name+"_greedy_cog", nres=nrestraints)

        ####Analysis:
        data = analysis(res=res, data=data, file_name=file_name, approach="greedy_cog", duration=duration)

        #### Write out:
        write_out(res, outpath=approach_path+"/"+file_name+"_greedy_cog.disres")

        ### shortest
        startt = datetime.datetime.now()
        opt = Optimizer.TreeHeuristicOptimizer(filtered_atoms)
        opt.get_args(lambda x: (nrestraints, distance_treshold, 'shortest', 'convex_hull'))
        res = opt.make_restraints()
        endt = datetime.datetime.now()
        duration = (endt - startt).seconds

        ####Vis:
        vis(res=res, path_prefix=approach_path+"/"+file_name+"_greedy_shortest", nres=nrestraints)

        ####Analysis:
        data = analysis(res=res, data=data, file_name=file_name, approach="greedy_shortest", duration=duration)

        continue

        ### biased avg.
        startt = datetime.datetime.now()
        opt = Optimizer.TreeHeuristicOptimizer(filtered_atoms)
        opt.get_args(lambda x: (nrestraints, distance_treshold, 'biased_avg', None))
        res = opt.make_restraints()
        endt = datetime.datetime.now()
        duration = (endt - startt).seconds

        ####Vis:
        vis(res=res, path_prefix=approach_path + "/"+file_name+"_greedy_biasedavg")

        ####Analysis:
        data = analysis(res=res, data=data, file_name=file_name, approach="greedy_biasedavg", duration=duration)

        ### prim
        startt = datetime.datetime.now()
        opt = Optimizer.TreeHeuristicOptimizer(filtered_atoms)
        opt.get_args(lambda x: (nrestraints, distance_treshold, 'prim', None))
        res = opt.make_restraints()
        endt = datetime.datetime.now()
        duration = (endt - startt).seconds

        ####Vis:
        vis(res=res, path_prefix=approach_path + "/"+file_name+"_greedy_prim", nres=nrestraints)

        ####Analysis:
        data = analysis(res=res, data=data, file_name=file_name, approach="greedy_prim", duration=duration)


        ## Brute Force convexHull
        if(nrestraints>3):
            startt = datetime.datetime.now()
            opt = Optimizer.BruteForceRingOptimzer(filtered_atoms)
            opt.get_args(lambda x: (nrestraints, distance_treshold, 'convex_hull', 'None'))
            res = opt.make_restraints()
            endt = datetime.datetime.now()
            duration = (endt - startt).seconds

            ####Vis:
            vis(res=res, path_prefix=approach_path + "/"+file_name+"_brute_force_ch", nres=nrestraints)

            ####Analysis:
            data = analysis(res=res, data=data, file_name=file_name, approach="brute_force_ch", duration=duration)

        ## Brute Force dist.
        startt = datetime.datetime.now()
        opt = Optimizer.BruteForceRingOptimzer(filtered_atoms)
        opt.get_args(lambda x: (nrestraints, distance_treshold, 'dist', 'None'))
        res = opt.make_restraints()
        endt = datetime.datetime.now()
        duration = (endt - startt).seconds

        ####Vis:
        vis(res=res, path_prefix=approach_path + "/"+file_name+"_brute_force_dist", nres=nrestraints)

        ####Analysis:
        data = analysis(res=res, data=data, file_name=file_name, approach="brute_force_dist", duration=duration)



        ## MY brute force
        ## Brute Force - dist:
        startt = datetime.datetime.now()

        atoms = []
        all_possible_restraints = []
        all_res = {}
        t = []
        ind = 0
        for i,a1 in enumerate(filtered_atoms):
            for a2 in filtered_atoms[i:]:
                d =  np.sqrt(np.sum(np.square(np.array([a1.x-a2.x, a1.y-a2.y, a1.z-a2.z]))))
                if(a1 != a2 and d < distance_treshold and a1.resn != a2.resn):
                    cog_atom = np.array([a1.x + a2.x, a1.y + a2.y, a1.z + a2.z]) / 2

                    all_possible_restraints.append(cog_atom)
                    atoms_ind = [a1.id, a2.id]
                    all_res.update({ind: {"atoms":(a1, a2), "atomsind":atoms_ind, "cog": cog_atom}})
                    ind +=1


        print("possible restraints ", len(all_res))

        molecule_ids = []
        cmd.iterate("res1", "molecule_ids.append(ID)", space={"molecule_ids": molecule_ids})
        all_res_ind = list(all_res.keys())
        ind_combinations = permutations(all_res_ind, nrestraints)
        cleaned_combinations = {}
        bruteF_distance_d = {}
        selected = []
        iterator = tqdm.tqdm(ind_combinations, leave=False)
        for res_comb in iterator:
            unique_atoms = set(sorted(np.unique(np.concatenate(list(map(lambda x: all_res[x]['atomsind'], res_comb))))))
            if (len(unique_atoms) == nrestraints * 2 and not unique_atoms in selected):
                selected.append(unique_atoms)
                cogs = np.array(list(map(lambda x: all_res[x]['cog'], res_comb)))
                d = np.round(np.sum([np.sum(
                    [np.sqrt(np.sum(np.square([a1[0] - a2[0], a1[1] - a2[1], a1[2] - a2[2]]))) for a2 in cogs[ind:]])
                                     for ind, a1 in enumerate(cogs)]), 5)
                if(nrestraints >3):
                    cleaned_combinations.update({res_comb: round(ConvexHull(cogs).volume, 4)})
                else:
                    cleaned_combinations.update({res_comb: np.nan})
                bruteF_distance_d.update({res_comb: d})

        if(nrestraints>3):
            max_key = max(cleaned_combinations, key=lambda x: cleaned_combinations[x])
            max_vol = cleaned_combinations[max_key]
        else:
            max_key = None
            max_vol = np.nan


        max_key = max(bruteF_distance_d, key=lambda x: bruteF_distance_d[x])
        max_dist = bruteF_distance_d[max_key]


        max_atoms = list(map(lambda x: all_res[x]['atomsind'], max_key))


        endt = datetime.datetime.now()
        duration = (endt - startt).seconds


        vol_dist = list(cleaned_combinations.values())
        dist_dist = list(bruteF_distance_d.values())

        print("all 4 res_combinations: ", len(cleaned_combinations))
        print("maxV: ", max_key, max_vol, max_atoms)
        print("BRUTEF Duration: ", duration)
        data[file_name].update({"total_myBruteF_volume": max_vol, "total_myBruteF_distance": max_dist, "total_BruteF_t": duration,
                               "myBruteF_all_V":vol_dist , "myBruteF_all_d":dist_dist,
                               })

        # Write out file:
        pd.DataFrame(data).T.to_csv(approach_path + "/tmp_algorithm_ana2.tsv", sep="\t")

print("fini")
exit()

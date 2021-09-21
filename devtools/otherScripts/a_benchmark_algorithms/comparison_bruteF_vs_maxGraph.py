"""
This is a script, that starts a restraint maker session.

"""
import datetime
import glob, os, time
import pickle
import numpy as np
import pymol
from pymol import cmd
from scipy.spatial import ConvexHull

#Distance restraints
from restraintmaker.algorithm import Optimizer
from restraintmaker.interface_Pymol.pymol_utilities import pymol_utitlities as up
from build_particleClouds import build_particle_clouds
from pygromos.utils import bash

#PARAMS:
distance_treshold = 100
visualization = False
iterations = 5
random_seeds = [7, 42, 64, 91, 128]
max_particle = 11

#PATHS
approach_path = "/home/bschroed/Documents/projects/restraintmaker/devtools/otherScripts/a_benchmark_algorithms"
data_dir = approach_path+"/data"
out_dir = approach_path + "/out"
out_stat_path = out_dir+"/out_state.obj"


#FUNCS
def analysis(res,pair_name, approach, data, duration)->dict:
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

    data[pair_name].update({ approach:{"volume": v_greed, "distance": d_greed, "t": duration}})

    return data

def vis( res, path_prefix, c=["firebrick", "forest", 'purple', 'salmon', "gold", "red"]):
    cmd.color("bluewhite", "resi 1")
    cmd.color("wheat", "resi 2")
    cmd.show("spheres")
    time.sleep(2)

    for ind, r in enumerate(res):
        for a in r.atoms:
            cmd.show("spheres", "id " + str(a.id)+" and resi "+str(a.resi))
            cmd.color(c[ind], "id " + str(a.id)+" and resi "+str(a.resi))
        cmd.distance("d" + str(ind), "id " + str(r.atoms[0].id)+" and resi "+str(r.atoms[0].resi), "id " + str(r.atoms[1].id)+" and resi "+str(r.atoms[1].resi))
        cmd.hide("labels")
    cmd.center()
    #cmd.move("zoom", -10)

    cmd.ray(1200)
    time.sleep(2)
    cmd.set("grid_mode", 0)
    cmd.png(path_prefix + "_restraints.png")
    cmd.delete("d*")

    time.sleep(2)


pymol.finish_launching()
result_data = {}
for it in range(iterations):
    data = {}
    out_it_dir = bash.make_folder(out_dir+'/iteration'+str(it))
    pairs = build_particle_clouds(out_dir=data_dir+"/iteration"+str(it), seedInt=random_seeds[it], max_particle=max_particle)  #Randomize point Clouds

    for nrestraints in [4]: #3,6
        for ind, pair_path in enumerate(pairs):
            cmd.reinitialize()
            time.sleep(2)
            cmd.bg_colour("white")
            cmd.set('ray_opaque_background', 'off')
            cmd.set("label_color", "black")
            cmd.set("sphere_scale", "0.2")

            pair_name = os.path.basename(pair_path).replace(".pdb", "")
            molA_name = "L01"
            molB_name = "L02"
            data.update({pair_name:{}})

            print(pair_name)

            ###############################
            #BUILD DISRES
            ##load data
            cmd.load(pair_path)
            time.sleep(1)

            cmd.center()
            cmd.zoom()
            cmd.show("spheres")

            # set nice scene
            obj_list = cmd.get_object_list()
            cmd.create("res1", "resn "+molA_name)
            cmd.create("res2", "resn "+molB_name)

            cmd.delete(obj_list[0])

            obj_list = cmd.get_object_list()
            cmd.set("pdb_retain_ids", 0)

            cmd.sync()


            ## GET ATOMS
            atom_list = filtered_atoms = up.pymol_selection_to_atom_list("all")
            cmd.color("bluewhite", "resi 1")
            cmd.color("wheat", "resi 2")

            ### shortest
            print("shortest")
            startt = datetime.datetime.now()

            opt = Optimizer.TreeHeuristicOptimizer(filtered_atoms)
            opt.get_args(lambda x: (nrestraints, distance_treshold, 'shortest', None))
            res = opt.make_restraints()

            endt = datetime.datetime.now()
            duration = (endt - startt)
            duration = duration.seconds+duration.microseconds/1000000

            ####Vis:
            if(visualization): vis(res=res, path_prefix=out_it_dir+"/"+pair_name+"_greedy_shortest")
            ####Analysis:
            data = analysis(res=res, data=data, pair_name=pair_name, approach="greedy_shortest", duration=duration)


            ## Brute Force convexHull
            print("Brute")
    
            if(nrestraints>3):
                startt = datetime.datetime.now()
                opt = Optimizer.BruteForceRingOptimzer(filtered_atoms)
                opt.get_args(lambda x: (nrestraints, distance_treshold, 'convex_hull', 'None'))
                res = opt.make_restraints()
                endt = datetime.datetime.now()
                duration = (endt - startt)
                duration = duration.seconds + duration.microseconds / 1000000
                ####Vis:
                if(visualization): vis(res=res, path_prefix=out_dir + "/"+pair_name+"_brute_force_ch")
    
                ####Analysis:
                data = analysis(res=res, data=data, pair_name=pair_name, approach="brute_force_ch", duration=duration)
    
            ## Brute Force dist.
            startt = datetime.datetime.now()
            opt = Optimizer.BruteForceRingOptimzer(filtered_atoms)
            opt.get_args(lambda x: (nrestraints, distance_treshold, 'dist', 'None'))
            res = opt.make_restraints()
            endt = datetime.datetime.now()
            duration = (endt - startt)
            duration = duration.seconds+duration.microseconds/1000000
    
            ####Vis:
            if(visualization): vis(res=res, path_prefix=out_dir + "/"+pair_name+"_brute_force_dist")
    
            ####Analysis:
            data = analysis(res=res, data=data, pair_name=pair_name, approach="brute_force_dist", duration=duration)    


            ## Random
            print("Random")
            import random
            from restraintmaker.utils import Types

            tmp_data = {pair_name:{}}
            sample_random = 100
            duration=0
            for x in range(sample_random):

                atoms_res1 = [a for a in filtered_atoms if(a.resi == '1')]
                atoms_res2 = [a for a in filtered_atoms if(a.resi == '2')]

                selected_atoms_res1 = random.sample(atoms_res1, k=4)
                selected_atoms_res2 = random.sample(atoms_res2, k=4)

                res = [Types._Restraint([a,b]) for a,b in  zip(selected_atoms_res1, selected_atoms_res2)]

                ####Analysis:
                tmp_data = analysis(res=res, data=tmp_data, pair_name=pair_name, approach="random"+"_"+str(x), duration=duration)

            ####Vis:
            if(visualization):  vis(res=res, path_prefix=out_it_dir+"/"+pair_name+"_random")

            data[pair_name].update({"random": {
                "t": 0,
                "volume":np.mean([d["volume"] for k, d in tmp_data[pair_name].items()]),
                "distance":np.mean([d["distance"] for k, d in tmp_data[pair_name].items()]),
                "volume_std": np.std([d["volume"] for k, d in tmp_data[pair_name].items()]),
                "distance_std": np.std([d["distance"] for k, d in tmp_data[pair_name].items()]),
                    }
               })

            # Write out file:
            row_keys = list(data.keys())
            column_keys = list(sorted(data[row_keys[0]].keys()))
            header1_keys = "approach\t"+"".join([key+"\t\t\t" for key in column_keys])+"\n"
            header2_keys = "points\t"+"".join(["time[s]\tvolume[A^3]\tdistance[A]\t" for key in column_keys])+"vol_std\tdistance_std\n"

            content = []
            for row in row_keys:
                print(row)
                line = str(row.split("_")[-1])+"\t"
                for approach in column_keys:
                    d = data[row][approach]
                    line += str(d["t"])+"\t"+str(d["volume"])+"\t"+str(d["distance"])+"\t"
                    if(approach == "random"):
                        line+=str(d["volume_std"])+"\t"+str(d["distance_std"])+"\t"
                line += "\n"
                content.append(line)

            #write:
            msg = header1_keys+header2_keys+"".join(content)
            print(msg)
            out_file = open(out_it_dir + "/tmp_algorithm_ana.tsv", 'w')
            out_file.write(msg)
            out_file.close()

    #Write Out
    result_data.update({"it"+str(it): data})
    f = open(out_stat_path, "wb")
    pickle.dump(result_data, f)
    f.close()


print("fini")
exit()



"""
ALTERNATIVE ALGORITHMS< THAT WERE NOT USED!
# Optimizers
## Greedy
### COG
print("COG")
startt = datetime.datetime.now()
opt = Optimizer.TreeHeuristicOptimizer(filtered_atoms)
opt.get_args(lambda x: (nrestraints, distance_treshold, 'cog', None))
res = opt.make_restraints()
endt = datetime.datetime.now()
duration = (endt - startt)
duration = duration.seconds+duration.microseconds/1000000

####Vis:
vis(res=res, path_prefix=out_dir+"/"+pair_name+"_greedy_cog")

####Analysis:
data = analysis(res=res, data=data, pair_name=pair_name, approach="greedy_cog", duration=duration)
"""

### biased avg.
"""
print("bias Avg")

startt = datetime.datetime.now()
opt = Optimizer.TreeHeuristicOptimizer(filtered_atoms)
opt.get_args(lambda x: (nrestraints, distance_treshold, 'biased_avg', None))
res = opt.make_restraints()
endt = datetime.datetime.now()
duration = (endt - startt)
duration = duration.seconds+duration.microseconds/1000000

####Vis:
#vis(res=res, path_prefix=out_dir + "/"+pair_name+"_greedy_biasedavg")

####Analysis:
data = analysis(res=res, data=data, pair_name=pair_name, approach="greedy_biasedavg", duration=duration)

### prim
print("Prim")

startt = datetime.datetime.now()
opt = Optimizer.TreeHeuristicOptimizer(filtered_atoms)
opt.get_args(lambda x: (nrestraints, distance_treshold, 'prim', None))
res = opt.make_restraints()
endt = datetime.datetime.now()
duration = (endt - startt)
duration = duration.seconds+duration.microseconds/1000000

####Vis:
#vis(res=res, path_prefix=out_dir + "/"+pair_name+"_greedy_prim")

####Analysis:
data = analysis(res=res, data=data, pair_name=pair_name, approach="greedy_prim", duration=duration)
"""
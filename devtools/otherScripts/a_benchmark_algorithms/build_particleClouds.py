import time
from random import seed
from random import random
from pymol import cmd
from pygromos.utils import bash
import numpy as np


def build_particle_clouds(out_dir:str, max_particle:int=15, min_particles:int= 6, particle_steps:int=1, seedInt=42)->str:
    seed(seedInt)

    particles = np.arange(min_particles, max_particle, particle_steps)
    out_dir = bash.make_folder(out_dir)

    coordsA = []
    coordsB = []

    out_pdbs = []
    step = 0

    time.sleep(2)
    cmd.reinitialize()
    for i in particles:
        densF = i//2 if(i<12) else 12//2
        print("particles: ", i)

        #cmd.reinitialize()
        for j in range(step, i):
            cmd.pseudoatom(resn="A")
            f=0
            while True:
                dx = -5 + random() * densF
                dy = -5 + random() * densF
                dz = -5 + random() * densF

                if(len(coordsA) == 0 or all([np.sqrt(np.sum(np.square([c[0]-dx, c[1]-dy, c[2]-dz]))) > 1 for c in coordsA])):
                    coordsA.append([dx, dy, dz])
                    break
                f+=1


            cmd.translate([dx, dy, dz], 'resn A')
            cmd.alter('resn A', 'resn="L01"')
            cmd.alter('resn L01', 'resi=1')

        for j in range(step, i):
            cmd.pseudoatom(resn="B")

            while True:
                dx = -5 + random() * densF
                dy = -5 + random() * densF
                dz = -5 + random() * densF
                if (len(coordsB) == 0 or all(
                        [np.sqrt(np.sum(np.square([c[0] - dx, c[1] - dy, c[2] - dz]))) > 1 for c in coordsB])):
                    coordsB.append([dx, dy, dz])
                    break

            cmd.translate([dx, dy, dz], 'resn B')
            cmd.alter('resn B', 'resn="L02"')
            cmd.alter('resn L02', 'resi=2')

        out_pdb_path = out_dir+"/punktwolke_" + str(i) + ".pdb"
        cmd.save(out_pdb_path)
        out_pdbs.append(out_pdb_path)
        step = i

    return out_pdbs
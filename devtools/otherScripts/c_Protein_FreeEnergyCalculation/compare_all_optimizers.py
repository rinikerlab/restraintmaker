'''
module compare_all_optimziers
Based on demo. Loop over all sets of Ligands to compare the Optimizers using the Optimizer module;s compare_pair_optimizers function
'''
import time
from restraintmaker.algorithm import Optimizer

if __name__ == "__main__":
    import sys
    import os
    import src.Interface_Pymol.RestraintMaker_PyMOL as diswiz


def start_pymol():
    #checkout if PyMol Package is present
    try:
        try:
            #start Interface_Pymol
            import pymol
            pymol.finish_launching()
            from pymol import cmd
        except Exception as err:
            print("WARNING: could not find PyMOL in enviroment. Try installing via anaconda3")
            if("conda" in sys.modules):
                import conda.cli as cli
                cli.main('conda', 'install',  '-y', '-c schrodinger' , 'Interface_Pymol')
            else:
                raise Exception("Could not find conda enviroment.")
    except Exception as err:
        print("\nERROR: Failed to start PyMol!\n\tIncluding PyMOL did not work and installing did also not work, as anaconda enviroment was not found!\n\tPlease Use an Anaconda3 enviroment to run the plugin.")
        exit(1)

    return cmd

def load_ligands(ligand_file):
    print(ligand_file)
    cmd.load(os.path.dirname(__file__) +'/'+ ligand_file)
    cmd.center()
    cmd.zoom()
    pymol_module = diswiz.Restraints_Wizard()
    cmd.set_wizard(pymol_module)
    return pymol_module

ligands_dir = 'data/ligand_system/'
ligand_files =['CHK1_5Ligands.pdb','BRD4_7Ligs.pdb','PNMT_9ligs.pdb','21_large_Ligands.pdb']#['PNMT_9ligs.pdb','BRD4_7Ligs.pdb',


ligands_path = [ligands_dir+lig_file for lig_file in ligand_files]
cmd=start_pymol()

all_results_dir = 'tmp_results'
if(not os.path.exists(all_results_dir)):
    os.mkdir(all_results_dir)


for lig in ligands_path:
    cmd.delete('all')
    pymol_module=load_ligands(lig)


    print('COMPARING '+lig)
    print('len(all_atoms)',len(pymol_module.logic_handler.all_atoms))
    print('len(selected_atoms)',len(pymol_module.logic_handler.selected_atoms))
    print('len(selected_restraints)',len(pymol_module.logic_handler.selected_restraints))

    start_time = time.time()
    #Note; Changed criterion from convex_hull back to pca unscaled on Fr 24.05
    Optimizer.compare_pair_optimizers(criterion=Optimizer._calculate_value_unscaled_pca_2d, \
                                  atoms=pymol_module.logic_handler.all_atoms, \
                                  opt_types=[Optimizer.TreeHeuristicOptimizer, Optimizer.TreeHeuristicOptimizer,
                                              Optimizer.TreeHeuristicOptimizer, Optimizer.TreeHeuristicOptimizer,Optimizer.BestMoleculeRingOptimizer,Optimizer.BruteForceRingOptimzer,
                                              Optimizer.BruteForceRingOptimzer], \
                                   opt_args=[(4, 1.2, 'prim', None), (4, 1.2, 'cog', None),
                                             (4, 1.2, 'shortest', None), (4, 1.2, 'biased_avg', None),(4, 1.2, 'pca','pca_2d'),(4, 1.2, 'pca', 'None'),
                                             (4, 1.2, 'convex_hull', 'None')], \
                                  out_path=all_results_dir,\
                                  new_dir_name=lig[len(ligands_dir):-4]) #Name of the Ligands without path and .pdb

    end_time=time.time()
    passed_time = end_time-start_time
    print('COMPARED '+lig+' IN '+str(passed_time)+'s')
    print('')
    print('______________________________________________________________________________________________________________')
    print('______________________________________________________________________________________________________________')
    print('______________________________________________________________________________________________________________')
    print('')
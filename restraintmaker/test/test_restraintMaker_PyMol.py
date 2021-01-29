import time, typing as t
from pymol import cmd

from restraintmaker.interface_Pymol import pymol_utilities as pu
from restraintmaker.utils.Utilities import print
from restraintmaker.algorithm import Optimizer

# ----------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------    Testing       -----------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------
# TODO CREATE AN OWN CLASS OPTIMIZER TESTER
# Would be nice if we had just one object my_tester_that we can deal with

def _check_optimzer_results_pairwise_6(self, show_next: bool = True, exit_test_mode=False):
    '''Creates groups for all pairs which are connected. Then allow to cycle thorugh these this groups as in test mode 4
    '''

    # CASE 0: EXIT TEST MODE
    if exit_test_mode:
        for pair in self.check_result_dict['pairs']:
            self.cmd.delete('pair_' + pair[0] + "_" + pair[1])

        self.cmd.set('grid_mode', 0)
        self.cmd.enable('all')
        self.check_result_dict.clear()
        self.check_results_mode = 0
        cmd.refresh_wizard()

        return

    # CASE 1: TEST MODE STARTED FRESHLY

    if not 'pairs' in self.check_result_dict.keys():
        self.check_results_mode = 6
        cmd.disable('all')

        # 1) Check which molecule pairs have connections

        pair_exists = {}
        for i_m1 in range(len(self.pymol_molecule_objects)):
            for i_m2 in range(len(self.pymol_molecule_objects)):
                pair_exists.update({str(i_m1 + 1) + str(i_m2 + 1): False})

        for r in self.logic_handler.selected_restraints:
            pair_exists.update({str(r.atoms[0].resi) + str(r.atoms[1].resi): True})

        mol_pairs = [(str(i_m1 + 1), str(i_m2 + 1)) for i_m1 in range(len(self.pymol_molecule_objects) - 1) for i_m2
                     in
                     range(i_m1, len(self.pymol_molecule_objects)) if
                     pair_exists[str(i_m1 + 1) + str(i_m2 + 1)] or pair_exists[str(i_m2 + 1) + str(i_m1 + 1)]]

        # Create   pairwise restraints. Distance restraints are added to the specifed group upon creation
        # Distance objects cant be copied, so I have to vreate two sepearte ones here
        for r in self.logic_handler.selected_restraints:
            __mol1 = r.atoms[0].resi
            __mol2 = r.atoms[1].resi

            if int(__mol1) > int(__mol2):
                __mol1, __mol2 = __mol2, __mol1

            restr_name = 'rest_' + str(__mol1) + '_' + str(__mol2)

            cmd.distance(name=restr_name, \
                         selection1=pu._atom_list_to_pymol_selection([r.atoms[0]]), \
                         selection2=pu._atom_list_to_pymol_selection([r.atoms[1]]))
            cmd.disable(restr_name)  # Will be enabled on demand

            cmd.distance(name='copy_' + restr_name, \
                         selection1=pu._atom_list_to_pymol_selection([r.atoms[0]]), \
                         selection2=pu._atom_list_to_pymol_selection([r.atoms[1]]))
            cmd.enable('copy_' + restr_name)  # Must be enabled, so it is shown when the grpout is enabled

        group_names = []
        group_expressions = []
        for p in mol_pairs:
            m1, m2 = p[0], p[1]
            pair = (m1, m2)

            # Create copy for group shown separatelt
            name_copy_m1 = 'mol_' + m1 + 'copy' + m1 + "_" + m2
            name_copy_m2 = 'mol_' + m2 + 'copy' + m1 + "_" + m2
            cmd.copy(name_copy_m1, 'mol_' + m1)
            cmd.copy(name_copy_m2, 'mol_' + m2)

            # Create copy of restraints

            group_name = 'pair_' + m1 + '_' + m2
            cmd.group(group_name, name_copy_m1 + ' ' + name_copy_m2 + ' ' + 'copy_rest_' + m1 + '_' + m2)

            # Important: An object is only shown when it AND its group are enabled
            # => Enable all objects, but disable the groups.
            # cmd.enable(group_expression)
            cmd.disable(group_name)

            group_names.append(group_name)

        cmd.set('grid_mode', 1)
        cmd.refresh_wizard()
        cmd.reset()

        self.check_result_dict.update({'pairs': mol_pairs})
        self.check_result_dict.update({'i_pair': -1})  # Will be adjusted in next block

    # CASE 2: Show next pair also called in Case 1
    # A) Disable old displays
    i_pair = self.check_result_dict['i_pair']
    m1, m2 = self.check_result_dict['pairs'][i_pair]
    cmd.disable('mol_' + m1)
    cmd.disable('mol_' + m2)
    cmd.disable('pair_' + m1 + "_" + m2)
    cmd.disable('rest_' + m1 + "_" + m2)

    # B) Enable next pair
    i_pair += 1 if show_next else -1
    i_pair %= len(self.check_result_dict['pairs'])
    self.check_result_dict.update({'i_pair': i_pair})
    m1, m2 = self.check_result_dict['pairs'][i_pair]

    cmd.enable('mol_' + m1)
    cmd.set('grid_slot', '1', 'mol_' + m1)
    cmd.enable('mol_' + m2)
    cmd.set('grid_slot', '2', 'mol_' + m2)
    cmd.enable('pair_' + m1 + "_" + m2)
    cmd.set('grid_slot', '3', 'pair_' + m1 + "_" + m2)
    cmd.enable('rest_' + m1 + "_" + m2)
    cmd.set('grid_slot', '4', 'rest_' + m1 + "_" + m2)

    self.redraw()
    # self.cmd.reset()


def _check_optimzer_results_pairwise_5(self, exit_test_mode=False):
    '''Creates groups for all pairs which are connected. And show them all in grid mode
        :warning: only designed for pairwise restraints
    '''

    self.check_results_mode = 5
    # 1) Check which molecule pairs have connections

    pair_exists = {}
    for i_m1 in range(len(self.pymol_molecule_objects)):
        for i_m2 in range(len(self.pymol_molecule_objects)):
            pair_exists.update({str(i_m1 + 1) + str(i_m2 + 1): False})

    for r in self.logic_handler.selected_restraints:
        pair_exists.update({str(r.atoms[0].resi) + str(r.atoms[1].resi): True})

    mol_pairs = [(str(i_m1 + 1), str(i_m2 + 1)) for i_m1 in range(len(self.pymol_molecule_objects) - 1) for i_m2 in
                 range(i_m1, len(self.pymol_molecule_objects)) if
                 pair_exists[str(i_m1 + 1) + str(i_m2 + 1)] or pair_exists[str(i_m2 + 1) + str(i_m1 + 1)]]
    print(mol_pairs, mv=1)
    cmd.disable('all')
    for p in mol_pairs:
        m1, m2 = p[0], p[1]
        m1_name = m1 + '-' + m2 + 'mol_' + m1
        m2_name = m1 + '-' + m2 + 'mol_' + m2
        cmd.copy(m1 + '-' + m2 + 'mol_' + m1, 'mol_' + m1)
        cmd.copy(m1 + '-' + m2 + 'mol_' + m2, 'mol_' + m2)
        group_expression = m1_name + ' ' + m2_name
        cmd.group('pair_' + m1 + '_' + m2, group_expression)
        cmd.enable(group_expression)

    cmd.set('grid_mode', 1)
    cmd.reset()

    if (exit_test_mode):
        self.check_results_mode = 0

        for p in mol_pairs:
            m1_name = m1 + '-' + m2 + 'mol_' + m1
            m2_name = m1 + '-' + m2 + 'mol_' + m2
            group_expression = m1_name + ' ' + m2_name
            cmd.delete(m1_name)
            cmd.delete(m2_name)
            cmd.ungroup('group_expression')

        cmd.set('grid_mode', 0)
        cmd.enable('all')


def _check_optimzer_results_pairwise_4(self, mol1: str = None, mol2: str = None, exit_test_mode=False):
    # TODO: Make the group of copies instead of the objects themselves. Then show the two moleucles and the restraints between them in grid mode

    # Look at 2 molecules at once, use interface_Pymol group command to select relevant molecules and restraints
    # Stores in cro: mol1, mol2. The mmoleucles we are currently looking at

    # -1) DELETE OLD GROUPS, if they exits
    if 'mol1' in self.check_result_dict.keys() and 'mol2' in self.check_result_dict.keys():
        self.cmd.ungroup('rest_' + self.check_result_dict['mol1'])
        self.cmd.ungroup('rest_' + self.check_result_dict['mol2'])
        self.cmd.delete(self.check_result_dict['mol1'] + '_copy')
        self.cmd.delete(self.check_result_dict['mol2'] + '_copy')
        self.cmd.delete(self.check_result_dict['group_name'])

    # 0) CHECK IF WE SHOUDL EXIST TEST MODE
    if exit_test_mode:
        self.check_results_mode = 0
        self.check_result_dict = {}
        self.cmd.enable('all')
        return

    # 1) CHECK WHAT MOLECULES SHOULD BE DISPLAYED
    # 1 - Check that the user either provided two inputs, or decide what moeules to display from dcurretn input

    # 1a) Case 1: The user speciefied both molecules => Check that they are valid molecules. Leave them as they are
    self.check_results_mode = 4
    if mol1 != None and mol2 != None:
        if not self.check_objects_exists(mol1):
            print('There is no molecule called ' + mol1, mv=4)
            return
        elif not self.check_objects_exists(mol2):
            print('There is no molecule called ' + mol1, mv=4)
            return

    # 1b) Case 2: The user speciefied one molecule -> Check that that Molecule is valid, and let the other one unchanged
    elif (mol1 == None) ^ (mol2 == None):  # ^ i s xor
        print(str(mol1 == None), str(mol2 == None), str((mol1 == None) ^ (mol2 == None)), mv=4)
        if mol1 != None:
            if not self.check_objects_exists(mol1):
                print('There is no molecule called ' + mol1, mv=4)
                return
            else:
                mol2 = self.check_result_dict['mol2']
        elif mol2 != None:
            if not self.check_objects_exists(mol2):
                print('There is no molecule called ' + mol2, mv=4)
                return
            else:
                mol1 = self.check_result_dict['mol1']


    # 1c) Case 3: The user specified no moleule => Move thorugh the list by one molecule
    else:
        mol1 = self.check_result_dict['mol2']
        all_mols = self.pymol_molecule_objects
        mol2 = all_mols[(all_mols.index(self.check_result_dict['mol2']) + 1) % len(all_mols)]

    # 2) Update the currently displayed moleucles
    if mol1 == None or mol2 == None:
        print('Please specify which molecules you want to look at using test mode 4', mv=4)
        return

    self.cmd.copy(mol1 + '_copy', mol1)
    self.cmd.copy(mol2 + '_copy', mol2)

    group_name = mol1 + '+' + mol2
    self.check_result_dict = {'mol1': mol1, 'mol2': mol2, 'group_name': group_name}

    group_expression = mol1 + '_copy or ' + mol2 + '_copy or rest_' + mol1 + ' and rest_' + mol2

    self.cmd.set('grid_mode', 1)
    self.cmd.disable('all')
    self.cmd.group(name=group_name, members=group_expression)
    #        self.cmd.show('sticks',group_name)
    self.cmd.enable(group_name)
    self.cmd.enable(group_expression)
    self.cmd.enable(mol1)
    self.cmd.enable(mol2)
    self.cmd.reset()

    self.cmd.refresh_wizard()

    # TODO: Draw a fourth grop containing only the restraints, not overlayed on the molecule. Not so easz becasue we can not copy restraints
    # => Need to create new objects


def _test_optimizer_1(self, sele, Optimizer):
    """
   Only used during programming to test some things about the optimizer.
    :return:
    :rtype:
    """
    # all_atoms = pymol_selection_to_atom_list("all")
    cons = Optimizer.get_all_short_connections(atoms=pu.pymol_selection_to_atom_list('all'), max_dis=1.2)

    # for a in all_atoms:
    if True:
        a = pu.pymol_selection_to_atom_list(sele)[0]
        neighbours = []
        # TODO: If I ever use this in efficient context I can speed it up with the knowledge, that the atoms are ordered: the lower of the two atoms is always the first one in the list
        for c in cons:
            if c.Atom1 == a:
                neighbours.append(c.Atom2)
            elif c.Atom2 == a:
                neighbours.append(c.Atom1)

        self.cmd.color("green", pu._atom_list_to_pymol_selection(neighbours))
        self.cmd.color("red", pu._atom_list_to_pymol_selection([a]))

        self.cmd.delete('mypseudo')
        self.cmd.pseudoatom('mypseudo', pos='[' + str(a.x) + ',' + str(a.y) + ',' + str(a.z) + ']')
        self.cmd.show('spheres', 'mypseudo')
        self.cmd.set('sphere_scale', 1.2)
        self.cmd.set('sphere_transparency', 0.5, selection='SphericalSelection')

        # time.sleep(0.1)
    #  pass #Set a breakpoint here to watch result


def _test_optimizer_2(self):
    # TODO: Can not create selection in constructor,because molecules ar not yet loeaded. => Find reasonable order

    # Need to define some global variables, so the show next molecule fkt can see them. Mark them with _to2_ for test_optimizer_2
    self._to2_i_m = 0
    self._to2_n_restraints = 4

    a = pu.pymol_selection_to_atom_list('all') if len(
        self.logic_handler.selected_atoms) == 0 else self.logic_handler.selected_atoms
    my_Optimizer = Optimizer.Volume_Optimzier(a)
    my_Optimizer.get_args(lambda x: self._to2_n_restraints, lambda x: 'avg')
    self._to2_restraints = my_Optimizer.make_restraints()

    self.cmd.set('grid_mode', 1)
    self._optimzer_show_next_molecule_pair_2()

    # Display the restraints one after the other


def _test_optimizer_3(self, timer_runs: int = 0):
    '''_test_optimizer_3 can be used to compare the results of the different optimizers. if timer_runs is specified the optimizers will be timed by running their make_restraint functions timer_runs time.
    The results will be printed to the log.
    '''
    self._to3_i_m = 0
    self._to3_n_restraints = 3

    a = pu.pymol_selection_to_atom_list('all') if len(
        self.logic_handler.selected_atoms) == 0 else self.logic_handler.selected_atoms

    prim_Optimizer = pu.Optimizer.MoleculeRingOptimizer_1_0(a)
    shortest_Optimizer = pu.Optimizer.MoleculeRingOptimizer_1_0(a)
    avg_Optimizer = pu.Optimizer.MoleculeRingOptimizer_1_0(a)

    prim_Optimizer.get_args(lambda x: self._to3_n_restraints, lambda x: 'prim')
    shortest_Optimizer.get_args(lambda x: self._to3_n_restraints, lambda x: 'shortest')
    avg_Optimizer.get_args(lambda x: self._to3_n_restraints, lambda x: 'avg')

    if timer_runs > 0:
        print('PRIM')
        self.timer(func=prim_Optimizer.make_restraints, runs=timer_runs)
        print('SHORTEST')
        self.timer(func=shortest_Optimizer.make_restraints, runs=timer_runs)
        print('AVG')
        self.timer(func=shortest_Optimizer.make_restraints, runs=timer_runs)

    self._prim_restraints = prim_Optimizer.make_restraints()
    self._shortest_restraints = shortest_Optimizer.make_restraints()
    self._avg_restraints = avg_Optimizer.make_restraints()

    print([r.atoms for r in self._prim_restraints])
    print([r.atoms for r in self._shortest_restraints])
    print([r.atoms for r in self._avg_restraints])

    for m in self.cmd.get_object_list():
        self.cmd.create(m + '_prim', m)
        self.cmd.create(m + '_shortest', m)
        self.cmd.create(m + '_avg', m)
        self.cmd.delete(m)

    self.cmd.set('grid_mode', 1)
    self._optimzer_show_next_molecule_pair_3()


def _optimzer_show_next_molecule_pair_2(self):
    '''
    Only for debugging and testing.
    :return:
    :rtype:
    '''

    molecules = self.cmd.get_object_list()

    color = ['red', 'green', 'magenta', 'cyan']
    self.redraw()  # TODO: There also seems to be a self.cmd.redraw attribute. => Change name of my function to avoid confusion
    self.cmd.hide('spheres', 'all')
    self.cmd.disable('all')
    self.cmd.enable(molecules[self._to2_i_m])
    self.cmd.enable(molecules[(self._to2_i_m + 1) % len(molecules)])
    self.cmd.reset()
    self.cmd.refresh()
    for r in self._to2_restraints[
             self._to2_n_restraints * self._to2_i_m: self._to2_n_restraints * self._to2_i_m + self._to2_n_restraints]:
        self.cmd.color(color[self._to2_restraints.index(r) % self._to2_n_restraints],
                       pu._atom_list_to_pymol_selection(r.atoms))
        self.cmd.show('spheres', pu._atom_list_to_pymol_selection(r.atoms))
    self.cmd.reset()
    self.cmd.refresh()


def _optimzer_show_next_molecule_pair_3(self):
    '''
    Only for debugging and testing.
    :return:
    :rtype:
    '''

    molecules = (self.cmd.get_object_list())

    color = ['orange', 'green', 'magenta', 'brown']
    self.redraw()  # TODO: There also seems to be a self.cmd.redraw attribute. => Change name of my function to avoid confusion
    self.cmd.hide('spheres', 'all')
    self.cmd.disable('all')
    m1 = self._to3_i_m + 1  # convert ot counting from 1 to n, instead of 0 to n-1
    m2 = m1 + 1 if int(m1) != len(molecules) / 3 else '1'

    self.cmd.enable(str(m1) + '_prim')
    self.cmd.enable(str(m1) + '_shortest')
    self.cmd.enable(str(m1) + '_avg')
    self.cmd.enable(str(m2) + '_prim')
    self.cmd.enable(str(m2) + '_shortest')
    self.cmd.enable(str(m2) + '_avg')

    my_slice = slice(self._to3_n_restraints * self._to3_i_m, self._to3_n_restraints * (self._to3_i_m + 1))

    for r_prim, r_shortest, r_avg in zip(self._prim_restraints[my_slice], self._shortest_restraints[my_slice],
                                         self._avg_restraints[my_slice]):
        col = color[self._prim_restraints.index(r_prim) % self._to3_n_restraints]
        self.cmd.color(col, pu._atom_list_to_pymol_selection(r_prim.atoms) + ' & ' + str(
            m1) + '_prim')  # After using the create command the corresponding atoms have the same id => need to specify molecule
        self.cmd.color(col, pu._atom_list_to_pymol_selection(r_prim.atoms) + ' & ' + str(m2) + '_prim')

        self.cmd.color(col, pu._atom_list_to_pymol_selection(r_shortest.atoms) + ' & ' + str(m1) + '_shortest')
        self.cmd.color(col, pu._atom_list_to_pymol_selection(r_shortest.atoms) + ' & ' + str(m2) + '_shortest')

        self.cmd.color(col, pu._atom_list_to_pymol_selection(r_avg.atoms) + ' & ' + str(m1) + '_avg')
        self.cmd.color(col, pu._atom_list_to_pymol_selection(r_avg.atoms) + ' & ' + str(m2) + '_avg')

        self.cmd.show('spheres', pu._atom_list_to_pymol_selection(r_prim.atoms + r_shortest.atoms + r_avg.atoms))

    print(pu._atom_list_to_pymol_selection(r_avg.atoms) + ' & m1_avg')
    self.cmd.reset()
    self.cmd.refresh()
    # for r in self._to2_restraints[
    #          self._to2_n_restraints * self._to2_i_m: self._to2_n_restraints * self._to2_i_m + self._to2_n_restraints]:
    #     self.cmd.color(color[self._to2_restraints.index(r) % self._to2_n_restraints],
    #                    atom_list_to_pymol_selection(r.atoms))
    #     self.cmd.show('spheres', atom_list_to_pymol_selection(r.atoms))


def timer(self, runs: int, func: t.Callable):
    '''Timer: Can be used'''
    t_start = time.time()
    for i in range(runs):
        func()
    t_end = time.time()

    # Correct for time of loop execution etc
    t_corr_start = time.time()
    for i in range(runs):
        pass
    t_corr_end = time.time()

    # TODO: Round to sign digits. No direct built in function in python
    upper = (t_end - t_start) / runs
    lower = ((t_end - t_start) - (t_corr_end - t_corr_start)) / runs
    print('Timer estimates:', min_verbosity=0)
    print('  Lower: ' + '%10.3e' % lower + ' s', min_verbosity=0)
    print('  Upper: ' + '%10.3e' % upper + ' s', min_verbosity=0)

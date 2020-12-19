"""
    This module plugs in the restraintmaker logic into PyMOL.
    It basically contains the GUI-realisation and the linkeage to the logic

"""
import time
import typing as t

import pymol
from pymol import cmd
from pymol.wizard import Wizard

import restraintmaker.interface_Pymol.Pymol_Utitlities as pu
from restraintmaker.algorithm import Optimizer
from restraintmaker.algorithm.Selection import SphericalSelection  # Necessary for special case checks
from restraintmaker.interface_Pymol import RestraintMaker_Logic as Logic
from restraintmaker.utils import Utilities as u
from restraintmaker.utils.Utilities import print


# from interface_Pymol.exporting import _resn_to_aa as one_letter
#            d['onelettercode'] = one_letter.get(d['resn'], d['resn'])

class Restraints_Wizard(Wizard):
    """
    This class helps generating dis res.
    Sources: label(Wizard)
    """

    def __init__(self, _self=cmd):
        self.logic_handler = Logic.Logic_Handler(pu.pymol_selection_to_atom_list('all'))

        Wizard.__init__(self, _self)
        self.set_wizard_style()
        # reset gui selections:
        self.cmd.unpick()
        self.cmd.deselect()
        self.cmd.set('mouse_selection_mode', 0)

        # TODO: Try to directly pass the type to set_*_type() instead of a string. Should work with f.__name__. Problem; When function is called the Filter Module is out of scope
        ##generate Menu drop down

        # Pymols python interpreter will not be able to understand a reference to self. But we can find the Way back here via cmd.get_wizard.
        # We can not use the __class__  attribute: The function will be passed as a string to the interface_Pymol Python interpreter. So its ars have to be strings themselves
        self.menu = {

            'Import': [[2, "Importer", '']] +
                      [[1, i.__name__, 'cmd.get_wizard()._action_button_pressed(\"' + i.__name__ + '\")'] for i in
                       self.logic_handler.available_importer_types],
            'Restraint': [[2, 'Restraint', '']] +
                         [[1, r.__name__, 'cmd.get_wizard()._action_button_pressed(\"' + r.__name__ + '\")'] for r in
                          self.logic_handler.available_restraint_types],
            'Selection': [[2, "Selection", '']] +
                         [[1, s.__name__, 'cmd.get_wizard()._action_button_pressed(\"' + str(s.__name__) + '\")'] for s
                          in self.logic_handler.available_selection_types],

            'Filter': [[2, "Filter", '']] +
                      [[1, f.__name__, 'cmd.get_wizard()._action_button_pressed(\"' + str(f.__name__) + '\")'] for f in
                       self.logic_handler.available_filter_types],
            'Optimizer': [[2, 'Optimizer', '']] +
                         [[1, o.__name__, 'cmd.get_wizard()._action_button_pressed(\"' + str(o.__name__) + '\")'] for o
                          in self.logic_handler.available_optimizer_types],
            # 'Save': [[1,'Save','']],
            'Exporter': [[2, 'Exporter', '']] +
                        [[1, e.__name__, 'cmd.get_wizard()._action_button_pressed(\"' + e.__name__ + '\")'] for e in
                         self.logic_handler.available_exporter_types],
            'Done': [[1, 'Done', 'cmd.set_wizard()']]

        }
        self.show_message = True
        self.cmd.set('retain_order')  # NEVER CHANGE IDS OF ATOMS
        self.pymol_molecule_objects = pu.create_pymol_objects_for_molecules()  # Must be done vefore logic handler, because it changes atom ids!
        # hashes to keep work of interface_Pymol: self.redraw only redraw things that actually changed\
        self.last_hashes = {'atoms': 0, 'selection': 0, 'restraints': 0, 'test_mode': 0}
        # crd (check_results_dict) is an dict that can beused by the test mode to store information
        self.check_results_mode = 0
        self.crd = {}


        if(len(self.cmd.get_object_list())>0):
            self.cmd.delete(self.cmd.get_object_list()[0])  # Delete the object containing all molecules. If not all Molecules are doubled

        self.redraw()
        self.cmd.refresh_wizard()
        self.cmd.viewing.reset()
        self.tmp_cnt = 0
        self.supress_prompt = False
        self.quicksnapshot_nr = 1

        # Check results mode indicates if we are using a special test mode at the moment.
        # The different check_opzimier_results functions indicate in whihc mode they are called

    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    # ----------------------------------------------------------Important TODOs----------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    # TODO: Find a better way to recoulor automatically after each event. Eithre allow backward communcation between logic_unit and GUI_unit or, make an event function of the GUI,that will 1) call the appropriate logic_handler.react_t-_8_event and then update the colour

    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    # -----------------------------------------------------------Pymol Event handling-----------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    # The trigger-event method is used to translate the interface_Pymol events into events the selections can deal with, by calling the corresponding Logic_handler.react_to_*_event function
    # TODO STUC: Decide wether to redraw automatically after each event, or wether to trigger the redrawing from the Logic module.

    def get_event_mask(self):
        return Wizard.event_mask_pick + Wizard.event_mask_select + Wizard.event_mask_special + Wizard.event_mask_key \
               + Wizard.event_mask_dirty  # Event_maks_special = Arrow keys

    # TODO: Might be cleaner if I pass key word args as dict (**) instead of normal args as tuple (*) => unpack dict for the _Selection.update_x function and check arguments there (TypeError)
    #   Problem with that: When the triggere event function is called, the key words must be given. ACTUALLY. THAT IS NOT A BAD THING. WHEN PROGRAMMING IT YOu have to know anyway what args to use. By using **kwargs I can enforce that and check that it is correctly done
    #   TO avaid confusion I can just catch the Error and explain what is going on i the ERROR message

    def trigger_event(self, event_type: u.EventType, **kw_args: t.Any):
        '''
        trigger_event will call the logic_handler's appropriate react_to_*_event function, depending on event_type
        and cause the wizard and the view to apply the changes caused by the event.

        :warning: TypeError will be thrown if trigger_event is called without explicitly stating the keywords of the target update function
        :param event_type: what event has been caused? Acceptable values: 'select', 'move', 'size', 'confirm'
        :type event_type: u.EventType
        :param kw_args: Arguments the relevant update functon will need. KEYWORDS OF THE CORRESPONDING _Selection._update_* function MUST BE GIVEN. (** => dict)
        :type args: t.Any
        :raises: ValueError if event_type can not be resolved to a logic.handler.react_to_*_event function
        :raises: u.BadArgumentException if args to not match the the corresponding logic_handler.react_to_*_event_function
        :return: -
        :rtype: None
        '''

        # Note: The ** Operator does 2 different things here depending on context:
        #       In the function definiton **kw_args means: Accept arbitrarily many keyword arguments and pack them into a dict called kw_args
        #       In function calls it means: Unpack kw_args and pass the content as single keyword arguments to the next function

        self.logic_handler.react_to_event(event_type, **kw_args)
        self.redraw()
        self.cmd.refresh_wizard()

    def do_special(self, k, x, y, mod):
        print('special')

        # TODO CLEAN: Do not acess selection directly
        if self.logic_handler.my_selection == None or self.logic_handler.my_selection.has_finished:
            return
        if k == 102:  # right
            self.trigger_event(u.EventType.SIZE, increase=True)  # increase = True
            if hasattr(self.logic_handler.my_selection, 'radius'):  # Can not test directlz if it is Spherical Selection
                self.cmd.set('sphere_scale', self.logic_handler.my_selection.radius,
                             selection='SphericalSelection')  # TODO CLEAN: GUI should not acess Selection directly.. Integrate that into display function or something
        elif k == 100:
            self.trigger_event(u.EventType.SIZE, increase=False)  # increase = False
            if hasattr(self.logic_handler.my_selection, 'radius'):  # Can not test directlz if it is Spherical Selection
                self.cmd.set('sphere_scale', self.logic_handler.my_selection.radius, selection='SphericalSelection')

        # TODO CLEAN: Find a more intuitive way for chanign the mode
        elif k == 101:  # +
            self.cmd.controlling.button('Left', 'None', 'MovA')
        elif k == 103:  # down
            self.cmd.controlling.button('Left', 'None', 'Rota')

        # TODO CLEAN: SHOULD NOT ACESS SELECTION BY ITSELF. => Inroduce a mousclicked event / finish input event or something, Maybe with enter? DoubleClick?
        if k == 101 and mod == 1:  # Shift + up
            self.trigger_event(u.EventType.CONFIRM)
            self.cmd.controlling.button('Left', 'None', 'Rota')
            self.cmd.controlling.unmask('all')
            # self.cmd.hide('spheres','SphericalSelection')
            self.cmd.delete('SphericalSelection')
            self.logic_handler._set_selection_type(None)

    def do_dirty(self):
        print('dirty', mv=0)

        self.redraw()  # only for debugging
        # TODO: Find a better way to detect movement of pseudoatom: Do_position? / Save last positon and only call do move if it actually moved
        # TODO STRUC: Move this out to a provider function. The GUI handling program should not deal with special cases for certain selction types.
        if isinstance(self.logic_handler.my_selection, SphericalSelection):
            if not self.check_objects_exists('SphericalSelection'):
                self.provide_services_to_spherical_selection()

            # Find the coordinates of the Spherical selection
            coords = [None]
            self.cmd.iterate_state(-1, selection='SphericalSelection', expression='coords[0]=(x,y,z)', space=locals())
            x, y, z = coords[0]
            self.trigger_event(u.EventType.MOVE, new_x=x, new_y=y, new_z=z)

    def do_key(self, k, x, y, mod):

        # Does not deal with events passed down to Selection
        # TODO: Introduce a flag for testing, or simply remove the event_mask_key and only set i, once one of the test functions is called
        # Check what the current test mode is and do the appropriate action at each key. Problem: FBefore test mode 4 i did not use crd and self.test_mode but just checked if attributes set by the test functions are present

        print(k, x, y, mod, mv=0)

        # TODO: Create function check_reults, which will then call the fitting check_results_x function instead of testing them all here
        if k == 13 and mod == 0:  # Enter
            if self.check_results_mode == 4:
                self._check_optimzer_results_pairwise_4()
            elif self.check_results_mode == 6:
                self._check_optimzer_results_pairwise_6()  # Show_next+True by default
            elif hasattr(self, '_to2_i_m'):
                self._to2_i_m = (self._to2_i_m + 1) % len(self.cmd.get_object_list())
                self._optimzer_show_next_molecule_pair_2()
            elif hasattr(self, '_to3_i_m'):
                self._to3_i_m = int((self._to3_i_m + 1) % (len(self.cmd.get_object_list()) / 3))
                self._optimzer_show_next_molecule_pair_3()

        if k == 13 and mod == 1:  # shift Enter
            # Start or end test mode 5
            if self.check_results_mode == 4:
                self._check_optimzer_results_pairwise_4(exit_test_mode=True)
            elif self.check_results_mode == 5:
                self._check_optimzer_results_pairwise_5(exit_test_mode=True)
            elif self.check_results_mode == 6:
                self._check_optimzer_results_pairwise_6(exit_test_mode=True)
            elif len(self.logic_handler.selected_restraints) != 0:
                self._check_optimzer_results_pairwise_6()


        elif k == 8:  # Backspace
            if self.check_results_mode == 6:
                self._check_optimzer_results_pairwise_6(show_next=False)
            if hasattr(self, '_to2_i_m'):
                self._to2_i_m = (self._to2_i_m - 1) % len(self.cmd.get_object_list())
                self._optimzer_show_next_molecule_pair_2()
            elif hasattr(self, 'to_3_i_m'):
                self._to3_i_m = int((self._to3_i_m - 1) % (len(self.cmd.get_object_list()) / 3))
                self._optimzer_show_next_molecule_pair_3()
        elif k == 32:  # Space
            self.cmd.set('grid_mode', (int(self.cmd.get('grid_mode')) + 1) % 2)
        elif k == 81 or k == 113:  # q or Q
            self.cmd.enable('all')

        elif k in range(48, 58):  # Number
            # TODO: ALLOW SHIFT + NUMBER TO BE USED TO SWITCH INTO THE CORRESPONDING TEST VIEW MODE

            # Show the corresponidng moleucle in test mode
            # TODO; Find a way that works for more than 10 ligands
            num = k - 48
            print(str(num), mv=1)
            if self.check_results_mode == 4:
                if num > 0 and num <= len(self.pymol_molecule_objects):
                    self._check_optimzer_results_pairwise_4(mol1=self.pymol_molecule_objects[(num - 1)],
                                                            mol2=self.crd['mol1'])

        # KEYBOARD SHORTCUTS FOR USEFULL FUNCTIONS. Warning: Do not use alt or Ctrl for shortcuts. Many keys are already occupied
        # TODO: Some of those could be executed from Logic module
        # h: Hide instructions
        if k == 104:
            self.show_message = False  # The only instruction will be how to show instructions
            self.cmd.refresh_wizard()
        # s: show instructions
        if k == 115:
            self.show_message = True
            self.cmd.refresh_wizard()

        # Shf R: Delete all restraints
        if k == 82 and mod == 1:
            self.logic_handler.selected_restraints = []
            self.logic_handler.set_action_states()
            self.get_panel()
            self.redraw()
            self.cmd.refresh_wizard()

        # Shf Q: Choose all atoms, start brute force optimizer, start test mode 4
        if k == 81 and mod == 1:
            self._action_button_pressed('UniversalSelection')
            self.trigger_event(u.EventType.CONFIRM)
            self._action_button_pressed('TreeHeuristicOptimizer')

            # Make sure Optimizatyion was sucessfull
            if (len(self.logic_handler.selected_restraints) > 0):
                self._check_optimzer_results_pairwise_6()

        # Shf T: Test Optimizers
        if k == 84 and mod == 1:
            print('Testing Optimizers...', mv=3)
            # TODO: Move away from interface_Pymol modules.  I do not like to import everything here
            import Optimizer
            from PyQt5.QtWidgets import QFileDialog, QInputDialog

            Optimizer.compare_pair_optimizers(criterion=Optimizer._calculate_value_convex_hull, \
                                              atoms=self.logic_handler.selected_atoms, \
                                              opt_types=[Optimizer.TreeHeuristicOptimizer,
                                                         Optimizer.TreeHeuristicOptimizer,
                                                         Optimizer.TreeHeuristicOptimizer,
                                                         Optimizer.BruteForceRingOptimzer,
                                                         Optimizer.BruteForceRingOptimzer], \
                                              opt_args=[(4, 1.2, 'prim', None), (4, 1.2, 'cog', None),
                                                        (4, 1.2, 'shortest', None), (4, 1.2, 'pca', None),
                                                        (4, 1.2, 'convex_hull', None)], \
                                              out_path=QFileDialog.getExistingDirectory(
                                                  caption='Where should all restraints be saved?'),
                                              new_dir_name=
                                              QInputDialog.getText(cmd.gui.get_qtwindow(), 'Directory name',
                                                                   'What should the new directory be calledd?')[0])

            print('... tested optimizers', mv=3)

        # Shft C: Readjust color
        if k == 67 and mod == 1:
            pu.colorize(self.logic_handler.all_atoms)

        # Shft S: Snapshot
        if k == 83 and mod == 1:
            # TODO: refresh does not redraw instatnly, but 'as soon as the operatoing system allows it'. So the prompt might not get hidden for the picture
            # self.supress_prompt = True
            self.cmd.refresh()
            # file=u.create_file_save_dialog()('Where should the file be saved?')
            file = u.create_file_save_dialog()(
                "Where should the snapshot be saved?")  # '~/Desktop/'+str(self.quicksnapshot_nr)+'.png'
            self.quicksnapshot_nr += 1
            self.cmd.png(file)
            # self.supress_prompt=False
            self.cmd.refresh()
            self.cmd.refresh_wizard()
            print('Saved picture at ' + file, mv=3)

    file_nr = 1

    def do_pick(self, bondFlag):
        self.do_select('pk1')
        self.cmd.unpick()

    def do_select(self, sele):
        """
        ..function: do_generate_restraints
            todo: reimplement
            todo: move to Selections and out_source the interface_Pymol commands here into extra functions
        :return:
        """

        newly_selected_atoms: t.List[u.Atom] = pu.pymol_selection_to_atom_list(
            sele)  # sele is deleted at the end of the function. So it now only contains atoms selected by mouse since  the last call of do select

        # TODO: Intoduce swtich table to specify action dpending on the chosen selection =>

        self.trigger_event(u.EventType.SELECT, new_atoms=newly_selected_atoms)
        self.cmd.delete(
            'sele')  # Necessary, to assure that in the next call of do_select we do not get all atoms selected up to now as well

        # TODO: Change program flow:
        # Mode1: IF we have an active restraint type, the selection will be used to create a new restraint and add this to a list of restraints
        # Mode2: If there is no active restraint type, we select atoms and do not create restraints and just wildly select and filter atoms
        # At the moment we only have mode 1

    #  ---------------------Buttton functions-------------------------

    # TODO STRUC: MIght be cleaner to do it in one function:
    def _action_button_pressed(self, x_name: str):
        # TODO STRUC: ALSO ROUTE THIS THROUGH THE EVENT HANDLING< SO THE logic_MODUEL updates after any presed buttion!!!

        '''
        button_pressed should be called when any of the action buttons is pressed.

            It will call the releveant logic_handler.set_*_type method, redraw the molecule and enable/disabnle buttons if necessary.
            I think it is easier to pack all buttons into one function, and go through a switch table, than to define seperate functions for importer, Selection etc,
            plus a function called in all cases

        :param x_name: The __name__ attribute of a type Importer, Restraint, Selection, Filter, Optimizer or Exporter
        :type x_name: str
        :raises: TypeError if x_type is not a type of Importer, Restraint, Selection, Filter, Optimizer or Exporter
        :return -
        :rtype: None
        :warning: The x_type attribute MUST be of type string, not type, because of technical issues:\
         To be able to pass this function to the interface_Pymol pyhton interpreter all its arguments need to be of type string\
         So the interface_Pymol module will work with __name__ instead of __class__ attributes

        '''

        # Because we get __name__ and not __class__ attributes for x_type (->see warning in docstring), we have to tediously convert it back to  a type
        # Using types instead of names is much more sensible, so I want to do the conversion here and not in the logic_module, as this is a interface_Pymol specific problem.

        # SWITCH-TABLE to convert from name to type:
        x_type = None

        # Importer
        if x_type == None:
            for i_type in self.logic_handler.available_importer_types:
                if i_type.__name__ == x_name:
                    x_type = i_type

        # Restraint
        if x_type == None:
            for r_type in self.logic_handler.available_restraint_types:
                if r_type.__name__ == x_name:
                    x_type = r_type

        # Selection
        if x_type == None:
            for s_type in self.logic_handler.available_selection_types:
                if s_type.__name__ == x_name:
                    x_type = s_type

        # Filter
        if x_type == None:
            for f_type in self.logic_handler.available_filter_types:
                if f_type.__name__ == x_name:
                    x_type = f_type

        # Optimizer
        if x_type == None:
            for o_type in self.logic_handler.available_optimizer_types:
                if o_type.__name__ == x_name:
                    x_type = o_type

        # Exporter
        if x_type == None:
            for e_type in self.logic_handler.available_exporter_types:
                if e_type.__name__ == x_name:
                    x_type = e_type
        # Not Found
        if x_type == None:
            raise ValueError(
                'The x_name argument of button_pressed must be the __name__ of an Importer, Restraint, Selection, Filter, Optimizer or Exporter. Not: ' + str(
                    x_name))

        self.logic_handler.set_action_type(x_type)
        self.redraw()  # Actually only necessary after Optimizer
        self.cmd.refresh_wizard()

    def _toggle_button_pressed(self, button):
        '''
        _mode_button pressed is called when one of the 2 buttons defining the user mode is pressed.

        :param button: Which button was pressed
        :type button: str
        :return: -
        :rtype: -
        '''
        if button == 'toggle_delete':
            self.logic_handler.set_mode(select_delete=not self.logic_handler.select_or_delete_mode,
                                        atom_restraint=self.logic_handler.atom_or_restraint_mode)
        elif button == 'toggle_atom':
            self.logic_handler.set_mode(select_delete=self.logic_handler.select_or_delete_mode,
                                        atom_restraint=not self.logic_handler.atom_or_restraint_mode)
        else:
            raise ValueError(
                '_toggle_button_pressed needs to be called with \'button\' = toggle_delte or toggle_atom. Not ' + button)
        self.get_panel()
        self.cmd.refresh_wizard()

    # -------Utility functions--------------------------

    def get_prompt(self):
        '''called automatically upon refresh wizard'''
        # TODO: It must be possible to use a smaller font
        if self.supress_prompt:
            print('SUPRESSING PROMPT')
            return []

        prompt = []
        if not self.show_message:
            prompt.append('show instructions: s')
        else:
            prompt.append('hide instructions: h')
            prompt.append('')

            # Current mode
            mode, instruction = '', ''

            if self.logic_handler.atom_or_restraint_mode:
                if self.logic_handler.select_or_delete_mode:
                    mode = 'SELECT ATOMS MODE'
                    instruction = 'select atoms using selections'
                else:
                    mode = 'DELETE ATOMS MODE'
                    instruction = 'delete atoms using selections'
            else:
                if self.logic_handler.select_or_delete_mode:
                    mode = 'SELECT RESTRAINTS MODE'
                    instruction = 'create (distance) restraints by clicking atoms'
                    # TODO: Generalize to all restraints.
                else:
                    mode = 'DELETE RESTRAINTS MODE'
                    instruction = 'delete restraints by clicking atoms'

            prompt += [mode, instruction]

            # How to use a selection
            # TODO: Would be nice to only show the events that the current selection supports.
            if self.logic_handler.current_selection_type != None:
                prompt += ['', 'SELECTION', \
                           'select: lft clck', \
                           'cnfrm: shft+up']
                if 'SphericalSelection' in self.cmd.get_object_list():
                    prompt += [
                        'size: left / right' \
                        'move/rotate: mouse, left', \
                        'activate rota: up',
                        'activate move: down']
                    # TODO: Set a bool indicating the current mode, so we can check here

            # How to get into test mode
            # TODO: Implement old test modes
            if len(self.logic_handler.selected_restraints) > 0 and self.check_results_mode == 0:
                prompt += ['', 'TEST MODES', \
                           'Show restraints pairwise: Shift + Enter']

            # If in test mode
            if self.check_results_mode == 4 or self.check_results_mode == 6:
                prompt += ['', 'TEST MODE: Pairwise', \
                           'leave check mode: Shft+Enter', \
                           'show next pair: Enter']
                if self.check_results_mode == 4:
                    prompt += ['show molecule #: number-keys']
                if self.check_results_mode == 6:
                    prompt += ['show previous pair #: backspace']

            # Shortcuts
            prompt += ['', 'SHORTCUTS', \
                       'Shft r: del. all rest.', \
                       'Shft q: instant tree', \
                       'Shft s: snapshot', \
                       'spcace: grid view']

        return prompt

    def get_panel(self):
        # TODO: Do not name the buttons after the Objects they call, but after what the user can do (select, instead of selection. select all instead of universal selection)
        empty_panel = [[0, '', '']]  # Empty Space

        # shorthand to get the name of an  class
        def get_name(x_type) -> str:
            return '-' if x_type == None else x_type.__name__

        empty_panel = [[0, '', '']]

        text_select_delete = 'select (x) | delete    ( )' if self.logic_handler.select_or_delete_mode \
            else 'select ( ) | delete    (x)'

        text_atom_restraint_mode = 'atoms  (x) | restraints ( )' if self.logic_handler.atom_or_restraint_mode \
            else 'atoms  ( ) | restraints (x)'

        menu = []
        menu += empty_panel
        menu += [[str(self.logic_handler._action_states['toggle_select_delete']), text_select_delete,
                  'cmd.get_wizard()._toggle_button_pressed(\"toggle_delete\")']]
        menu += [[str(self.logic_handler._action_states['toggle_select_delete']), text_atom_restraint_mode,
                  'cmd.get_wizard()._toggle_button_pressed(\"toggle_atom\")']]

        menu += empty_panel
        menu += [[str(self.logic_handler._action_states['Importer']),
                  'Import:' + get_name(self.logic_handler.current_importer_type), 'Import']]
        menu += [[str(self.logic_handler._action_states['Restraint']),
                  'Restraint:' + get_name(self.logic_handler.current_restraint_type), 'Restraint']]
        menu += [[str(self.logic_handler._action_states['Selection']),
                  'Selection: ' + get_name(self.logic_handler.current_selection_type), 'Selection']]
        menu += [[str(self.logic_handler._action_states['Filter']),
                  'Filter: ' + get_name(self.logic_handler.current_filter_type), 'Filter']]
        menu += [[str(self.logic_handler._action_states['Optimizer']),
                  'Optimizer: ' + get_name(self.logic_handler.current_optimizer_type), 'Optimizer']]
        menu += [[str(self.logic_handler._action_states['Exporter']),
                  'Exporter: ' + get_name(self.logic_handler.current_exporter_type), 'Exporter']]
        menu += [[str(self.logic_handler._action_states['Done']), 'Done', 'cmd.set_wizard()']]

        return menu

    def set_wizard_style(self):
        cmd.color("grey", "elem C")
        cmd.do("set sphere_scale, 0.2")

    def start(self):
        """
        Opens a new Pymol Window
        :return: -
        :rtype: -
        """
        # checkout if PyMol Package is present
        # TODO: Make it work
        # TODO: Remove hardcoded loading of a certain Protein. Give a FileChooserDialog

        pymol.finish_launching()

        cmd.load("/home/rhinerc/project/code/restraintmaker/test/test_files/ligand_system/9_similar_Ligands.pdb")

        cmd.show("sticks", "resn HIS")
        cmd.center()
        cmd.zoom()
        cmd.set_wizard(self)

    # TODO STRUCTURE SPHSELE 0: The current situation for SphericalSelection is a dirty, tedious, convoluted solution.
    # Find a proper way to generalize:\
    # aybe: every Selection define as 'static' class attributes in which events in needs to be notified. Then rm_pymol sets flags, that tells each do function, if it should call a provider function
    # if I have to define a provider function for every selection and every do_function it is very tedious. But maybe I can just use update. YES: I CAN GIVE UPDATE AN EXTRSA ARGUMENT INDICATING WHICH do_function called it.

    def provide_services_to_spherical_selection(self):
        # 1st call: Draw Pseudoatom, set mouse function. Register with dirty listener, mask atoms

        # TODO STRUC: We also have to properly stop  the selection if it is discarded
        self.cmd.pseudoatom(
            "SphericalSelection")  # TODO STRUCTURE SELE: We could store the name of the Pseudoatom IN the selection
        self.cmd.mask('not SphericalSelection')

        self.cmd.show('sphere', 'SphericalSelection')
        self.cmd.set('sphere_scale', 5, 'SphericalSelection')
        self.cmd.set('sphere_transparency', 0.5, selection='SphericalSelection')

        # TODO STRUC: Find mor intuituve solution, that still allows to move selection, and that is not overridden automatically
        self.cmd.controlling.button('Left', 'None', 'MovA')

    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    # -----------------------------------------------------------------Utility functions that need to be part of the Distance_Restraints class-----------------------------------------------------------------
    # -----------------------------------------------------------------------------------------------------------------------------------------------------

    def create_pymol_objects_for_restraints(self):
        # TODO: Move to utilities. We can acess the wizrd via cmd.get_wizard instead of self

        # DELTE ALL PYMOL RESTRAINT OBJECTS (They are not on the object list => Need to reconstruct their names here)
        for m in self.cmd.get_object_list('all'):
            self.cmd.delete('rest_' + m)

        # CREATE NEW RESTRAINT OBJECTS FOR EACH MOLECULE SEPARATELY. (Not every Molecule pair becasue it has to be generalizeable to other Restraints than PairRestraint)
        # => the same restraint will appear in several lists

        def _connect_atoms(atom_list: t.List[u.Atom], restraint_name: str):
            for i_a1, a1 in enumerate(atom_list[:-1]):
                for a2 in atom_list[i_a1 + 1:]:
                    self.cmd.distance(name=restraint_name, \
                                      selection1=pu._atom_list_to_pymol_selection([a1]), \
                                      selection2=pu._atom_list_to_pymol_selection([a2]))

        for mol in self.pymol_molecule_objects:
            for res in self.logic_handler.selected_restraints:
                for atm in res.atoms:
                    # TODO: FIND A WAY TO GET THE NAME OF THE MOLCULE INTO THE ATOM
                    if 'mol_' + atm.resi == mol:
                        _connect_atoms(atom_list=res.atoms, restraint_name='rest_' + mol)
                        break

    def check_objects_exists(self, obj: str) -> bool:
        """
        check_if_onject_exists checks if a interface_Pymol object of a certain name existst.
        :warning: Pymol does not count Selections as objects.

        :param sele: Name of a interface_Pymol selection
        :type sele: str
        :return: Selection contains at least one object: True. Selection does not exist or is empty: False
        :rtype: bool
        """
        return obj in self.cmd.get_object_list()

    def redraw(self):
        """
        redraw paints all selected atoms, and the atoms in the active selection in different colours. Sets all other atoms back to standard colours

        Different colours for 1) Atoms in the curently active selections (not yet confirmed)
                              2) Selected Atoms
                              3) Atoms in restraints
                              4) Atoms in restraints we are looking at at the moment (in tst mode))
        :return: -
        :rtype: -
        """

        # Has anything that influences the coloring of atoms changed?

        atoms_hash = hash(str(self.logic_handler.selected_atoms))
        selection_hash = hash(
            str(self.logic_handler.my_selection.atoms)) if self.logic_handler.my_selection != None else hash(
            str([None]))
        restraints_hash = hash(str(self.logic_handler.selected_restraints))
        test_mode_hash = hash(str(self.check_results_mode) + str(self.crd))

        if atoms_hash != self.last_hashes['atoms'] or selection_hash != self.last_hashes['selection'] \
                or restraints_hash != self.last_hashes['restraints'] or test_mode_hash != self.last_hashes['test_mode']:

            # 0) Recolour all atoms to standard colour
            try:
                self.cmd.util.cba(
                    'vanadium')  # The colour for carbon atoms needs to be given explicitly. WARNING; The python cba command does not work properly. It does not work at all with gray. It raises an exception with vanadium, but it does paint it corectly.
            except:
                pass

            # 1) Colour selected atoms
            pu.help_pymol_with_big_atom_list(self.cmd.color, self.logic_handler.selected_atoms, color="orange")

            # 2) Colour atoms in current selection
            if self.logic_handler.my_selection != None:
                pu.help_pymol_with_big_atom_list(self.cmd.color, self.logic_handler.my_selection.atoms, color='yellow')

            # 3) Colour atoms in restraints

            # Have the restraints changed?
            if restraints_hash != self.last_hashes['restraints']:
                self.last_hashes.update(restraints=restraints_hash)
                self.create_pymol_objects_for_restraints()

            restrained_atoms = []
            for r in self.logic_handler.selected_restraints:
                for a in r.atoms: restrained_atoms.append(a)

            # 3a) Colour all restrained atoms
            pu.help_pymol_with_big_atom_list(self.cmd.color, restrained_atoms, color='marine')

            # 3b) In test mode: Colour atoms in restraints we are looking at at the moment
            if self.check_results_mode in [4, 6] and test_mode_hash != self.last_hashes['test_mode']:
                restrained_atoms = []

                if self.check_results_mode == 4:
                    for r in self.logic_handler.selected_restraints:
                        if (r.atoms[0].resi == self.crd['mol1'][4:] and r.atoms[1].resi == self.crd['mol2'][4:]) \
                                or (
                                r.atoms[0].resi == self.crd['mol2'][4:] and r.atoms[1].resi == self.crd['mol1'][4:]):
                            for a in r.atoms: restrained_atoms.append(a)
                if self.check_results_mode == 6:
                    for r in self.logic_handler.selected_restraints:
                        m1, m2 = self.crd['pairs'][self.crd['i_pair']]
                        if (r.atoms[0].resi == m1 and r.atoms[1].resi == m2) \
                                or (r.atoms[0].resi == m2 and r.atoms[1].resi == m1):
                            for a in r.atoms: restrained_atoms.append(a)

                pu.help_pymol_with_big_atom_list(self.cmd.color, restrained_atoms, color='cyan')

            self.last_hashes.update(atoms=atoms_hash, selection=selection_hash, restraints=restraints_hash,
                                    test_mode=test_mode_hash)
        # -----------END OF if has has changed clause-------------------------------------

        # SPECIAL CASE FOR SPHERICAL SELECTION; CHECK IF SPHERICAL SELECTION IS STILL ACTIVE
        if (not isinstance(self.logic_handler.my_selection, SphericalSelection)) and self.check_objects_exists(
                'SphericalSelection'):
            self.cmd.unmask('all')
            self.cmd.delete('SphericalSelection')
            self.cmd.button('Left', 'None', "Rota")

    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    # -----------------------------------------------------------------    Testing       -----------------------------------------------------------------
    # -----------------------------------------------------------------------------------------------------------------------------------------------------
    # TODO CREATE AN OWN CLASS OPTIMIZER TESTER

    # Would be nice if we had just one object my_tester_that we can deal with

    # TODO: put this into the old format making it a fucntion test_optimizer_4 and show_next_molecule_pair_4

    def _check_optimzer_results_pairwise_6(self, show_next: bool = True, exit_test_mode=False):
        '''Creates groups for all pairs which are connected. Then allow to cycle thorugh these this groups as in test mode 4
        '''

        # CASE 0: EXIT TEST MODE
        if exit_test_mode:
            for pair in self.crd['pairs']:
                self.cmd.delete('pair_' + pair[0] + "_" + pair[1])

            self.cmd.set('grid_mode', 0)
            self.cmd.enable('all')
            self.crd.clear()
            self.check_results_mode = 0
            cmd.refresh_wizard()

            return

        # CASE 1: TEST MODE STARTED FRESHLY

        if not 'pairs' in self.crd.keys():
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

            self.crd.update({'pairs': mol_pairs})
            self.crd.update({'i_pair': -1})  # Will be adjusted in next block

        # CASE 2: Show next pair also called in Case 1
        # A) Disable old displays
        i_pair = self.crd['i_pair']
        m1, m2 = self.crd['pairs'][i_pair]
        cmd.disable('mol_' + m1)
        cmd.disable('mol_' + m2)
        cmd.disable('pair_' + m1 + "_" + m2)
        cmd.disable('rest_' + m1 + "_" + m2)

        # B) Enable next pair
        i_pair += 1 if show_next else -1
        i_pair %= len(self.crd['pairs'])
        self.crd.update({'i_pair': i_pair})
        m1, m2 = self.crd['pairs'][i_pair]

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
        if 'mol1' in self.crd.keys() and 'mol2' in self.crd.keys():
            self.cmd.ungroup('rest_' + self.crd['mol1'])
            self.cmd.ungroup('rest_' + self.crd['mol2'])
            self.cmd.delete(self.crd['mol1'] + '_copy')
            self.cmd.delete(self.crd['mol2'] + '_copy')
            self.cmd.delete(self.crd['group_name'])

        # 0) CHECK IF WE SHOUDL EXIST TEST MODE
        if exit_test_mode:
            self.check_results_mode = 0
            self.crd = {}
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
                    mol2 = self.crd['mol2']
            elif mol2 != None:
                if not self.check_objects_exists(mol2):
                    print('There is no molecule called ' + mol2, mv=4)
                    return
                else:
                    mol1 = self.crd['mol1']


        # 1c) Case 3: The user specified no moleule => Move thorugh the list by one molecule
        else:
            mol1 = self.crd['mol2']
            all_mols = self.pymol_molecule_objects
            mol2 = all_mols[(all_mols.index(self.crd['mol2']) + 1) % len(all_mols)]

        # 2) Update the currently displayed moleucles
        if mol1 == None or mol2 == None:
            print('Please specify which molecules you want to look at using test mode 4', mv=4)
            return

        self.cmd.copy(mol1 + '_copy', mol1)
        self.cmd.copy(mol2 + '_copy', mol2)

        group_name = mol1 + '+' + mol2
        self.crd = {'mol1': mol1, 'mol2': mol2, 'group_name': group_name}

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

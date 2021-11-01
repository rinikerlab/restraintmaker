"""
    This module plugs in the restraintmaker logic into PyMOL.
    It basically contains the GUI-realisation and the linkeage to the logic.

"""
import typing as t

from pymol import cmd
from pymol.wizard import Wizard


from restraintmaker.algorithm.Selection import SphericalSelection  # Necessary for special case checks
from restraintmaker.interface_Pymol import RestraintMaker_Logic as Logic
from restraintmaker.interface_Pymol.pymol_utilities import program_states, qt_dialogs
import restraintmaker.interface_Pymol.pymol_utilities.pymol_utitlities as pu

from restraintmaker.utils import Utilities as u
from restraintmaker.utils.Utilities import print


class Restraints_Wizard(Wizard):
    """
    This class helps generating dis res.
    Sources: label(Wizard)
    """
    file_nr = 1

    def __init__(self, _self=cmd):
        """
            Initialize PyMol to the wizard!

        Parameters
        ----------
        _self
        """
        self.logic_handler = Logic.Logic_Handler(pu.pymol_selection_to_atom_list('all'))

        Wizard.__init__(self, _self)
        self.set_wizard_style()

        # reset gui selections:
        self.cmd.unpick()
        self.cmd.deselect()
        self.cmd.set('mouse_selection_mode', 0)

        ##generate Menu drop down
        ### Pymols python interpreter will not be able to understand a reference to self. But we can find the Way back here via cmd.get_wizard.
        ### We can not use the __class__  attribute: The function will be passed as a string to the interface_Pymol Python interpreter. So its ars have to be strings themselves
        self.update_menu()

        ##Settings
        self.show_message = True
        self.supress_prompt = False
        self.quicksnapshot_nr = 1

        self.last_hashes = {'atoms': 0, 'selection': 0, 'restraints': 0,
                            'test_mode': 0}  # hashes to keep work of interface_Pymol: self.redraw only redraw things that actually changed\

        self.check_result_dict = {}  # crd (check_results_dict) is an dict that can beused by the test mode to store information
        self.check_results_mode = 0  # Check results mode indicates if we are using a special test mode at the moment.
        self.tmp_cnt = 0  # The different check_opzimier_results functions indicate in whihc mode they are called

        ## Set up objs
        self.cmd.set('retain_order')  # NEVER CHANGE IDS OF ATOMS
        self.pymol_molecule_objects = pu.create_pymol_objects_for_molecules()  # Must be done before logic handler, because it changes atom ids!
        if (len(self.cmd.get_object_list()) > 0):
            self.cmd.delete(self.cmd.get_object_list()[
                                0])  # Delete the object containing all molecules. If not all Molecules are doubled

        ## Update
        self.redraw()
        self.cmd.refresh_wizard()
        self.cmd.viewing.reset()

    """
        Pymol event handling
    """

    def update_menu(self):
        self.menu = {
            'Import': [[2, "Importer", '']] +
                      [[1, i.__name__.replace("import_", ""), 'cmd.get_wizard()._action_button_pressed(\"' + i.__name__ + '\")'] for i in
                       self.logic_handler.available_importer_types],
            'Restraint': [[2, 'Restraint', '']] +
                         [[1, r.__name__.replace("_Type", ""), 'cmd.get_wizard()._action_button_pressed(\"' + r.__name__.replace("_Type", "") + '\")'] for r in
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

            'Exporter': [[2, 'Exporter', '']] +
                        [[1, e.__name__.replace("export_", ""), 'cmd.get_wizard()._action_button_pressed(\"' + e.__name__ + '\")'] for e in
                         self.logic_handler.available_exporter_types],
            'Done': [[1, 'Done', 'cmd.set_wizard()']]
        }

    # The trigger-event method is used to translate the interface_Pymol events into events the selections can deal with, by calling the corresponding Logic_handler.react_to_*_event function
    # TODO STUC: Decide wether to redraw automatically after each event, or wether to trigger the redrawing from the Logic module.

    def get_event_mask(self):
        return Wizard.event_mask_pick + Wizard.event_mask_select + Wizard.event_mask_special + Wizard.event_mask_key \
               + Wizard.event_mask_dirty  # Event_maks_special = Arrow keys

    def trigger_event(self, event_type: program_states.EventType, **kw_args: t.Any):
        """
            trigger_event will call the logic_handler's appropriate react_to_*_event function, depending on event_type
            and cause the wizard and the view to apply the changes caused by the event.

            TypeError will be thrown if trigger_event is called without explicitly stating the keywords of the target update function

        Parameters
        ----------
        event_type :program_states.EventType
            what event has been caused? Acceptable values: 'select', 'move', 'size', 'confirm'
        kw_args : t.Any
            Arguments the relevant update functon will need. KEYWORDS OF THE CORRESPONDING _Selection._update_* function MUST BE GIVEN. (** => dict)

        Returns
        -------
        NoReturn

        Raises
        ------
        ValueError
            if event_type can not be resolved to a logic.handler.react_to_*_event function
        u.BadArgumentException
            if args to not match the the corresponding logic_handler.react_to_*_event_function
        """

        self.logic_handler.react_to_event(event_type, **kw_args)
        self.redraw()
        self.cmd.refresh_wizard()

    def do_special(self, k, x, y, mod):
        """
            control special actions of the keyboard

        Parameters
        ----------
        k : int
            pressed key
        x : int
            coordinate x
        y : float
            coordinate y
        mod : int
            pressed shift?


        Returns
        -------
        None

        """
        print('special')

        if self.logic_handler.my_selection == None or self.logic_handler.my_selection.has_finished:
            return

        if k == 102:  # right - arrow key
            self.trigger_event(program_states.EventType.SIZE, increase=True)  # increase = True
            if hasattr(self.logic_handler.my_selection, 'radius'):  # if it is Spherical Selection increase sphere size
                self.cmd.set('sphere_scale', self.logic_handler.my_selection.radius, selection='SphericalSelection')
        elif k == 100:  # left - arrow key
            self.trigger_event(program_states.EventType.SIZE, increase=False)  # increase = False
            if hasattr(self.logic_handler.my_selection, 'radius'):  # if it is Spherical Selection reduce sphere size
                self.cmd.set('sphere_scale', self.logic_handler.my_selection.radius, selection='SphericalSelection')
        elif k == 101:  # + -key
            self.cmd.controlling.button('Left', 'None', 'MovA')
        elif k == 103:  # down - key
            self.cmd.controlling.button('Left', 'None', 'Rota')

        if k == 101 and mod == 1:  # up + Shift - keys
            self.trigger_event(program_states.EventType.CONFIRM)
            self.cmd.controlling.button('Left', 'None', 'Rota')
            self.cmd.controlling.unmask('all')
            self.cmd.delete('SphericalSelection')
            self.logic_handler._set_selection_type(None)

    def do_dirty(self):
        """
            TODO: Find a better way to detect movement of pseudoatom: Do_position? / Save last positon and only call do move if it actually moved
            TODO STRUC: Move this out to a provider function. The GUI handling program should not deal with special cases for certain selction types.

        Returns
        -------

        """
        print('dirty', mv=0)

        self.redraw()  # only for debugging

        if isinstance(self.logic_handler.my_selection, SphericalSelection):
            if not self.check_objects_exists('SphericalSelection'):
                self.provide_services_to_spherical_selection()

            # Find the coordinates of the Spherical selection
            coords = [None]
            self.cmd.iterate_state(-1, selection='SphericalSelection', expression='coords[0]=(x,y,z)', space=locals())
            x, y, z = coords[0]
            self.trigger_event(program_states.EventType.MOVE, new_x=x, new_y=y, new_z=z)

    def do_key(self, k, x, y, mod):
        """

        Parameters
        ----------
        k : int
            pressed key
        x :  float
            coordinate x
        y : float
            coordinate y
        mod :
            shift pressed

        Returns
        -------

        """
        # Does not deal with events passed down to Selection
        # Check what the current test mode is and do the appropriate action at each key. Problem: FBefore test mode 4 i did not use crd and self.test_mode but just checked if attributes set by the test functions are present

        print(k, x, y, mod, mv=0)

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

        elif k == 13 and mod == 1:  # shift Enter
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

        elif k in range(48, 58):  # Number - #TODO: seems not to work
            # Show the corresponidng moleucle in test mode
            num = k - 48
            print(str(num), mv=4)
            print(k, x, y, mod, mv=4)

            if self.check_results_mode == 4:
                if num > 0 and num <= len(self.pymol_molecule_objects):
                    self._check_optimzer_results_pairwise_4(mol1=self.pymol_molecule_objects[(num - 1)],
                                                            mol2=self.check_result_dict['mol1'])

        # KEYBOARD SHORTCUTS FOR USEFULL FUNCTIONS. Warning: Do not use alt or Ctrl for shortcuts. Many keys are already occupied
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
            self.trigger_event(program_states.EventType.CONFIRM)
            self._action_button_pressed('TreeHeuristicOptimizer')

            # Make sure Optimizatyion was sucessfull
            if (len(self.logic_handler.selected_restraints) > 0):
                self._check_optimzer_results_pairwise_6()

        # Shf T: Test Optimizers
        if k == 84 and mod == 1:
            print('Testing Optimizers...', mv=3)
            from restraintmaker.algorithm import Optimizer
            from PyQt5.QtWidgets import QFileDialog, QInputDialog

            Optimizer.compare_pair_optimizers(criterion=Optimizer._calculate_value_convex_hull, \
                                              atoms=self.logic_handler.selected_atoms, \
                                              opt_types=[Optimizer.GreedyGraphOptimizer,
                                                         Optimizer.GreedyGraphOptimizer,
                                                         Optimizer.GreedyGraphOptimizer,
                                                         Optimizer.BruteForceRingOptimzer,
                                                         Optimizer.BruteForceRingOptimzer], \
                                              opt_args=[(4, 1.0, 'prim', 'convex_hull'), (4, 1.0, 'cog', 'convex_hull'),
                                                        (4, 1.0, 'shortest', 'convex_hull'), (4, 1.0, 'pca', 'convex_hull'),
                                                        (4, 1.0, 'convex_hull', 'convex_hull')], \
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
            self.cmd.refresh()
            file = qt_dialogs.create_file_save_dialog()(
                "Where should the snapshot be saved?")  # '~/Desktop/'+str(self.quicksnapshot_nr)+'.png'
            self.quicksnapshot_nr += 1
            self.cmd.png(file)
            self.cmd.refresh()
            self.cmd.refresh_wizard()
            print('Saved picture at ' + file, mv=3)

    def do_pick(self, bondFlag):
        self.do_select('pk1')
        self.cmd.unpick()

    def do_select(self, sele):
        """
        ..function: do_generate_restraints

        """

        newly_selected_atoms: t.List[u.Atom] = pu.pymol_selection_to_atom_list(
            sele)  # sele is deleted at the end of the function. So it now only contains atoms selected by mouse since  the last call of do select

        # TODO: Introduce switch table to specify action depending on the chosen selection =>

        self.trigger_event(program_states.EventType.SELECT, new_atoms=newly_selected_atoms)
        self.cmd.delete(
            'sele')  # Necessary, to assure that in the next call of do_select we do not get all atoms selected up to now as well

        # TODO: Change program flow:
        # Mode1: IF we have an active restraint type, the selection will be used to create a new restraint and add this to a list of restraints
        # Mode2: If there is no active restraint type, we select atoms and do not create restraints and just wildly select and filter atoms
        # At the moment we only have mode 1

    """
        Buttton functions
    """

    def _action_button_pressed(self, x_name: str):

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
                if r_type.__name__.replace("_Type", "")  == x_name:
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
        self.update_menu()
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
        self.update_menu()
        self.get_panel()
        self.cmd.refresh_wizard()

    def _reset(self):
        self.cmd.unpick()
        self.cmd.deselect()
        self.cmd.set('mouse_selection_mode', 0)

        self.check_result_dict = {}  # crd (check_results_dict) is an dict that can beused by the test mode to store information
        self.check_results_mode = 0  # Check results mode indicates if we are using a special test mode at the moment.
        self.tmp_cnt = 0  # The different check_opzimier_results functions indicate in whihc mode they are called

        self.logic_handler.selected_atoms = []
        self.logic_handler.selected_restraints = []

        self.logic_handler.set_action_states()
        self.redraw()
        cmd.hide("spheres") #this is a bit hacky

        self.cmd.refresh_wizard()
        pass


    """
        Utility functions
    """

    def get_prompt(self):
        '''called automatically upon refresh wizard'''
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
                else:
                    mode = 'DELETE RESTRAINTS MODE'
                    instruction = 'delete restraints by clicking atoms'

            prompt += [mode, instruction]

            # How to use a selection
            if self.logic_handler.current_selection_type != None:
                prompt += ['', 'SELECTION', \
                           'select: lft clck', \
                           'cnfrm: shft+up']
                if 'SphericalSelection' in self.cmd.get_object_list():
                    prompt += [
                        'size: left / right', \
                        'move/rotate: mouse, left', \
                        'activate rota: up',
                        'activate move: down']

            # How to get into test mode
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
        menu += [[str(program_states.ActionState.CURRENTLY_DISABLED),  "RestraintMaker", ""]]  #Should actually be Always_Disabled, but not nice colored!
        menu += [[0, '--------------', '']]
        menu += [[str(program_states.ActionState.CURRENTLY_DISABLED),  "Selection Options: ", ""]]  #Should actually be Always_Disabled, but not nice colored!
        menu += [[str(self.logic_handler._action_states['toggle_select_delete']), text_atom_restraint_mode,
                  'cmd.get_wizard()._toggle_button_pressed(\"toggle_atom\")']]
        menu += [[str(self.logic_handler._action_states['toggle_select_delete']), text_select_delete,
                  'cmd.get_wizard()._toggle_button_pressed(\"toggle_delete\")']]
        menu += [[str(self.logic_handler._action_states['Reset']), 'Reset',  'cmd.get_wizard()._reset()' ]]
        menu += empty_panel
        menu += [[str(self.logic_handler._action_states['Restraint']),
                  'Type:' + get_name(self.logic_handler.current_restraint_type).replace("_Type", ""), 'Restraint']]
        menu += [[str(self.logic_handler._action_states['Importer']),
                  'Import:' + get_name(self.logic_handler.current_importer_type), 'Import']]
        menu += [[str(self.logic_handler._action_states['Selection']),
                  'Mode: ' + get_name(self.logic_handler.current_selection_type), 'Selection']]
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

    def provide_services_to_spherical_selection(self):
        # 1st call: Draw Pseudoatom, set mouse function. Register with dirty listener, mask atoms
        self.cmd.pseudoatom("SphericalSelection")
        self.cmd.mask('not SphericalSelection')

        self.cmd.show('sphere', 'SphericalSelection')
        self.cmd.set('sphere_scale', self.logic_handler.my_selection.radius, 'SphericalSelection')
        self.cmd.set('sphere_transparency', 0.5, selection='SphericalSelection')

        self.cmd.controlling.button('Left', 'None', 'MovA')

    """
        Utility functions that need to be part of the Distance_Restraints class
    """

    def create_pymol_objects_for_restraints(self):
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
        test_mode_hash = hash(str(self.check_results_mode) + str(self.check_result_dict))

        if atoms_hash != self.last_hashes['atoms'] or selection_hash != self.last_hashes['selection'] \
                or restraints_hash != self.last_hashes['restraints'] or test_mode_hash != self.last_hashes['test_mode']:

            # 0) Recolour all atoms to standard colour
            try:
                self.cmd.util.cba('vanadium')  # The colour for carbon atoms needs to be given explicitly. WARNING; The python cba command does not work properly. It does not work at all with gray. It raises an exception with vanadium, but it does paint it corectly.
            except:
                pass

            # 1) Colour selected atoms
            pu.help_pymol_with_big_atom_list(self.cmd.color, self.logic_handler.selected_atoms, color='copper')

            # 2) Colour atoms in current selection
            if self.logic_handler.my_selection != None:
                pu.help_pymol_with_big_atom_list(self.cmd.color, self.logic_handler.my_selection.atoms, color='hassium')

            # 3) Colour atoms in restraints
            # Have the restraints changed?
            if restraints_hash != self.last_hashes['restraints']:
                self.last_hashes.update(restraints=restraints_hash)
                self.create_pymol_objects_for_restraints()

            already_colored=False
            restrained_atoms = []
            for i, r in enumerate(self.logic_handler.selected_restraints):
                pair = []
                for a in r.atoms:
                    if(a in restrained_atoms):
                        already_colored = True
                    restrained_atoms.append(a)
                    pair.append(a)
                if(already_colored):
                    pu.help_pymol_with_big_atom_list(self.cmd.color, pair, color=i%54 )
                else:
                    pu.help_pymol_with_big_atom_list(self.cmd.set, pair, name="sphere_color", value=i % 54)

                already_colored=False

            # 3a) Colour all restrained atoms
            cmd.hide("spheres", "all")
            pu.help_pymol_with_big_atom_list(self.cmd.show, restrained_atoms, representation='spheres')


            # 3b) In test mode: Colour atoms in restraints we are looking at at the moment
            if self.check_results_mode in [4, 6] and test_mode_hash != self.last_hashes['test_mode']:
                restrained_atoms = []

                if self.check_results_mode == 4:
                    for r in self.logic_handler.selected_restraints:
                        if (r.atoms[0].resi == self.check_result_dict['mol1'][4:] and r.atoms[1].resi ==
                            self.check_result_dict['mol2'][4:]) \
                                or (
                                r.atoms[0].resi == self.check_result_dict['mol2'][4:] and r.atoms[1].resi ==
                                self.check_result_dict['mol1'][4:]):
                            for a in r.atoms: restrained_atoms.append(a)
                if self.check_results_mode == 6:
                    for r in self.logic_handler.selected_restraints:
                        m1, m2 = self.check_result_dict['pairs'][self.check_result_dict['i_pair']]
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

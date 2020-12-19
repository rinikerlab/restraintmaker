"""
.. automodule::  restraintmaker_Logic
    :members:
    TODO:write DOCU!

Contains the class Logic_Handler, which defines the program flow, decides which user activities are allowed and manages Importers, Selections,
Restraints, Filters, Optimzers and Exporters
"""
import time
import typing as t

from restraintmaker.algorithm import Filter
from restraintmaker.algorithm import Optimizer
from restraintmaker.algorithm import Selection
from restraintmaker.interface_Pymol.Pymol_Utitlities import \
    _convert_molecules_to_pdb  # TODO: Get independent of interface_Pymol +. wirte own conversion function
from restraintmaker.io import Exporter
from restraintmaker.io import Importer
from restraintmaker.types import Restraint
from restraintmaker.utils import Utilities as u


# Override print functions
# Program verbosity can be set in utilities
def print(*args, mv=0):
    u.print(*args, mv=mv)


class Logic_Handler:
    """class Logic_Handler: Handles the program flow. By defining it inside a class we get proper encapsulation and prevent problems with multiple imports

        When a GUI Program Launches it should create one instance of a Logic_Handler and use it to control the program flow.
        It should pass GUI-events to it (by calling the react_to_event method.). The GUI Module should NOT need to directly access any other modules. (Except utilities).
        Warning: To keep things clean there is not backwards communicatopn from the Logic_Handler to the GUI module. It is the GUI Modules responsibility to check on the current state of the program, after creating events.
        Warning: I really do NOT want the GUI-Program to inherit form the Logic_Handler. It will inherit from a GUI Element (e.g. interface_Pymol.Wizard). Python suppports multiple inheritance, but I don't.

        RECOMMENDATIONS TO WRITE A NEW GUI MODULE:
        The GUI Module should have the following functions/Elements

        1) Buttons/Lists... for the different actions that can be performed (Import, RestraintType, Selection ...). These should call the function set_importer_type, set_restraint_type ...
                The GUI Program should use the action_states dict, to check which actions are currently available.
        2) 2 Toggle buttons / radio buttons... to switch between select / delete mode and atom / restraint mode

        2) Events: The GUI Module is responsible for translating whatever events it understands into events relevant to the selections. It can do this by calling the function react_to_event(event_type, kw_args)
           Recommended translations:    SELECT: When the user clicks a certain Atom
                                        MOVE: Mouse movement / keyboard arrows ...
                                        SIZE: Mouse wheel / keyboard arrows
                                        CONFIRM: Double Click, Enter ...

             The GUI module should check the state of the program (event_state dict and selected_atoms, selected-restraints) every time after calling react_to_event.
             I recommend wrapping the react_to_event function into a function of the GUI module to do that automatically.

        3) Some Selections (SphericalSelection & subclasses) are more intuitive if the GUI program draws a sphere at the corresponding position. => If you want that you have to introduce special for the GUI module to check the current selection
    """

    def __init__(self, all_atoms: t.List[u.Atom]):
        '''
        :param all_atoms: List of all atoms in the Molecules, which should be connected
        :type all_atoms: t.List[u.Atom]
        '''

        self.all_atoms: t.List[u.Atom] = all_atoms
        self.selected_atoms: t.List[u.Atom] = []
        self.selected_restraints: t.List[Restraint._Restraint] = []

        # Define Restraint-Types
        self.available_restraint_types = u.get_all_subclasses(Restraint._Restraint)
        self.current_restraint_type = None
        self.my_restraint = None

        # Define Importer-Types
        self.available_importer_types = u.get_all_subclasses(Importer._Importer)
        self.current_importer_type = None
        self.my_importer = None

        # Define Selection-Types
        self.available_selection_types: t.list[t.type] = u.get_all_subclasses(Selection._Selection)
        self.current_selection_type = None
        self.my_selection = None  # Unlike my_filter my_selection needs to be accessed by several methods: do_pick, set_selection

        # Define Filter-Types
        self.available_filter_types: t.list[t.type] = u.get_all_subclasses(Filter._Filter)
        self.current_filter_type = None
        self.my_filter = None

        # Define-Optimizer Types
        self.available_optimizer_types: t.List[t.type] = u.get_all_subclasses(Optimizer._Optimizer)
        self.current_optimizer_type = None
        self.my_optimizer = None

        # Define-Exporter-Types
        self.available_exporter_types = u.get_all_subclasses(Exporter._Exporter)
        self.current_exporter_type = None
        self.my_exporter = None

        # Modes & Action_States (Modes: What do to atoms, when they are selected. Action states: Which actions are possible right now)
        # TODO CLEAN (1): If we ever use more porgram modes it mioght be cleaner to make the enum u.ActionStates a class,
        # which contains one attribute indicating actions AND one for each mode
        self.select_or_delete_mode = True  # True = select, False = delte
        self.atom_or_restraint_mode = True  # True = atom, False = Restraint

        self._action_states: dict = {
            'toggle_select_delete': u.ActionState.ALWAYS_ENABLED,
            'toggle_atom_restraint': u.ActionState.ALWAYS_ENABLED,
            'Importer': u.ActionState.ALWAYS_DISABLED,
            'Restraint': u.ActionState.ALWAYS_DISABLED,
            'Selection': u.ActionState.ALWAYS_DISABLED,
            'Filter': u.ActionState.ALWAYS_DISABLED,
            'Optimizer': u.ActionState.ALWAYS_DISABLED,
            'Exporter': u.ActionState.ALWAYS_DISABLED,
            'Done': u.ActionState.ALWAYS_ENABLED}
        self.set_action_states()

    # END OF __init___

    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    # -----------------------------------------------------------------event_functions--------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    # These event functions will be called by the GUI-Program to indicate an action. It is the GUI program that defines what exactly consitutues what event. (e.g: a move event might be triggered by a mouse movement or keyboatrd arrows etc.depending on the GUI-Program)
    # CAREFULL: Never call the Slections _update*functions alone. Always check their state (finished or not & selected_atoms) directly afterwards in teh GUI program

    def react_to_event(self, event_type: u.EventType, **kw_args: t.Any):
        '''
        react_to_event should be called by the GUI-module to report an event

        react_to_event will then call the relevant update function of the active selection and refresh all atom lists etc.

        :warning: react_to_event will change the state of the program (i.e. selected_atoms, selected_restraints, mode and actions states)
        It is strongly recommended to wrap react_to_event into a function of the GUI-module, which will check these and redraw the scene and de?activate buttons, if necessary/.

        :param event_type: What kind of event was triggered?
        :type event_type: u.EventType
        :param kw_args: Arguments the relevant update function will need
        :type kw_args: t.Any
        :raises u.BadArgumentException if args to not match the the corresponding logic_handler.react_to_*_event_function
        :raises NotImplementedError if a new event_type is used, which the function does not know about.
        :return: None
        :rtype: None
        '''

        # Note: Wierd python Syntax:  The * operator means two different things here: in the function definition **kw_args indicates, that I want to take an arbitrary number of keyword arguments and pack them all INTO a dict called args.
        #                                                in the function calls ** kw_args means, that I want to UNPACK a dict called kwargs s and pass each value as a separate keyword-argument.

        if self.my_selection == None or self.my_selection.has_finished:
            return

        try:
            if event_type == u.EventType.SELECT:
                self.my_selection._update_select(**kw_args)
            elif event_type == u.EventType.MOVE:
                self.my_selection._update_move(**kw_args)
            elif event_type == u.EventType.SIZE:
                self.my_selection._update_size(**kw_args)
            elif event_type == u.EventType.CONFIRM:
                self.my_selection._update_confirm(**kw_args)
            else:
                raise NotImplementedError('No action defined for event_type ' + str(event_type))

        except TypeError as err:
            raise u.BadArgumentException(
                'The provided arguments do not fit the function logic_handler.react to_ ' + str(
                    event_type) + ' _event') from err

        self._check_selection_status()
        self.set_action_states()

    def set_action_type(self, x_type):
        '''set_action should be called by the GUI module when the user selects an Importer, Restraint, Selection, Filter, Optimizer or Exporter

        It will then call the relevant _set_X_type function.

        :warning: set_action_type will change the state of the program (i.e. selected_atoms, selected_restraints, mode and actions states)
        It is strongly recommended to wrap react_to_event into a function of the GUI-module, which will check these and redraw the scene and de?activate buttons, if necessary/.
        :param type: A type of Importer, Selection, Restraint, Filter, Optimizer or Exporter
        :type type: type
        :raise: TypeError if the specified type is not recognized
        '''

        # SWITCH-TABLE
        if x_type in self.available_importer_types:
            self._set_importer_type(x_type)
        elif x_type in self.available_restraint_types:
            self._set_restraint_type(x_type)
        elif x_type in self.available_selection_types:
            self._set_selection_type(x_type)
        elif x_type in self.available_filter_types:
            self._set_filter_type(x_type)
        elif x_type in self.available_optimizer_types:
            self._set_optimizer_type(x_type)
        elif x_type in self.available_exporter_types:
            self._set_exporter_type(x_type)
        # Not Found
        else:
            raise TypeError(str(
                x_type) + ' is not an available type of Importer, Restraint, Selection, Filter, Optimizer  or Exporter. ')

        self._check_selection_status()
        self.set_action_states()

    # TODO CLEAN (1): If we ever use more porgram modes it might be cleaner to combine action state and mode into one class
    # and set_action_state and set_mode into one function.
    def set_mode(self, select_delete: bool, atom_restraint: bool):
        '''
        set_mode should be called by the GUI module when the user toggles the buttons to change the program mode

        warning: set_action_type will change the state of the program (i.e. selected_atoms, selected_restraints, mode and actions states)
        It is strongly recommended to wrap react_to_event into a function of the GUI-module, which will check these and redraw the scene and de?activate buttons, if necessary/.
        :param select_delete: Select mode (True) or Delete mode (False)?
        :type select_delete: bool
        # :param atom_restraint: Atoms mode (True) or Delete mode (False)?
        :type atom_restraint: bool
        :return: None
        :rtype: None
        '''
        self.select_or_delete_mode = select_delete
        self.atom_or_restraint_mode = atom_restraint

        # 'Depending' on the mode some actions (e.g creating restraints) are not possible

        # Select or delte atoms mode
        if self.atom_or_restraint_mode:
            self._set_restraint_type(None)  # restraint Option not needed

        # Select restraint mode
        elif (not self.atom_or_restraint_mode) and self.select_or_delete_mode:
            self._set_restraint_type(
                None)  # Restraint should be chosen by hand. Reset to None in case an old Restraint type is still active
            self._set_selection_type(None)  # Selection Type will be set by the Restraint

        # Delte Restraint Mode
        else:
            self._set_restraint_type(None)  # To clerly show, that any restraint can be deleted
            self._set_selection_type(Selection.SingleAtomSelection)

        self.set_action_states()

    def _check_selection_status(self):
        """
        _check_selection_status needs to be called at the end of every react_to_*_event function, after updating the selection.

         It checks if the selection has finished and adds the selected atoms to the list, if necessary.
        :return: -
        :rtype: -
        """
        if self.my_selection == None:
            return

        if self.my_selection.has_finished:

            # CHECK ALL FOR POSSIBLE MODES

            # SELECT ATOM MODE
            if self.atom_or_restraint_mode and self.select_or_delete_mode:
                self.selected_atoms.extend(filter(lambda x: not x in self.selected_atoms,
                                                  self.my_selection.atoms))  # Do not allow addition of duplicates
                self.my_selection.reset()

            # DELETE ATOM MODE
            elif self.atom_or_restraint_mode and not self.select_or_delete_mode:
                self.selected_atoms = list(filter(lambda x: not (x in self.my_selection.atoms), self.selected_atoms))
                self.my_selection.reset()

            # SELECT RESTRAINT MODE
            elif not self.atom_or_restraint_mode and self.select_or_delete_mode:
                new_restraint = self.current_restraint_type(self.my_selection.atoms)
                self.selected_restraints.append(new_restraint)
                self.my_selection.reset()

            # DELETE RESTRAINT MODE
            else:
                # In delete Restraint Mode the Selection is a Single Atom Selection
                # TODO DELETE_RESTRAINT_MODE (2): Allow to use every kind of Selection in Delete-Restraint_mode. Upon confirmation: Delete everz restraint containing an atom in that selection
                atom = self.my_selection.atoms[0]
                # Remove all restraints containing the selected atom
                self.selected_restraints = [r for r in self.selected_restraints if not atom in r.atoms]
                # open a new SingleAtomSelection
                self.my_selection.reset()

    # END OF _check_selection_status

    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    # -----------------------------------------------------------------set_*_type()-----------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------------------------------------------

    # Are called by set_action_type

    # used to be called do_import
    def _set_importer_type(self, new_importer_type: t.Type[Importer._Importer]):
        '''
        _set_importer_type changes the current Importer type, and creates and apllies a new Importer of that type.

        :warning: Should only be called via set_action_type which handles type checking
        :param new_importer_type: An available type of Importer
        :type new_importer_type: type[Importer]
        :return: None
        :rtype: None
        :raises TypeError if new_importer_type is not in available_importer_types
        '''

        if new_importer_type == None:
            self.current_importer_type = None
            self.my_importer = None
            return

        if not new_importer_type in self.available_importer_types:
            raise TypeError(str(new_importer_type) + ' is not an available type of Importer')

        self.my_importer = new_importer_type(self.all_atoms)
        if try_to_get_args(self.my_importer, u.create_file_open_dialog()):
            self.selected_restraints = self.my_importer.import_restraints()

        else:
            self.my_importer = None
            self.current_importer_type = None

    def _set_restraint_type(self, new_restraint_type: t.Type[Restraint._Restraint]):
        '''
        _set_restraint_type changes the current type of Restraint.

        :warning: Should only be called via set_action_type which handles type checking
        :param new_restraint_type:  A type of restraint
        :type new_restraint_type: type[Restraint._Restraint}
        :return: None
        :rtype: None
        :raises: TypeError if new_restraint_type is not in available_restraint_types
        '''

        if new_restraint_type == None:
            self.current_restraint_type = None
            self.my_restraint == None
            return

        if not new_restraint_type in self.available_restraint_types + [None]:
            raise TypeError(str(new_restraint_type) + ' is not an available type of Restraint')

        print("SET TYPE: " + str(new_restraint_type))
        self.current_restraint_type = new_restraint_type

        # Set the Selection type to the preferred selection of this restraint type. set_mode will deactivate the selection, button, so it cant be changed
        # TODO: In theory a restraint could accept several selection types. Would be nice to give the user the possibilitz to choose between all of them
        # TODO: If the first selection in the list is e.g a Limited selection the user will be asked for input. At a moment I just check for all lists that that is not the case.
        #   Would be nicer if there was a way to pass an input function that does not require user input in this context
        self._set_selection_type(self.current_restraint_type.accepted_selection_types[0])

    def _set_selection_type(self, new_selection_type: type, ):
        '''
        _set_selection_type chagnes the current type of Selection and starts a new selection of that type.

        :warning: Should only be called via set_action_type which handles type checking
        :param new_selection_type: type of the new Selection
        :type new_selection_type: t.Type[Selection]
        :return: None
        :rtype: None
        '''

        if new_selection_type == None:
            self.current_selection_type = None
            self.my_selection = None
            return

        if not new_selection_type in self.available_selection_types:
            raise TypeError(str(new_selection_type) + ' is not an  available type of Selection')

        self.current_selection_type = new_selection_type
        self.my_selection = new_selection_type(self.all_atoms)

        # Switchtable for different provider funcs
        # TODO STRUC: Remove special cases, use a generic input function

        input_function = u.do_nothing
        if self.my_selection.__class__ == Selection.LimitedSelection:
            input_function = u.create_input_dialog(parent_window=None, title='Selection')
        elif isinstance(self.my_selection, Selection.SphericalSelection):
            input_function = lambda _: 5
        elif isinstance(self.my_selection, Selection.MCS_Selection):
            input_function = lambda dummy_arg: _convert_molecules_to_pdb()

        if not try_to_get_args(self.my_selection, input_function):
            self.my_selection = None
            self.current_selection_type = None

        # Unlike other classes Selection does not call its own method (update(...)) here. It is just started. update(...) is called at every pick

    # Actually: Also consider just keeping it simple, as it is
    def _set_filter_type(self, new_filter_type: type):
        '''
        _set_filter_type changes the current type of Filter and creates and applies a Filter of that type.

        :warning: Should only be called via set_action_type which handles type checking
        :param new_filter_type: type of the new Filter
        :type new_filter_type: t.Type[Filter]
        :return: None
        :rtype: None
        '''

        if new_filter_type == None:
            self.current_filter_type = None
            self.my_filter = None
            return

        # Check for erros the GUI module hsould already have dealt with
        if not new_filter_type in self.available_filter_types + [None]:
            raise TypeError(str(new_filter_type) + ' is not an available type of Filter')

        # Create list of all atoms in the current selection:

        self.my_filter = new_filter_type(self.selected_atoms)

        # Switchtable for different provider funcs
        # TODO STRUC: Remove special cases, use a generic input function

        input_function = u.do_nothing
        if isinstance(self.my_filter, Filter.PropertyFilter):
            input_function = u.create_input_dialog(parent_window=None, title='PropertyFilter')
        elif isinstance(self.my_filter, Filter.ElementFilter):
            input_function = u.create_input_dialog(parent_window=None, title='ElementFilter')
        elif isinstance(self.my_filter, Filter.RingFilter):
            # TODO: GET INDEPENDENT OF PYMOL HERE => Write an own pdb function
            input_function = lambda dummy_arg: _convert_molecules_to_pdb()

        if try_to_get_args(self.my_filter, input_function):
            self.selected_atoms = self.my_filter.filter()
        else:
            self.my_filter = None
            self.current_filter_type = None

    def _set_optimizer_type(self, new_optimizer_type: type):
        '''
        _set_optimizer_type changes the current type of Optimizer and creates and applies a Optimizer of that type.

        :warning: Should only be called via set_action_type which handles type checking.
        :param new_optimizer_type: type of the new Filter
        :type new_optimizer_type: t.Type[Filter]
        :return: None
        :rtype: None
        '''
        if new_optimizer_type == None:
            self.current_optimizer_type = None
            self.my_optimizer = None
            return

        if not new_optimizer_type in self.available_optimizer_types + [None]:
            raise TypeError(str(new_optimizer_type) + ' is not an available type of Optimizer')

        self.my_optimizer = new_optimizer_type(self.selected_atoms)

        # TODO: See TODO in u.creat_multi_dialog
        # SWITCHTABLE FOR input function
        input_functions = u.do_nothing()
        # TODO: More elegant than these hardcoded input functions: Give each Optimizer etc. an attribute specifzing which args it needs, which we can than pass to create_multi_dialog
        # TODO: Generate new field in pyqt win to select if you want a chain, ring, or all to all connection

        if isinstance(self.my_optimizer, Optimizer.TreeHeuristicOptimizer):
            input_function = u.create_multi_dialog(title='Parameters for TreeHeuristicOptimizer', \
                                                   inputs=['number of restraints', 'maximal distance of restraints',
                                                           'tree-algorithm', 'optimize molecules pairs by'], \
                                                   options={'tree-algorithm': ['shortest', 'cog', 'prim', "biased_avg"],
                                                            'optimize molecules pairs by': ['None', 'convex_hull',
                                                                                            'pca_2d']}, \
                                                   default={'number of restraints': '4',
                                                            'maximal distance of restraints': '1.2',
                                                            'algorithm': 'shortest',
                                                            'optimize molecules pairs by': 'pca_2d'})


        elif isinstance(self.my_optimizer, Optimizer.BruteForceRingOptimzer):
            input_function = u.create_multi_dialog(title='Parameters for BruteForceOptimizer', \
                                                   inputs=['number of restraints', 'maximal distance of restraints',
                                                           'algorithm', 'optimize molecules pairs by'], \
                                                   options={'algorithm': ['convex_hull', 'pca'],
                                                            'optimize molecules pairs by': ['None', 'convex_hull',
                                                                                            'pca_2d']}, \
                                                   default={'number of restraints': '4',
                                                            'maximal distance of restraints': '1.2', 'algorithm': 'pca',
                                                            'optimize molecules pairs by': 'pca_2d'})

        elif isinstance(self.my_optimizer, Optimizer.MetaMoleculeRingOptimizer):
            input_function = u.create_multi_dialog(title='Parameters for BestMoleculeRingOptimizer', \
                                                   inputs=['number of restraints', 'maximal distance of restraints',
                                                           'algorithm', 'optimize molecules pairs by'], \
                                                   options={'algorithm': ['convex_hull', 'pca'],
                                                            'optimize molecules pairs by': ['convex_hull',
                                                                                            'pca_2d']}, \
                                                   default={'number of restraints': '4',
                                                            'maximal distance of restraints': '1.2', 'algorithm': 'pca',
                                                            'optimize molecules pairs by': 'pca_2d'})

        if try_to_get_args(self.my_optimizer, input_function):
            time_start = time.time()

            try:
                self.selected_restraints = self.my_optimizer.make_restraints()
                time_stop = time.time()
                print('Optimized in ' + '{:0.1f}'.format(time_stop - time_start) + ' s', mv=3)
            except u.NoOptimalSolutionException as ex:
                print('Failed to find a optimal solution under the given conditions:', ex, mv=4)
        else:
            self.my_optimizer = None
            self.current_optimizer_type = None

    def _set_exporter_type(self, new_exporter_type: type):
        '''
            _set_exporter_type changes the current type of Exporter and creates and applies an Exporter of that type.

            :warning: Should only be called via set_action_type which handles type checking.
            :param new_exporter_type: type of the new Filter
            :type new_exporter_type: t.Type[Exporter]
            :return: None
            :rtype: None
            '''

        if new_exporter_type == None:
            self.current_exporter_type == None
            self.my_exporter == None
            return

        if not new_exporter_type in self.available_exporter_types + [None]:
            raise TypeError(str(new_exporter_type) + ' is not an available type of Exporter')

        self.my_exporter = new_exporter_type(self.selected_restraints)
        if try_to_get_args(self.my_exporter, u.create_file_save_dialog()):  # u.create_input_dialog(None, 'Exporter')
            self.my_exporter.export_restraints()
        else:
            self.my_exporter = None
            self.current_exporter_type = None

    def set_action_states(self):
        """set_action_statess checks the current state of the program and decides, which actions are allowed at the moment"""

        # ALL CHANGES TO THE ACTION STATE SHOULD BE DONE IN THIS FUNCTION. THE STATES SHOULD BE DETERMINED ONLY BY THE CURRENT STAte of THE PROGRAM
        # ORDER: Set types for each action one after the other
        # if the structure beomces more complex consider first checking the toggle modes (atom?restraint, delete/select) and go through each state in iach state

        self._action_states['toggle_select_delete']: u.ActionState.ALWAYS_ENABLED
        self._action_states['toggle_atom_restraint']: u.ActionState.ALWAYS_ENABLED
        self._action_states['Importer'] = u.ActionState.CURRENTLY_ENABLED if len(
            self.selected_restraints) == 0 else u.ActionState.CURRENTLY_DISABLED  # TODO Later: Check if there is any Molecules loaded.
        self._action_states['Restraint'] = u.ActionState.CURRENTLY_ENABLED if (
                                                                                  not self.atom_or_restraint_mode) and self.select_or_delete_mode \
            else u.ActionState.CURRENTLY_DISABLED  # TODO Later: Always Enabled, can be changed as much as we like
        self._action_states[
            'Selection'] = u.ActionState.CURRENTLY_ENABLED if self.atom_or_restraint_mode else u.ActionState.CURRENTLY_DISABLED  # TODO Later: Check if any atoms are loaded
        self._action_states['Filter'] = u.ActionState.CURRENTLY_ENABLED if len(self.selected_atoms) > 0 \
                                                                           and len(
            self.selected_restraints) == 0 else u.ActionState.CURRENTLY_DISABLED
        self._action_states['Optimizer'] = u.ActionState.CURRENTLY_ENABLED if len(self.selected_atoms) >= 2 \
                                                                              and len(
            self.selected_restraints) == 0 else u.ActionState.CURRENTLY_DISABLED
        # TODO: Check if we have Atoms from at least 2 Molecules
        # TODO: Check indivudiual Optimizers depending on Restraint type
        self._action_states['Exporter'] = u.ActionState.CURRENTLY_ENABLED if len(self.selected_restraints) > 0 \
            else u.ActionState.CURRENTLY_DISABLED
        self._action_states['Done'] = u.ActionState.ALWAYS_ENABLED


# ------------------END OF CLASS------------------------------------------------------
def try_to_get_args(my_x, input_function) -> bool:
    '''
    wrapper function for get_args. Takes care of the error handling

    :param my_x: An Object that should call get args
    :type my_x: _Importer,_Restraint,_Selection,_Filter,_Optimzer or _Exporter
    :param input_function: The input_function argument of get_args
    :type input_function: t.Callable[str]
    :return: Could get_args be executed without errors?
    :rtype: bool
    '''

    failed = False

    try:
        my_x.get_args(input_function)
    except u.BadArgumentException as err:
        print('Failed to create a ' + my_x.__class__.__name__ + ' with the given arguments: \n' + str(err), mv=4)
        failed = True
    except IndexError:
        print('Did not receive enough arguments to create a ' + my_x.__class__.__name__, mv=4)
        # TODO: Would be nice to know how many args are needed here. See also TODO before u.create_multi_dialog
        failed = True
    except NotImplementedError:
        print(my_x.__class__.__name__ + 'has not been implemented', mv=4)
        failed = True
    # Do not expect Argumetn error that will ocurr if my_x does not have a get_args function. That is a fault that really should stop the program . => Leave it to standard error handling
    return not failed

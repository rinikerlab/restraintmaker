"""
    Contains the class Logic_Handler, which defines the program flow, decides which user activities are allowed and manages Importers, Selections,
    Restraints, Filters, Optimzers and Exporters.

"""
import time
import typing as t


from restraintmaker.algorithm import Filter
from restraintmaker.algorithm import Optimizer
from restraintmaker.algorithm import Selection
from restraintmaker.interface_Pymol.pymol_utilities.pymol_utitlities import _convert_molecules_to_pdb, try_to_get_args

from restraintmaker.io import Exporter, Importer
from restraintmaker.utils import Utilities as u, Restraints, Restraint_Types
from restraintmaker.interface_Pymol.pymol_utilities import program_states, qt_dialogs
from restraintmaker.utils.Utilities import print


class Logic_Handler:

    def __init__(self, all_atoms: t.List[u.Atom]):
        """
            Handles the program flow. By defining it inside a class we get proper encapsulation and prevent problems with multiple imports

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

        Parameters
        ----------
        all_atoms : t.List[u.Atom]
             List of all atoms in the Molecules, which should be connected

        """

        self.all_atoms: t.List[u.Atom] = all_atoms
        self.selected_atoms: t.List[u.Atom] = []
        self.selected_restraints: t.List[Restraints._Restraint] = []

        # Define Restraint-Types
        self.all_restraint_types = u.get_all_subclasses(Restraint_Types._Restraint_Type)
        self.available_restraint_types = []
        self.current_restraint_type:Restraint_Types._Restraint_Type = Restraint_Types.Distance_Restraint_Type
        self.my_restraint = None

        # Define Importer-Types
        self.available_importer_types = u.get_all_subclasses(Importer._Importer)
        self.current_importer_type = None
        self.my_importer = None

        # Define Selection-Types
        self.available_selection_types: t.List[t.Type] = u.get_all_subclasses(Selection._Selection)
        self.current_selection_type = None
        self.my_selection = None  # Unlike my_filter my_selection needs to be accessed by several methods: do_pick, set_selection

        # Define Filter-Types
        self.available_filter_types: t.List[t.Type] = u.get_all_subclasses(Filter._Filter)
        self.current_filter_type = None
        self.my_filter = None

        # Define-Optimizer Types
        self.available_optimizer_types: t.List[t.Type] = [x for x in u.get_all_subclasses(Optimizer._Optimizer) if (not x == Optimizer.MetaMoleculeRingOptimizer)]
        self.current_optimizer_type = None
        self.my_optimizer = None

        # Define-Exporter-Types
        self.all_exporter_types = u.get_all_subclasses(Exporter._Exporter)
        self.available_exporter_types = u.get_all_subclasses(Exporter._Export_Distance_Restraints)
        self.current_exporter_type = None
        self.my_exporter = None

        # Modes & Action_States (Modes: What do to atoms, when they are selected. Action states: Which actions are possible right now)
        # which contains one attribute indicating actions AND one for each mode

        self._action_states: dict = {
            'toggle_select_delete': program_states.ActionState.ALWAYS_ENABLED,
            'toggle_atom_restraint': program_states.ActionState.ALWAYS_ENABLED,
            'Reset':  program_states.ActionState.ALWAYS_ENABLED,
            'Importer': program_states.ActionState.ALWAYS_DISABLED,
            'Restraint': program_states.ActionState.ALWAYS_DISABLED,
            'Selection': program_states.ActionState.ALWAYS_DISABLED,
            'Filter': program_states.ActionState.ALWAYS_DISABLED,
            'Optimizer': program_states.ActionState.ALWAYS_DISABLED,
            'Exporter': program_states.ActionState.ALWAYS_DISABLED,
            'Done': program_states.ActionState.ALWAYS_ENABLED}

        #initialize
        self.set_mode(select_delete=True, atom_restraint=True)
        self.set_action_states()
        self._set_restraint_type(self.current_restraint_type)


    """
        event_functions
    """

    # These event functions will be called by the GUI-Program to indicate an action. It is the GUI program that defines what exactly consitutues what event. (e.g: a move event might be triggered by a mouse movement or keyboatrd arrows etc.depending on the GUI-Program)
    # CAREFULL: Never call the Slections _update*functions alone. Always check their state (finished or not & selected_atoms) directly afterwards in teh GUI program

    def react_to_event(self, event_type: program_states.EventType, **kw_args: t.Any):
        """
                react_to_event should be called by the GUI-module to report an event

            react_to_event will then call the relevant update function of the active selection and refresh all atom lists etc.

         @Waring react_to_event will change the state of the program (i.e. selected_atoms, selected_restraints, mode and actions states)
         It is strongly recommended to wrap react_to_event into a function of the GUI-module, which will check these and redraw the scene and de?activate buttons, if necessary/.


        Parameters
        ----------
        event_type : program_states.EventType
             What kind of event was triggered?
        kw_args :  dict
             Arguments the relevant update function will need

        Returns
        -------
        NoReturn

        Raises
        ------
        BadArgumentException
            if args to not match the the corresponding logic_handler.react_to_*_event_function
        NotImplementedError
            if a new event_type is used, which the function does not know about.

        """
        if self.my_selection == None or self.my_selection.has_finished:
            return

        try:
            if event_type == program_states.EventType.SELECT:
                self.my_selection._update_select(**kw_args)
            elif event_type == program_states.EventType.MOVE:
                self.my_selection._update_move(**kw_args)
            elif event_type == program_states.EventType.SIZE:
                self.my_selection._update_size(**kw_args)
            elif event_type == program_states.EventType.CONFIRM:
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
        """set_action should be called by the GUI module when the user selects an Importer, Restraint, Selection, Filter, Optimizer or Exporter

            It will then call the relevant _set_X_type function.

            @Warnings set_action_type will change the state of the program (i.e. selected_atoms, selected_restraints, mode and actions states)
            It is strongly recommended to wrap react_to_event into a function of the GUI-module, which will check these and redraw the scene and de?activate buttons, if necessary/.

        Parameters
        ----------
        x_type : type
            A type of Importer, Selection, Restraint, Filter, Optimizer or Exporter

        Returns
        -------
        NoReturn

        Raises
        ------
        TypeError
            if the specified type is not recognized

        """

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

    def set_mode(self, select_delete: bool, atom_restraint: bool):
        """
            set_mode should be called by the GUI module when the user toggles the buttons to change the program mode

        warning: set_action_type will change the state of the program (i.e. selected_atoms, selected_restraints, mode and actions states)
        It is strongly recommended to wrap react_to_event into a function of the GUI-module, which will check these and redraw the scene and de?activate buttons, if necessary/.

        Parameters
        ----------
        select_delete : bool
            Select mode (True) or Delete mode (False)?
        atom_restraint : bool
            Atoms mode (True) or Delete mode (False)?

        Returns
        -------
        NoReturn

        """

        self.select_or_delete_mode = select_delete
        self.atom_or_restraint_mode = atom_restraint

        # 'Depending' on the mode some actions (e.g creating restraints) are not possible

        # Select or delete atoms mode
        if self.atom_or_restraint_mode:
            self.available_restraint_types = [x for x in self.all_restraint_types if(len(x.accepted_optimizer_types)>0)]
            set_type = self.current_restraint_type if(self.current_restraint_type in self.available_restraint_types) else self.available_restraint_types[0]
            self._set_restraint_type(set_type)  # restraint Option not needed

        # Select restraint mode
        elif (not self.atom_or_restraint_mode) and self.select_or_delete_mode:
            self.available_restraint_types = self.all_restraint_types
            set_type = self.current_restraint_type if(self.current_restraint_type in self.available_restraint_types) else self.available_restraint_types[0]
            self._set_restraint_type(set_type)

            # Delete Restraint Mode
        else:
            self._set_restraint_type(None)  # To clearly show, that any restraint can be deleted
            self._set_selection_type(Selection.SingleAtomSelection)

        self.set_action_states()

    def _check_selection_status(self):
        """
         _check_selection_status needs to be called at the end of every react_to_*_event function, after updating the selection.

         It checks if the selection has finished and adds the selected atoms to the list, if necessary.

        Returns
        -------
        NoReturn

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
                if(Restraint_Types.Position_Restraint_Type == self.current_restraint_type):
                    new_restraint = [self.current_restraint_type.restraint(atomA=x) for x in self.my_selection.atoms]
                    self.selected_restraints.extend(new_restraint)
                else:
                    new_restraint = self.current_restraint_type.restraint(*self.my_selection.atoms)
                    self.selected_restraints.append(new_restraint)
                self.my_selection.reset()

            # DELETE RESTRAINT MODE
            else:
                # In delete Restraint Mode the Selection is a Single Atom Selection
                atoms = self.my_selection.atoms
                # Remove all restraints containing the selected atom
                self.selected_restraints = [r for r in self.selected_restraints if not any([atom in r.atoms for atom in atoms])]
                # open a new SingleAtomSelection
                self.my_selection.reset()

    """
        set_*_type - functions:
            Are called by set_action_type
            used to be called do_import
    """

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
        if try_to_get_args(self.my_importer, qt_dialogs.create_file_open_dialog()):
            self.selected_restraints = self.my_importer.import_restraints()

        else:
            self.my_importer = None
            self.current_importer_type = None

    def _set_restraint_type(self, new_restraint_type: t.Type[Restraints._Restraint]):
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
            self.my_restraint = None
            return

        if not new_restraint_type in self.available_restraint_types + [None]:
            raise TypeError(str(new_restraint_type) + ' is not an available type of Restraint')

        print("SET TYPE: " + str(new_restraint_type))
        self.current_restraint_type = new_restraint_type

        # Set the Selection type to the preferred selection of this restraint type. set_mode will deactivate the selection, button, so it cant be changed
        # TODO: In theory a restraint could accept several selection types. Would be nice to give the user the possibilitz to choose between all of them
        # TODO: If the first selection in the list is e.g a Limited selection the user will be asked for input. At a moment I just check for all lists that that is not the case.
        #   Would be nicer if there was a way to pass an input function that does not require user input in this context
        self.set_action_states()
        self.available_selection_types = self.current_restraint_type.accepted_selection_types
        self._set_selection_type(self.available_selection_types[0])
        self.available_optimizer_types = self.current_restraint_type.accepted_optimizer_types
        self.available_exporter_types = self.current_restraint_type.accepted_exporter_types
        self.available_importer_types = self.current_restraint_type.accepted_importer_types


    def _set_selection_type(self, new_selection_type: type, ):
        """
             chagnes the current type of Selection and starts a new selection of that type.

            @warnings: Should only be called via set_action_type which handles type checking

        Parameters
        ----------
        new_selection_type : t.Type[Selection]
            type of the new Selection

        Returns
        -------
        NoReturn

        """

        if new_selection_type == None:
            self.current_selection_type = None
            self.my_selection = None
            return

        if not new_selection_type in self.available_selection_types and self.select_or_delete_mode:
            raise TypeError(str(new_selection_type) + ' is not an  available type of Selection')

        self.current_selection_type = new_selection_type
        self.my_selection = new_selection_type(self.all_atoms)

        # Switchtable for different provider funcs
        # TODO STRUC: Remove special cases, use a generic input function

        input_function = u.do_nothing
        if self.my_selection.__class__ == Selection.LimitedSelection:
            input_function = qt_dialogs.create_input_dialog(parent_window=None, title='Selection')
        elif isinstance(self.my_selection, Selection.SphericalSelection):
            input_function = lambda _: 5
        elif isinstance(self.my_selection, Selection.MCS_Selection):
            input_function = lambda dummy_arg: _convert_molecules_to_pdb()

        if not try_to_get_args(self.my_selection, input_function):
            self.my_selection = None
            self.current_selection_type = None

        # Unlike other classes Selection does not call its own method (update(...)) here. It is just started. update(...) is called at every pick

    def _set_filter_type(self, new_filter_type: type):
        """
            changes the current type of Filter and creates and applies a Filter of that type.

            @Warnings should only be called via set_action_type which handles type checking

        Parameters
        ----------
        new_filter_type :  t.Type[Filter]
            type of the new Filter

        Returns
        -------
        NoReturn

        """

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
            input_function = qt_dialogs.create_input_dialog(parent_window=None,
                                                            title='PropertyFilter')
        elif isinstance(self.my_filter, Filter.ElementFilter):
            input_function = qt_dialogs.create_input_dialog(parent_window=None,
                                                            title='ElementFilter')
        elif isinstance(self.my_filter, Filter.RingFilter):
            # TODO: GET INDEPENDENT OF PYMOL HERE => Write an own pdb function
            input_function = lambda dummy_arg: _convert_molecules_to_pdb()

        if try_to_get_args(self.my_filter, input_function):
            self.selected_atoms = self.my_filter.filter()
        else:
            self.my_filter = None
            self.current_filter_type = None

    def _set_optimizer_type(self, new_optimizer_type: type):
        """
            changes the current type of Optimizer and creates and applies a Optimizer of that type.
            @Warnings  Should only be called via set_action_type which handles type checking.

        Parameters
        ----------
        new_optimizer_type : t.Type[Filter]
            type of the new Filter

        Returns
        -------
        NoReturn

        """

        if new_optimizer_type == None:
            self.current_optimizer_type = None
            self.my_optimizer = None
            return

        if not new_optimizer_type in self.available_optimizer_types + [None]:
            raise TypeError(str(new_optimizer_type) + ' is not an available type of Optimizer')

        self.my_optimizer = new_optimizer_type(self.selected_atoms)

        # SWITCHTABLE FOR input function
        input_functions = u.do_nothing()

        if isinstance(self.my_optimizer, Optimizer.GreedyGraphOptimizer):
            input_function = qt_dialogs.create_multi_dialog(
                title='Parameters for GreedyGraphOptimizer', \
                inputs=['number of restraints', 'maximal distance of restraints',
                        'tree-algorithm', 'optimize molecules pairs by'], \
                options={'tree-algorithm': ['minmax', 'cog', 'prim', "biased_avg"],
                         'optimize molecules pairs by': ['None', 'convex_hull',
                                                         'pca_2d']}, \
                default={'number of restraints': '4',
                         'maximal distance of restraints': '1.0',
                         'algorithm': 'minmax',
                         'optimize molecules pairs by': 'convex_hull'})


        elif isinstance(self.my_optimizer, Optimizer.BruteForceRingOptimzer):
            input_function = qt_dialogs.create_multi_dialog(
                title='Parameters for BruteForceOptimizer', \
                inputs=['number of restraints', 'maximal distance of restraints',
                        'algorithm', 'optimize molecules pairs by'], \
                options={'algorithm': ['convex_hull', 'pca'],
                         'optimize molecules pairs by': ['None', 'convex_hull',
                                                         'pca_2d']}, \
                default={'number of restraints': '4',
                         'maximal distance of restraints': '1.2', 'algorithm': 'pca',
                         'optimize molecules pairs by': 'pca_2d'})

        elif isinstance(self.my_optimizer, Optimizer.MetaMoleculeRingOptimizer):
            input_function = qt_dialogs.create_multi_dialog(
                title='Parameters for BestMoleculeRingOptimizer', \
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
        """
            changes the current type of Exporter and creates and applies an Exporter of that type.
            @Warnings  Should only be called via set_action_type which handles type checking.

        Parameters
        ----------
        new_exporter_type :  t.Type[Exporter]
            type of the new Filter

        Returns
        -------
        NoReturn

        """

        if new_exporter_type == None:
            self.current_exporter_type == None
            self.my_exporter == None
            return

        if not new_exporter_type in self.available_exporter_types + [None]:
            raise TypeError(str(new_exporter_type) + ' is not an available type of Exporter')

        self.my_exporter = new_exporter_type(self.selected_restraints)
        if try_to_get_args(self.my_exporter,
                           qt_dialogs.create_file_save_dialog()):  # u.create_input_dialog(None, 'Exporter')
            self.my_exporter.export_restraints()
        else:
            self.my_exporter = None
            self.current_exporter_type = None

    def set_action_states(self):
        """
            checks the current state of the program and decides, which actions are allowed at the moment

            many TODOs
            
        Returns
        -------
        NoReturn

        """

        # ALL CHANGES TO THE ACTION STATE SHOULD BE DONE IN THIS FUNCTION. THE STATES SHOULD BE DETERMINED ONLY BY THE CURRENT STAte of THE PROGRAM
        # ORDER: Set types for each action one after the other
        # if the structure beomces more complex consider first checking the toggle modes (atom?restraint, delete/select) and go through each state in iach state
        self._action_states['toggle_select_delete']: program_states.ActionState.ALWAYS_ENABLED
        self._action_states['toggle_atom_restraint']: program_states.ActionState.ALWAYS_ENABLED

        self._action_states['Importer'] = program_states.ActionState.CURRENTLY_ENABLED if len(
            self.selected_restraints) == 0 else program_states.ActionState.CURRENTLY_DISABLED  # TODO Later: Check if there is any Molecules loaded.
        self._action_states['Restraint'] = program_states.ActionState.CURRENTLY_ENABLED

        self._action_states['Selection'] = program_states.ActionState.CURRENTLY_ENABLED if (self.atom_or_restraint_mode or (not self.atom_or_restraint_mode and self.current_restraint_type == Restraint_Types.Position_Restraint_Type)) else \
            program_states.ActionState.CURRENTLY_DISABLED  # TODO Later: Check if any atoms are loaded
        self._action_states['Filter'] = program_states.ActionState.CURRENTLY_ENABLED if len(
            self.selected_atoms) > 0 \
                                                                                        and len(
            self.selected_restraints) == 0 else program_states.ActionState.CURRENTLY_DISABLED
        self._action_states['Optimizer'] = program_states.ActionState.CURRENTLY_ENABLED if len(
            self.selected_atoms) >= 2 \
            and len(self.selected_restraints) == 0 else program_states.ActionState.CURRENTLY_DISABLED
        # TODO: Check if we have Atoms from at least 2 Molecules
        # TODO: Check indivudiual Optimizers depending on Restraint type
        self._action_states['Exporter'] = program_states.ActionState.CURRENTLY_ENABLED if len(
            self.selected_restraints) > 0 \
            else program_states.ActionState.CURRENTLY_DISABLED

        self._action_states['Done'] = program_states.ActionState.ALWAYS_ENABLED
        self._action_states['Reset']: program_states.ActionState.ALWAYS_ENABLED


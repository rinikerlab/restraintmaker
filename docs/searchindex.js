Search.setIndex({docnames:["_source/modules","_source/restraintmaker","_source/restraintmaker.algorithm","_source/restraintmaker.interface_Pymol","_source/restraintmaker.interface_Pymol.pymol_utilities","_source/restraintmaker.io","_source/restraintmaker.io.Files","_source/restraintmaker.test","_source/restraintmaker.test.test_files","_source/restraintmaker.tools_Rdkit","_source/restraintmaker.utils","examples/example_script_execution","examples/index","index","introduction"],envversion:{"sphinx.domains.c":2,"sphinx.domains.changeset":1,"sphinx.domains.citation":1,"sphinx.domains.cpp":3,"sphinx.domains.index":1,"sphinx.domains.javascript":2,"sphinx.domains.math":2,"sphinx.domains.python":2,"sphinx.domains.rst":2,"sphinx.domains.std":1,"sphinx.ext.intersphinx":1,"sphinx.ext.viewcode":1,nbsphinx:3,sphinx:56},filenames:["_source/modules.rst","_source/restraintmaker.rst","_source/restraintmaker.algorithm.rst","_source/restraintmaker.interface_Pymol.rst","_source/restraintmaker.interface_Pymol.pymol_utilities.rst","_source/restraintmaker.io.rst","_source/restraintmaker.io.Files.rst","_source/restraintmaker.test.rst","_source/restraintmaker.test.test_files.rst","_source/restraintmaker.tools_Rdkit.rst","_source/restraintmaker.utils.rst","examples/example_script_execution.ipynb","examples/index.rst","index.rst","introduction.rst"],objects:{"":{restraintmaker:[1,0,0,"-"]},"restraintmaker.algorithm":{Filter:[2,0,0,"-"],Optimizer:[2,0,0,"-"],Selection:[2,0,0,"-"]},"restraintmaker.algorithm.Filter":{ElementFilter:[2,1,1,""],PropertyFilter:[2,1,1,""],RingFilter:[2,1,1,""]},"restraintmaker.algorithm.Filter.ElementFilter":{__init__:[2,2,1,""],filter:[2,2,1,""],get_args:[2,2,1,""]},"restraintmaker.algorithm.Filter.PropertyFilter":{__init__:[2,2,1,""],filter:[2,2,1,""],get_args:[2,2,1,""]},"restraintmaker.algorithm.Filter.RingFilter":{__init__:[2,2,1,""],filter:[2,2,1,""],get_args:[2,2,1,""]},"restraintmaker.algorithm.Optimizer":{BruteForceRingOptimzer:[2,1,1,""],MetaMoleculeRingOptimizer:[2,1,1,""],TreeHeuristicOptimizer:[2,1,1,""],_calculate_value_convex_hull:[2,3,1,""],_calculate_value_scaled_pca_2d:[2,3,1,""],_calculate_value_unscaled_pca_2d:[2,3,1,""],calculate_value_pca_relative_to_pca_of_both_molecules:[2,3,1,""],cog_distance_restraint:[2,3,1,""],get_all_short_pair_restraints:[2,3,1,""],maximal_spanning_tree_greedy:[2,3,1,""],maximal_weight_ring:[2,3,1,""],priority_node:[2,1,1,""]},"restraintmaker.algorithm.Optimizer.BruteForceRingOptimzer":{__init__:[2,2,1,""],connect_two_molecules:[2,2,1,""],get_args:[2,2,1,""]},"restraintmaker.algorithm.Optimizer.MetaMoleculeRingOptimizer":{__init__:[2,2,1,""],connect_two_molecules:[2,2,1,""],get_args:[2,2,1,""]},"restraintmaker.algorithm.Optimizer.TreeHeuristicOptimizer":{__init__:[2,2,1,""],connect_two_molecules:[2,2,1,""],get_args:[2,2,1,""]},"restraintmaker.algorithm.Optimizer.priority_node":{node:[2,2,1,""],priority:[2,2,1,""]},"restraintmaker.algorithm.Selection":{LimitedSelection:[2,1,1,""],MCS_Selection:[2,1,1,""],PaintSelection:[2,1,1,""],PairSelection:[2,1,1,""],SingleAtomSelection:[2,1,1,""],SphericalSelection:[2,1,1,""],UniversalSelection:[2,1,1,""]},"restraintmaker.algorithm.Selection.LimitedSelection":{__init__:[2,2,1,""],_update_select:[2,2,1,""],get_args:[2,2,1,""]},"restraintmaker.algorithm.Selection.MCS_Selection":{__init__:[2,2,1,""],_update_confirm:[2,2,1,""],_update_move:[2,2,1,""],_update_select:[2,2,1,""],_update_size:[2,2,1,""],get_args:[2,2,1,""],update:[2,2,1,""]},"restraintmaker.algorithm.Selection.PaintSelection":{__init__:[2,2,1,""],select_atoms_within_sphere:[2,2,1,""]},"restraintmaker.algorithm.Selection.PairSelection":{__init__:[2,2,1,""],get_args:[2,2,1,""]},"restraintmaker.algorithm.Selection.SingleAtomSelection":{__init__:[2,2,1,""],get_args:[2,2,1,""]},"restraintmaker.algorithm.Selection.SphericalSelection":{__init__:[2,2,1,""],_update_confirm:[2,2,1,""],_update_move:[2,2,1,""],_update_size:[2,2,1,""],atom_within_sphere:[2,2,1,""],get_args:[2,2,1,""],select_atoms_within_sphere:[2,2,1,""]},"restraintmaker.algorithm.Selection.UniversalSelection":{__init__:[2,2,1,""],_update_confirm:[2,2,1,""],_update_move:[2,2,1,""],_update_select:[2,2,1,""],_update_size:[2,2,1,""],get_args:[2,2,1,""],update:[2,2,1,""]},"restraintmaker.interface_Pymol":{RestraintMaker_Logic:[3,0,0,"-"],RestraintMaker_Pymol:[3,0,0,"-"],pymol_utilities:[4,0,0,"-"]},"restraintmaker.interface_Pymol.RestraintMaker_Logic":{Logic_Handler:[3,1,1,""]},"restraintmaker.interface_Pymol.RestraintMaker_Logic.Logic_Handler":{__init__:[3,2,1,""],_check_selection_status:[3,2,1,""],_set_exporter_type:[3,2,1,""],_set_filter_type:[3,2,1,""],_set_importer_type:[3,2,1,""],_set_optimizer_type:[3,2,1,""],_set_restraint_type:[3,2,1,""],_set_selection_type:[3,2,1,""],react_to_event:[3,2,1,""],set_action_states:[3,2,1,""],set_action_type:[3,2,1,""],set_mode:[3,2,1,""]},"restraintmaker.interface_Pymol.RestraintMaker_Pymol":{Restraints_Wizard:[3,1,1,""]},"restraintmaker.interface_Pymol.RestraintMaker_Pymol.Restraints_Wizard":{__init__:[3,2,1,""],_action_button_pressed:[3,2,1,""],_toggle_button_pressed:[3,2,1,""],check_objects_exists:[3,2,1,""],create_pymol_objects_for_restraints:[3,2,1,""],do_dirty:[3,2,1,""],do_key:[3,2,1,""],do_pick:[3,2,1,""],do_select:[3,2,1,""],do_special:[3,2,1,""],file_nr:[3,4,1,""],get_event_mask:[3,2,1,""],get_panel:[3,2,1,""],get_prompt:[3,2,1,""],provide_services_to_spherical_selection:[3,2,1,""],redraw:[3,2,1,""],set_wizard_style:[3,2,1,""],trigger_event:[3,2,1,""]},"restraintmaker.interface_Pymol.pymol_utilities":{program_states:[4,0,0,"-"],pymol_utitlities:[4,0,0,"-"],qt_dialogs:[4,0,0,"-"]},"restraintmaker.interface_Pymol.pymol_utilities.program_states":{ActionState:[4,1,1,""],EventType:[4,1,1,""]},"restraintmaker.interface_Pymol.pymol_utilities.program_states.ActionState":{ALWAYS_DISABLED:[4,4,1,""],ALWAYS_ENABLED:[4,4,1,""],CURRENTLY_DISABLED:[4,4,1,""],CURRENTLY_ENABLED:[4,4,1,""]},"restraintmaker.interface_Pymol.pymol_utilities.program_states.EventType":{CONFIRM:[4,4,1,""],MOVE:[4,4,1,""],SELECT:[4,4,1,""],SIZE:[4,4,1,""]},"restraintmaker.interface_Pymol.pymol_utilities.pymol_utitlities":{_atom_list_to_pymol_selection:[4,3,1,""],_convert_molecules_to_pdb:[4,3,1,""],colorize:[4,3,1,""],create_pymol_objects_for_molecules:[4,3,1,""],help_pymol_with_big_atom_list:[4,3,1,""],pymol_selection_to_atom_list:[4,3,1,""],try_to_get_args:[4,3,1,""]},"restraintmaker.interface_Pymol.pymol_utilities.qt_dialogs":{create_file_open_dialog:[4,3,1,""],create_file_save_dialog:[4,3,1,""],create_input_dialog:[4,3,1,""],create_multi_dialog:[4,3,1,""]},"restraintmaker.io":{Exporter:[5,0,0,"-"],Files:[6,0,0,"-"],Importer:[5,0,0,"-"]},"restraintmaker.io.Exporter":{Gromos_Distance_Restraint_Exporter:[5,1,1,""],_Exporter:[5,1,1,""]},"restraintmaker.io.Exporter.Gromos_Distance_Restraint_Exporter":{__init__:[5,2,1,""],export_restraints:[5,2,1,""],get_args:[5,2,1,""]},"restraintmaker.io.Exporter._Exporter":{__init__:[5,2,1,""],export_restraints:[5,2,1,""],get_args:[5,2,1,""]},"restraintmaker.io.Files":{Gromacs_files:[6,0,0,"-"],Gromos_blocks:[6,0,0,"-"],Gromos_files:[6,0,0,"-"]},"restraintmaker.io.Files.Gromos_blocks":{atomP:[6,1,1,""],atomV:[6,1,1,""],atom_pair_distanceRes:[6,1,1,""],atom_pos_block:[6,1,1,""],atom_ref_pos_block:[6,1,1,""],atom_vel_block:[6,1,1,""],distance_res_spec_block:[6,1,1,""],distant_Restraint_Type:[6,1,1,""],genbox_block:[6,1,1,""],geometric_code:[6,1,1,""],lattice_shift:[6,1,1,""],lattice_shifts_block:[6,1,1,""],timestep_block:[6,1,1,""],title_block:[6,1,1,""]},"restraintmaker.io.Files.Gromos_blocks.atomP":{__init__:[6,2,1,""],to_string:[6,2,1,""]},"restraintmaker.io.Files.Gromos_blocks.atomV":{__init__:[6,2,1,""],to_string:[6,2,1,""]},"restraintmaker.io.Files.Gromos_blocks.atom_pair_distanceRes":{__init__:[6,2,1,""],to_string:[6,2,1,""]},"restraintmaker.io.Files.Gromos_blocks.atom_pos_block":{__init__:[6,2,1,""],block_to_string:[6,2,1,""]},"restraintmaker.io.Files.Gromos_blocks.atom_ref_pos_block":{__init__:[6,2,1,""],block_to_string:[6,2,1,""]},"restraintmaker.io.Files.Gromos_blocks.atom_vel_block":{__init__:[6,2,1,""],block_to_string:[6,2,1,""]},"restraintmaker.io.Files.Gromos_blocks.distance_res_spec_block":{__init__:[6,2,1,""],block_to_string:[6,2,1,""]},"restraintmaker.io.Files.Gromos_blocks.distant_Restraint_Type":{full_harmonic_distance_restraint:[6,4,1,""],half_harmonic_attractive:[6,4,1,""],half_harmonic_repulsive:[6,4,1,""]},"restraintmaker.io.Files.Gromos_blocks.genbox_block":{__init__:[6,2,1,""],block_to_string:[6,2,1,""]},"restraintmaker.io.Files.Gromos_blocks.geometric_code":{pseudo_H_atom_goc_H_atoms_CH3:[6,4,1,""],pseudo_H_atoms_goc_of_three_CH3:[6,4,1,""],pseudo_H_atoms_goc_of_two_CH3:[6,4,1,""],real_atom:[6,4,1,""],virtual_H_atom_aliphaticC:[6,4,1,""],virtual_H_atom_aliphaticC_strange:[6,4,1,""],virtual_H_atom_aromaticC:[6,4,1,""],virtual_H_atoms_goc_aliph:[6,4,1,""],virtual_atoms_cog:[6,4,1,""],virtual_atoms_com:[6,4,1,""]},"restraintmaker.io.Files.Gromos_blocks.lattice_shift":{__init__:[6,2,1,""],to_string:[6,2,1,""]},"restraintmaker.io.Files.Gromos_blocks.lattice_shifts_block":{__init__:[6,2,1,""],block_to_string:[6,2,1,""]},"restraintmaker.io.Files.Gromos_blocks.timestep_block":{__init__:[6,2,1,""],block_to_string:[6,2,1,""]},"restraintmaker.io.Files.Gromos_blocks.title_block":{__init__:[6,2,1,""],block_to_string:[6,2,1,""],order:[6,4,1,""]},"restraintmaker.io.Files.Gromos_files":{disres:[6,1,1,""],general_gromos_file:[6,1,1,""]},"restraintmaker.io.Files.Gromos_files.disres":{__init__:[6,2,1,""],parse_disres_file:[6,2,1,""],read_disres_file:[6,2,1,""],required_blocks:[6,4,1,""]},"restraintmaker.io.Files.Gromos_files.general_gromos_file":{__init__:[6,2,1,""],_read_gromos_block:[6,2,1,""],add_block:[6,2,1,""],blocksset:[6,4,1,""],content:[6,4,1,""],path:[6,4,1,""],required_blocks:[6,4,1,""],write:[6,2,1,""]},"restraintmaker.io.Importer":{Gromos_Distance_Restraint_Importer:[5,1,1,""]},"restraintmaker.io.Importer.Gromos_Distance_Restraint_Importer":{__init__:[5,2,1,""],get_args:[5,2,1,""],import_restraints:[5,2,1,""]},"restraintmaker.restraintMaker":{_check_importing_packages:[1,3,1,""],run_plugin:[1,3,1,""],run_plugin_gui:[1,3,1,""]},"restraintmaker.test":{pickle_input_test_file:[7,0,0,"-"],test_exporter:[7,0,0,"-"],test_files:[8,0,0,"-"],test_importer:[7,0,0,"-"],test_optimizer:[7,0,0,"-"]},"restraintmaker.test.pickle_input_test_file":{test_file_dir:[7,5,1,""]},"restraintmaker.test.test_exporter":{test_Exporter:[7,1,1,""]},"restraintmaker.test.test_exporter.test_Exporter":{in_disresObj_path1:[7,4,1,""],in_disresObj_path2:[7,4,1,""],in_solution_disres_path1:[7,4,1,""],in_solution_disres_path2:[7,4,1,""],out_gromos_disres_path1:[7,4,1,""],out_gromos_disres_path2:[7,4,1,""],setUp:[7,2,1,""],test_Exporter_Gromos_construct:[7,2,1,""],test_Exporter_Gromos_export_disresDat:[7,2,1,""],test_Exporter_Gromos_getargs:[7,2,1,""],test_file_dir:[7,4,1,""],test_files_IO:[7,4,1,""]},"restraintmaker.test.test_importer":{test_Importer:[7,1,1,""]},"restraintmaker.test.test_importer.test_Importer":{in_disres1:[7,4,1,""],in_disres2:[7,4,1,""],in_pdb1:[7,4,1,""],in_pdb2:[7,4,1,""],setUp:[7,2,1,""],test_Importer_Gromos_construct:[7,2,1,""],test_Importer_Gromos_import_disresDat:[7,2,1,""],test_Importer_Gromos_import_file_notFound:[7,2,1,""],test_file_dir:[7,4,1,""],test_files_io:[7,4,1,""],test_files_structures:[7,4,1,""]},"restraintmaker.test.test_optimizer":{test_Optimizer:[7,1,1,""]},"restraintmaker.test.test_optimizer.test_Optimizer":{check_restraint_results:[7,2,1,""],test_bruteForce_optimizer:[7,2,1,""],test_cog_optimizer:[7,2,1,""],test_prim_optimizer:[7,2,1,""],test_shortest_optimizer:[7,2,1,""]},"restraintmaker.tools_Rdkit":{Rdkit_Functions:[9,0,0,"-"]},"restraintmaker.tools_Rdkit.Rdkit_Functions":{PolyArea:[9,3,1,""],PolygonSort:[9,3,1,""],_calc_pca_without_scaling:[9,3,1,""],mcs_selection:[9,3,1,""],parse_pdb_blocks_to_rdkit:[9,3,1,""],ring_atom_filter:[9,3,1,""]},"restraintmaker.utils":{Types:[10,0,0,"-"],Utilities:[10,0,0,"-"]},"restraintmaker.utils.Types":{CoM_Restraint:[10,1,1,""],Distance_Restraint:[10,1,1,""],Position_restraint:[10,1,1,""]},"restraintmaker.utils.Types.CoM_Restraint":{__init__:[10,2,1,""]},"restraintmaker.utils.Types.Distance_Restraint":{__init__:[10,2,1,""],atomA:[10,2,1,""],atomB:[10,2,1,""],atom_limit:[10,4,1,""],distance:[10,2,1,""]},"restraintmaker.utils.Types.Position_restraint":{__init__:[10,2,1,""],atomA:[10,2,1,""],atom_limit:[10,4,1,""],distance_to_reference_position:[10,2,1,""],reference_atom:[10,2,1,""]},"restraintmaker.utils.Utilities":{Atom:[10,1,1,""],BadArgumentException:[10,6,1,""],NoOptimalSolutionException:[10,6,1,""],RestraintPair:[10,1,1,""],check_for_doubles:[10,3,1,""],check_or_convert_argument:[10,3,1,""],check_restraint_pairs_for_doubles:[10,3,1,""],convert_atoms_to_pdb_molecules:[10,3,1,""],do_nothing:[10,3,1,""],execute_at_different_verbosity:[10,3,1,""],find_atom_by_property:[10,3,1,""],get_all_subclasses:[10,3,1,""],order_atoms_by_molecule:[10,3,1,""],print:[10,3,1,""]},"restraintmaker.utils.Utilities.Atom":{alt:[10,2,1,""],b:[10,2,1,""],chain:[10,2,1,""],elem:[10,2,1,""],id:[10,2,1,""],label:[10,2,1,""],name:[10,2,1,""],resi:[10,2,1,""],resn:[10,2,1,""],x:[10,2,1,""],y:[10,2,1,""],z:[10,2,1,""]},"restraintmaker.utils.Utilities.RestraintPair":{distance:[10,2,1,""],r1:[10,2,1,""],r2:[10,2,1,""]},restraintmaker:{algorithm:[2,0,0,"-"],interface_Pymol:[3,0,0,"-"],io:[5,0,0,"-"],restraintMaker:[1,0,0,"-"],test:[7,0,0,"-"],tools_Rdkit:[9,0,0,"-"],utils:[10,0,0,"-"],wizard:[1,0,0,"-"]}},objnames:{"0":["py","module","Python module"],"1":["py","class","Python class"],"2":["py","method","Python method"],"3":["py","function","Python function"],"4":["py","attribute","Python attribute"],"5":["py","data","Python data"],"6":["py","exception","Python exception"]},objtypes:{"0":"py:module","1":"py:class","2":"py:method","3":"py:function","4":"py:attribute","5":"py:data","6":"py:exception"},terms:{"00000":11,"00100":11,"00200":11,"00300":11,"00400":11,"004582659776816626":11,"00500":11,"00600":11,"00700":11,"01100":11,"01157615128508411":11,"01200008392334":11,"012530026346328164":11,"015166066564349288":11,"016552840275103946":11,"017720267741949906":11,"018999099731445":11,"019000053405762":11,"019105510292366224":11,"019261684801060098":11,"020000457763672":11,"02299976348877":11,"025376894134306997":11,"027000427246094":11,"027676766479476103":11,"02842549067573466":11,"0318273766301751":11,"033660026883969185":11,"033970555741711594":11,"03921695674451211":11,"041641362066013494":11,"04263797390290723":11,"04675501269373408":11,"05418529920277415":11,"057000160217285":11,"06499988115794252":11,"07078801180767609":11,"07400":11,"076000213623047":11,"07700":11,"077000141143799":11,"079000473022461":11,"0829999446868896":11,"0880000591278076":11,"1000":4,"100000381469727":11,"104":11,"105":11,"106":11,"107":11,"108":11,"109":11,"109999656677246":11,"110":11,"11236545326164007":11,"119000434875488":11,"122":11,"123":11,"123000144958496":11,"124":11,"132":11,"133":11,"134":11,"135":11,"136":11,"137":11,"138":11,"139":11,"141":11,"142":11,"142000198364258":11,"153":11,"1540000438690186":11,"1570000648498535":11,"177000045776367":11,"182000160217285":11,"189000129699707":11,"197000026702881":11,"197999954223633":11,"200":4,"202000141143799":11,"203":11,"20399999618530273":11,"206999778747559":11,"211999893188477":11,"215999603271484":11,"2160000056028366":11,"217000007629395":11,"220000267028809":11,"229000091552734":11,"23000000417232513":11,"235000133514404":11,"236000061035156":11,"239999771118164":11,"240999937057495":11,"243000030517578":11,"25600004196167":11,"257999897003174":11,"259000062942505":11,"272000312805176":11,"275000095367432":11,"2779998779296875":11,"2789998054504395":11,"281999588012695":11,"284000396728516":11,"284999847412109":11,"286999702453613":11,"309000015258789":11,"309999465942383":11,"310999870300293":11,"313999891281128":11,"319000244140625":11,"323999881744385":11,"324000358581543":11,"324999809265137":11,"327000141143799":11,"335000038146973":11,"338000297546387":11,"342000007629395":11,"3429999351501465":11,"34499979019165":11,"345000267028809":11,"3460001945495605":11,"347000122070312":11,"347999572753906":11,"354999542236328":11,"355999946594238":11,"359999656677246":11,"36400032043457":11,"371999740600586":11,"371999979019165":11,"373000144958496":11,"3799999952316284":11,"380000114440918":11,"387999534606934":11,"3880000114440918":11,"3889999389648438":11,"3919999599456787":11,"395999908447266":11,"39799976348877":11,"4019999504089355":11,"404999732971191":11,"4060001373291":11,"408999919891357":11,"410999774932861":11,"415999889373779":11,"4230000972747803":11,"430000066757202":11,"432000160217285":11,"4420000314712524":11,"442999839782715":11,"444000244140625":11,"448999404907227":11,"4529999494552612":11,"458999633789062":11,"461000442504883":11,"461999893188477":11,"463000297546387":11,"4670000076293945":11,"47800064086914":11,"4790000915527344":11,"4800000190734863":11,"480999946594238":11,"494999885559082":11,"501999855041504":11,"502000093460083":11,"505000114440918":11,"5069999694824219":11,"5099999904632568":11,"515999794006348":11,"520999908447266":11,"527999877929688":11,"529999732971191":11,"534000396728516":11,"5369999408721924":11,"539999961853027":11,"545000076293945":11,"545999526977539":11,"54699993133545":11,"548999786376953":11,"550000190734863":11,"550999641418457":11,"557999610900879":11,"562000274658203":11,"562999725341797":11,"563000202178955":11,"567000389099121":11,"570000171661377":11,"574000358581543":11,"574999809265137":11,"57800006866455":11,"578000068664551":11,"581000328063965":11,"590000152587891":11,"592000007629395":11,"593999862670898":11,"5960001945495605":11,"5969998836517334":11,"599999904632568":11,"607999801635742":11,"611000061035156":11,"61299991607666":11,"6299999952316284":11,"631999969482422":11,"63599967956543":11,"642000198364258":11,"642999649047852":11,"647000312805176":11,"6610000133514404":11,"666999816894531":11,"675000190734863":11,"703000068664551":11,"730999946594238":11,"734000205993652":11,"736000061035156":11,"73900032043457":11,"7421280231109166":11,"744999885559082":11,"745999813079834":11,"746000289916992":11,"746999740600586":11,"753999710083008":11,"765164079974097":11,"770999908447266":11,"781999588012695":11,"791000366210938":11,"811000108718872":11,"811999797821045":11,"824999809265137":11,"830999851226807":11,"836999893188477":11,"8529999852180481":11,"8569999933242798":11,"8619999885559082":11,"880000114440918":11,"883999824523926":11,"8889999985694885":11,"8899999856948853":11,"902999997138977":11,"906000137329102":11,"951000213623047":11,"9670000076293945":11,"9679999351501465":11,"9690001010894775":11,"973999977111816":11,"980999946594238":11,"9839999675750732":11,"984000205993652":11,"985000133514404":11,"986999988555908":11,"996000051498413":11,"998000144958496":11,"case":[2,3,7,10],"catch":2,"class":[2,3,4,5,6,7,10],"default":[4,10],"enum":[4,6],"export":[0,1,3,7],"float":[2,3,6],"function":[1,2,3,4,5,6,9,10],"import":[0,1,3,7,11],"int":[2,3,6,9,10],"new":[2,3],"return":[2,3,4,5,9,10],"switch":3,"throw":4,"true":[2,3,5,9],"try":[5,10],"while":2,ADDING:11,AND:5,BUT:2,But:[2,10],For:[2,5,13,14],IDS:4,MCS:[2,9],NOT:[2,3],Not:2,Ries:[13,14],THE:[3,4,5],The:[2,3,4,5,9,10],Then:2,These:[2,3],Use:2,WAS:6,Will:4,With:4,__class__:3,__init__:[2,3,5,6,10],__name__:3,_action_button_press:3,_atom_list_to_pymol_select:4,_calc_pca:9,_calc_pca_without_sc:9,_calculate_value_convex_hul:2,_calculate_value_scaled_pca_2d:2,_calculate_value_unscaled_pca_2d:2,_check_importing_packag:1,_check_selection_statu:3,_convert_molecules_to_pdb:4,_event:3,_event_funct:3,_export:[4,5],_filter:[2,4],_generic_field:6,_generic_gromos_block:6,_import:[3,4,5],_mode_button:3,_moleculeringoptim:2,_optimz:4,_read_gromos_block:6,_restraint:[2,3,4,5,10],_restrainttyp:2,_select:[2,3,4],_self:3,_set_exporter_typ:3,_set_filter_typ:3,_set_importer_typ:3,_set_optimizer_typ:3,_set_restraint_typ:3,_set_selection_typ:3,_set_x_typ:3,_toggle_button_press:3,_type:3,_update_:3,_update_confirm:2,_update_mov:2,_update_priority_:2,_update_s:2,_update_select:2,abl:3,abolish:4,about:[2,3],accept:[2,3,4,10],acceptable_valu:10,access:[2,3],accord:2,action:[3,4],action_st:3,actionst:4,activ:3,actual:[2,3],add:3,add_block:6,added:2,addit:[2,5],after:[2,3],again:2,algorithm:[0,1,11,13,14],alia:[2,10],all:[1,2,3,4,5,10,11],all_atom:[2,3,5],allerror:5,allow:[2,3,4],alreadi:[2,5],also:4,alt:[10,11],always_dis:4,always_en:4,anaconda3:3,anaconda:[13,14],analog:9,angl:6,ani:[2,3,4,5,6,9,10],apart:2,aplli:3,appli:[3,10],appropri:3,area:2,arg:[2,3,4,5,6,10],argument:[2,3,4,5,10],around:4,arrow:3,ask:4,assign:[2,5],atom:[2,3,4,5,6,9,10],atom_limit:10,atom_list:[4,11],atom_list_to_pymol_select:4,atom_pair_distancer:6,atom_pos_block:6,atom_ref_pos_block:6,atom_restraint:3,atom_vel_block:6,atom_within_spher:2,atoma:10,atomb:10,atomid:6,atomp:6,atomtyp:6,atomv:6,attribut:[2,3,5],autoclass:7,automat:3,avail:3,available_importer_typ:3,available_restraint_typ:3,avaiod:2,averag:2,back:3,backward:3,badargu:5,badargument_except:10,badargumentexcept:[2,3,5,10],base:[2,3,4,5,6,7,10],bash:[13,14],basic:3,becaus:[2,3,10],becom:2,been:[2,3],befor:[2,4,7,9],belong:2,benjamin:[13,14],best:2,bestmoleculesringoptim:2,better:[2,3,10],between:[2,3],big:2,biggest:2,bit:2,block:10,block_to_str:6,blocksset:6,blocktitl:6,bondflag:3,bool:[2,3,4,5,6,9,10],borrow:6,brd4:11,brd4_7lig:7,brute:2,bruteforceringoptimz:2,bschro:[3,7,11],build:2,bulk:1,button:[3,4],button_press:3,bymolecul:2,c10:11,c11:11,c12:11,c14:11,c16:11,c18:11,c20:11,c21:11,c22:11,c23:11,c24:11,c25:11,c26:11,c27:11,c28:11,c29:11,calcul:2,calculate_value_pca_relative_to_pca_of_both_molecul:2,call:[2,3,4,5,9,10],callabl:[2,4,5,10],can:[2,3,4,5,10,13,14],care:4,careful:2,carefulli:2,caus:3,center:4,certain:[2,3],chagn:[2,3,10],chain:[10,11],chang:[2,3,10],check:[1,2,3,4,6,7,9,10],check_for_doubl:10,check_if_onject_exist:3,check_objects_exist:3,check_or_convert_argu:10,check_restraint_pairs_for_doubl:10,check_restraint_result:7,chem:9,chk1_5lig:[7,11],choos:[2,4],chosen:2,clean:3,cleaner:2,clemen:[13,14],click:[2,3],clone:[13,14],close:4,closeest:2,closest:2,cmd:[3,11],code:[5,6,7,11],cog:2,cog_distance_restraint:2,color:4,colorize_atom:4,colour:[3,4],com:[10,13,14],com_restraint:10,combin:2,comment:6,commun:2,communicatopn:3,compar:2,compat:[2,5],complementari:4,complet:2,conda:[13,14],confirm:[2,3,4],confus:2,connect:[2,3,11],connect_two_molecul:2,connectivii:2,consid:[2,4],construct:[2,7],constructor:[4,6],contain:[2,3,4,5,10],content:0,control:3,conv:2,convers:[2,10],convert:[4,5,9,10],convert_atoms_to_pdb_molecul:[10,11],convert_molecules_to_pdb:4,convex:2,coord:9,coordin:[3,9],core:[13,14],corner:9,correct:[4,7],correspond:3,cost:2,could:[2,4,5],count:[3,10],cours:2,covari:2,creat:[2,3,4,13,14],create_file_dialog:4,create_file_open_dialog:4,create_file_save_dialog:4,create_input_dialog:4,create_multi_dialog:4,create_pymol_objects_for_molecul:4,create_pymol_objects_for_restraint:3,creation:2,criteria:[2,10],criterion:2,curent:3,current:[2,3,4],currently_dis:4,currently_en:4,cutoff:2,dat:[7,11],data:11,deal:[2,3,4],debbug:10,debug:[4,10],decid:3,decreas:2,defin:[2,3,4,10],definit:10,delet:3,deliv:5,depend:[1,2,3],desir:10,desired_typ:10,detect:3,develop:10,deviat:2,devtool:[13,14],dialog:4,dict:[3,4,5,6,9,10],dictionari:10,differ:[2,3,4],dim:9,dimens:2,diplai:2,directli:[2,3,5],dirti:2,dis:3,disabnl:3,discard:11,disinterest:10,displai:2,disr:[5,6,11],disres_text:11,distanc:[1,2,5,10],distance_res_spec_block:6,distance_restraint:[2,4,5,10,11],distance_to_reference_posit:10,distancerespec:6,distanceresspec:11,distant_restraint_typ:6,do_dirti:3,do_generate_restraint:3,do_kei:3,do_noth:10,do_pick:3,do_posit:3,do_select:3,do_speci:3,docu:6,document:[7,11],doe:[2,3,4,5,9,10],don:3,done:[2,9],doubl:[3,10],down:2,draw:3,dstanc:2,due:2,dummi:[2,10],dummy_str:10,dure:[2,10],duribng:10,each:[2,4,5],easier:3,edg:2,edge_list:2,eigenv:2,eigenvalu:2,either:4,elem:[10,11],element:[2,3,10,11],element_filt:2,elementfilt:[2,11],empti:[2,3],enabl:3,encapsul:3,end:[3,11],enmpti:4,enough:4,enter:3,entri:10,enumer:6,env:[3,13,14],enviorn:[13,14],environ:[13,14],environment_unix:[13,14],equal:2,error:[2,4,5,10],estim:1,etc:[3,4,5,9,10],euler:6,even:[2,10],event:[3,4],event_st:3,event_typ:3,eventtyp:[3,4],everi:[2,3,5],everyth:10,everz:2,exactli:[4,9],example_gui:[13,14],except:[2,3,10],exclud:4,excpet:[5,10],execut:[4,10,12,13,14],execute_at_different_verbos:10,exercis:7,exist:3,existst:3,expect:9,expected_restraint:7,explicitli:[2,3],export_restraint:[5,11],express:4,extra:2,fail:2,fals:[2,3,5,6,10],far:2,farthest:2,field:[2,10],field_seper:11,file:[1,4,5,7,11,13,14],file_nr:3,filter:[0,1,3,4,9,10],filtered_atom:11,find:[2,3,10],find_atom:10,find_atom_by_properti:10,finish:[3,4],finsih:2,fire:1,first:[1,2,6,10],fix:2,fixtur:7,fkt:4,flow:3,follow:3,fool:4,forbidden:2,forc:2,form:[2,3],format:[2,4,5,7,10],forward:10,found:2,found_restraint:[7,11],from:[1,2,3,4,5,6,9,10,11,13,14],fulfil:10,full:2,full_harmonic_distance_restraint:6,fulli:2,func:10,funciton:4,functiom:4,function_lib:11,functon:3,futur:2,genbox_block:6,gener:[2,3,11,13,14],general_gromos_fil:6,geometr:6,geometric_cod:6,get:[2,3,4,5,9,10,11],get_all_short_connect:2,get_all_short_pair_restraint:2,get_all_subclass:10,get_arg:[2,4,5,10,11],get_event_mask:3,get_panel:3,get_prior:2,get_prompt:3,git:[13,14],github:[13,14],give:[2,4,5],given:[1,2,3],going:9,grai:4,grandchild:10,graph:2,greedli:2,gromacs_fil:[1,5],gromo:[5,6,10],gromos_block:[1,5],gromos_distance_restraint_export:[5,11],gromos_distance_restraint_import:5,gromos_export:5,gromos_fil:[1,5],gromos_import:5,group:11,guarante:2,gui:[1,2,3],half_harmonic_attract:6,half_harmonic_repuls:6,handl:[2,3,4],handler:[3,4],happen:5,happpi:4,hardcod:2,has:[2,3,4,5],hate:4,have:[2,3,4,10],heapq:2,help:[3,9],help_pymol_with_big_atom_list:4,help_pymol_with_big_atoms_list:4,here:[1,2,5],heurist:2,heurstic:2,home:[3,7,11],hook:7,how:[2,10],http:[13,14],hull:2,idea:[4,5],ids:5,implement:[2,10],import_restraint:5,in_5ligs_disr:7,in_7ligs_disr:7,in_disres1:7,in_disres2:7,in_disresobj_path1:7,in_disresobj_path2:7,in_path:[5,6],in_pdb1:7,in_pdb2:7,in_solution_disres_path1:7,in_solution_disres_path2:7,includ:[2,4,10],include_private_sublcass:10,increas:2,ind:2,index:13,indic:[2,9,10],info:2,inherit:[2,3],initi:[2,3,4],inner:2,input:[2,4,5,7,10],input_funct:[2,4,5],insid:[2,3,9],instanc:[2,3],instead:[2,3],interconnect:1,interest:10,interfac:[5,13,14],interface_pymol:[0,1,5,11],interpret:3,introduc:3,intuit:3,involv:2,issu:3,its:[2,3,4,10],itself:2,join:11,just:[2,4,10],karg:10,kdisc:[6,11],kdish:[6,11],keep:[2,3],kei:[3,4,10],keyboard:3,keyword:[3,4],kind:[2,3],kl17:11,kl19:11,kl1:11,kl20:11,kl21:11,know:[3,4],kw_arg:[3,4],kwarg:10,label:[3,10,11],lambda:[2,10,11],last:[2,3],lattice_shift:6,lattice_shifts_block:6,launch:3,lead:2,least:[2,3],len:11,length:[4,6],lib:[3,11],like:2,limited_filt:2,limited_select:2,limitedselect:2,line:6,line_seper:11,link:9,linkeag:3,list:[2,3,4,5,6,9,10],load:11,logic:[2,3,4],logic_handl:3,look:[3,10],lost:2,lot:2,loud:9,made:2,main:1,make:[2,4,13,14],make_restraint:11,manag:[3,13,14],mani:[2,3,9],map:11,mark:2,match:3,max_atom:2,max_di:2,max_siz:2,maxim:2,maximal_restraint_dist:11,maximal_spanning_tree_greedi:2,maximal_spanning_tree_prim:2,maximal_weight_r:2,mcs_select:[2,9],mere:2,messag:[2,4,10],meta:2,metamoleculeringoptim:2,method:[2,3,4,7,11],methodnam:7,middl:2,might:[2,10],mimiz:2,min_mcs_siz:9,minim:[9,10],mod:3,mode:3,modul:[0,13],mol:9,molecul:[1,2,3,4,6,9,10,11],moment:3,more:[3,4,10],mous:3,move:[2,3,4],movement:3,multipl:[2,3],must:[3,5,10],mutabl:2,my_x:4,n13:11,n17:11,n_node:2,name:[3,4,10,11],necessari:[2,3,5,10],necessart:2,necessesari:2,need:[1,2,3,4,5,10],new_atom:2,new_exporter_typ:3,new_filter_typ:3,new_i:2,new_importer_typ:3,new_optimizer_typ:3,new_restraint_typ:3,new_selection_typ:3,new_verbosity_threshold:10,new_x:2,new_z:2,newli:2,nicer:10,node:2,noisi:9,non:[2,9],none:[3,6,7,10,11],nooptimalsolutionexcept:10,noreturn:[2,3,5,10],noth:10,notifi:2,notimplementederror:[2,3,5],now:[2,4],nrestraint:11,number:[2,4,10],numpi:1,o15:11,o24:11,o26:11,o27:11,obj:[3,9],object:[3,4,5,6],ofr:2,often:2,old:2,onc:2,one:[2,3,4,10],ones:2,onli:[2,3,4,5,9,10],opeat:2,open:4,oper:2,opt:10,optim:[0,1,3,4,7,10],optimz:3,optimzi:2,option:[2,4,6,10],order:[2,6,10],order_atoms_by_molecul:10,origin:6,other:[2,3],othr:2,out:[2,3,6,11],out_5ligs_disr:7,out_7ligs_disr:7,out_disr:11,out_gromos_disres_path1:7,out_gromos_disres_path2:7,out_path:[5,11],output:7,over:2,overkil:2,overload:2,overrid:2,overridden:[2,5,10],overriden:2,pack:3,packag:[0,13],page:13,paint:3,paintselect:2,pair:[2,10],pair_restraint:2,pairselect:2,pairwis:[2,9,10],pairwise_filt:2,param:4,paramet:[2,3,4,5,9,10],paramt:11,parent:[4,5,10],parent_window:4,parentn:10,parse_disres_fil:6,parse_pdb_blocks_to_rdkit:9,part:[1,2,4,13,14],pass:[2,3,4,10],path:[6,11],pbc:6,pca:2,pca_2d:11,pdb:[2,4,7,10,11],pdb_block:11,pdb_mol:9,pdb_path:11,pdblock:9,perfect:2,perform:3,pick:2,pickle_input_test_fil:[0,1],pipelin:12,place:2,placement:2,plu:3,plug:3,plugin:[1,13,14],point:1,polyarea:9,polygonsort:9,posit:[2,3],position_restraint:10,positon:3,possibl:[2,4,10],possilb:2,postion:2,potential_restraint:2,predefin:10,predefini:4,prefer:2,prepar:4,preselect:2,present:[2,4],press:[3,4],prevent:3,previou:10,prim:2,print:[5,10,11],prioriti:2,priority_algo:2,priority_nod:2,privat:[5,10],problem:3,produc:10,product:2,program:[2,3,10,11,13,14],program_st:[1,3],programm:1,progress:5,project:[10,11],proper:3,properti:[2,10],property_filt:2,property_nam:10,property_valu:10,propertyfilt:2,protein:4,provid:[2,3,4,5],provide_services_to_spherical_select:3,pseudo_h_atom_goc_h_atoms_ch3:6,pseudo_h_atoms_goc_of_three_ch3:6,pseudo_h_atoms_goc_of_two_ch3:6,pseudoatom:3,purpos:[13,14],pygromostool:6,pyhton:3,pymol:[1,3,4,5,11,13,14],pymol_funct:4,pymol_selection_to_atom_list:[4,11],pymol_util:[1,3,11],pymol_utitl:[1,3,11],pymool:4,pyqt5:4,python3:3,python:[1,3,11,13,14],pyumol:4,qt_dialog:[1,3],quantiti:2,quickli:2,quiet:4,radio:3,radiu:2,radius_0:6,rah:[6,11],rais:[2,3,4,5,10],rdchem:9,rdkit:[1,9],rdkit_funct:[0,1],react_to_:3,react_to_ev:3,read:[2,4,5],read_disres_fil:6,readabl:[2,10],real_atom:6,realis:[3,5],realli:[2,3],reason:[2,9],receiv:10,recogn:3,recommend:3,reconstruct:2,redraw:3,redund:2,reference_atom:10,refresh:[2,3],reiniti:11,rel:2,relev:3,remov:[2,9],report:3,repositori:[13,14],repres:2,requir:[1,2,4,5,10],required_arg:4,required_block:6,res:3,rescal:9,reset:10,resi:[10,11],resid:6,resn:[10,11],resnam:6,resolv:3,respons:3,restraint:[1,2,3,5,6,10,11,13,14],restraint_po:2,restraint_typ:6,restrainthead:6,restraintmak:12,restraintmaker_log:[0,1,4],restraintmaker_pymol:[0,1,2,4],restraintmakerd:3,restraintpair:[2,10],restraints_wizard:3,restrainttyp:[2,3],restraintyp:6,result:[9,11],retriev:[13,14],reus:2,rhiner:[13,14],ring:[2,9,11],ring_atom_filt:9,ringfilt:[2,11],ringoptim:2,rinik:11,rinikerlab:[13,14],rtype:[3,4],rule:10,run:1,run_plugin:1,run_plugin_gui:1,runtest:7,sacrif:2,same:2,samemolecul:2,satisfi:10,save:[3,5,9],scale:2,scene:3,scheme:4,scipi:1,search:[9,13],second:6,seem:9,selction:3,sele:[3,4],selecet:4,select:[0,1,3,4,5,9,10],select_atoms_within_spher:2,select_delet:3,selected_atom:3,selected_mol:9,selected_restraint:3,separ:4,seper:[3,4],set:[2,3,4,5,7,10,11],set_:3,set_act:3,set_action_st:3,set_action_typ:3,set_importer_typ:3,set_mod:3,set_restraint_typ:3,set_wizard_styl:3,set_x_typ:4,setup:[2,7],sever:4,shall:10,shift:3,shortcut:2,shorter:2,shortest:11,should:[2,3,4,5,9,10],show:2,similar:2,sinc:2,singl:2,singleatomselect:2,site:3,size:[2,3,4,9],slection:10,slow:2,small:4,solut:[2,4,10],solv:4,some:[2,3,4],somehow:4,sort:2,sourc:[1,2,3,4,5,6,7,9,10],span:2,special:3,specief:4,specifi:[3,4,10],spectrum:4,speed:2,speher:2,spher:2,sphere:[2,3],spheric:2,spherical_select:2,sphericalselect:[2,3],split:[13,14],squar:2,stabl:4,standalon:[13,14],standard:3,start:[1,2,3],state:[3,4],statement:10,step1:11,step2:11,step:6,stil:9,still:[2,10],stop:2,str:[2,3,4,5,6,9,10,11],string:[2,3,4,10],strongli:3,struc:3,structur:2,stuff:2,sub:4,subclass:[2,3,5,10],subcont:6,submethod:5,submodul:0,subpackag:0,suppport:3,system:11,systema:[7,11],systemb:7,tabl:3,tak:2,take:4,target:3,technic:3,tediouli:2,test:[0,1,2,4,11,13,14],test_bruteforce_optim:7,test_cog_optim:7,test_export:[0,1],test_exporter_gromos_construct:7,test_exporter_gromos_export_disresdat:7,test_exporter_gromos_getarg:7,test_fil:[1,7,11],test_file_dir:7,test_files_io:7,test_files_structur:7,test_import:[0,1],test_importer_gromos_construct:7,test_importer_gromos_import_disresdat:7,test_importer_gromos_import_file_notfound:7,test_optim:[0,1],test_prim_optim:7,test_restraintmaker_pymol:[0,1],test_shortest_optim:7,test_system1_pdb:11,test_system:[7,11],testcas:7,text:5,than:[2,3,4,10],thedialog:4,thefil:5,thei:[2,4],them:[1,2,4],thi:[1,2,3,4,5,6,7,9,10,13,14],thing:3,think:[2,3],thiswai:2,those:[2,4],though:2,three:2,threshold:10,through:3,thrown:3,time:[2,3,4,9],timestep_block:6,titl:[4,6,9,11],title_block:6,to_str:6,todo:[2,3,4,5,6,9,10],toggl:3,too:2,tools_rdkit:[0,1,2,13,14],total:2,touch:2,transfrom:9,translat:[3,5],tree:2,treeheuristicoptim:[2,11],trigger:[1,3],trigger_ev:3,try_to_get_arg:4,tst:3,tupl:[2,10],two:[2,13,14],twod:2,twomolecul:2,type1:6,type2:6,type:[0,1,2,3,4,5,9,11],typeerror:3,understand:[3,4],unittest:7,univers:2,universalselect:2,universi:2,updat:[2,3],upon:3,url:[13,14],use:[2,3],used:[2,3,10],useful:[2,4,10],user:[2,3,4,5,10],uses:2,using:[2,4,5],usr:[13,14],util:[0,1,2,3,4,5,11,13,14],valu:[2,3,4,6,10,11],value0:4,value_list:2,valueerror:[3,4,10],varaibl:5,variabl:[2,10],varianc:2,verbos:[5,6,9,10],verbosity_thresheold:10,verobs:10,version:2,via:[3,13,14],view:3,virtual_atoms_cog:6,virtual_atoms_com:6,virtual_h_atom_aliphaticc:6,virtual_h_atom_aliphaticc_strang:6,virtual_h_atom_aromaticc:6,virtual_h_atoms_goc_aliph:6,volum:2,vortex:2,wai:[3,4],walk:2,want:[2,3,4],ware:3,warn:[2,3,4,9,10],weight:[2,6],well:2,what:[2,3,4,9],whatev:3,wheel:3,when:[2,3,4],where:2,which:[2,3,4,5,10],who:[2,4],whole:[2,5],willbeassign:4,window:4,wip:[13,14],within:[2,4],withing:2,without:[3,4],wizard:[0,3,11],work:[2,3],would:2,wrap:[3,4],wrapper:[4,10],write:[3,4,5,6,11],wrote:11,x_name:3,x_type:3,yaml:[13,14],yet:[2,3,9],you:[3,4,13,14]},titles:["restraintmaker","restraintmaker package","restraintmaker.algorithm package","restraintmaker.interface_Pymol package","restraintmaker.interface_Pymol.pymol_utilities package","restraintmaker.io package","restraintmaker.io.Files package","restraintmaker.test package","restraintmaker.test.test_files package","restraintmaker.tools_Rdkit package","restraintmaker.utils package","Example Execution for RestraintMaker","Examples","Welcome to restraintMaker\u2019s documentation!","Welcome RestraintMaker"],titleterms:{"export":[5,11],"import":5,Using:[13,14],acknowledg:[13,14],algorithm:2,atom:11,author:[13,14],content:[1,2,3,4,5,6,7,8,9,10,13,14],develop:[13,14],document:13,exampl:[11,12,13,14],execut:11,file:6,filter:[2,11],gromacs_fil:6,gromos_block:6,gromos_fil:6,indic:13,instal:[13,14],interface_pymol:[3,4],introduct:[13,14],modul:[1,2,3,4,5,6,7,8,9,10],optim:[2,11],packag:[1,2,3,4,5,6,7,8,9,10],pickle_input_test_fil:7,pipelin:11,program_st:4,pymol_util:4,pymol_utitl:4,qt_dialog:4,rdkit_funct:9,restraintmak:[0,1,2,3,4,5,6,7,8,9,10,11,13,14],restraintmaker_log:3,restraintmaker_pymol:3,select:2,submodul:[1,2,3,4,5,6,7,9,10],subpackag:[1,3,5,7],tabl:13,test:[7,8],test_export:7,test_fil:8,test_import:7,test_optim:7,test_restraintmaker_pymol:7,tools_rdkit:9,type:10,util:10,welcom:[13,14],wizard:1}})
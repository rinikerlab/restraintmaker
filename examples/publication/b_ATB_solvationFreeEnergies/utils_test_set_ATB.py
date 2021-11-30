import pandas as pd
translate={
    647: {"orig_name":"F313", "short":"F31"},
    824: {"orig_name":"G078", "short":"G07"},
    984: {"orig_name":"G209", "short":"G20"},
    1045:  {"orig_name":"G277", "short":"G27"},
    1181: {"orig_name":"M097", "short":"M09"},
    1922: {"orig_name":"8018", "short":"801"},
    2905: {"orig_name":"S002", "short":"S00"},
    3348: {"orig_name":"M030", "short":"M03"},
    3596: {"orig_name":"M218", "short":"M21"},
    8021: {"orig_name":"_O6T", "short":"TO6"},
    8028: {"orig_name":"_O70", "short":"TO7"},
    8029: {"orig_name":"_O71", "short":"O71"},
    9378: {"orig_name":"_P8I", "short":"TP8"},
    26135: {"orig_name":"6J29", "short":"6J2"},
    30212: {"orig_name":"TVVS", "short":"TVV"},
    30253: {"orig_name":"E1VB", "short":"E1V"},
    30491: {"orig_name":"6KET", "short":"6KE"},
}

ATB_LIGANDS = {
    647: {"orig_name":"F313",  },
    824: {"orig_name":"G078",  },
    984: {"orig_name":"G209",  },
    1045:  {"orig_name":"G277",},
    1181: {"orig_name":"M097", },
    1922: {"orig_name":"8018", },
    2905: {"orig_name":"S002", },
    3348: {"orig_name":"M030", },
    3596: {"orig_name":"M218", },
    8021: {"orig_name":"_O6T", },
    8028: {"orig_name":"_O70", },
    8029: {"orig_name":"_O71", },
    9378: {"orig_name":"_P8I", },
    26135: {"orig_name":"6J29",},
    30212: {"orig_name":"TVVS",},
    30253: {"orig_name":"E1VB",},
    30491: {"orig_name":"6KET",},
}

#TI == ATB Absolute Hydr. FE
#
solvation_free_energy_kJ = pd.DataFrame({
    647: {"exp": -31.30, "TI": -30.38, "exp_std": 1.255, "TI_std": 0.44},
    824: {"exp": -11.80, "TI": -6.65, "exp_std": 1.255, "TI_std": 0.465},
    984: {"exp": -16.65, "TI": -9.98, "exp_std": 0, "TI_std": 0.43},
    1045: {"exp": -24.20, "TI": -17.18, "exp_std": 0, "TI_std":0.465},
    1181: {"exp": -20.54, "TI": -18.8, "exp_std": 1.255, "TI_std": 0.525},
    1922: {"exp": -14.40, "TI": -23.03, "exp_std": 0.200, "TI_std": 0.36},
    2905: {"exp": -9.96, "TI": -15.14, "exp_std": 0.0420, "TI_std": 0.615},
    3348: {"exp": -3.77, "TI": -0.2, "exp_std": 1.255, "TI_std": 0.34},
    3596: {"exp": -19.62, "TI": -25.14, "exp_std": 1.255, "TI_std": 0.45},
    8021: {"exp": -22.30, "TI": -26.42, "exp_std": 0.210, "TI_std": 0.74},
    8028: {"exp": -15.70, "TI": -14.63, "exp_std": 0.450, "TI_std": 0.58},
    8029: {"exp": -18.60, "TI": -22.58, "exp_std": 0.900, "TI_std": 0.615},
    9378: {"exp": -19.66, "TI": -15.54, "exp_std": 1.255, "TI_std": 0.455},
    26135: {"exp": -39.90, "TI": -48.43, "exp_std": 0.600, "TI_std": 0.285},
    30212: {"exp": -29.29, "TI": -26.33, "exp_std": 1.255, "TI_std": 0.44},
    30253: {"exp": -5.40, "TI": -5.4, "exp_std": 1.255, "TI_std": 0.605},
    30491: {"exp": -32.05, "TI": -36.29, "exp_std": 1.255, "TI_std": 0.54},
}).T

multistate_ligand_sets = {"all": ['_O6T', '_O71', 'G277', 'S002', '8018', 'M030', '6KET', 'F313', '_P8I', 'M097', 'G078', 'M218', '6J29', 'G209', 'E1VB', 'TVVS', '_O70'],
                          "singles": ['_O6T', '_O71', 'G277', 'S002', 'M030', '6KET', 'F313', 'M097', 'M218', 'TVVS', ],
                          "flat": ['_O6T', '_O71', 'S002', 'M030', '6KET', 'F313', 'M097', 'G078', 'M218', 'TVVS'],
                          "easy": ['_O6T', 'G277', "M030", 'M097', '6KET', 'F313']
                          }
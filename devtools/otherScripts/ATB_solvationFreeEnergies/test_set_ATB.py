import pandas as pd

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
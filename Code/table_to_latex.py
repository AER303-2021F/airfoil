

import pandas as pd
import numpy as np


# -----------------------------------------------------------------------------
# Manometer airfoil pressures.
# -----------------------------------------------------------------------------


data = np.loadtxt("./Data/calibrated_manometer_airfoil.csv", delimiter=",")
columns = data[:, 0]
data = data[:, 1:].T
data = np.round(data, decimals=1)
data = data.astype("<U24")

# Import uncertainties.
error = np.loadtxt("./Data/calibrated_manometer_airfoil_error.csv", delimiter=",")
error = error[:, 1:].T
error = np.round(error, decimals=1)
error = error.astype("<U24")

combined_data = np.char.add(np.char.add(np.char.add(np.char.add("$", data), "\pm"), error), "$")
data_df = pd.DataFrame(combined_data, index=range(1, 20), columns=columns)
print(data_df.to_latex(escape=False))


# -----------------------------------------------------------------------------
# Manometer wake pressures a.
# -----------------------------------------------------------------------------


# Import data.
data = np.loadtxt("./Data/calibrated_manometer_rake_a.csv", delimiter=",", skiprows=1)
columns = data[:, 0]
data = data[:, 1:].T
data = np.round(data, decimals=1)
data = data.astype("<U24")

# Import uncertainties.
error = np.loadtxt("./Data/calibrated_manometer_rake_a_uncertainty.csv", delimiter=",", skiprows=1)
error = error[:, 1:].T
error = np.round(error, decimals=1)
error = error.astype("<U24")

combined_data = np.char.add(np.char.add(np.char.add(np.char.add("$", data), "\pm"), error), "$")
data_df = pd.DataFrame(combined_data, index=range(1, 18), columns=columns)
print(data_df.to_latex(escape=False))


# -----------------------------------------------------------------------------
# Manometer wake pressures a.
# -----------------------------------------------------------------------------


# Import data.
data = np.loadtxt("./Data/calibrated_manometer_rake_b.csv", delimiter=",", skiprows=0)
columns = data[:, 0]
data = data[:, 1:].T
data = np.round(data, decimals=1)
data = data.astype("<U24")

# Import uncertainties.
error = np.loadtxt("./Data/calibrated_manometer_rake_b_uncertainty.csv", delimiter=",", skiprows=0)
error = error[:, 1:].T
error = np.round(error, decimals=1)
error = error.astype("<U24")

combined_data = np.char.add(np.char.add(np.char.add(np.char.add("$", data), "\pm"), error), "$")
data_df = pd.DataFrame(combined_data, index=range(1, 18), columns=columns)
print(data_df.to_latex(escape=False))


# -----------------------------------------------------------------------------
# Scanivalve airfoil pressures.
# -----------------------------------------------------------------------------


data = np.loadtxt("./Data/scanivalve_airfoil_pressure_data.csv", delimiter=",", skiprows=1)
columns = data[:, 0]
data = data[:, 1:].T
data = np.round(data, decimals=2)
data = data.astype("<U24")

# Import uncertainties.
error = np.loadtxt("./Data/scanivalve_error.csv", delimiter=",") * 115
print(error)
error = np.round(error, decimals=2)
error = error.astype("<U24")

combined_data = np.char.add(np.char.add(np.char.add(np.char.add("$", data), "\pm"), error), "$")
data_df = pd.DataFrame(combined_data, index=range(1, 20), columns=columns)
print(data_df.to_latex(escape=False))


# -----------------------------------------------------------------------------
# Scanivalve wake pressures a.
# -----------------------------------------------------------------------------


# Import data.
data = np.loadtxt("./Data/calibrated_scanivalve_rake_a.csv", delimiter=",", skiprows=1)
columns = data[:, 0]
data = data[:, 1:].T
data = np.round(data, decimals=2)
data = data.astype("<U24")

# Import uncertainties.
error = np.loadtxt("./Data/scanivalve_error.csv", delimiter=",") * 115
error = np.round(error, decimals=2)
error = error.astype("<U24")

combined_data = np.char.add(np.char.add(np.char.add(np.char.add("$", data), "\pm"), error), "$")
data_df = pd.DataFrame(combined_data, index=range(1, 18), columns=columns)
print(data_df.to_latex(escape=False))


# -----------------------------------------------------------------------------
# Scanivalve wake pressures a.
# -----------------------------------------------------------------------------


# Import data.
data = np.loadtxt("./Data/calibrated_scanivalve_rake_b.csv", delimiter=",", skiprows=1)
columns = data[:, 0]
data = data[:, 1:].T
data = np.round(data, decimals=2)
data = data.astype("<U24")

# Import uncertainties.
error = np.loadtxt("./Data/scanivalve_error.csv", delimiter=",") * 115
error = np.round(error, decimals=2)
error = error.astype("<U24")

combined_data = np.char.add(np.char.add(np.char.add(np.char.add("$", data), "\pm"), error), "$")
data_df = pd.DataFrame(combined_data, index=range(1, 18), columns=columns)
print(data_df.to_latex(escape=False))


# -----------------------------------------------------------------------------
# Create drag latex table.
# -----------------------------------------------------------------------------


mano_total_drag = np.loadtxt("code/wake_analysis/total_drag_manometer.csv", delimiter=",")
mano_total_error_drag = np.loadtxt("code/wake_analysis/total_drag_error_manometer.csv", delimiter=",")
mano_pressure_drag = np.loadtxt("Data/mano_pressure_drag.csv", delimiter=",")
mano_pressure_drag_error = np.loadtxt("Data/mano_pressure_drag_error.csv", delimiter=",")

scani_total_drag = np.loadtxt("code/wake_analysis/total_drag_scanivalve.csv", delimiter=",")
scani_total_error_drag = np.loadtxt("code/wake_analysis/total_drag_error_scanivalve.csv", delimiter=",")
scani_pressure_drag = np.loadtxt("Data/scani_pressure_drag.csv", delimiter=",")
scani_pressure_drag_error = np.loadtxt("Data/scani_pressure_drag_error.csv", delimiter=",")

mano_total_drag = mano_total_drag.reshape((-1, 1))
mano_total_error_drag = mano_total_error_drag.reshape((-1, 1))
mano_pressure_drag = mano_pressure_drag.reshape((-1, 1))
mano_pressure_drag_error = mano_pressure_drag_error.reshape((-1, 1))
scani_total_drag = scani_total_drag.reshape((-1, 1))
scani_total_error_drag = scani_total_error_drag.reshape((-1, 1))
scani_pressure_drag = scani_pressure_drag.reshape((-1, 1))
scani_pressure_drag_error = scani_pressure_drag_error.reshape((-1, 1))

aoa = np.array([0, 3, 6, 8, 10, 11, 13, 15, 16, 17, 20]).reshape((-1, 1))
aoa_error = np.repeat(0.25, aoa.shape[0]).reshape((-1, 1))

combined_data = np.concatenate((aoa, scani_pressure_drag, mano_pressure_drag, scani_total_drag, mano_total_drag), axis=1)
combined_error_data = np.concatenate((aoa_error, scani_pressure_drag_error, mano_pressure_drag_error, scani_total_error_drag, mano_total_error_drag), axis=1)

combined_data = np.round(combined_data, decimals=3)
combined_data = combined_data.astype("<U24")

combined_error_data = np.round(combined_error_data, decimals=3)
combined_error_data = combined_error_data.astype("<U24")

combined_data = np.char.add(np.char.add(np.char.add(np.char.add("$", combined_data), "\pm"), combined_error_data), "$")
print(combined_data.shape)
columns=["$\alpha^\circ$", "$D'_p$ (Scanivalve)", "$D'_p$ (Manometer)", "$D'$ (Scanivalve)", "$D'$ (Manometer)"]
data_df = pd.DataFrame(combined_data, index=range(1, 12), columns=columns)
print(data_df.to_latex(escape=False, index=False))


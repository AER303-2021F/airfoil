

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

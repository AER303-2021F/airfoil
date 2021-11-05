

import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import pandas as pd


# -----------------------------------------------------------------------------
# Define program functions.
# -----------------------------------------------------------------------------


AOA_UNCERTAINTY = 0.25  # Degrees.


# -----------------------------------------------------------------------------
# Define program functions.
# -----------------------------------------------------------------------------


def airfoil_theta(xl, yl, xu, yu):
    """Return the angle between the pressure and airfoil normals.

    Parameters
    ----------
    xl : np.ndarray
        Array of shape (n,) containing tap x positions on the lower wing
        surface.
    yl : np.ndarray
        Array of shape (n,) containing tap y positions on the lower wing
        surface.
    xu : np.ndarray
        Array of shape (n,) containing tap x positions on the upper wing
        surface.
    yu : np.ndarray
        Array of shape (n,) containing tap y positions on the upper wing
        surface.

    Returns
    -------
    np.ndarray, np.ndarray
        A tuple of arrays, each of shape (n,), containing the angle between the
        pressure and airfoil normals in radians.
    """

    num = - np.diff(yl)
    denom = np.diff(xl)
    theta_l = np.arctan(num / denom)

    num = - np.diff(yu)
    denom = np.diff(xu)
    theta_u = np.arctan(num / denom)

    return theta_l, theta_u


def normal_force(xl, yl, xu, yu, pl, pu, thetal, thetau):
    """Return normal force of the airfoil.

    Parameters
    ----------
    xl : np.ndarray
        Array of shape (n,) containing tap x positions on the lower wing
        surface.
    yl : np.ndarray
        Array of shape (n,) containing tap y positions on the lower wing
        surface.
    xu : np.ndarray
        Array of shape (n,) containing tap x positions on the upper wing
        surface.
    yu : np.ndarray
        Array of shape (n,) containing tap y positions on the upper wing
        surface.
    pl : [type]
        Array of shape (n,) containing the pressure at the lower wing taps.
    pu : [type]
        Array of shape (n,) containing the pressure at the upper wing taps.
    thetal : [type]
        Array of shape (n,) containing the angle between the pressure
        and lower airfoil surface normals.
    thetau : [type]
        Array of shape (n,) containing the angle between the pressure
        and upper airfoil surface normals.

    Returns
    -------
    float
        The normal force of the airfoil.
    """

    ds_lower = np.sqrt(np.diff(xl) ** 2 + np.diff(yl) ** 2)
    ds_upper = np.sqrt(np.diff(xu) ** 2 + np.diff(yu) ** 2)

    N = 0

    for n in range(len(ds_upper)):
        N -= 0.5 * (pu[n] + pu[n + 1]) * np.cos(thetau[n]) * ds_upper[n]
    
    for n in range(len(ds_lower)):
        N += 0.5 * (pl[n] + pl[n + 1]) * np.cos(thetal[n]) * ds_lower[n]

    return N


def axial_force(xl, yl, xu, yu, pl, pu, thetal, thetau):
    """Return axial force of the airfoil.

    Parameters
    ----------
    xl : np.ndarray
        Array of shape (n,) containing tap x positions on the lower wing
        surface.
    yl : np.ndarray
        Array of shape (n,) containing tap y positions on the lower wing
        surface.
    xu : np.ndarray
        Array of shape (n,) containing tap x positions on the upper wing
        surface.
    yu : np.ndarray
        Array of shape (n,) containing tap y positions on the upper wing
        surface.
    pl : [type]
        Array of shape (n,) containing the pressure at the lower wing taps.
    pu : [type]
        Array of shape (n,) containing the pressure at the upper wing taps.
    thetal : [type]
        Array of shape (n,) containing the angle between the pressure
        and lower airfoil surface normals.
    thetau : [type]
        Array of shape (n,) containing the angle between the pressure
        and upper airfoil surface normals.

    Returns
    -------
    float
        The axial force in Newtons.
    """

    ds_lower = np.sqrt(np.diff(xl) ** 2 + np.diff(yl) ** 2)
    ds_upper = np.sqrt(np.diff(xu) ** 2 + np.diff(yu) ** 2)

    A = 0

    for n in range(len(ds_upper)):
        A -= 0.5 * (pu[n] + pu[n + 1]) * np.sin(thetau[n]) * ds_upper[n]
    
    for n in range(len(ds_lower)):
        A += 0.5 * (pl[n] + pl[n + 1]) * np.sin(thetal[n]) * ds_lower[n]

    return A


def moment(xl, yl, xu, yu, pl, pu, thetal, thetau):
    """Return the moment of the airfoil in Newtons.

    Parameters
    ----------
    xl : np.ndarray
        Array of shape (n,) containing tap x positions on the lower wing
        surface.
    yl : np.ndarray
        Array of shape (n,) containing tap y positions on the lower wing
        surface.
    xu : np.ndarray
        Array of shape (n,) containing tap x positions on the upper wing
        surface.
    yu : np.ndarray
        Array of shape (n,) containing tap y positions on the upper wing
        surface.
    pl : [type]
        Array of shape (n,) containing the pressure at the lower wing taps.
    pu : [type]
        Array of shape (n,) containing the pressure at the upper wing taps.
    thetal : [type]
        Array of shape (n,) containing the angle between the pressure
        and lower airfoil surface normals.
    thetau : [type]
        Array of shape (n,) containing the angle between the pressure
        and upper airfoil surface normals.

    Returns
    -------
    float
        The moment of the airfoil in Newtons.
    """

    ds_lower = np.sqrt(np.diff(xl) ** 2 + np.diff(yl) ** 2)
    ds_upper = np.sqrt(np.diff(xu) ** 2 + np.diff(yu) ** 2)

    M = 0

    for n in range(len(ds_upper)):
        M += 0.5 * (pu[n] + pu[n + 1]) * np.cos(thetau[n]) * xu[n] * ds_upper[n]
        M -= 0.5 * (pu[n] + pu[n + 1]) * np.sin(thetau[n]) * yu[n] * ds_upper[n]
    
    for n in range(len(ds_lower)):
        M += 0.5 * (pl[n] + pl[n + 1]) * np.cos(thetal[n]) * xl[n] * ds_lower[n]
        M += 0.5 * (pl[n] + pl[n + 1]) * np.sin(thetal[n]) * yl[n] * ds_lower[n]

    return M


def lift(N, A, aoa):
    """Return the lift force in Newtons.

    Parameters
    ----------
    N : float
        The airfoils normal force.
    A : float
        The airfoils axial force.
    aoa : int
        The angle of attack of the airfoil.

    Returns
    -------
    float
        The lift force in Newtons.
    """

    aoa = np.radians(aoa)
    return N * np.cos(aoa) - A * np.sin(aoa)


def pressure_drag(N, A, aoa):
    """Return the pressure drag force in Newtons.

    Parameters
    ----------
    N : float
        The airfoils normal force.
    A : float
        The airfoils axial force.
    aoa : int
        The angle of attack of the airfoil.

    Returns
    -------
    float
        The pressure drag force in Newtons.
    """
    
    aoa = np.radians(aoa)
    return N * np.sin(aoa) + A * np.cos(aoa)


def coefficients(L, D, M, velocity):
    """Return the lift, drag, and moment coefficients.

    Parameters
    ----------
    L : float
        The lift force in Newtons.
    D : float
        The drag force in Newtons.
    M : float
        The moment force in Newtons.
    velocity : np.ndarray
        Array of shape (n,) containing the free stream velocity at a given
        anngle of attack.

    Returns
    -------
    float, float, float
        A tuple containing the lift coefficient, drag coefficient, and
        moment coefficient.
    """

    q_infty = 0.5 * 1.225 * velocity ** 2

    C_L = L / q_infty
    C_D = D / (q_infty * 0.1)
    C_M = M / q_infty

    return C_L, C_D, C_M


# -----------------------------------------------------------------------------
# Import free stream velocity data.
# -----------------------------------------------------------------------------


AoA = np.array([0, 3, 6, 8, 10, 11, 13, 15, 16, 17, 20])

mano_velocity_data = np.loadtxt("./Data/velocities_manometer.csv", delimiter=",")
mano_velocity = mano_velocity_data[0, :]
mano_velocity_uncertainty = mano_velocity_data[1, :]

scani_velocity_data = np.loadtxt("./Data/velocities_scanivalve.csv", delimiter=",")
scani_velocity = scani_velocity_data[0, :]
scani_velocity_uncertainty = scani_velocity_data[1, :]


# -----------------------------------------------------------------------------
# Import and interpolate airfoil data at tap locations.
# -----------------------------------------------------------------------------


airfoil_data = np.loadtxt("./XFOIL Analysis/clark_y.dat", skiprows=1)
n = airfoil_data.shape[0] // 2

airfoil_x_upper = airfoil_data[:n + 1, 0]
airfoil_x_lower = airfoil_data[n:, 0]
airfoil_y_upper = airfoil_data[:n + 1, 1]
airfoil_y_lower = airfoil_data[n:, 1]

locations_x_upper = [0, 0.03, 0.06, 0.10, 0.15, 0.20, 0.30, 0.40, 0.55, 0.70, 0.85, 1.00]
locations_x_lower = [0.05, 0.10, 0.20, 0.30, 0.40, 0.60, 0.90]
locations_y_upper = []
locations_y_lower = []

for loc in locations_x_upper:
    locations_y_upper.append(interp1d(airfoil_x_upper, airfoil_y_upper, kind='cubic')(loc))

for loc in locations_x_lower:
    locations_y_lower.append(interp1d(airfoil_x_lower, airfoil_y_lower, kind='cubic')(loc))

# plt.plot(airfoil_x_upper, airfoil_y_upper, label="Actual Upper")
# plt.scatter(locations_x_upper, locations_y_upper, label="Interpolated Upper")
# plt.plot(airfoil_x_lower, airfoil_y_lower, label="Actual Lower")
# plt.scatter(locations_x_lower, locations_y_lower, label="Interpolated Lower")
# plt.title("Airfoil Interpolation")
# plt.legend()
# plt.grid()
# plt.show()


# -----------------------------------------------------------------------------
# Import and segregate the upper and lower pressure measurements.
# -----------------------------------------------------------------------------


scanivalve_data_file = "./Data/scanivalve_airfoil_pressure_data.csv"
manometer_data_file = "./Data/calibrated_manometer_airfoil.csv"

pressure_data_s = np.loadtxt(scanivalve_data_file, delimiter=",", skiprows=1)
pressure_upper_s = pressure_data_s[:, 1:13]
pressure_lower_s = np.flip(pressure_data_s[:, 13:], axis=1)

pressure_data_m = np.loadtxt(manometer_data_file, delimiter=",", skiprows=0)
pressure_upper_m = pressure_data_m[:, 1:13]
pressure_lower_m = np.flip(pressure_data_m[:, 13:], axis=1)

AoA = pressure_data_s[:, 0].reshape((-1,))


# -----------------------------------------------------------------------------
# Import and segregate the upper and lower pressure uncertainties.
# -----------------------------------------------------------------------------


scanivalve_error_file = "./Data/scanivalve_error.csv"
manometer_data_file = "./Data/calibrated_manometer_airfoil_error.csv"

pressure_error_s = np.loadtxt(scanivalve_error_file, delimiter=",") * 115
pressure_error_m = np.loadtxt(manometer_data_file, delimiter=",", skiprows=0)
pressure_error_m = pressure_error_m[:, 1:]


# -----------------------------------------------------------------------------
# Get theta angles.
# -----------------------------------------------------------------------------


thetal, thetau = airfoil_theta(locations_x_lower, locations_y_lower, locations_x_upper, locations_y_upper)


# -----------------------------------------------------------------------------
# Get coefficients for scanivalve data.
# -----------------------------------------------------------------------------


print("SCANIVALVE DATA")

L_vs_aoa_s = np.zeros(AoA.size)
D_vs_aoa_s = np.zeros(AoA.size)
M_vs_aoa_s = np.zeros(AoA.size)
cL_vs_aoa_s = np.zeros(AoA.size)
cD_vs_aoa_s = np.zeros(AoA.size)
cM_vs_aoa_s = np.zeros(AoA.size)

delta_N = np.zeros(AoA.size)
delta_A = np.zeros(AoA.size)
delta_M = np.zeros(AoA.size)

for i in range(AoA.size):

    aoa = int(AoA[i])

    print(f"Calculations for AOA={aoa}")

    params = [
        locations_x_lower,
        locations_y_lower,
        locations_x_upper,
        locations_y_upper,
        pressure_lower_s[i, :],
        pressure_upper_s[i, :],
        thetal,
        thetau
    ]

    # Calculate values.
    N, A, M = normal_force(*params), axial_force(*params), moment(*params)
    L, D = lift(N, A, aoa), pressure_drag(N, A, aoa)
    cL, cD, cM = coefficients(L, D/100, M, scani_velocity[i])

    L_vs_aoa_s[i], D_vs_aoa_s[i], M_vs_aoa_s[i] = L, D, M
    cL_vs_aoa_s[i], cD_vs_aoa_s[i], cM_vs_aoa_s[i] = cL, cD, cM

    print(f"Lift = {L}, Drag = {D}, Moment = {M}")
    print(f"C_L = {cL}, C_D = {cD}, C_M = {cM}")

    # Calculate uncertainties (assuming delta_rho and delta_c are zero).

    delta_N_term_1 = np.sum((-np.cos(thetau) * np.sqrt(np.diff(locations_x_upper) ** 2 + np.diff(locations_y_upper) ** 2) * pressure_error_s[i]) ** 2)
    delta_N_term_2 = np.sum((np.cos(thetal) * np.sqrt(np.diff(locations_x_lower) ** 2 + np.diff(locations_y_lower) ** 2) * pressure_error_s[i]) ** 2)
    delta_N[i] = np.sqrt(delta_N_term_1 + delta_N_term_2)

    delta_A_term_1 = np.sum((-np.sin(thetau) * np.sqrt(np.diff(locations_x_upper) ** 2 + np.diff(locations_y_upper) ** 2) * pressure_error_s[i]) ** 2)
    delta_A_term_2 = np.sum((np.cos(thetal) * np.sqrt(np.diff(locations_x_lower) ** 2 + np.diff(locations_y_lower) ** 2) * pressure_error_s[i]) ** 2)
    delta_A[i] = np.sqrt(delta_A_term_1 + delta_A_term_2)

    delta_M_term_1 = np.sum(((np.cos(thetau) * locations_x_upper[1:] - np.sin(thetau) * locations_y_upper[1:]) * np.sqrt(np.diff(locations_x_upper) ** 2 + np.diff(locations_y_upper) ** 2) * pressure_error_s[i]) ** 2)
    delta_M_term_2 = np.sum(((np.cos(thetal) * locations_x_lower[1:] + np.sin(thetal) * locations_y_lower[1:]) * np.sqrt(np.diff(locations_x_lower) ** 2 + np.diff(locations_y_lower) ** 2) * pressure_error_s[i]) ** 2)
    delta_M[i] = np.sqrt(delta_M_term_1 + delta_M_term_2)

delta_aoa = AOA_UNCERTAINTY

delta_L = np.sqrt((np.cos(AoA) * delta_N) ** 2 + (np.sin(AoA) * delta_A) ** 2 + ((N * np.sin(AoA) + A * np.cos(AoA)) * delta_aoa) ** 2)
cL_uncertainty_s = abs(cL_vs_aoa_s) * np.sqrt((delta_L / L_vs_aoa_s) ** 2 + (2 * scani_velocity_uncertainty / scani_velocity) ** 2)

delta_D = np.sqrt((np.sin(AoA) * delta_N) ** 2 + (np.cos(AoA) * delta_A) ** 2 + ((N * np.cos(AoA) - A * np.sin(AoA)) * delta_aoa) ** 2)
cD_uncertainty_s = abs(cD_vs_aoa_s) * np.sqrt((delta_D / D_vs_aoa_s) ** 2 + (2 * scani_velocity_uncertainty / scani_velocity) ** 2)

np.savetxt("./Data/scani_pressure_drag.csv", D_vs_aoa_s / 100, delimiter=",")
np.savetxt("./Data/scani_pressure_drag_error.csv", delta_D / 100, delimiter=",")

cM_uncertainty_s = abs(cM_vs_aoa_s) * np.sqrt((delta_M / M_vs_aoa_s) ** 2 + (2 * scani_velocity_uncertainty / scani_velocity) ** 2)


# -----------------------------------------------------------------------------
# Get coefficients for manometer data.
# -----------------------------------------------------------------------------

    
print("MANOMETER DATA")

L_vs_aoa_m = np.zeros(AoA.size)
D_vs_aoa_m = np.zeros(AoA.size)
M_vs_aoa_m = np.zeros(AoA.size)
cL_vs_aoa_m = np.zeros(AoA.size)
cD_vs_aoa_m = np.zeros(AoA.size)
cM_vs_aoa_m = np.zeros(AoA.size)

delta_N = np.zeros(AoA.size)
delta_A = np.zeros(AoA.size)
delta_M = np.zeros(AoA.size)

for i in range(AoA.size):

    aoa = int(AoA[i])

    print(f"Calculations for AOA={aoa}")

    params = [
        locations_x_lower,
        locations_y_lower,
        locations_x_upper,
        locations_y_upper,
        pressure_lower_m[i, :],
        pressure_upper_m[i, :],
        thetal,
        thetau
    ]

    N, A, M = normal_force(*params), axial_force(*params), moment(*params)
    L, D = lift(N, A, aoa), pressure_drag(N, A, aoa)
    cL, cD, cM = coefficients(L, D / 100, M, mano_velocity[i])

    L_vs_aoa_m[i], D_vs_aoa_m[i], M_vs_aoa_m[i] = L, D, M
    cL_vs_aoa_m[i], cD_vs_aoa_m[i], cM_vs_aoa_m[i] = cL, cD, cM

    print(f"Lift = {L}, Drag = {D}, Moment = {M}")
    print(f"C_L = {cL}, C_D = {cD}, C_M = {cM}")

    # Calculate uncertainties (assuming delta_rho and delta_c are zero).

    delta_N_term_1 = np.sum((-np.cos(thetau) * np.sqrt(np.diff(locations_x_upper) ** 2 + np.diff(locations_y_upper) ** 2) * pressure_error_m[i, :11]) ** 2)
    delta_N_term_2 = np.sum((np.cos(thetal) * np.sqrt(np.diff(locations_x_lower) ** 2 + np.diff(locations_y_lower) ** 2) * pressure_error_m[i, 13:]) ** 2)
    delta_N[i] = np.sqrt(delta_N_term_1 + delta_N_term_2)

    delta_A_term_1 = np.sum((-np.sin(thetau) * np.sqrt(np.diff(locations_x_upper) ** 2 + np.diff(locations_y_upper) ** 2) * pressure_error_m[i, :11]) ** 2)
    delta_A_term_2 = np.sum((np.cos(thetal) * np.sqrt(np.diff(locations_x_lower) ** 2 + np.diff(locations_y_lower) ** 2) * pressure_error_m[i, 13:]) ** 2)
    delta_A[i] = np.sqrt(delta_A_term_1 + delta_A_term_2)

    delta_M_term_1 = np.sum(((np.cos(thetau) * locations_x_upper[1:] - np.sin(thetau) * locations_y_upper[1:]) * np.sqrt(np.diff(locations_x_upper) ** 2 + np.diff(locations_y_upper) ** 2) * pressure_error_s[i]) ** 2)
    delta_M_term_2 = np.sum(((np.cos(thetal) * locations_x_lower[1:] + np.sin(thetal) * locations_y_lower[1:]) * np.sqrt(np.diff(locations_x_lower) ** 2 + np.diff(locations_y_lower) ** 2) * pressure_error_s[i]) ** 2)
    delta_M[i] = np.sqrt(delta_M_term_1 + delta_M_term_2)

delta_aoa = AOA_UNCERTAINTY

delta_L = np.sqrt((np.cos(AoA) * delta_N) ** 2 + (np.sin(AoA) * delta_A) ** 2 + ((N * np.sin(AoA) + A * np.cos(AoA)) * delta_aoa) ** 2)
cL_uncertainty_m = abs(cL_vs_aoa_s) * np.sqrt((delta_L / L_vs_aoa_s) ** 2 + (2 * scani_velocity_uncertainty / scani_velocity) ** 2)

delta_D = np.sqrt((np.sin(AoA) * delta_N) ** 2 + (np.cos(AoA) * delta_A) ** 2 + ((N * np.cos(AoA) - A * np.sin(AoA)) * delta_aoa) ** 2)
cD_uncertainty_m = abs(cD_vs_aoa_s) * np.sqrt((delta_D / D_vs_aoa_s) ** 2 + (2 * scani_velocity_uncertainty / scani_velocity) ** 2)

np.savetxt("./Data/mano_pressure_drag.csv", D_vs_aoa_m / 100, delimiter=",")
np.savetxt("./Data/mano_pressure_drag_error.csv", delta_D / 100, delimiter=",")

cM_uncertainty_m = abs(cM_vs_aoa_s) * np.sqrt((delta_M / M_vs_aoa_s) ** 2 + (2 * scani_velocity_uncertainty / scani_velocity) ** 2)


# -----------------------------------------------------------------------------
# Plot scanivalve vs validation data.
# -----------------------------------------------------------------------------


val_data = np.loadtxt("./Data/validation_coefficients.txt", skiprows=12, usecols=(0, 1, 2, 4))
plt.plot(val_data[:, 0], val_data[:, 1], label="cL xfoil", color="red", linestyle="--", linewidth=1)
plt.plot(val_data[:, 0], val_data[:, 2], label="cD xfoil", color="blue", linestyle="--", linewidth=1)
plt.plot(val_data[:, 0], val_data[:, 3], label="cM xfoil", color="green", linestyle="--", linewidth=1)

plt.plot(AoA, cL_vs_aoa_s, label="cL, scani", color="red", linewidth=1)
plt.errorbar(AoA, cL_vs_aoa_s, yerr=cL_uncertainty_s, fmt="none", capsize=2, color="red", linewidth=1)
plt.plot(AoA, cD_vs_aoa_s, label="cD, scani", color="blue", linewidth=1)
plt.errorbar(AoA, cD_vs_aoa_s, yerr=cD_uncertainty_s, fmt="none", capsize=2, color="blue", linewidth=1)
plt.plot(AoA, cM_vs_aoa_s, label="cM, scani", color="green", linewidth=1)
plt.errorbar(AoA, cM_vs_aoa_s, yerr=cM_uncertainty_s, fmt="none", capsize=2, color="green", linewidth=1)

plt.title("$C_L$, $C_D$, and $C_M$ vs $\\alpha$, Scanivalve Data")
plt.xlabel("$\\alpha^\circ$")
plt.ylabel("Coefficient")
plt.legend()
plt.grid()
plt.show()


# -----------------------------------------------------------------------------
# Plot manometer vs validation data.
# -----------------------------------------------------------------------------


val_data = np.loadtxt("./Data/validation_coefficients.txt", skiprows=12, usecols=(0, 1, 3, 4))
plt.plot(val_data[:, 0], val_data[:, 1], label="cL xfoil", color="red", linestyle="--", linewidth=1)
plt.plot(val_data[:, 0], val_data[:, 2], label="cD xfoil", color="blue", linestyle="--", linewidth=1)
plt.plot(val_data[:, 0], val_data[:, 3], label="cM xfoil", color="green", linestyle="--", linewidth=1)

plt.plot(AoA, cL_vs_aoa_m, label="cL, mano", color="red", linewidth=1)
plt.errorbar(AoA, cL_vs_aoa_m, yerr=cL_uncertainty_m, fmt="none", capsize=2, color="red", linewidth=1)
plt.plot(AoA, cD_vs_aoa_m, label="cD, mano", color="blue", linewidth=1)
plt.errorbar(AoA, cD_vs_aoa_m, yerr=cD_uncertainty_m, fmt="none", capsize=2, color="blue", linewidth=1)
plt.plot(AoA, cM_vs_aoa_m, label="cM, mano", color="green", linewidth=1)
plt.errorbar(AoA, cM_vs_aoa_m, yerr=cM_uncertainty_m, fmt="none", capsize=2, color="green", linewidth=1)

plt.title("$C_L$, $C_D$, and $C_M$ vs $\\alpha$, Manometer Data")
plt.xlabel("$\\alpha^\circ$")
plt.ylabel("Coefficient")
plt.legend()
plt.grid()
plt.show()


# -----------------------------------------------------------------------------
# Save scanivalve and manometer coefficients and uncertainties to csv.
# -----------------------------------------------------------------------------


AoA = AoA.reshape((AoA.size, 1))

cL_vs_aoa_s = cL_vs_aoa_s.reshape((cL_vs_aoa_s.size, 1))
cD_vs_aoa_s = cD_vs_aoa_s.reshape((cD_vs_aoa_s.size, 1))
cM_vs_aoa_s = cM_vs_aoa_s.reshape((cM_vs_aoa_s.size, 1))
scani_coeff_data = np.concatenate((AoA, cL_vs_aoa_s, cD_vs_aoa_s, cM_vs_aoa_s), axis=1)

np.savetxt("./Data/scanivalve_coefficients.csv", scani_coeff_data, delimiter=",", header="AoA, cL, cDp, cM")

cL_uncertainty_s = cL_uncertainty_s.reshape((cL_uncertainty_s.size, 1))
cD_uncertainty_s = cD_uncertainty_s.reshape((cD_uncertainty_s.size, 1))
cM_uncertainty_s = cM_uncertainty_s.reshape((cM_uncertainty_s.size, 1))
scani_uncertainty_data = np.concatenate((np.repeat(0.25, AoA.size).reshape((-1, 1)), cL_uncertainty_s, cD_uncertainty_s, cM_uncertainty_s), axis=1)

np.savetxt("./Data/scanivalve_coefficient_uncertainty.csv", scani_uncertainty_data, delimiter=",", header="AoA, cL Uncertainty, cDp Uncertainty, cM Uncertainty")

scani_coeff_data = np.round(scani_coeff_data, decimals=7).astype("<U24")
scani_uncertainty_data = np.round(scani_uncertainty_data, decimals=7).astype("<U24")
scani_out = np.char.add(np.char.add(np.char.add(np.char.add("$", scani_coeff_data), "\pm"), scani_uncertainty_data), "$")
columns = ["$\alpha^\circ$", "$C_L$", "$C_{D,p}$", "$C_M$"]
df = pd.DataFrame(scani_out, columns=columns)
print(df.to_latex(escape=False, index=False))

cL_vs_aoa_m = cL_vs_aoa_m.reshape((cL_vs_aoa_m.size, 1))
cD_vs_aoa_m = cD_vs_aoa_m.reshape((cD_vs_aoa_m.size, 1))
cM_vs_aoa_m = cM_vs_aoa_m.reshape((cM_vs_aoa_m.size, 1))
mano_coeff_data = np.concatenate((AoA, cL_vs_aoa_m, cD_vs_aoa_m, cM_vs_aoa_m), axis=1)

np.savetxt("./Data/manometer_coefficients.csv", mano_coeff_data, delimiter=",", header="AoA, cL, cDp, cM")

cL_uncertainty_m = cL_uncertainty_m.reshape((cL_uncertainty_m.size, 1))
cD_uncertainty_m = cD_uncertainty_m.reshape((cD_uncertainty_m.size, 1))
cM_uncertainty_m = cM_uncertainty_m.reshape((cM_uncertainty_m.size, 1))
mano_uncertainty_data = np.concatenate((np.repeat(0.25, AoA.size).reshape((-1, 1)), cL_uncertainty_m, cD_uncertainty_m, cM_uncertainty_m), axis=1)

np.savetxt("./Data/manometer_coefficient_uncertainty.csv", mano_uncertainty_data, delimiter=",", header="AoA, cL Uncertainty, cDp Uncertainty, cM Uncertainty")

mano_coeff_data = np.round(mano_coeff_data, decimals=7).astype("<U24")
mano_uncertainty_data = np.round(mano_uncertainty_data, decimals=7).astype("<U24")
mano_out = np.char.add(np.char.add(np.char.add(np.char.add("$", mano_coeff_data), "\pm"), mano_uncertainty_data), "$")
columns = ["$\alpha^\circ$", "$C_L$", "$C_{D,p}$", "$C_M$"]
df = pd.DataFrame(mano_out, columns=columns)
print(df.to_latex(escape=False, index=False))

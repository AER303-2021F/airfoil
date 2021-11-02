

import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from typing import Tuple


# -----------------------------------------------------------------------------
# Interpolate airfoil data at tap locations.
# -----------------------------------------------------------------------------


airfoil_data = np.loadtxt("./XFOIL Analysis/clark_y.dat", skiprows=1)
n = airfoil_data.shape[0] // 2

airfoil_x_upper = airfoil_data[:n + 1, 0]
airfoil_x_lower = airfoil_data[n:, 0]
airfoil_y_upper = airfoil_data[:n + 1, 1]
airfoil_y_lower = airfoil_data[n:, 1]

locations = [0, 0.03, 0.06, 0.10, 0.15, 0.20, 0.30, 0.40, 0.55, 0.70, 0.85, 1.00, 0.90, 0.60, 0.40, 0.30, 0.20, 0.10, 0.05]

locations_x_upper = [0, 0.03, 0.06, 0.10, 0.15, 0.20, 0.30, 0.40, 0.55, 0.70, 0.85, 1.00]
locations_x_lower = [0.90, 0.60, 0.40, 0.30, 0.20, 0.10, 0.05]
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
# Segregate the upper and lower pressure measurements.
# -----------------------------------------------------------------------------


pressure_data = np.loadtxt("Data/pressure_data_csv.csv", delimiter=",", skiprows=1)
AoA = pressure_data[:, 0].reshape((-1,))

pressure_upper = pressure_data[:, 1:13]
pressure_lower = pressure_data[:, 13:]


# -----------------------------------------------------------------------------
# Define functions.
# -----------------------------------------------------------------------------


def airfoil_theta(xl, yl, xu, yu):

    num = - np.diff(yl)
    denom = np.diff(xl)
    theta_l = np.arctan(num / denom)

    num = np.diff(yu)
    denom = np.diff(xu)
    theta_u = np.arctan(num / denom)

    return theta_l, theta_u

def normal_force(xl, yl, xu, yu, pl, pu, thetal, thetau):

    ds_upper = np.sqrt(np.diff(xl) ** 2 + np.diff(yl) ** 2)
    ds_lower = np.sqrt(np.diff(xu) ** 2 + np.diff(yu) ** 2)

    N = 0

    for n in range(len(ds_upper) - 1):
        N -= pu[n] * np.cos(thetau[n]) * ds_upper[n]
    
    for n in range(len(ds_lower) - 1):
        N += pl[n] * np.cos(thetal[n]) * ds_lower[n]

    return N


def axial_force(xl, yl, xu, yu, pl, pu, thetal, thetau):

    ds_upper = np.sqrt(np.diff(xl) ** 2 + np.diff(yl) ** 2)
    ds_lower = np.sqrt(np.diff(xu) ** 2 + np.diff(yu) ** 2)

    A = 0

    for n in range(len(ds_upper) - 1):
        A -= pu[n] * np.sin(thetau[n]) * ds_upper[n]
    
    for n in range(len(ds_lower) - 1):
        A += pl[n] * np.sin(thetal[n]) * ds_lower[n]

    return A


def moment(xl, yl, xu, yu, pl, pu, thetal, thetau):

    ds_upper = np.sqrt(np.diff(xl) ** 2 + np.diff(yl) ** 2)
    ds_lower = np.sqrt(np.diff(xu) ** 2 + np.diff(yu) ** 2)

    M = 0

    for n in range(len(ds_upper) - 1):
        M += pu[n] * np.cos(thetau[n]) * ds_upper[n] * xu[n]
        M -= pu[n] * np.sin(thetau[n]) * ds_upper[n] * yu[n]
    
    for n in range(len(ds_lower) - 1):
        M += pl[n] * np.sin(thetal[n]) * ds_lower[n]

        M -= pl[n] * np.cos(thetal[n]) * ds_lower[n] * xl[n]
        M += pl[n] * np.sin(thetal[n]) * ds_lower[n] * yl[n]

    return M


def lift(N, A, aoa):
    
    return N * np.cos(aoa) - A * np.sin(aoa)


def pressure_drag(N, A, aoa):
    
    return N * np.sin(aoa) + A * np.cos(aoa)


def coefficients(L, D, M):
    # Assuming chord equals one.

    q_infty = 0.5 * 1.225 * 29.214 ** 2

    C_L = L / q_infty
    C_D = D / q_infty
    C_M = M / q_infty

    return C_L, C_D, C_M


def main():

    L_vs_aoa = np.zeros(AoA.size)
    D_vs_aoa = np.zeros(AoA.size)
    M_vs_aoa = np.zeros(AoA.size)
    cL_vs_aoa = np.zeros(AoA.size)
    cD_vs_aoa = np.zeros(AoA.size)
    cM_vs_aoa = np.zeros(AoA.size)

    thetal, thetau = airfoil_theta(locations_x_lower, locations_y_lower, locations_x_upper, locations_y_upper)

    for i in range(AoA.size):

        aoa = int(AoA[i])

        print(f"Calculations for AOA={aoa}")

        N = normal_force(locations_x_lower, locations_y_lower, locations_x_upper, locations_y_upper, pressure_lower[aoa, :], pressure_upper[aoa, :], thetal, thetau)
        A = axial_force(locations_x_lower, locations_y_lower, locations_x_upper, locations_y_upper, pressure_lower[aoa, :], pressure_upper[aoa, :], thetal, thetau)
        M = moment(locations_x_lower, locations_y_lower, locations_x_upper, locations_y_upper, pressure_lower[aoa, :], pressure_upper[aoa, :], thetal, thetau)

        L = lift(N, A, AoA[i])
        D = pressure_drag(N, A, aoa)

        cL, cD, cM = coefficients(L, D, M)

        L_vs_aoa[i] = L
        D_vs_aoa[i] = D
        M_vs_aoa[i] = M
        cL_vs_aoa[i] = cL
        cD_vs_aoa[i] = cD
        cM_vs_aoa[i] = cM

        print(f"Lift = {L}, Drag = {D}, Moment = {M}")
        print(f"C_L = {cL}, C_D = {cD}, C_M = {cM}")
    
    plt.plot(AoA, L_vs_aoa, label="Lift")
    plt.plot(AoA, D_vs_aoa, label="Drag")
    plt.plot(AoA, M_vs_aoa, label="Moment")
    plt.legend()
    plt.grid()
    plt.show()

    plt.plot(AoA, cL_vs_aoa, label="cL")
    plt.plot(AoA, cD_vs_aoa, label="cD")
    plt.plot(AoA, cM_vs_aoa, label="cM")
    plt.legend()
    plt.grid()
    plt.show()


if __name__ == "__main__":
    main()

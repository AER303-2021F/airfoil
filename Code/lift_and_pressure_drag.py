

import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from typing import Tuple


locations = [0, 0.03, 0.06, 0.10, 0.15, 0.20, 0.30, 0.40, 0.55, 0.70, 0.85, 1.00, 0.90, 0.60, 0.40, 0.30, 0.20, 0.10, 0.05]
LOCATION_TOP = 12


def airfoil_theta(airfoil_data: np.ndarray, locations: np.ndarray) -> float:
    """Return the angle between the chord line's normal and the aifoil normal.

    Parameters
    ----------
    airfoil_data : np.ndarray
        Array of shape (n, 2) containing the x and y coordinates of the airfoil.

    Returns
    -------
    float
        Array of shape (n,) containing the angle between the chord line's
        normal and the aifoil normal.
    """

    # plt.plot(airfoil_data[:, 0], airfoil_data[:, 1])
    num = np.diff(airfoil_data[:, 1])
    denom = np.diff(airfoil_data[:, 0])
    theta = np.arctan(num / denom)
    # plt.plot(airfoil_data[1:, 0], theta)

    min_airfoil_x_val = np.min(airfoil_data[:, 0])
    min_airfoil_x = np.where(airfoil_data[:, 0] == min_airfoil_x_val)[0][0] - 1

    up_x, up_y = airfoil_data[:min_airfoil_x + 2, 0], theta[:min_airfoil_x + 2]
    low_x, low_y = airfoil_data[min_airfoil_x + 1: -1, 0], theta[min_airfoil_x + 1:]
    upper_interp = interp1d(up_x, up_y, kind='cubic')
    lower_interp = interp1d(low_x, low_y, kind='cubic')

    up_interp = upper_interp(locations[:LOCATION_TOP])
    low_interp = lower_interp(locations[LOCATION_TOP:])
    theta = np.concatenate((up_interp, low_interp))

    # plt.scatter(locations, theta)
    # plt.grid()
    # plt.show()

    return theta


def normal_force(pressure_data: np.ndarray, airfoil_data: np.ndarray, locations: np.ndarray) -> float:
    
    theta = airfoil_theta(airfoil_data, locations)

    vertical_pressure = pressure_data * np.cos(theta)

    sign = np.ones(len(vertical_pressure))
    sign[:len(vertical_pressure) // 2] = -1
    signed_vertical_pressure = sign * vertical_pressure

    normal_force = np.trapz(signed_vertical_pressure, locations, dx=0.01)

    return normal_force


def axial_force(pressure_data: np.ndarray, airfoil_data: np.ndarray, locations: np.ndarray) -> float:
    
    theta = airfoil_theta(airfoil_data, locations)

    horizontal_pressure = pressure_data * np.sin(theta)
    signed_horizontal_pressure = np.sign(theta) * horizontal_pressure

    axial_force = np.trapz(signed_horizontal_pressure, locations, dx=0.01)

    return axial_force


def moment(pressure_data: np.ndarray, airfoil_data: np.ndarray, locations: np.ndarray) -> float:
    
    theta = airfoil_theta(airfoil_data, locations)

    vertical_pressure = pressure_data * np.cos(theta)
    signed_vertical_pressure = np.sign(theta) * vertical_pressure

    horizontal_pressure = pressure_data * np.sin(theta)
    signed_horizontal_pressure = np.sign(theta) * horizontal_pressure

    moment = np.trapz(signed_vertical_pressure, locations, dx=0.01)

    return moment


def lift(N: float, A: float, aoa: float) -> float:
    
    return N * np.cos(aoa) - A * np.sin(aoa)


def pressure_drag(N: float, A: float, aoa: float) -> float:
    
    return N * np.sin(aoa) + A * np.cos(aoa)


def coefficients(L: float, D: float, M: float) -> Tuple[float, float, float]:
    # Assuming chord equals one.

    q_infty = 0.5 * 1.225 * 29.214 ** 2

    C_L = L / q_infty
    C_D = D / q_infty
    C_M = M / q_infty

    return C_L, C_D, C_M


def main():

    airfoil_data = np.loadtxt("./XFOIL Analysis/clark_y.dat", skiprows=1)
    
    pressure_data = np.loadtxt("Data/pressure_data_csv.csv", delimiter=",", skiprows=1)
    AoA = pressure_data[:, 0].reshape((-1,))

    L_vs_aoa = np.zeros(AoA.size)
    D_vs_aoa = np.zeros(AoA.size)
    M_vs_aoa = np.zeros(AoA.size)
    cL_vs_aoa = np.zeros(AoA.size)
    cD_vs_aoa = np.zeros(AoA.size)
    cM_vs_aoa = np.zeros(AoA.size)

    for i in range(AoA.size):

        print(f"Calculations for AOA={AoA[i]}")

        N = normal_force(pressure_data[i, 1:], airfoil_data, locations)
        A = axial_force(pressure_data[i, 1:], airfoil_data, locations)
        M = moment(pressure_data[i, 1:], airfoil_data, locations)

        L = lift(N, A, AoA[i])
        D = pressure_drag(N, A, AoA[i])

        cL, cD, cM = coefficients(L, D, 1)

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

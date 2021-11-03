

import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt


# -----------------------------------------------------------------------------
# Constants.
# -----------------------------------------------------------------------------


FREE_STREAM_VELOCITY = {
    0: 30.1224,
    3: 30.2637,
    6: 30.4689,
    8: 30.7127,
    10: 30.6284,
    11: 30.7055,
    13: 30.7579,
    15: 30.9922,
    16: 31.1288,
    17: 31.2287,
    20: 31.4522
}


# -----------------------------------------------------------------------------
# Interpolate airfoil data at tap locations.
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
# Segregate the upper and lower pressure measurements.
# -----------------------------------------------------------------------------


scanivalve_data_file = "./Data/pressure_data_csv.csv"
manometer_data_file = "./Data/calibrated_manometer_airfoil.csv"

scani = False

if scani:
    data_file = scanivalve_data_file
    skiprows = 1
    title = "$C_L$, $C_D$, and $C_M$ vs $\\alpha$, Scanivalve Data"
else:
    data_file = manometer_data_file
    skiprows = 0
    title = "$C_L$, $C_D$, and $C_M$ vs $\\alpha$, Manometer Data"

pressure_data_s = np.loadtxt(scanivalve_data_file, delimiter=",", skiprows=1)
pressure_upper_s = pressure_data_s[:, 1:13]
pressure_lower_s = np.flip(pressure_data_s[:, 13:], axis=1)

pressure_data_m = np.loadtxt(manometer_data_file, delimiter=",", skiprows=0)
pressure_upper_m = pressure_data_m[:, 1:13]
pressure_lower_m = np.flip(pressure_data_m[:, 13:], axis=1)

AoA = pressure_data_s[:, 0].reshape((-1,))


# -----------------------------------------------------------------------------
# Define functions.
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

    # for n in range(len(ds_upper)):
    #     M += 0.5 * (pu[n] + pu[n + 1]) * np.cos(thetau[n]) * xu[n] * ds_upper[n]
    #     M -= 0.5 * (pu[n] + pu[n + 1]) * np.sin(thetau[n]) * yu[n] * ds_upper[n]
    
    # for n in range(len(ds_lower)):
    #     M -= 0.5 * (pl[n] + pl[n + 1]) * np.cos(thetal[n]) * xl[n] * ds_lower[n]
    #     M += 0.5 * (pl[n] + pl[n + 1]) * np.sin(thetal[n]) * yl[n] * ds_lower[n]

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


def coefficients(L, D, M, aoa):
    """Return the lift, drag, and moment coefficients.

    Parameters
    ----------
    L : float
        The lift force in Newtons.
    D : float
        The drag force in Newtons.
    M : float
        The moment force in Newtons.
    aoa : int
        The angle of attack of the airfoil.

    Returns
    -------
    float, float, float
        A tuple containing the lift coefficient, drag coefficient, and
        moment coefficient.
    """

    q_infty = 0.5 * 1.225 * FREE_STREAM_VELOCITY[aoa] ** 2

    C_L = L / q_infty
    C_D = D / q_infty
    C_M = M / q_infty

    return C_L, C_D, C_M


def main():

    thetal, thetau = airfoil_theta(locations_x_lower, locations_y_lower, locations_x_upper, locations_y_upper)

    print("SCANIVALVE DATA")

    L_vs_aoa_s = np.zeros(AoA.size)
    D_vs_aoa_s = np.zeros(AoA.size)
    M_vs_aoa_s = np.zeros(AoA.size)
    cL_vs_aoa_s = np.zeros(AoA.size)
    cD_vs_aoa_s = np.zeros(AoA.size)
    cM_vs_aoa_s = np.zeros(AoA.size)

    for i in range(AoA.size):

        aoa = int(AoA[i])

        print(f"Calculations for AOA={aoa}")

        N = normal_force(locations_x_lower, locations_y_lower, locations_x_upper, locations_y_upper, pressure_lower_s[i, :], pressure_upper_s[i, :], thetal, thetau)
        A = axial_force(locations_x_lower, locations_y_lower, locations_x_upper, locations_y_upper, pressure_lower_s[i, :], pressure_upper_s[i, :], thetal, thetau)
        M = moment(locations_x_lower, locations_y_lower, locations_x_upper, locations_y_upper, pressure_lower_s[i, :], pressure_upper_s[i, :], thetal, thetau)

        L = lift(N, A, aoa)
        D = pressure_drag(N, A, aoa)

        cL, cD, cM = coefficients(L, D, M, aoa)

        L_vs_aoa_s[i] = L
        D_vs_aoa_s[i] = D
        M_vs_aoa_s[i] = M
        cL_vs_aoa_s[i] = cL
        cD_vs_aoa_s[i] = cD
        cM_vs_aoa_s[i] = cM

        print(f"Lift = {L}, Drag = {D}, Moment = {M}")
        print(f"C_L = {cL}, C_D = {cD}, C_M = {cM}")
    
    print("MANOMETER DATA")

    L_vs_aoa_m = np.zeros(AoA.size)
    D_vs_aoa_m = np.zeros(AoA.size)
    M_vs_aoa_m = np.zeros(AoA.size)
    cL_vs_aoa_m = np.zeros(AoA.size)
    cD_vs_aoa_m = np.zeros(AoA.size)
    cM_vs_aoa_m = np.zeros(AoA.size)

    for i in range(AoA.size):

        aoa = int(AoA[i])

        print(f"Calculations for AOA={aoa}")

        N = normal_force(locations_x_lower, locations_y_lower, locations_x_upper, locations_y_upper, pressure_lower_m[i, :], pressure_upper_m[i, :], thetal, thetau)
        A = axial_force(locations_x_lower, locations_y_lower, locations_x_upper, locations_y_upper, pressure_lower_m[i, :], pressure_upper_m[i, :], thetal, thetau)
        M = moment(locations_x_lower, locations_y_lower, locations_x_upper, locations_y_upper, pressure_lower_m[i, :], pressure_upper_m[i, :], thetal, thetau)

        L = lift(N, A, aoa)
        D = pressure_drag(N, A, aoa)

        cL, cD, cM = coefficients(L, D, M, aoa)

        L_vs_aoa_m[i] = L
        D_vs_aoa_m[i] = D
        M_vs_aoa_m[i] = M
        cL_vs_aoa_m[i] = cL
        cD_vs_aoa_m[i] = cD
        cM_vs_aoa_m[i] = cM

        print(f"Lift = {L}, Drag = {D}, Moment = {M}")
        print(f"C_L = {cL}, C_D = {cD}, C_M = {cM}")
    
    # Plot scanivalve and validation data.

    val_data = np.loadtxt("./Data/validation_coefficients.txt", skiprows=12, usecols=(0, 1, 2, 4))
    plt.plot(val_data[:, 0], val_data[:, 1], label="cL xfoil")
    plt.plot(val_data[:, 0], val_data[:, 2], label="cD xfoil")
    plt.plot(val_data[:, 0], val_data[:, 3], label="cM xfoil")

    plt.plot(AoA, cL_vs_aoa_s, label="cL, scani")
    plt.plot(AoA, cD_vs_aoa_s, label="cD, scani")
    plt.plot(AoA, cM_vs_aoa_s, label="cM, scani")

    plt.title("$C_L$, $C_D$, and $C_M$ vs $\\alpha$, Scanivalve Data")
    plt.xlabel("$\\alpha^\circ$")
    plt.ylabel("Coefficient")
    plt.legend()
    plt.grid()
    plt.show()

    # Plot manometer and validation data.

    val_data = np.loadtxt("./Data/validation_coefficients.txt", skiprows=12, usecols=(0, 1, 2, 4))
    plt.plot(val_data[:, 0], val_data[:, 1], label="cL xfoil")
    plt.plot(val_data[:, 0], val_data[:, 2], label="cD xfoil")
    plt.plot(val_data[:, 0], val_data[:, 3], label="cM xfoil")

    plt.plot(AoA, cL_vs_aoa_m, label="cL, mano")
    plt.plot(AoA, cD_vs_aoa_m, label="cD, mano")
    plt.plot(AoA, cM_vs_aoa_m, label="cM, mano")

    plt.title("$C_L$, $C_D$, and $C_M$ vs $\\alpha$, Manometer Data")
    plt.xlabel("$\\alpha^\circ$")
    plt.ylabel("Coefficient")
    plt.legend()
    plt.grid()
    plt.show()

    # Plot scanivalve and manometer data.

    plt.plot(AoA, cL_vs_aoa_s, label="cL, scani")
    plt.plot(AoA, cD_vs_aoa_s, label="cD, scani")
    plt.plot(AoA, cM_vs_aoa_s, label="cM, scani")

    plt.plot(AoA, cL_vs_aoa_m, label="cL, mano")
    plt.plot(AoA, cD_vs_aoa_m, label="cD, mano")
    plt.plot(AoA, cM_vs_aoa_m, label="cM, mano")

    plt.title("$C_L$, $C_D$, and $C_M$ vs $\\alpha$")
    plt.xlabel("$\\alpha^\circ$")
    plt.ylabel("Coefficient")
    plt.legend()
    plt.grid()
    plt.show()


if __name__ == "__main__":
    main()

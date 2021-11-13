# AER303 - Airfoil Lab (Group 6) - INSTRUCTIONS FOR RUNNING CODE
# Prerequisites
* MATLAB2019b+
* Python 3.7+ with `numpy` and `matplotlib`

# Pressure Tap Locations
* A MATLAB variable `location` is saved at `/Code/location.mat`. It contains the pressure tab locations on the airfoil, stored as a row vector. It it not necessary to load this manually, as the code will do so automatically.

# Pre-Processing
1. Open MATLAB to the root folder of the repository.
2. Run `/Code/Interpolation.mlx`.
    * This generates the pressure distributions from the scanivalve data.
    * The wing pressure data (`p0, p3, ... , p20`) are stored in `/Code/Calibrated_Airfoil_Pressure_Pa.mat` as a MATLAB variable.
    * The wake pressure data (`p_rakea_0, p_rakeb_0, p_rakea_3, p_rakeb_3, …`) are stored in `/Code/Calibrated_Rake_Pa.mat` as a MATLAB variable.
3. Change directory in MATLAB to `/Code`
4. Run `Run manometer_conversion.m`
    * This generates the pressure distributions from the manometer data.
    * This data is stored as the following files:
        * `/Data/calibrated_manometer_data_Pa.mat`; containing data from all pressure taps and their uncertainties
        * `/Data/calibrated_manometer_airfoil.csv` and `/Data/calibrated_manometer_airfoil_error.csv`; containing the same data but in CSV form for Python analysis
5. Run `freestream_velocity.m`
    * This generates the freestream velocities for scanivalve & manometer data from the rake velocities.
    * This data is stored as `/Data/velocities_manometer.csv` and `/Data/velocities_scanivalve.csv`; containing freestream velocities as well as the corresponding uncertainty.

# Autocorrelation
1. Open MATLAB to the root directory of the repository.
2. Run `/Code/Autocorrelation.mlx`. 
    * This livescript calculates the autocorrelation of the Scanivalve data using port 36 on all 11 AOAs. 
    * Subsequent variable `err` is saved in `/Data/autocorrelation_error_deltaP.mat`.

# C_p Plots
1. Open MATLAB to the root directory of the repository.
2. Run `/Code/Data_vis_forloop.mlx`.
    * This file is used to generate the calibrated C_p plots from Scanivalve, Manometer, and XFOIL for all AOAs.
    * The uncertainty calculations for Scanivalve and Manometer are also computed in the script.
    * Images are saved at `Data/C_p Graphs/AOA%d.png`, where %d is the angle of attack of the corresponding graph. 
3. Note: the filename of this script is listed as `C_p_Calc.mlx` in Appendix C.5 of the report. However, the name has been reverted to `Data_vis_forloop.mlx` for unknown reasons after the final commit. Please ignore the label mistake and follow this latest readme. 

# Aerodynamic Forces
1. In a Python IDE or similar, open the root of the repository. Then, run `/Code/airfoil_force_analysis/lift_and_pressure_drag.py` file. This file is used to calculate and save the lift, drag, and moment coefficients along with their uncertainties. This will update data files within the `/Data` directory (8 in total).
2. In MATLAB, change directory to `/Code`.
3. Run `coefficient_plots.m`. This file will call upon the files generated in the first step and automatically save the generated plots to the code directory.

# Wake and Drag
1. In MATLAB, change directory to `/Code/wake_analysis`.
2. Run `wake_analysis_scanivalve.m`. This generates the c_D and total drag data for the scanivalve measurements, along with the plots for wake velocity and c_D as `.csv` files in the `wake_analysis` folder.
3. Run `wake_analysis_manometer.m`. This generates the c_D and total drag data for the manometer measurements, along with the plots for wake velocity and c_Das `.csv` files in the `wake_analysis` folder.

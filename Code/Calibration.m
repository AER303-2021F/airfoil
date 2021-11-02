clear
format long
close all
clc


coeffs        = [115, -5];              % Calibration coefficients for the linear fit 
Scani_V       = linspace(-10, 10, 100);   % Scanivalve output 
P_mmH20       = polyval(coeffs, Scani_V); % Betz manometer output in mmH20. Convert these to pascals

% Use Scani_V, P_mmH20 values as a look-up table using 'interp1' function
% to convert the actual scanivalve voltages from the 36 ports to
% corresponding pressure readings. 
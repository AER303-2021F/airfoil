clear
format long
close all
clc

% For AOA: 0, 7, 10, 12, 15, 17, 19
coeffs        = [115, 38];              % Calibration coefficients for the linear fit 
Scani_V       = linspace(-10, 10, 100);   % Scanivalve output 
P_mmH20       = polyval(coeffs, Scani_V); % Betz manometer output in mmH20. Convert these to pascals

% Use Scani_V, P_mmH20 values as a look-up table using 'interp1' function
% to convert the actual scanivalve voltages from the 36 ports to
% corresponding pressure readings. 


% For AOA: 3,9,11,16. Somehow the calibration got changed for these angles
% of attack. 

coeffs        = [110, -5];              % Calibration coefficients for the linear fit 
Scani_V       = linspace(-10, 10, 100);   % Scanivalve output 
P_mmH20       = polyval(coeffs, Scani_V); % Betz manometer output in mmH20. Convert these to pascals

% For AOA: 14

coeffs        = [110, 5];              % Calibration coefficients for the linear fit 
Scani_V       = linspace(-10, 10, 100);   % Scanivalve output 
P_mmH20       = polyval(coeffs, Scani_V); % Betz manometer output in mmH20. Convert these to pascals

% Discard AOA = 9 deg data. Does not look right. Something went wrong
% here.....................................................................
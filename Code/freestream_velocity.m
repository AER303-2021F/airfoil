% Calculation of free-stream velocity for each AoA

load("../Data/calibrated_manometer_data_Pa.mat");
% Load manometer pressures
% Manometer error, Pa
load("../Data/Calibrated_Rake_Pa.mat");
% Load scanivalve pressures for rake

% Load manometer and scanivalve pressure data
RHO = 1.225;
d_p_scani = 115 * readmatrix("../Data/scanivalve_error.csv");
% Error from all scanivalve time series, Pa

alphas = [0 3 6 8 10 11 13 15 16 17 20];
v_inf_scanivalve = zeros(size(alphas));
d_v_inf_scanivalve = zeros(size(alphas));
v_inf_manometer = zeros(size(alphas));
d_v_inf_manometer = zeros(size(alphas));

for i = 1:11
    alpha = alphas(i);
    combined_rakes_scanivalve = [eval(sprintf("p_rakeb_%d", alpha)); eval(sprintf("p_rakea_%d", alpha))];
    combined_rakes_scanivalve = combined_rakes_scanivalve(:);
    v_inf_scanivalve(i) = 1/2 * (sqrt(2 * combined_rakes_scanivalve(1) / RHO) + sqrt(2 * combined_rakes_scanivalve(end) / RHO));
    d_v_inf_scanivalve(i) = sqrt((d_p_scani(i)/(2*v_inf_scanivalve(i)*RHO))^2 + (d_p_scani(i)/(2*v_inf_scanivalve(i)*RHO))^2);
    
    combined_rakes_manometer = [rake_b_manometer(i,:); rake_a_manometer(i,:)];
    combined_rakes_manometer = combined_rakes_manometer(:);
    d_combined_rakes_manometer = [d_rake_b_manometer(i,:); d_rake_a_manometer(i,:)];
    d_combined_rakes_manometer = d_combined_rakes_manometer(:);
    
    
    v_inf_manometer(i) = 1/2 * (sqrt(2 * combined_rakes_manometer(1) / RHO) ...
                                + sqrt(2 * combined_rakes_manometer(end) / RHO));
    d_v_inf_manometer(i) = sqrt((d_combined_rakes_manometer(1)/(2*v_inf_manometer(i)*RHO))^2 ...
                                + (d_combined_rakes_manometer(end)/(2*v_inf_manometer(i)*RHO))^2);
    
    writematrix([v_inf_manometer; d_v_inf_manometer], "../Data/velocities_manometer.csv");
    writematrix([v_inf_scanivalve; d_v_inf_scanivalve], "../Data/velocities_scanivalve.csv");
end
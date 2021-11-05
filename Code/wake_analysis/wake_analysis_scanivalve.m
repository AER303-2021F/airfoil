% http://brennen.caltech.edu/fluidbook/externalflows/drag/dragNwake.pdf
load("../../Data/Calibrated_Rake_Pa.mat");
XFOIL_data = readmatrix("../../XFOIL Analysis/clark_y_coefficients", "numheaderlines", 12);
alphas = [0 3 6 8 10 11 13 15 16 17 20];
c_D = zeros(size(alphas));
d_c_D = zeros(size(alphas));
drag_force = zeros(size(alphas));
d_drag_force = zeros(size(alphas));
chord = 0.1; % m

d_p_scani = 9.80665 * readmatrix("../../Data/scanivalve_error.csv");
% Load scanivalve error

% Read freestream velocities from data folder
freestream = readmatrix("../../Data/velocities_scanivalve.csv");
d_freestream = freestream(2, :);
freestream = freestream(1, :);

figure
hold on

title("Wake Velocity Profile for varying AoA [Scanivalve]")
xlabel("Velocity Deficit (m/s)", "interpreter", "latex")
ylabel("Vertical Position (m)", "interpreter", "latex")

for i = 1:11
    alpha = alphas(i);
    RHO = 1.225; % kg/m^3
    locations_a = [0, 1.67 3.33 5 6 7 8 9 10 11 12 13 14 15 16.67 18.33, 20] + 0.5; % Raised position
    locations_b = [0, 1.67 3.33 5 6 7 8 9 10 11 12 13 14 15 16.67 18.33, 20]; % Lowered position

    % Merge the rake data together
    combined_locations = [locations_b; locations_a]; % indices in ascending order
    combined_locations = combined_locations(:) / 100; % locations in m
    combined_rakes = [eval(sprintf("p_rakeb_%d", alpha)); eval(sprintf("p_rakea_%d", alpha))];
    combined_rakes = combined_rakes(:);

    velocities = sqrt(2 * combined_rakes / RHO);
    d_velocities = d_p_scani(i) ./ (2 * velocities * RHO);
    
    vel_diff_err = sqrt(d_freestream(i)^2 + d_velocities.^2);

    errorbar(freestream(i)- velocities, combined_locations, vel_diff_err, "horizontal", "DisplayName", sprintf("%d", alpha));
    
    drag_force(i) = RHO * trapz(combined_locations,  velocities .* (freestream(i) - velocities));
    
    % Calculate drag error
    integrand_error_squared = abs(velocities .* d_freestream(i) +...
                (freestream(i)-2.*velocities) .* d_velocities);
    d_drag_force(i) = sqrt(sum(integrand_error_squared));
    
    q = 1/2 * RHO * freestream(i).^2;
    d_q = d_freestream(i) * RHO * freestream(i);
    
    c_D(i) = drag_force(i) / (q * chord);
    d_c_D = sqrt((d_drag_force / (q * chord)).^2 + (c_D(i) / q * d_q).^2);
end
legend
grid
saveas(gcf,'../../source_latex/figures/scanivalve_wake_velocities.png')

figure
errorbar(alphas, c_D, d_c_D, "DisplayName", "Experiment")
hold on
plot(XFOIL_data(:, 1), XFOIL_data(:,3), "DisplayName", "XFOIL")
xlabel("$\alpha$ (degrees)", "interpreter", "latex")
ylabel("$c_D$ (dimensionless)", "interpreter", "latex")
title("Comparison of XFOIL with Actual c_D variation [Scanivalve]")
legend
grid
saveas(gcf,'../../source_latex/figures/scanivalve_cd.png')

figure
errorbar(alphas, drag_force, d_drag_force)
xlabel("$\alpha$ (degrees)", "interpreter", "latex")
ylabel("Drag Force per Unit Span (N/m) [Scanivalve]")
title("Drag Force")
grid
saveas(gcf,'../../source_latex/figures/scanivalve_drag_force.png')

% Export CSV
writematrix(c_D, "cd_scanivalve.csv")
writematrix(d_c_D, "cd_error_scanivalve.csv")
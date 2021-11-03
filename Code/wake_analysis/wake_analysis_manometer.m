% http://brennen.caltech.edu/fluidbook/externalflows/drag/dragNwake.pdf
load("../../Data/calibrated_manometer_data_Pa.mat");
XFOIL_data = readmatrix("../../XFOIL Analysis/clark_y_coefficients", "numheaderlines", 12);
alphas = [0 3 6 8 10 11 13 15 16 17 20];
c_D = zeros(size(alphas));
drag_force = zeros(size(alphas));
chord = 0.1; % m

figure
hold on

title("Wake Velocity Profile for varying AoA")
xlabel("Velocity Deficit (m/s)", "interpreter", "latex")
ylabel("Vertical Position (m)", "interpreter", "latex")

for i = 1:11
    alpha = alphas(i);
    RHO = 1.225; % kg/m^3
    locations_a = [0, 1.67 3.33 5 6 7 8 9 10 11 12 13 14 15 16.67 18.33, 20] + 0.5; % Raised position
    locations_b = [0, 1.67 3.33 5 6 7 8 9 10 11 12 13 14 15 16.67 18.33, 20]; % Lowered position

    combined_locations = [locations_b; locations_a]; % indices in ascending order
    combined_locations = combined_locations(:) / 100; % locations in m
    combined_rakes = [rake_b_manometer(i,:); rake_a_manometer(i,:)];
    combined_rakes = combined_rakes(:);
    velocities = sqrt(2 * combined_rakes / RHO);
    
    V_inf = 1/2 * (sqrt(2 * combined_rakes(1) / RHO) + sqrt(2 * combined_rakes(end) / RHO));

    plot(V_inf- velocities, combined_locations, "DisplayName", sprintf("%d", alpha));

    drag_force(i) = RHO * trapz(combined_locations,  velocities .* (V_inf - velocities));
    c_D(i) = drag_force(i) / (1/2 * RHO * V_inf^2 * chord);
    %fprintf("D = %f N, CD = %f for AoA = %d\n", drag_force, c_D(i), alpha);
end
legend
grid
saveas(gcf,'manometer_wake_velocities.png')

figure
plot(alphas, c_D, "DisplayName", "Experiment")
hold on
plot(XFOIL_data(:, 1), XFOIL_data(:,3), "DisplayName", "XFOIL")
xlabel("$\alpha$ (degrees)", "interpreter", "latex")
ylabel("$c_D$ (dimensionless)", "interpreter", "latex")
title("Comparison of XFOIL with Actual c_D variation")
legend
grid
saveas(gcf,'manometer_cd.png')

figure
plot(alphas, drag_force)
xlabel("$\alpha$ (degrees)", "interpreter", "latex")
ylabel("Drag Force (N)")
title("Drag Force")
grid
saveas(gcf,'manomater_drag_force.png')

% http://brennen.caltech.edu/fluidbook/externalflows/drag/dragNwake.pdf
load("alpha_0_wake_pressure_test.mat")
load("../../Data/Calibrated_Rake_Pa.mat")
alphas = [0 3 6 8 10 11 13 15 16 17 20];
for i = 1:11
    alpha = alphas(i);
    RHO = 1.225; % kg/m^3
    locations_a = [0, 1.67 3.33 5 6 7 8 9 10 11 12 13 14 15 16.67 18.33, 20] + 0.5; % Raised position
    locations_b = [0, 1.67 3.33 5 6 7 8 9 10 11 12 13 14 15 16.67 18.33, 20]; % Lowered position

    combined_locations = [locations_b; locations_a]; % indices in ascending order
    combined_locations = combined_locations(:) / 100; % locations in m
    combined_rakes = [eval(sprintf("p_rakeb_%d", alpha))'; eval(sprintf("p_rakea_%d", alpha))'];
    combined_rakes = combined_rakes(:);
    velocities = sqrt(2 * combined_rakes / RHO);

%     plot(V_mps - velocities, combined_locations)
%     title("Wake Velocity Profile for $\alpha$ = " + sprintf("%d", alpha) + "$^\circ$", "interpreter", "latex")
%     xlabel("Velocity Deficit (m/s)", "interpreter", "latex")
%     ylabel("Vertical Position (m)", "interpreter", "latex")

    drag = RHO * trapz(combined_locations,  velocities .* (V_mps - velocities));
    c_D = drag / (1/2 * RHO * V_mps^2 * 0.1);
    fprintf("D = %f N, CD = %f for AoA = %d\n", drag, c_D, alpha);
end

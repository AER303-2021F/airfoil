

% Import scanivalve data.
scani_coeff_data = importdata("../Data/scanivalve_coefficients.csv");
scani_coeff_uncertainty = importdata("../Data/scanivalve_coefficient_uncertainty.csv");

% Get scanivalve data.
aoa = scani_coeff_data.data(:, 1);
cl = scani_coeff_data.data(:, 2);
cdp = scani_coeff_data.data(:, 3);
cm = scani_coeff_data.data(:, 4);

% Get scanivalve uncertainty data.
cl_err = scani_coeff_uncertainty.data(:, 2);
cdp_err = scani_coeff_uncertainty.data(:, 3);
cm_err = scani_coeff_uncertainty.data(:, 4);

% Plot coefficients.
figure
errorbar(aoa, cl, cl_err)
hold on
errorbar(aoa, cdp, cdp_err)
hold on
errorbar(aoa, cm, cm_err)
xlabel('$\alpha^\circ$', "interpreter", "latex")
ylabel('Coefficients')
title('$C_L$, $C_{D,p}$, and $C_M$ vs $\alpha$ for Scanivalve Data', 'interpreter', 'latex')
legend('$C_L$', '$C_{D,p}$', '$C_M$', 'interpreter', 'latex')
grid on
set(gca, 'FontSize', 15)
saveas(gcf,'scanivalve_cl_cd_cm.png')

% Plot cL vs cD.
figure
errorbar(cdp, cl, cl_err, cl_err, cdp_err, cdp_err)
xlabel('$C_{D,p}$', "interpreter", "latex")
ylabel('$C_L$', "interpreter", "latex")
title('$C_L$ vs $C_{D,p}$ for Scanivalve Data', 'interpreter', 'latex')
grid on
set(gca, 'FontSize', 15)
saveas(gcf,'scanivalve_cl_vs_cd.png')

% Import manometer data.
mano_coeff_data = importdata("../Data/manometer_coefficients.csv");
mano_coeff_uncertainty = importdata("../Data/manometer_coefficient_uncertainty.csv");

% Get manometer data.
aoa = mano_coeff_data.data(:, 1);
cl = mano_coeff_data.data(:, 2);
cdp = mano_coeff_data.data(:, 3);
cm = mano_coeff_data.data(:, 4);

% Get manometer uncertainty data.
cl_err = mano_coeff_uncertainty.data(:, 2);
cdp_err = mano_coeff_uncertainty.data(:, 3);
cm_err = mano_coeff_uncertainty.data(:, 4);

% Plot coefficients.
figure
errorbar(aoa, cl, cl_err)
hold on
errorbar(aoa, cdp, cdp_err)
hold on
errorbar(aoa, cm, cm_err)
xlabel('$\alpha^\circ$', "interpreter", "latex")
ylabel('Coefficients')
title('$C_L$, $C_{D,p}$, and $C_M$ vs $\alpha$ for Manometer Data', 'interpreter', 'latex')
legend('$C_L$', '$C_{D,p}$', '$C_M$', 'interpreter', 'latex')
grid on
set(gca, 'FontSize', 15)
saveas(gcf,'manometer_cl_cd_cm.png')

% Plot cL vs cD.
figure
errorbar(cdp, cl, cl_err, cl_err, cdp_err, cdp_err)
xlabel('$C_{D,p}$', "interpreter", "latex")
ylabel('$C_L$', "interpreter", "latex")
title('$C_L$ vs $C_{D,p}$ for Manometer Data', 'interpreter', 'latex')
grid on
set(gca, 'FontSize', 15)
saveas(gcf,'manometer_cl_vs_cd.png')

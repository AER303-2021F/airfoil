

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

total_drag_scani = importdata("wake_analysis/cd_scanivalve.csv");
total_drag_scani_error = importdata("wake_analysis/cd_error_scanivalve.csv");

figure
errorbar(total_drag_scani, cl, cl_err, cl_err, total_drag_scani_error, total_drag_scani_error, 'linewidth', 0.5)
xlabel('$C_{D}$', "interpreter", "latex")
ylabel('$C_L$', "interpreter", "latex")
title('$C_L$ vs $C_{D}$ for Scanivalve Data', 'interpreter', 'latex')
grid on
set(gca, 'FontSize', 15)
saveas(gcf,'scanivalve_cl__vs_cd_tot.png')

% Plot coefficients.
figure
yyaxis left
errorbar(aoa, cl, cl_err, 'linewidth', 0.5)
ylabel('$C_L$ Scale', 'interpreter', 'latex')
hold on
yyaxis right
errorbar(aoa, cdp, cdp_err, 'linewidth', 0.5)
hold on
errorbar(aoa, cm, cm_err, 'linewidth', 0.5)
xlabel('$\alpha^\circ$', "interpreter", "latex")
ylabel('$C_{D,p}$ and $C_M$ Scale', 'interpreter', 'latex')
title('$C_L$, $C_{D,p}$, and $C_M$ vs $\alpha$ for Scanivalve Data', 'interpreter', 'latex')
legend('$C_L$', '$C_{D,p}$', '$C_M$', 'interpreter', 'latex')
grid on
set(gca, 'FontSize', 15)
saveas(gcf,'scanivalve_cl_cd_cm.png')

% Plot cL vs cD.
figure
errorbar(cdp, cl, cl_err, cl_err, cdp_err, cdp_err, 'linewidth', 0.5)
xlabel('$C_{D,p}$', "interpreter", "latex")
ylabel('$C_L$', "interpreter", "latex")
title('$C_L$ vs $C_{D,p}$ for Scanivalve Data', 'interpreter', 'latex')
grid on
set(gca, 'FontSize', 15)
saveas(gcf,'scanivalve_cl_vs_cd.png')

val_data = importdata("../Data/validation_coefficients.csv");
val_aoa = val_data(:, 1);
val_cl = val_data(:, 2);
val_cdp = val_data(:, 3);
val_cm = val_data(:, 4);

figure
yyaxis left
plot(val_aoa, val_cl)
hold on
errorbar(aoa, cl, cl_err)
yyaxis right
hold on
plot(val_aoa, val_cdp)
hold on
errorbar(aoa, cdp, cdp_err)
hold on
plot(val_aoa, val_cm)
hold on
errorbar(aoa, cm, cm_err)
xlabel('$\alpha^\circ$', "interpreter", "latex")
ylabel('$C_{D,p}$ and $C_M$ Scale', 'interpreter', 'latex')
title('$C_L$, $C_{D,p}$, and $C_M$ vs $\alpha$ for Scanivalve and XFOIL Data', 'interpreter', 'latex')
legend('$C_L$', '$C_L$, XFOIL', '$C_{D,p}$', '$C_{D,p}$, XFOIL', '$C_M$', '$C_M$, XFOIL', 'interpreter', 'latex')
grid on
set(gca, 'FontSize', 15)
saveas(gcf,'scanivalve_XFOIL_cl_vs_cd.png')


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

total_drag_mano = importdata("wake_analysis/cd_manometer.csv");
total_drag_mano_error = importdata("wake_analysis/cd_error_manometer.csv");

figure
errorbar(total_drag_mano, cl, cl_err, cl_err, total_drag_mano_error, total_drag_mano_error, 'linewidth', 0.5)
xlabel('$C_{D}$', "interpreter", "latex")
ylabel('$C_L$', "interpreter", "latex")
title('$C_L$ vs $C_{D}$ for Manometer Data', 'interpreter', 'latex')
grid on
set(gca, 'FontSize', 15)
saveas(gcf,'manometer_cl__vs_cd_tot.png')

% Plot coefficients.
figure
yyaxis left
errorbar(aoa, cl, cl_err, 'linewidth', 0.5)
ylabel('$C_L$ Scale', 'interpreter', 'latex')
hold on
yyaxis right
errorbar(aoa, cdp, cdp_err, 'linewidth', 0.5)
hold on
errorbar(aoa, cm, cm_err, 'linewidth', 0.5)
xlabel('$\alpha^\circ$', "interpreter", "latex")
ylabel('$C_{D,p}$ and $C_M$ Scale', 'interpreter', 'latex')
title('$C_L$, $C_{D,p}$, and $C_M$ vs $\alpha$ for Manometer Data', 'interpreter', 'latex')
legend('$C_L$', '$C_{D,p}$', '$C_M$', 'interpreter', 'latex', 'location', 'southeast')
grid on
set(gca, 'FontSize', 15)
saveas(gcf,'manometer_cl_cd_cm.png')

% Plot cL vs cD.
figure
errorbar(cdp, cl, cl_err, cl_err, cdp_err, cdp_err, 'linewidth', 0.5)
xlabel('$C_{D,p}$', "interpreter", "latex")
ylabel('$C_L$', "interpreter", "latex")
title('$C_L$ vs $C_{D,p}$ for Manometer Data', 'interpreter', 'latex')
grid on
set(gca, 'FontSize', 15)
saveas(gcf,'manometer_cl_vs_cd.png')

figure
yyaxis left
plot(val_aoa, val_cl)
hold on
errorbar(aoa, cl, cl_err)
yyaxis right
hold on
plot(val_aoa, val_cdp)
hold on
errorbar(aoa, cdp, cdp_err)
hold on
plot(val_aoa, val_cm)
hold on
errorbar(aoa, cm, cm_err)
xlabel('$\alpha^\circ$', "interpreter", "latex")
ylabel('$C_{D,p}$ and $C_M$ Scale', 'interpreter', 'latex')
title('$C_L$, $C_{D,p}$, and $C_M$ vs $\alpha$ for Manometer and XFOIL Data', 'interpreter', 'latex')
legend('$C_L$', '$C_L$, XFOIL', '$C_{D,p}$', '$C_{D,p}$, XFOIL', '$C_M$', '$C_M$, XFOIL', 'interpreter', 'latex', 'location', 'northwest')
grid on
set(gca, 'FontSize', 15)
saveas(gcf, 'manometer_XFOIL_cl_vs_cd.png')

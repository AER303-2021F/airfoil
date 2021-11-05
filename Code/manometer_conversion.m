alphas = [0 3 6 8 10 11 13 15 16 17 20];

theta = 61 * pi/180; % degrees
d_theta = 0.5 * pi/180; % degrees uncertainty

d_cm = 0.1; % cm uncertainty on manometer

cmH20_to_Pa = 98.0665; % Pa / (cmH20)

REF_HEIGHT = 18.2; % cm, initial guess 15.6, tuned manually to get c_p

% Arrays for exporting
airfoil_export = zeros(11, 20);
d_airfoil_export = zeros(11, 20);
airfoil_manometer = zeros(11, 19);
d_airfoil_manometer = zeros(11, 19);
rake_a_manometer = zeros(11, 17);
d_rake_a_manometer = zeros(11, 17);
rake_b_manometer = zeros(11, 17);
d_rake_b_manometer = zeros(11, 17);

for i = 1:11
    alpha = alphas(i);
    airfoil_export(i, 1) = alpha;
    d_airfoil_export(i, 1) = alpha;
    manometer_heights_a = readmatrix(sprintf(...
        "../Data/jai_enumerated_data/data_%da.csv", alpha));
    manometer_heights_b = readmatrix(sprintf(...
        "../Data/jai_enumerated_data/data_%db.csv", alpha));
    
    p_airfoil = -(manometer_heights_a(1:19, 2)'- REF_HEIGHT)...
        * sin(theta) * cmH20_to_Pa;
    airfoil_export(i, 2:end) = p_airfoil;
    airfoil_manometer(i,:) = p_airfoil;
    d_airfoil_manometer(i, :) = sqrt((-cmH20_to_Pa * sin(theta) * d_cm).^2 ...
    + (p_airfoil / tan(theta) * d_theta).^2);
    d_airfoil_export(i, 2:end) = d_airfoil_manometer(i, :);
    
    p_rake_1 = -(manometer_heights_a(20:end, 2)' - REF_HEIGHT) * sin(theta) ...
        * cmH20_to_Pa;
    rake_a_manometer(i, :) = p_rake_1;
    d_rake_a_manometer(i, :) = sqrt((-cmH20_to_Pa * sin(theta) * d_cm).^2 ...
    + (p_rake_1 / tan(theta) * d_theta).^2);
    
    p_rake_2 = - (manometer_heights_b(20:end, 2)' - REF_HEIGHT) * sin(theta) ...
        * cmH20_to_Pa;
    rake_b_manometer(i, :) = p_rake_2;
    d_rake_b_manometer(i, :) = sqrt((-cmH20_to_Pa * sin(theta) * d_cm).^2 ...
        + (p_rake_2 / tan(theta) * d_theta).^2);
    
end
writematrix(airfoil_export, "../Data/calibrated_manometer_airfoil.csv");
writematrix(d_airfoil_export, "../Data/calibrated_manometer_airfoil_error.csv");
save("../Data/calibrated_manometer_data_Pa.mat", ...
    "airfoil_manometer", "rake_a_manometer", "rake_b_manometer", ...
    "d_airfoil_manometer", "d_rake_a_manometer", "d_rake_b_manometer");


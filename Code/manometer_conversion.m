alphas = [0 3 6 8 10 11 13 15 16 17 20];

theta = 61 * pi/180; % degrees
d_theta = 0.1; % degrees uncertainty

cmH20_to_Pa = 98.0665; % Pa / (cmH20)

REF_HEIGHT = 16; % cm, initial guess 15.6

airfoil_export = zeros(11, 20);
airfoil_manometer = zeros(11, 19);
rake_a_manometer = zeros(11, 17);
rake_b_manometer = zeros(11, 17);

for i = 1:11
    alpha = alphas(i);
    airfoil_export(i, 1) = alpha;
    manometer_heights_a = readmatrix(sprintf("../Data/jai_enumerated_data/data_%da.csv", alpha));
    manometer_heights_b = readmatrix(sprintf("../Data/jai_enumerated_data/data_%db.csv", alpha));
    
    p_airfoil = -(manometer_heights_a(1:19, 2)'- REF_HEIGHT) * sin(theta) * cmH20_to_Pa;
    airfoil_export(i, 2:end) = p_airfoil;
    airfoil_manometer(i,:) = p_airfoil;
    p_rake_1 = -(manometer_heights_a(20:end, 2)' - REF_HEIGHT) * sin(theta) * cmH20_to_Pa;
    rake_a_manometer(i, :) = p_rake_1;
    p_rake_2 = - (manometer_heights_b(20:end, 2)' - REF_HEIGHT) * sin(theta) * cmH20_to_Pa;
    rake_b_manometer(i, :) = p_rake_2;
end
writematrix(airfoil_export, "../Data/calibrated_manometer_airfoil.csv");
save("../Data/calibrated_manometer_data_Pa.mat", "airfoil_manometer", "rake_a_manometer", "rake_b_manometer");


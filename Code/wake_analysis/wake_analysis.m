% http://brennen.caltech.edu/fluidbook/externalflows/drag/dragNwake.pdf
locations_a = [0, 1.67 3.33 5 6 7 8 9 10 11 12 13 14 15 16.67 18.33, 20] + 1;
locations_b = [0, 1.67 3.33 5 6 7 8 9 10 11 12 13 14 15 16.67 18.33, 20] + 0.5;

combined_locations = [locations_b; locations_a]; % indices in ascending order
combined_locations = combined_locations(:)

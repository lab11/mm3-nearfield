header_distance_all = [header_distance_0, header_distance_1, header_distance_2, header_distance_3, header_distance_4];
seed_one_16 = [16, 32, 3, 6, 12, 24, 48, 35, 5, 10, 20, 40, 19, 38, 15];
seed_one_8 = [16, 32, 3, 6, 12, 24, 48];
seed_16_all = [seed_one_16, seed_one_16, seed_one_16, seed_one_16, seed_one_16];
seed_8_all = [seed_one_8, seed_one_8, seed_one_8, seed_one_8, seed_one_8];
one_array = ones(1 ,35);
X = [seed_8_all; one_array]';
Y = header_distance_all';

B = inv((X' * X)) * X' * Y;
estimated_Y = X * B;
error = abs(Y - estimated_Y);

for n = 1:15
    error_percent(n) = error(n)/Y(n);
end
clc;
clear all;
close all;
gr_read_complex_binary('../iq_data_v8_chip_18_FAFA.dat');
raw = abs(ans);
find_peak(raw);
pos = ans;
header_pos = pos(1:32);
header_distance = cal_dis(header_pos);
[error, error_percent, B] = calibrate_fn(header_distance);
for n = 1: length(pos) - 1
    all_distance(n) = pos(n+1) - pos(n); 
end
even_dis = (all_distance(18) ...
            + all_distance(20) ...
            + all_distance(22) ...
            + all_distance(24) ...
            + all_distance(26) ...
            + all_distance(28) ...
            + all_distance(30))/7;
odd_dis = (all_distance(17) ...
            + all_distance(19) ...
            + all_distance(21) ...
            + all_distance(23) ...
            + all_distance(25) ...
            + all_distance(27) ...
            + all_distance(29) ...
            + all_distance(31))/8;   
all_dis = (odd_dis * 8 +even_dis * 7)/15;


for k = 1:15
    if(error_percent(k) > 0.01) 
        error('unable to calibrate, check input\n');
    end
end
fprintf('calibration succeeded\n');

fileID = fopen('calibrate.txt','w');
offset = B(1)/(542.4 * 2);
value = abs(1 - offset) / 0.0015 + 20;
fprintf(fileID, 'offset is: %f\n', floor(value) + 0.5);
fprintf('offset is: %f\n', floor(value) + 0.5);

            

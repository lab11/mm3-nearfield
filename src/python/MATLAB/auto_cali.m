clc;
clear all;
close all;
gr_read_complex_binary('../calibrate.dat');

%gr_read_complex_binary('../iq_data_v8_chip_3_FAFA.dat');

raw = abs(ans);
find_peak(raw);
pos = ans;
num = 1;
for n = 1: length(pos) - 1
    if(pos(n+1) - pos(n) < 40000*2 && pos(n+1) - pos(n) > 3000*2)
        all_distance_raw(n) = pos(n+1) - pos(n);
    else
        all_distance_raw(n) = 0;
    end
end
all_distance(1) = 0;
start = 0;
for n = 1: length(all_distance_raw) - 1
    if(all_distance_raw(n) > 14000*2 && all_distance_raw(n) < 17000*2 && start == 0)
        start = 1;
    end
    if(all_distance_raw(n) > 0 && start == 1)
        all_distance(num) = all_distance_raw(n);
        num = num + 1;
    end
end
if(length(all_distance)>15) 
    header_distance = all_distance(1:15);
    [error, error_percent, B] = calibrate_fn(header_distance);
else
    error('unable to calibrate, check input\n');
end

% even_dis = (all_distance(18) ...
%             + all_distance(20) ...
%             + all_distance(22) ...
%             + all_distance(24) ...
%             + all_distance(26) ...
%             + all_distance(28) ...
%             + all_distance(30))/7;
% odd_dis = (all_distance(17) ...
%             + all_distance(19) ...
%             + all_distance(21) ...
%             + all_distance(23) ...
%             + all_distance(25) ...
%             + all_distance(27) ...
%             + all_distance(29) ...
%             + all_distance(31))/8;   
% all_dis = (odd_dis * 8 +even_dis * 7)/15;


for k = 1:15
    if(error_percent(k) > 0.03) 
        error('unable to calibrate, check input\n');
    end
end
fprintf('calibration succeeded\n');

fileID = fopen('calibrate.txt','w');
offset = B(1)/(542.4 * 2);
%offset = B(1)/(542.4);
if(abs(1 - offset) < 0) 
    value = abs(1 - offset) / 0.0015 + 20;
else
    value = (1 - offset)/0.0015 + 20;
end
fprintf(fileID, 'offset is: %f\n', floor(value) + 0.5);
fprintf('offset is: %f\n', floor(value) + 0.5);

            

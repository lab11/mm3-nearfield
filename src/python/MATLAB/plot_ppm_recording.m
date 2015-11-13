clear all;

raw = gr_read_complex_binary('../calibrate.dat');
data = raw;

figure(1);
plot((1:length(data)),abs(data));
%plot((1:length(data))./12.5e6,data);
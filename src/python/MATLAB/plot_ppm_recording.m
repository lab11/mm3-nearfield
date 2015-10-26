clear all;

raw = gr_read_complex_binary('../s_4_h_8_d_8_c_30_2.dat');
data = raw;

figure(1);
plot((1:length(data)),data);
%plot((1:length(data))./12.5e6,data);
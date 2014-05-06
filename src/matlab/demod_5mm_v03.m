% mm3 demod
% Nathan E. Roberts
% 04.23.2011

% usrp_out_5mm: data point 1 -> 25e6

% refresh
clear;
clc;

% variables
sample_rate = 12.5e6;       % samples/sec
seg_start = 12.5e6;         % segment beginning point
seg_end = 22.5e6;
% .mat file was generated using break_up_file_XXmm.m
% matload = load('usrp_out_5mm_1_25e6.mat');
% matload = load('usrp_out_5mm_25e6_50e6.mat');
% matload = load('usrp_out_5mm_50e6_75e6.mat');
% matload = load('usrp_out_5mm_75e6_100e6.mat');
% matload = load('usrp_out_5mm_100e6_125e6.mat');
matload = load('usrp_out_5mm_125e6_150e6.mat');
% matload = load('usrp_out_5mm_150e6_175e6.mat');
% matload = load('usrp_out_5mm_175e6_200e6.mat');
% matload = load('usrp_out_5mm_200e6_225e6.mat');
% matload = load('usrp_out_5mm_225e6_250e6.mat');
% matload = load('usrp_out_5mm_250e6_275e6.mat');
% matload = load('usrp_out_5mm_275e6_end.mat');
data = matload.data_iq;

% calculations
period = 1/sample_rate;
time = linspace(0,period*length(data)-period,length(data));

% calculate magnitude of signal and signal power
data_i = real(data(seg_start:seg_end));                      % I data
data_q = imag(data(seg_start:seg_end));                      % Q data
data_volt = sqrt((data_i.^2)+(data_q.^2));      % units = volts?
data_pwr = data_volt.^2;                        % signal power

%% plotting

figure(1);
plot(time(seg_start:seg_end),data_pwr,'b','LineWidth',1.5);
grid on;
xlabel('Time [s]','FontSize',14);
ylim([0 2]);
set(1,'Position',[300,300,1000,500]);
set(gca,'FontSize',14);

%% demodulation algorithm
% Data is OOK-modulated where a 1 is a pulse of energy and a 0 is no pulse.
% Consecutive 1's will look like multiple pulses at a consisent rate (not a
% single really long pulse).
% Steps:
% 1. Set threshold to clip data
% 2. Find a sync header, which is defined as 4 consecutive bits (1111)
%    Determine pulse width and pulse repetition frequency (PRF) from sync
%    header
% 3. Using that information demodulate the remaining N bits, where N is
%    provided.

% variables
N = 5;                      % provided information
threshold = 0.6;            % threshold set after observing data
SDR_sample_rate = 12.5e6;   % sample rate of software defined radio receiver
pulse_max = 550e-9;         % maximum pulse width
pulse_min = 450e-9;         % minimum pulse width
prf_max = 12e-3;            % maximum pulse repetition frequency
prf_min = 10e-3;            % minimum pulse repetition freuqnecy

% setup computations
sample_period = 1/SDR_sample_rate;

% Treating this like it's streaming data so I won't use some of the 
% efficient Matlab functions that would speed this up.
last_data = 0;      % previous received data point
pulse_count = 0;    % length of current pulse (1s)
prf_count = 0;      % length of current prf (0s)
sync = 0;           % current consecutive sync pulses
header = 4;         % # sync pulses needed
sync_prf = 0;       % prf counter in sync part of loop
sync_prf2 = 0;      % prf counter in sync part of loop (prf_window)
sync_pulse = 0;     % pulse counter in sync part of loop
prf_win_cnt = 0;    % counter for window of pulses between prf_min and prf_max in sync routine
valid_pulse = 0;    % flag for pulse detection
n = 0;              % counter for N

pulse_vec = [];
prf_vec = [];
demod_data = [];

clc

for x=1:length(data_pwr)
    % STEP 1 --------------------------------------------------------------
    if data_pwr(x) > threshold
        rx_data = 1;
    else
        rx_data = 0;
    end
    % STEP 2 --------------------------------------------------------------
    transition = rx_data - last_data;     % see if we are at an edge of a pulse
    % NOT SYNCHRONIZED YET
    if sync < header;       
        if transition > 0           % rising pulse
            prf_val = prf_count*sample_period;
            disp('prf_count');
            disp(prf_val);
            if prf_val > prf_min && prf_val < prf_max
                prf_vec = [prf_vec prf_val];
            else
                sync = 0;                     % reset the sync counter
                prf_vec = [];                 % reset prf vector
                pulse_vec = [];               % reset pulse vector
            end
            pulse_count = 1;     
            prf_count = 0;         
        elseif transition < 0       % falling pulse
            pulse_val = pulse_count*sample_period;
            disp('pulse_count');
            disp(pulse_val);
            if pulse_val > pulse_min && pulse_val < pulse_max
                sync = sync+1;              % valid pulse
                pulse_vec = [pulse_vec pulse_val];
                disp('valid pulse');
            else
                prf_count = pulse_count;    % it was noise...not pulse - count it towards the prf
            end
            pulse_count = 0;
            prf_count = prf_count+1;
        else    % no transition
            if rx_data == 1
                pulse_count = pulse_count+1;
            else
                prf_count = prf_count+1;
            end
        end
        last_data = rx_data;
    end
    % SYNCHRONIZED
    if sync == header;
       disp('x');
       disp(x);
       % find average pulse width and prf
       pulse_length = round(mean(pulse_vec)/sample_period);
       prf_length = round(mean(prf_vec)/sample_period);
       % using this, demodulate the next N bits
       % right after last sync pulse is a prf window
       pulse_length_min = pulse_length - 1;
       pulse_length_max = pulse_length + 1;
       prf_length_min = round(prf_length - prf_length*0.1);
       prf_length_max = round(prf_length + prf_length*0.1);
       prf_window = prf_length_max - prf_length_min;
       % advance to the edge of the prf_min window and look for a pulse
       % anywhere in the prf_max-prf_min frame.
       % then return to prf_length and repeat.
       if sync_prf < prf_length_min         % don't care until we get to the point where a pulse might come
           disp('sync_prf');
           disp(sync_prf);
           sync_prf = sync_prf + 1;
       else                                 % start looking for a pulse (we are in the prf window)
           disp('sync_prf2');
           disp(sync_prf2);
           sync_prf2 = sync_prf2 + 1;       % do this for prf window duration
           if rx_data == 1
               sync_pulse = sync_pulse + 1; % might be a pulse
           else
               if sync_pulse >= pulse_length_min & sync_pulse <= pulse_length_max     % valid pulse
                   valid_pulse = 1;
                   sav_pulse = sync_pulse;
               end
               sync_pulse = 0;              % reset
           end
           if valid_pulse == 1              % we found a pulse
               demod_data = [demod_data 1];
               n = n + 1;
               sync_prf = sav_pulse;        % prf is defined as rising edge of pulse to the next rising edge of pulse
               sync_prf2 = 0;
               valid_pulse = 0;             % reset
           elseif sync_prf2 == prf_window   % ran through whole prf window w/o finding pulse
               demod_data = [demod_data 0];
               n = n + 1;
               sync_prf = prf_length_max - prf_length;      % don't set to 0...reset to middle of window for timing
               sync_prf2 = 0;
           end
       end
       if n == N                            % we've looked for all the data
           sync = 0;                        % start over looking for sync
           prf_vec = [];
           pulse_vec = [];
           prf_count = 0;
           pulse_count = 0;
           n = 0;
           sync_prf = 0;
           sync_prf2 = 0;
           sync_pulse = 0;
           valid_pulse = 0;
%            break;
       end
    end
end
disp('done');


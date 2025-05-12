% ASK, FSK, BPSK, QPSK for a given binary input sequence.
clc;
clear;
close all;

% Parameters
bit_seq = [1 0 1 1 0 0 1 0]; % Example binary sequence
bit_rate = 1;                % Bit rate in bits/sec
fc = 10;                    % Carrier frequency (Hz)
fs = 1000;                  % Sampling frequency (samples per sec)
Tb = 1/bit_rate;            % Bit duration (seconds)
t_bit = 0:1/fs:Tb-1/fs;     % Time vector for one bit

% Time vectors for each modulation (aligned to total signal length)
t_full = [];
for k = 1:length(bit_seq)
    t_full = [t_full t_bit + (k-1)*Tb];
end

% Total time for QPSK (each symbol modulates 2 bits)
num_symbols = ceil(length(bit_seq)/2);
Tb_symbol = 2*Tb;            % QPSK symbol duration (2 bits)
t_symbol = 0:1/fs:Tb_symbol - 1/fs;
t_qpsk_full = [];
for k = 1:num_symbols
    t_qpsk_full = [t_qpsk_full t_symbol + (k-1)*Tb_symbol];
end

%% Amplitude Shift Keying (ASK)
% Bit 1 -> Amplitude 1; Bit 0 -> Amplitude 0
A1 = 1;  
A0 = 0;  

ask_signal = [];
for k = 1:length(bit_seq)
    if bit_seq(k) == 1
        A = A1;
    else
        A = A0;
    end
    ask_signal = [ask_signal A*cos(2*pi*fc*t_bit)];
end

%% Frequency Shift Keying (FSK)
% Bit 1 -> freq f1; Bit 0 -> freq f0
f1 = 15; 
f0 = 5;

fsk_signal = [];
for k = 1:length(bit_seq)
    if bit_seq(k) == 1
        freq = f1;
    else
        freq = f0;
    end
    fsk_signal = [fsk_signal cos(2*pi*freq*t_bit)];
end

%% Binary Phase Shift Keying (BPSK)
% Bit 1 -> phase 0; Bit 0 -> phase pi
bpsk_signal = [];
for k = 1:length(bit_seq)
    if bit_seq(k) == 1
        phase = 0;
    else
        phase = pi;
    end
    bpsk_signal = [bpsk_signal cos(2*pi*fc*t_bit + phase)];
end

%% Quadrature Phase Shift Keying (QPSK)
% Map 2 bits to one symbol:
% 00 -> 0, 01 -> pi/2, 11 -> pi, 10 -> 3pi/2

% Pad with zero if odd number of bits
if mod(length(bit_seq), 2) ~= 0
    bit_seq = [bit_seq 0];
end

qpsk_signal = [];
for k = 1:2:length(bit_seq)
    bits_pair = bit_seq(k:k+1);
    
    if isequal(bits_pair, [0 0])
        phase = 0;
    elseif isequal(bits_pair, [0 1])
        phase = pi/2;
    elseif isequal(bits_pair, [1 1])
        phase = pi;
    else % [1 0]
        phase = 3*pi/2;
    end
    
    qpsk_signal = [qpsk_signal cos(2*pi*fc*t_symbol + phase)];
end

%% Plotting

figure('Name', 'Digital Modulation Techniques', 'NumberTitle', 'off', 'Position', [100 100 1000 900]);

% Plot input binary data
subplot(5,1,1)
stairs([0:length(bit_seq)]*Tb, [bit_seq 0], 'LineWidth', 2);
ylim([-0.5 1.5]);
grid on;
title('Input Binary Data');
xlabel('Time (s)');
ylabel('Amplitude');

% ASK plot
subplot(5,1,2)
plot(t_full, ask_signal, 'b', 'LineWidth', 1.5);
grid on;
title('Amplitude Shift Keying (ASK)');
xlabel('Time (s)');
ylabel('Amplitude');

% FSK plot
subplot(5,1,3)
plot(t_full, fsk_signal, 'r', 'LineWidth', 1.5);
grid on;
title('Frequency Shift Keying (FSK)');
xlabel('Time (s)');
ylabel('Amplitude');

% BPSK plot
subplot(5,1,4)
plot(t_full, bpsk_signal, 'm', 'LineWidth', 1.5);
grid on;
title('Binary Phase Shift Keying (BPSK)');
xlabel('Time (s)');
ylabel('Amplitude');

% QPSK plot
subplot(5,1,5)
plot(t_qpsk_full, qpsk_signal, 'k', 'LineWidth', 1.5);
grid on;
title('Quadrature Phase Shift Keying (QPSK)');
xlabel('Time (s)');
ylabel('Amplitude');

sgtitle('Digital Modulation Techniques Demonstration', 'FontSize', 16);



% End of script

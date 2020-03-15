%{
    Tugas Besar EB3102 - Pengolahan Sinyal Biomedika
    Irfan Tito Kurniawan
    NIM 18317019

    Topik 1 - Analisis Sinyal

    Built-in functions used:
        - fft
        - fftshift
        - abs
        - log10
%}

% Clear the screen, clean all the variables, close all windows
clear; clc; close all;

% Open the file
filename_1 = "./ecg_data/psb_ecg1.dat";
filename_2 = "./ecg_data/psb_ecg2.dat";

x_1 = load(filename_1);
x_2 = load(filename_2);

% Signal properties
sample_frequency = 1000;
power_supply_frequency = 60;
max_frequency = sample_frequency / 2;

signal_length_1 = length(x_1);
n_1 = [0 : signal_length_1 - 1];

signal_length_2 = length(x_2);
n_2 = [0 : signal_length_2 - 1];

% Generate the FFT and the PSD (in dB) of the signal
x_1_fft = fftshift(fft(x_1));
x_1_psd = 20 * log10(abs(x_1_fft) .^ 2);

x_2_fft = fftshift(fft(x_2));
x_2_psd = 20 * log10(abs(x_2_fft) .^ 2);

% Plot Stuffs
% ECG Sample 1
figure;
subplot(2, 1, 1);
plot(n_1, x_1);
title('ECG Sample 1 in Time Domain');
xlabel('n (sample)')
ylabel('x(n)')
xlim([0, signal_length_1])

subplot(2, 1, 2);
stem((n_1 / signal_length_1 * sample_frequency) - (sample_frequency / 2), x_1_psd);
title('ECG Sample 1 PSD');
xlabel('f (Hz)');
ylabel('Power (dB)')

% ECG Sample 2
figure;
subplot(2, 1, 1);
plot(n_2, x_2);
title('ECG Sample 1 in Time Domain');
xlabel('n (sample)')
ylabel('x(n)')
xlim([0, signal_length_2])

subplot(2, 1, 2);
stem((n_2 / signal_length_2 * sample_frequency) - (sample_frequency / 2), x_2_psd);
title('ECG Sample 1 PSD');
xlabel('f (Hz)');
ylabel('Power (dB)')
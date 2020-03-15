%{
    Tugas Besar EB3102 - Pengolahan Sinyal Biomedika
    Irfan Tito Kurniawan
    NIM 18317019

    Topik 2 - Sistem sebagai Filter

    Built-in functions:
        - fft
        - fftshift
        - abs
        - unwrap
        - filter
        - freqz
        - conv
        - zplane (pkg signal)
%}

% Clear the screen, clean all the variables, close all windows
clear; clc; close all;

% Load signal for zplane function
pkg load signal

%{
    Function definitions:
        - design_notch: to design notch at a certain frequency
        - design_comb: to design comb at a certain frequency multiples
        - design_butterworth_lowpass: to design a nth order butterworth LPF
        - design_butterworth_highpass: to design a nth order butterworth HPF
        - design_butterworth_bandpass: to design a nth order butterworth BPF
%}

% Second order notch filter difference equation generator
function [b, a] = design_notch(notch_frequency, max_frequency, pole_distance_from_origin)
    notch_frequency_ratio = notch_frequency / max_frequency;
    cos_omega = cos(notch_frequency_ratio * pi);

    % Difference equation
    a_0 = 1;
    a_1 = -2 * pole_distance_from_origin * cos_omega;
    a_2 = pole_distance_from_origin .^ 2;
    b_0 = 1;
    b_1 = -2 * cos_omega;
    b_2 = 1;

    % Multiply with this so the gain is 1 at the passband
    b = ((pole_distance_from_origin .^ 2 + (2 * pole_distance_from_origin) + 1) / 4) .* [b_0, b_1, b_2];    
    a = [a_0, a_1, a_2];
end

function [b, a] = design_comb(natural_frequency, sample_frequency, pole_distance_from_origin)
    b = [1];
    a = [1];

    % Calculate the needed notch amount
    max_frequency = sample_frequency / 2;
    frequency_ratio = natural_frequency / sample_frequency;
    notch_amount = floor(1 / frequency_ratio);

    % Convolve all the notch filter difference equation
    for i = 1:notch_amount
        [dummy_b, dummy_a] = design_notch(i * natural_frequency, max_frequency, pole_distance_from_origin);
        b = conv(b, dummy_b);
        a = conv(a, dummy_a);
    end
end

% Butterworth lowpass filter difference equation generator
function [b, a] = design_butterworth_lowpass(n, cutoff_frequency, sample_frequency)
    a = [1];
    b = [1];

    % Constant for discretization
    gamma_const = cot(pi * cutoff_frequency / sample_frequency);

    if mod(n, 2) == 0
        for k = 0:(n / 2) - 1
            % Constant for each pole
            alpha_const = 2 * cos(2 * pi * (2 * k + n + 1) / (4 * n));
            dummy_a = [gamma_const .^ 2 - (alpha_const * gamma_const) + 1, -2 * (gamma_const .^ 2) + 2, ...
                gamma_const .^ 2 + (alpha_const * gamma_const) + 1];
            a = conv(a, dummy_a);

            dummy_b = [1, 2, 1];
            b = conv(b, dummy_b);
        end
    else
        for k = 0:((n - 1) / 2) - 1
            alpha_const = 2 * cos(2 * pi * (2 * k + n + 1) / (4 * n));
            dummy_a = [gamma_const .^ 2 - (alpha_const * gamma_const) + 1, -2 * (gamma_const .^ 2) + 2, ...
                gamma_const .^ 2 + (alpha_const * gamma_const) + 1];
            a = conv(a, dummy_a);

            dummy_b = [1, 2, 1];
            b = conv(b, dummy_b);
        end
        
        dummy_a = [gamma_const + 1, 1 - gamma_const];
        a = conv(a, dummy_a);

        dummy_b = [1, 1];
        b = conv(b, dummy_b);
    end
end

% Butterworth highpass filter difference equation generator
function [b, a] = design_butterworth_highpass(n, cutoff_frequency, sample_frequency)
    % Use lowpass filter generator to get the poles
    [b, a] = design_butterworth_lowpass(n, cutoff_frequency, sample_frequency);

    gamma_const = cot(pi * cutoff_frequency / sample_frequency);
    
    % Redefine the zeros
    b = [1];

    for k = 1:n
        dummy_b = gamma_const .* [1, -1];
        b = conv(b, dummy_b);
    end
end

% Butterworth bandpass filter difference equation generator
function [b, a] = design_butterworth_bandpass(n_lower, n_upper, cutoff_frequency_lower, cutoff_frequency_upper, sample_frequency)
    % Generate lowpass and highpass filter difference equation
    [b_low, a_low] = design_butterworth_lowpass(n_upper, cutoff_frequency_upper, sample_frequency);
    [b_high, a_high] = design_butterworth_highpass(n_lower, cutoff_frequency_lower, sample_frequency);

    % Combine both filters into one
    b = conv(b_low, b_high);
    a = conv(a_low, a_high);    
end

% Signal processing
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

% Comb Filter
pole_distance_from_origin = 0.995;
[b_comb, a_comb] = design_comb(power_supply_frequency, sample_frequency, pole_distance_from_origin);

figure;
zplane(b_comb, a_comb)

% Plot the comb filter's frequency response
[h_comb, w_comb] = freqz(b_comb, a_comb, sample_frequency, 'whole', sample_frequency);

h_comb = fftshift(h_comb);
w_comb = fftshift(w_comb);

figure;
subplot(2, 1, 1); 
plot(-max_frequency:max_frequency - 1, 20 * log10(abs(h_comb)));
xlabel('f (Hz)');
ylabel('Magnitude (dB)');
title('Comb (Notch) Magnitude Response');
subplot(2, 1, 2); 
plot(-max_frequency:max_frequency - 1, unwrap(w_comb));
xlabel('f (Hz)');
ylabel('Phase (rad)');
title('Comb (Notch) Phase Response');

% Bandpass Filter
high_pass_frequency = 0.5;
low_pass_frequency = 20;
bandpass_order = 4;

[b_bandpass, a_bandpass] = design_butterworth_bandpass(bandpass_order, bandpass_order, ...
    high_pass_frequency, low_pass_frequency, sample_frequency);

figure;
zplane(b_bandpass, a_bandpass)

% Plot the bandpass filter's frequency response
[h_bandpass, w_bandpass] = freqz(b_bandpass, a_bandpass, sample_frequency, 'whole', sample_frequency);

h_bandpass = fftshift(h_bandpass);
w_bandpass = fftshift(w_bandpass);

figure;
subplot(2, 1, 1); 
plot(-max_frequency:max_frequency - 1, 20 * log10(abs(h_bandpass)));
xlabel('f (Hz)');
ylabel('Magnitude (dB)');
title('Bandpass Magnitude Response');
subplot(2, 1, 2); 
plot(-max_frequency:max_frequency - 1, unwrap(w_bandpass));
xlabel('f (Hz)');
ylabel('Phase (rad)');
title('Bandpass Phase Response');

% Filter out the power line harmonic noises
y_1 = filter(b_comb, a_comb, x_1);
y_2 = filter(b_comb, a_comb, x_2);

% Filter out high and DC frequencies
y_1 = filter(b_bandpass, a_bandpass, y_1);
y_2 = filter(b_bandpass, a_bandpass, y_2);

% Check the resulting signal's PSD
y_1_fft = fftshift(fft(y_1));
y_1_psd = 20 * log10(abs(y_1_fft) .^ 2);

y_2_fft = fftshift(fft(y_2));
y_2_psd = 20 * log10(abs(y_2_fft) .^ 2);

% Plot time-domain signals
% ECG Sample 1
figure;
subplot(2, 1, 1); 
plot(n_1, x_1);
xlabel('n (sample)');
title('Unfiltered Signal');
xlim([0, signal_length_1]);

subplot(2, 1, 2); 
plot(n_1, y_1);
xlabel('n (sample)');
title('Filtered Signal');
xlim([0, signal_length_1]);

figure;
subplot(2, 1, 1);
stem((n_1 / signal_length_1 * sample_frequency) - (sample_frequency / 2), x_1_psd);
xlabel('f (Hz)');
ylabel('Power (dB)')
title('Unfiltered Signal PSD');

subplot(2, 1, 2); 
stem((n_1 / signal_length_1 * sample_frequency) - (sample_frequency / 2), y_1_psd);
xlabel('f (Hz)');
ylabel('Power (dB)')
title('Filtered Signal PSD');

% ECG Sample 2
figure;
subplot(2, 1, 1); 
plot(n_2, x_2);
xlabel('n (sample)');
title('Unfiltered Signal');
xlim([0, signal_length_2]);

subplot(2, 1, 2); 
plot(n_2, y_2);
xlabel('n (sample)');
title('Filtered Signal');
xlim([0, signal_length_2]);

figure;
subplot(2, 1, 1);
stem((n_2 / signal_length_2 * sample_frequency) - (sample_frequency / 2), x_2_psd);
xlabel('f (Hz)');
ylabel('Power (dB)')
title('Unfiltered Signal PSD');

subplot(2, 1, 2); 
stem((n_2 / signal_length_2 * sample_frequency) - (sample_frequency / 2), y_2_psd);
xlabel('f (Hz)');
ylabel('Power (dB)')
title('Filtered Signal PSD');
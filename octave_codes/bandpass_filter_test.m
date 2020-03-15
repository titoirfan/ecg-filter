%{
    Tugas Besar EB3102 - Pengolahan Sinyal Biomedika
    Irfan Tito Kurniawan
    NIM 18317019

    Butterworth bandpass filter performance proof on cosine function
%}

% Clear the screen, clean all the variables, close all windows
clear; clc; close all;

pkg load signal

% Signal properties
sample_frequency = 1000;
power_supply_frequency = 60;
max_frequency = sample_frequency / 2;
duration = 3;
signal_length = duration * sample_frequency;
tstep = duration / signal_length;

freq_1 = 400 * 2 * pi;
freq_2 = 0.001 * 2 * pi;
freq_3 = 10 * 2 * pi;

n = [0:tstep:duration];
n_2 = [0:signal_length];

x = cos(freq_1 * n) + cos(freq_2 * n) + cos(freq_3 * n);

x_fft = fftshift(fft(x));

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

% Bandpass Filter
high_pass_frequency = 1;
low_pass_frequency = 20;
bandpass_order = 4;

[b, a] = design_butterworth_bandpass(bandpass_order, bandpass_order, ...
    high_pass_frequency, low_pass_frequency, sample_frequency);

y = filter(b, a, x);

y_fft = fftshift(fft(y));

figure;
subplot(2, 1, 1); 
plot(n, x);
xlabel('t (s)');
ylabel('x(n)');
title('Original Signal');
subplot(2, 1, 2); 
plot(n, y);
xlabel('t (s)');
ylabel('y(n)');
title('Filtered Signal');

figure;
subplot(2, 1, 1); 
stem((n_2 / signal_length * sample_frequency) - (sample_frequency / 2), x_fft);
xlabel('f (Hz)');
ylabel('X(\omega)');
title('FFT of Original Signal');
subplot(2, 1, 2); 
stem((n_2 / signal_length * sample_frequency) - (sample_frequency / 2), y_fft);
xlabel('f (Hz)');
ylabel('X(\omega)');
title('FFT of Filtered Signal');
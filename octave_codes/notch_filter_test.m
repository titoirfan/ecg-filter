%{
    Tugas Besar EB3102 - Pengolahan Sinyal Biomedika
    Irfan Tito Kurniawan
    NIM 18317019

    Comb filter performance proof on cosine function
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

freq_1 = 60 * 2 * pi;
freq_2 = 120 * 2 * pi;
freq_3 = 10 * 2 * pi;

n = [0:tstep:duration];
n_2 = [0:signal_length];

x = (0.3 * cos(freq_1 * n)) + (0.3 * cos(freq_2 * n)) + cos(freq_3 * n);

x_fft = fftshift(fft(x));

% Second order notch filter difference equation generator
function [b, a] = design_notch(notch_frequency, max_frequency, pole_distance_from_origin)
    notch_frequency_ratio = notch_frequency / max_frequency;
    cos_omega = cos(notch_frequency_ratio * pi);

    a_0 = 1;
    a_1 = -2 * pole_distance_from_origin * cos_omega;
    a_2 = pole_distance_from_origin .^ 2;
    b_0 = 1;
    b_1 = -2 * cos_omega;
    b_2 = 1;

    b = ((pole_distance_from_origin .^ 2 + (2 * pole_distance_from_origin) + 1) / 4) .* [b_0, b_1, b_2];    
    a = [a_0, a_1, a_2];
end

function [b, a] = design_comb(natural_frequency, sample_frequency, pole_distance_from_origin)
    b = [1];
    a = [1];

    max_frequency = sample_frequency / 2;
    frequency_ratio = natural_frequency / sample_frequency;
    notch_amount = floor(1 / frequency_ratio);

    for i = 1:notch_amount
        [dummy_b, dummy_a] = design_notch(i * natural_frequency, max_frequency, pole_distance_from_origin);
        b = conv(b, dummy_b);
        a = conv(a, dummy_a);
    end
end

% Comb Filter
pole_distance_from_origin = 0.995;
[b, a] = design_comb(power_supply_frequency, sample_frequency, pole_distance_from_origin);

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
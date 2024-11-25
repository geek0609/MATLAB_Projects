% MATLAB Program to Verify the Sampling Theorem (Nyquist-Shannon)
% Author: [Your Name]
% Date: [Current Date]

clear all; close all; clc;

%% 1. Define the Original Signal
% Parameters
f1 = 50;    % Frequency component 1 (Hz)
f2 = 120;   % Frequency component 2 (Hz)
f3 = 200;   % Frequency component 3 (Hz)
f_max = max([f1, f2, f3]);  % Maximum frequency component

A1 = 1; A2 = 0.5; A3 = 0.3;  % Amplitudes

% Continuous time parameters
T = 1;              % Total time duration (seconds)
fs_cont = 100000;   % High sampling rate to simulate continuous signal
t_cont = 0:1/fs_cont:T;  % Continuous time vector

% Original signal: sum of three sine waves
x_cont = A1*sin(2*pi*f1*t_cont) + A2*sin(2*pi*f2*t_cont) + A3*sin(2*pi*f3*t_cont);

%% 2. Sampling the Signal

% Define two sampling rates
fs1 = 2.5 * f_max;  % Sampling rate above Nyquist rate
fs2 = 1.5 * f_max;  % Sampling rate below Nyquist rate

% Ensure sampling rates are integers
fs1 = ceil(fs1);
fs2 = ceil(fs2);

% Sampling time vectors
t_sample1 = 0:1/fs1:T;
t_sample2 = 0:1/fs2:T;

% Sampled signals
x_sample1 = A1*sin(2*pi*f1*t_sample1) + A2*sin(2*pi*f2*t_sample1) + A3*sin(2*pi*f3*t_sample1);
x_sample2 = A1*sin(2*pi*f1*t_sample2) + A2*sin(2*pi*f2*t_sample2) + A3*sin(2*pi*f3*t_sample2);

%% 3. Reconstruction using Sinc Interpolation

% Reconstruction for fs1 (Above Nyquist)
x_recon1 = sinc_interp(t_sample1, x_sample1, t_cont);

% Reconstruction for fs2 (Below Nyquist)
x_recon2 = sinc_interp(t_sample2, x_sample2, t_cont);

%% 4. Visualization

figure;
subplot(3,1,1);
plot(t_cont, x_cont, 'b', 'LineWidth', 1.5);
title('Original Continuous-Time Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

subplot(3,1,2);
stem(t_sample1, x_sample1, 'r', 'Marker', 'o');
hold on;
plot(t_cont, x_cont, 'b', 'LineWidth', 1);
title(['Sampling Above Nyquist Rate: fs = ', num2str(fs1), ' Hz']);
xlabel('Time (s)');
ylabel('Amplitude');
legend('Sampled Points', 'Original Signal');
grid on;
hold off;

subplot(3,1,3);
stem(t_sample2, x_sample2, 'r', 'Marker', 'o');
hold on;
plot(t_cont, x_cont, 'b', 'LineWidth', 1);
title(['Sampling Below Nyquist Rate: fs = ', num2str(fs2), ' Hz']);
xlabel('Time (s)');
ylabel('Amplitude');
legend('Sampled Points', 'Original Signal');
grid on;
hold off;

figure;
subplot(2,1,1);
plot(t_cont, x_cont, 'b', 'LineWidth', 1.5);
hold on;
plot(t_cont, x_recon1, 'g--', 'LineWidth', 1.5);
title('Reconstruction with Sampling Rate Above Nyquist');
xlabel('Time (s)');
ylabel('Amplitude');
legend('Original Signal', 'Reconstructed Signal');
grid on;
hold off;

subplot(2,1,2);
plot(t_cont, x_cont, 'b', 'LineWidth', 1.5);
hold on;
plot(t_cont, x_recon2, 'm--', 'LineWidth', 1.5);
title('Reconstruction with Sampling Rate Below Nyquist');
xlabel('Time (s)');
ylabel('Amplitude');
legend('Original Signal', 'Reconstructed Signal');
grid on;
hold off;

%% 5. Function for Sinc Interpolation
function x_recon = sinc_interp(t_sample, x_sample, t_recon)
    % Sinc interpolation to reconstruct the signal
    % t_sample: Sampling time instances
    % x_sample: Sampled signal values
    % t_recon: Time instances for reconstruction
    
    % Initialize reconstructed signal
    x_recon = zeros(size(t_recon));
    
    % Perform interpolation
    for k = 1:length(x_sample)
        x_recon = x_recon + x_sample(k) * sinc((t_recon - t_sample(k)) * fs_estimate(t_sample));
    end
end

function fs = fs_estimate(t_sample)
    % Estimate the sampling frequency based on sampling time vector
    % Assumes uniform sampling
    if length(t_sample) < 2
        fs = 1;
    else
        fs = 1 / (t_sample(2) - t_sample(1));
    end
end

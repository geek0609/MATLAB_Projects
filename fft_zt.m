% MATLAB Program to Calculate and Plot Fourier Transform and Z-Transform
% Author: [Your Name]
% Date: [Current Date]

clear all; close all; clc;

%% 1. Define the Discrete-Time Signal
% Example Signal: A finite-length signal
n = 0:9; % Sample indices (n = 0 to 9)
x = [1, 2, 3, 4, 5, 4, 3, 2, 1, 0]; % Example signal

% Plot the original signal
figure;
stem(n, x, 'filled');
title('Original Discrete-Time Signal x[n]');
xlabel('n');
ylabel('x[n]');
grid on;

%% 2. Calculate the Fourier Transform (DFT using FFT)
N = 512; % Number of FFT points (zero-padding for better frequency resolution)
X_fft = fft(x, N);
X_fft_shifted = fftshift(X_fft); % Shift zero frequency component to center
f = (-N/2:N/2-1)*(1/N); % Normalized frequency (cycles/sample)
omega = 2*pi*f; % Angular frequency

% Plot the magnitude of the Fourier Transform
figure;
subplot(2,1,1);
plot(f, abs(X_fft_shifted));
title('Magnitude of Fourier Transform');
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('|X(e^{j\omega})|');
grid on;

% Plot the phase of the Fourier Transform
subplot(2,1,2);
plot(f, angle(X_fft_shifted));
title('Phase of Fourier Transform');
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('Phase(X(e^{j\omega})) [radians]');
grid on;

%% 3. Calculate the Z-Transform
% Using Symbolic Toolbox for Z-Transform
syms z n_sym;
x_sym = sym(x); % Convert x to symbolic

% Compute Z-Transform: X(z) = sum_{n=0}^{N-1} x[n] z^{-n}
X_z = 0;
for k = 1:length(x)
    X_z = X_z + x(k)*z^(-k+1); % Adjusting index since MATLAB indexing starts at 1
end
X_z = simplify(X_z);

% Display the symbolic Z-transform
disp('Symbolic Z-Transform X(z):');
disp(X_z);

% Evaluate Z-Transform numerically over a grid
% Define grid for z = r * e^{j\theta}
r = linspace(0.5, 1.5, 100); % Radius from 0.5 to 1.5
theta = linspace(0, 2*pi, 200); % Angle from 0 to 2*pi
[R, Theta] = meshgrid(r, theta);
Z = R .* exp(1j*Theta); % Complex grid

% Substitute z with the grid values
% For efficiency, compute X(z) numerically
% X(z) = sum x[n] z^{-n}
X_z_num = zeros(size(Z));
for k = 1:length(x)
    X_z_num = X_z_num + x(k) .* Z.^(-(k-1));
end

% Compute magnitude and phase
X_z_mag = abs(X_z_num);
X_z_phase = angle(X_z_num);

%% 4. Visualization of Z-Transform
% Plot the magnitude of Z-Transform
figure;
surf(R .* cos(Theta), R .* sin(Theta), X_z_mag, 'EdgeColor', 'none');
title('Magnitude of Z-Transform |X(z)|');
xlabel('Real(z)');
ylabel('Imag(z)');
zlabel('|X(z)|');
colorbar;
view(2); % Top view
grid on;

% Plot the phase of Z-Transform
figure;
surf(R .* cos(Theta), R .* sin(Theta), X_z_phase, 'EdgeColor', 'none');
title('Phase of Z-Transform \angleX(z)');
xlabel('Real(z)');
ylabel('Imag(z)');
zlabel('Phase(X(z)) [radians]');
colorbar;
view(2); % Top view
grid on;

%% 5. Comparison Between Fourier Transform and Z-Transform
% Fourier Transform is the Z-Transform evaluated on the unit circle (r = 1)
% Extract X(z) on the unit circle
Z_unit = exp(1j*Theta(1, :)); % Only one theta value since r=1
X_z_unit = zeros(size(Z_unit));
for k = 1:length(x)
    X_z_unit = X_z_unit + x(k) .* Z_unit.^(-(k-1));
end

% Fourier Transform from Z-Transform on unit circle
figure;
subplot(2,1,1);
plot(theta(1, :), abs(X_z_unit));
title('Magnitude of Z-Transform on Unit Circle |X(e^{j\omega})|');
xlabel('Frequency \omega [rad/sample]');
ylabel('|X(e^{j\omega})|');
grid on;

subplot(2,1,2);
plot(theta(1, :), angle(X_z_unit));
title('Phase of Z-Transform on Unit Circle \angleX(e^{j\omega})');
xlabel('Frequency \omega [rad/sample]');
ylabel('Phase(X(e^{j\omega})) [radians]');
grid on;

% Compare with FFT results
% Since FFT was shifted, adjust theta accordingly
figure;
subplot(2,1,1);
plot(theta(1, :), abs(X_z_unit));
hold on;
plot(omega, abs(X_fft_shifted), '--');
title('Comparison of Fourier Transform and Z-Transform on Unit Circle');
xlabel('Frequency \omega [rad/sample]');
ylabel('|X(e^{j\omega})|');
legend('Z-Transform on Unit Circle', 'FFT Magnitude');
grid on;
hold off;

subplot(2,1,2);
plot(theta(1, :), angle(X_z_unit));
hold on;
plot(omega, angle(X_fft_shifted), '--');
title('Comparison of Fourier Transform and Z-Transform on Unit Circle');
xlabel('Frequency \omega [rad/sample]');
ylabel('Phase(X(e^{j\omega})) [radians]');
legend('Z-Transform on Unit Circle', 'FFT Phase');
grid on;
hold off;

%% 6. Function to Compute Z-Transform Numerically (Optional)
% If you prefer to compute Z-Transform without symbolic toolbox, use this function.
% However, symbolic approach is more straightforward for finite-length signals.

% Uncomment the following block to use the numerical Z-Transform function.

%{
function X_z = compute_Z_transform(x, z)
    % Computes the Z-Transform of a finite-length signal x at points z
    % x: input signal vector
    % z: complex points where Z-transform is evaluated
    % X_z: Z-transform evaluated at z
    
    N = length(x);
    X_z = zeros(size(z));
    for k = 1:N
        X_z = X_z + x(k) * z.^(-(k-1));
    end
end
%}


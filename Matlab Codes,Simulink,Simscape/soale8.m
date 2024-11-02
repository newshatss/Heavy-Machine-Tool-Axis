% clc;
% clear all;
% 
% % Define the given constants
% J_m = 0.0001;
% J_l = 0.001;
% K = 1;
% N = 5;
% C_m = 0.01;
% C_l = 0.01;
% 
% % Define the numerator and denominator coefficients
% numerator = K / N;
% denominator = [J_l * J_m, (J_l * C_m + C_l * J_m), ...
%     (J_l * K / N^2 + C_l * C_m + K * J_m), ...
%     (C_l * K / N^2 + K * C_m), 0];
% 
% % Create the open-loop transfer function
% G = tf(numerator, denominator);
% 
% % Define the closed-loop transfer function with unity feedback
% H = 1;
% L = G * H;
% 
% % Sensitivity function S = 1 / (1 + L)
% S = 1 / (1 + L);
% 
% % Complementary sensitivity function T = L / (1 + L)
% T = L / (1 + L);
% 
% % Bode plot data
% [mag_S, phase_S, w_S] = bode(S);
% [mag_T, phase_T, w_T] = bode(T);
% 
% % Convert magnitudes from dB to absolute
% mag_S = squeeze(mag_S);
% mag_T = squeeze(mag_T);
% 
% % Find frequency ranges
% freq_range_T = w_T(mag_T < 1);
% freq_range_S = w_S(mag_S < 1);
% 
% % Display the results
% fprintf('Frequency range where |T(jw)| < 1 (resistant to model changes):\n');
% fprintf('%.2f rad/s to %.2f rad/s\n', min(freq_range_T), max(freq_range_T));
% 
% fprintf('\nFrequency range where |S(jw)| < 1 (resistant to input disturbances):\n');
% fprintf('%.2f rad/s to %.2f rad/s\n', min(freq_range_S), max(freq_range_S));
% 
% % Plot the Bode plot for Sensitivity function S
% figure;
% bode(S);
% title('Bode Plot of Sensitivity Function S');
% grid on;
% 
% % Plot the Bode plot for Complementary Sensitivity function T
% figure;
% bode(T);
% title('Bode Plot of Complementary Sensitivity Function T');
% grid on;
%%
clc;
clear all;

% Define the given constants
J_m = 0.0001;
J_l = 0.001;
K = 1;
N = 5;
C_m = 0.01;
C_l = 0.01;

% Define the numerator and denominator coefficients
numerator = K / N;
denominator = [J_l * J_m, (J_l * C_m + C_l * J_m), ...
    (J_l * K / N^2 + C_l * C_m + K * J_m), ...
    (C_l * K / N^2 + K * C_m), 0];

% Create the open-loop transfer function
G = tf(numerator, denominator);

% Define the closed-loop transfer function with unity feedback
H = 1;
L = G * H;

% Sensitivity function S = 1 / (1 + L)
S = 1 / (1 + L);

% Complementary sensitivity function T = L / (1 + L)
T = L / (1 + L);

% Bode plot data
[mag_S, phase_S, w_S] = bode(S);
[mag_T, phase_T, w_T] = bode(T);

% Convert magnitudes from dB to absolute
mag_S = squeeze(mag_S);
mag_T = squeeze(mag_T);

% Find frequency ranges
freq_range_T = w_T(mag_T < 1);
freq_range_S = w_S(mag_S < 1);

% Determine the overlapping range
min_freq = max(min(freq_range_T), min(freq_range_S));
max_freq = min(max(freq_range_T), max(freq_range_S));

% Display the results
fprintf('Frequency range where |T(jw)| < 1 (resistant to model changes):\n');
fprintf('%.2f rad/s to %.2f rad/s\n', min(freq_range_T), max(freq_range_T));

fprintf('\nFrequency range where |S(jw)| < 1 (resistant to input disturbances):\n');
fprintf('%.2f rad/s to %.2f rad/s\n', min(freq_range_S), max(freq_range_S));

fprintf('\nOverlapping frequency range where both |S(jw)| and |T(jw)| < 1 (robustness range):\n');
fprintf('%.2f rad/s to %.2f rad/s\n', min_freq, max_freq);

% Plot the Bode plot for Sensitivity function S
figure;
bode(S);
title('Bode Plot of Sensitivity Function S');
grid on;

% Plot the Bode plot for Complementary Sensitivity function T
figure;
bode(T);
title('Bode Plot of Complementary Sensitivity Function T');
grid on;

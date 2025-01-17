clc;
clear all;
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

% Create the sensitivity function
S1 = -1 / (1 + G);

% Plot Bode diagram of G and S in one plot
figure;
bode(G);
hold on;
bode(S1);
legend('G', 'S1');
grid on;
title('bode for disturbance rejection')

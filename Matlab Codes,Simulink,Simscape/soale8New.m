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

% Create the closed-loop transfer function
T = feedback(G, 1);

% Plot Bode diagram of G and T in one plot
figure;
bode(G);
hold on;
bode(T);
legend('G', 'T');
grid on;
title ('bode for command tracking')
%%


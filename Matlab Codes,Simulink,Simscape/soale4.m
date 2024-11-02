% % Define the transfer function
% numerator = 0.2;
% denominator = [0.0000001, 0.000011, 0.00024, 0.0104, 0];
% G = tf(numerator, denominator);
% sisotool(G)
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

% Use pidtune to design a PD controller
desired_tf = 2.5;  % Desired natural frequency
damping_ratio = 0.8;  % Desired damping ratio

% Design the PD controller using pidtune
C = pidtune(G, 'PD', desired_tf / damping_ratio);

% Display the designed PD controller
disp('The designed PD controller is:');
disp(C);

% Plot the step response of the closed-loop system
closed_loop_sys = feedback(C * G, 1);
figure;
step(closed_loop_sys);
title('Step Response of the Closed-Loop System with PID Controller');
grid on;

% Calculate and display performance metrics
info = stepinfo(closed_loop_sys);
rise_time = info.RiseTime;
settling_time = info.SettlingTime;
overshoot = info.Overshoot;

fprintf('Rise Time (t_r): %.4f seconds\n', rise_time);
fprintf('Settling Time (t_s): %.4f seconds\n', settling_time);
fprintf('Overshoot: %.2f%%\n', overshoot);

% Define a sinusoidal input
t = 0:0.01:10;  % time vector from 0 to 10 seconds
f = 1;  % frequency of the sinusoidal input
input = sin(2 * pi * f * t);  % sinusoidal input

% Simulate the closed-loop response to the sinusoidal input
output = lsim(closed_loop_sys, input, t);

% Plot the sinusoidal input and the system's response
figure;
plot(t, input, 'r--', 'DisplayName', 'Sinusoidal Input');
hold on;
plot(t, output, 'b', 'DisplayName', 'System Response');
xlabel('Time (seconds)');
ylabel('Amplitude');
title('Closed-Loop Response to Sinusoidal Input');
legend;
grid on;

% Convert the closed-loop transfer function to a state-space model
closed_loop_sys_ss = ss(closed_loop_sys);

% Define initial conditions
initial_conditions = [1; 0; 0; 0];  % Example initial conditions (adjust as needed)

% Simulate the closed-loop response to initial conditions
figure;
initial(closed_loop_sys_ss, initial_conditions);
title('Closed-Loop Response to Initial Conditions');
grid on;

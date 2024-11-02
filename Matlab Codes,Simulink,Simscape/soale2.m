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

% Create the transfer function
open_loop_tf = tf(numerator, denominator);

% Define the time vector for a smoother response
t = 0:0.0001:10; % Adjust the time vector for better resolution and longer duration

% Plot the unit step response
figure;
step(open_loop_tf, t);
title('Unit Step Response of the Open-Loop System');
xlabel('Time (seconds)');
ylabel('Response');
grid on;

% Plot the unit impulse response
figure;
impulse(open_loop_tf, t);
title('Unit Impulse Response of the Open-Loop System');
xlabel('Time (seconds)');
ylabel('Response');
grid on;

% Get the open-loop system specifications
poles = pole(open_loop_tf);
time_constants = -1 ./ real(poles);


% Display the poles, time constants, and settling time
disp('Open-Loop System Poles:');
disp(poles);
disp('Open-Loop System Time Constants:');
disp(time_constants);

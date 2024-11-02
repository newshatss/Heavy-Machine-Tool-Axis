% Given Transfer Function
num = 0.2;
den = [1e-07 1.1e-05 0.00024 0.0104 0];
G = tf(num, den);

% Desired closed-loop poles for a 4th-order ITAE polynomial
n = 4; % order of the ITAE filter

% Predefined ITAE polynomial coefficients for n = 4 (standard form)
% The coefficients are typically given for normalized ITAE polynomials
itae_poly_coeffs = [1 3.8 6.1 4.6 1]; % ITAE polynomial coefficients for order 4

% Scaling the polynomial to the desired cut-off frequency
wc = 1; % cut-off frequency (can be adjusted as per requirements)
scaled_itae_poly_coeffs = itae_poly_coeffs * (wc^(n));

% Define the desired closed-loop characteristic polynomial
char_poly = scaled_itae_poly_coeffs(end:-1:1);

% Create the desired closed-loop transfer function
num_d = char_poly(end); % numerator of the desired closed-loop TF
den_d = char_poly; % denominator of the desired closed-loop TF
H_desired = tf(num_d, den_d);

% Determine the controller transfer function
% Open-loop transfer function: L(s) = G(s) * C(s)
% Closed-loop transfer function: T(s) = L(s) / (1 + L(s))

% We need to find C(s) such that T(s) = H_desired

% Let L(s) = G(s) * C(s) = H_desired / (1 - H_desired)
L_desired = H_desired / (1 - H_desired);

% Since L(s) = G(s) * C(s), we can solve for C(s)
C = L_desired / G;

% Display the controller transfer function
disp('The designed controller transfer function C(s) is:');
C

% Verify the closed-loop system with the designed controller
T = feedback(G * C, 1);
disp('The closed-loop transfer function T(s) with the designed controller is:');
T

% Plot the step response of the closed-loop system
figure;
step(T);
title('Step Response of the Closed-Loop System with Designed ITAE Controller');
grid on;

% Obtain step response information
info = stepinfo(T);

% Display overshoot, settling time, and rise time
disp('Step Response Characteristics:');
fprintf('Overshoot (O.S.): %.2f%%\n', info.Overshoot);
fprintf('Settling Time (Ts): %.2f seconds\n', info.SettlingTime);
fprintf('Rise Time (Tr): %.2f seconds\n', info.RiseTime);

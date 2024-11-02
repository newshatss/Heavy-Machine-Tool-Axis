% 
% %%
% clc;
% clear all;
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
% % Initialize an array to store the controllers for each w_cd
% controllers = [];
% legends = cell(1, 10); % Cell array to store legend entries
% 
% % Loop through w_cd from 1 to 10
% for w_cd = 1:1:10
%     % Define Td and Sd
%     s = tf('s');
%     Td = w_cd^4 / (s + w_cd)^4;
%     Sd = 1 - Td;
%     
%     % Calculate the controller C
%     C = Td / (Sd * G);
%     
%     % Store the controller in the array
%     controllers = [controllers; C];
%     
%     % Store legend entry
%     legends{w_cd} = ['w_{cd} = ', num2str(w_cd)];
% end
% 
% % Plot the closed-loop step responses
% figure;
% hold on;
% 
% for w_cd = 1:1:10
%     % Get the controller for the current w_cd
%     C = controllers(w_cd);
%     
%     % Define the closed-loop transfer function with unit feedback
%     closed_loop_tf = feedback(C * G, 1);
%     
%     % Plot the step response
%     step(closed_loop_tf);
% end
% 
% % Add legend to the plot
% legend(legends);
% 
% % Add title and labels
% title('Closed-Loop Step Responses for Different w_{cd} Values');
% xlabel('Time (seconds)');
% ylabel('Amplitude');
% 
% 
% 
% %%
% clc;
% clear all;
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
% % Initialize arrays to store the controllers and performance metrics
% controllers = [];
% legends = cell(1, 10); % Cell array to store legend entries
% overshoots = zeros(1, 10);
% rise_times = zeros(1, 10);
% settling_times = zeros(1, 10);
% 
% % Loop through w_cd from 1 to 10
% for w_cd = 1:1:10
%     % Define Td and Sd
%     s = tf('s');
%     Td = w_cd^4 / (s + w_cd)^4;
%     Sd = 1 - Td;
%     
%     % Calculate the controller C
%     C = Td / (Sd * G);
%     
%     % Store the controller in the array
%     controllers = [controllers; C];
%     
%     % Store legend entry
%     legends{w_cd} = ['w_{cd} = ', num2str(w_cd)];
% end
% 
% % Plot the closed-loop step responses
% figure;
% hold on;
% 
% for w_cd = 1:1:10
%     % Get the controller for the current w_cd
%     C = controllers(w_cd);
%     
%     % Define the closed-loop transfer function with unit feedback
%     closed_loop_tf = feedback(C * G, 1);
%     
%     % Plot the step response
%     [y, t] = step(closed_loop_tf);
%     plot(t, y);
%     
%     % Calculate overshoot, rise time, and settling time
%     S = stepinfo(closed_loop_tf);
%     overshoots(w_cd) = S.Overshoot;
%     rise_times(w_cd) = S.RiseTime;
%     settling_times(w_cd) = S.SettlingTime;
% end
% 
% % Add legend to the plot
% legend(legends);
% 
% % Add title and labels
% title('Closed-Loop Step Responses for Different w_{cd} Values');
% xlabel('Time (seconds)');
% ylabel('Amplitude');
% 
% % Display overshoot, rise time, and settling time for each plot
% disp('Performance Metrics for each w_{cd}:');
% for w_cd = 1:1:10
%     fprintf('w_{cd} = %d: Overshoot = %.2f%%, Rise Time = %.2f s, Settling Time = %.2f s\n', ...
%         w_cd, overshoots(w_cd), rise_times(w_cd), settling_times(w_cd));
% end

%
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
 s = tf('s');
    Td = 5^4 / (s + 5)^4;
    Sd = 1 - Td;
   
    % Calculate the controller C
    C = Td / (Sd * G);
    % Plot the step response of the closed-loop system
closed_loop_sys = feedback(C * G, 1);
figure;
step(closed_loop_sys);
title('Step Response of the Closed-Loop System with PD Controller');
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

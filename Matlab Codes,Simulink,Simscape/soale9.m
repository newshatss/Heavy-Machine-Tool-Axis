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

% Define the transfer function of the system
sys = tf(numerator, denominator);

% Design Butterworth filters
cutoff_frequency = 0.2;  % Normalized cutoff frequency (0 to 1)
[b1, a1] = butter(1, cutoff_frequency, 'low');  % 1st order Butterworth filter
[b2, a2] = butter(2, cutoff_frequency, 'low');  % 2nd order Butterworth filter

% Define transfer functions for the Butterworth filters
filter1 = tf(b1, a1);  % 1st Order Butterworth Filter
filter2 = tf(b2, a2);  % 2nd Order Butterworth Filter

% Define the closed-loop system with no filter (for comparison)
closed_loop_no_filter = feedback(sys, 1);

% Define the closed-loop system with 1st Order Butterworth Filter in the feedback path
sys_feedback_1 = feedback(sys * filter1, 1);

% Define the closed-loop system with 2nd Order Butterworth Filter in the feedback path
sys_feedback_2 = feedback(sys * filter2, 1);

% Generate a sample input signal
Fs = 1000;  % Sampling frequency
t = 0:1/Fs:1-1/Fs;  % Time vector
input_signal = sin(2*pi*7*t) + sin(2*pi*13*t) + 0.5*randn(size(t));  % Noisy input signal

% Simulate the system response to the input signal
output_no_filter = lsim(closed_loop_no_filter, input_signal, t);
output_filter1 = lsim(sys_feedback_1, input_signal, t);
output_filter2 = lsim(sys_feedback_2, input_signal, t);

% Plot the original and filtered output signals
figure;

subplot(4, 1, 1);
plot(t, input_signal, 'k');
title('Input Signal (Noisy)');
xlabel('Time (seconds)');
ylabel('Amplitude');
grid on;

subplot(4, 1, 2);
plot(t, output_no_filter, 'b');
title('Output Signal with No Filter in Feedback Path');
xlabel('Time (seconds)');
ylabel('Amplitude');
grid on;

subplot(4, 1, 3);
plot(t, output_filter1, 'r');
title('Output Signal with 1st Order Butterworth Filter in Feedback Path');
xlabel('Time (seconds)');
ylabel('Amplitude');
grid on;

subplot(4, 1, 4);
plot(t, output_filter2, 'g');
title('Output Signal with 2nd Order Butterworth Filter in Feedback Path');
xlabel('Time (seconds)');
ylabel('Amplitude');
grid on;

% Display Bode plot of the original system for reference
figure;
bode(sys);
title('Bode Plot of the Original System');
grid on;

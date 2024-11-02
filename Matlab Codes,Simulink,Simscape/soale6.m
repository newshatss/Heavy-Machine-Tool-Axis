clc;
clear all;
% Define the open-loop transfer function
numerator = 0.2;
denominator = [1e-07, 1.1e-05, 0.00024, 0.0104, 0];
G = tf(numerator, denominator);

% Define different gain values to analyze
gains = [0.05 , 0.08 , 0.1, 0.2, 0.5];

% Plot Nichols diagrams for different gains
figure;
hold on;
for k = gains
    G_k = k * G;
    nichols(G_k);
end
hold off;
title('Nichols Diagram for Different Gains');
legend(arrayfun(@(k) sprintf('Gain = %d', k), gains, 'UniformOutput', false));
grid on;

% Analyze phase margin and gain margin for different gains
for k = gains
    G_k = k * G;
    [Gm, Pm, Wcg, Wcp] = margin(G_k);
    fprintf('Gain = %d:\n', k);
    fprintf('  Gain Margin (dB): %f\n', 20*log10(Gm));
    fprintf('  Phase Margin (deg): %f\n', Pm);
    fprintf('  Gain Crossover Frequency (rad/s): %f\n', Wcg);
    fprintf('  Phase Crossover Frequency (rad/s): %f\n', Wcp);
    fprintf('\n');
end

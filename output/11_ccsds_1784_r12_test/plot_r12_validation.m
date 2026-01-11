% Plot R=1/2 Validation Results vs NASA Baseline
% Comparing: Original (8-bit tail) vs Modified (16-bit tail)
clear; clc; close all;

% 1. Load Simulation Data - Original (8-bit punctured tail)
data_old = readtable('ber_Turbo_LogMAP_K1784_R1-2_20260111_143738_merged.csv');
snr_old = data_old.Eb_N0;
fer_old = data_old.FER;

% 2. Load Simulation Data - Modified (16-bit full tail)
data_new = readtable('ber_Turbo_LogMAP_K1784_R1-2_20260111_165513_merged.csv');
snr_new = data_new.Eb_N0;
fer_new = data_new.FER;

% 3. Define NASA Baseline Data (Table A-3, 1784 bits, Rate 1/2)
nasa_snr = [0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2];
nasa_fer = [7.5000e-01, 3.8931e-01, 7.5529e-02, 7.9605e-03, ...
    3.2503e-04, 1.1620e-05, 3.7500e-06, 2.2500e-06, 7.5000e-07];

% 4. Plotting
figure('Name', 'CCSDS Turbo R=1/2 Validation', 'NumberTitle', 'off', ...
    'Position', [100, 100, 900, 600]);

% Original implementation (8-bit tail)
semilogy(snr_old, fer_old, 'b-o', 'LineWidth', 2, 'MarkerSize', 6, ...
    'DisplayName', 'Original (8-bit tail)');
hold on;

% Modified implementation (16-bit tail)
semilogy(snr_new, fer_new, 'g-^', 'LineWidth', 2, 'MarkerSize', 8, ...
    'DisplayName', 'Modified (16-bit tail)');

% NASA Baseline
semilogy(nasa_snr, nasa_fer, 'r--s', 'LineWidth', 2, 'MarkerSize', 8, ...
    'DisplayName', 'NASA Baseline (Table A-3)');

% 5. Formatting
grid on;
xlabel('E_b/N_0 (dB)', 'FontSize', 12, 'Interpreter', 'tex');
ylabel('Frame Error Rate (FER)', 'FontSize', 12);
title('CCSDS Turbo Code R=1/2 Performance Comparison (K=1784)', 'FontSize', 14);
legend('Location', 'SouthWest', 'FontSize', 10);
axis([-0.5 2.5 1e-5 1]);
set(gca, 'FontSize', 10);

% 6. Add annotation for the improvement
text(1.5, 0.1, sprintf('Modified: Full 16-bit tail\n(like R=1/3)'), ...
    'FontSize', 9, 'Color', [0 0.5 0], 'BackgroundColor', 'white');

% 7. Save
saveas(gcf, 'plot_r12_comparison.png');
disp('Plot saved as plot_r12_comparison.png');

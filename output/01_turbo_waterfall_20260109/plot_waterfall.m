%% Turbo Code Waterfall Region Plot
% 瀑布区数据（-1.0 ~ 4.0 dB）
% 数据来源：ber_turbo_20260109_213911.csv

clear; close all; clc;

%% === 读取 CSV 数据 ===
data_file = 'ber_turbo_20260109_213911.csv';
fprintf('正在读取: %s\n', data_file);
data = readtable(data_file, 'CommentStyle', '#');

snr_db = data.SNR_dB;
ber = data.BER;

fprintf('数据范围: SNR = [%.1f, %.1f] dB, 共 %d 个点\n', min(snr_db), max(snr_db), length(snr_db));

%% === 理论曲线 ===
EbN0_dB_theory = -2:0.2:6;
BER_uncoded = 0.5 * erfc(sqrt(10.^(EbN0_dB_theory/10)));

%% === 绘图 ===
figure('Position', [100 100 900 650], 'Color', 'white');

semilogy(EbN0_dB_theory, BER_uncoded, '-.', 'Color', [0.6 0.6 0.6], ...
    'LineWidth', 1.5, 'DisplayName', 'Uncoded BPSK (Theory)');
hold on;

semilogy(snr_db, ber, '-m', 'LineWidth', 2, 'MarkerSize', 7, ...
    'Marker', 'd', 'MarkerFaceColor', 'm', ...
    'DisplayName', 'Turbo Code (R=1/3, 8 iter)');

xlabel('$E_b/N_0$ (dB)', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Bit Error Rate (BER)', 'FontSize', 14, 'Interpreter', 'latex');
title({'Turbo Code Waterfall Region', ...
    'PCCC $R=1/3$, $k=1024$ bits, 100k frames/point'}, ...
    'FontSize', 14, 'Interpreter', 'latex');
legend('Location', 'southwest', 'FontSize', 11);
grid on;
xlim([-2 5]);
ylim([1e-8 1]);
set(gca, 'FontSize', 12);

saveas(gcf, 'plot_waterfall.png');
fprintf('已保存: plot_waterfall.png\n');

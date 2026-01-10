%% Turbo Code Error Floor Region Plot
% 错误平层数据（4.5 ~ 10.0 dB）
% 数据来源：ber_Turbo_LogMAP_20260110_121358_merged.csv

clear; close all; clc;

%% === 读取 CSV 数据 ===
data_file = 'ber_Turbo_LogMAP_20260110_121358_merged.csv';
fprintf('正在读取: %s\n', data_file);
data = readtable(data_file, 'CommentStyle', '#');

snr_db = data.SNR_dB;
ber = data.BER;

% 过滤 BER=0 的点
valid_idx = ber > 0;
snr_db = snr_db(valid_idx);
ber = ber(valid_idx);

fprintf('有效数据: SNR = [%.1f, %.1f] dB, 共 %d 个点\n', min(snr_db), max(snr_db), length(snr_db));

%% === 绘图 ===
figure('Position', [100 100 900 650], 'Color', 'white');

semilogy(snr_db, ber, '-b', 'LineWidth', 2, 'MarkerSize', 8, ...
    'Marker', 's', 'MarkerFaceColor', 'b', ...
    'DisplayName', 'Turbo Code Error Floor');
hold on;

yline(1e-7, '--r', 'BER = 10^{-7}', 'LineWidth', 1.5);

xlabel('$E_b/N_0$ (dB)', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Bit Error Rate (BER)', 'FontSize', 14, 'Interpreter', 'latex');
title({'Turbo Code Error Floor Region', ...
    'PCCC $R=1/3$, $k=1024$ bits, 1M frames/point'}, ...
    'FontSize', 14, 'Interpreter', 'latex');
legend('Location', 'southwest', 'FontSize', 11);
grid on;
xlim([4 6]);
ylim([1e-9 1e-6]);
set(gca, 'FontSize', 12);

saveas(gcf, 'plot_errorfloor.png');
fprintf('已保存: plot_errorfloor.png\n');

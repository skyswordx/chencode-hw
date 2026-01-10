%% Turbo Code Complete BER Performance
% 完整曲线（瀑布区 + 错误平层）
% 数据来源：ber_turbo_complete.csv

clear; close all; clc;

%% === 读取 CSV 数据 ===
data_file = 'ber_turbo_complete.csv';
fprintf('正在读取: %s\n', data_file);
data = readtable(data_file, 'CommentStyle', '#');

snr_db = data.SNR_dB;
ber = data.BER;

valid_idx = ber > 0;
snr_db = snr_db(valid_idx);
ber = ber(valid_idx);

fprintf('数据范围: SNR = [%.1f, %.1f] dB, 共 %d 个点\n', min(snr_db), max(snr_db), length(snr_db));

%% === 理论曲线 ===
EbN0_dB_theory = -2:0.2:8;
BER_uncoded = 0.5 * erfc(sqrt(10.^(EbN0_dB_theory/10)));

%% === 绘图 ===
figure('Position', [100 100 1000 700], 'Color', 'white');

semilogy(EbN0_dB_theory, BER_uncoded, '-.', 'Color', [0.6 0.6 0.6], ...
    'LineWidth', 1.5, 'DisplayName', 'Uncoded BPSK (Theory)');
hold on;

semilogy(snr_db, ber, '-', 'LineWidth', 2.5, 'Color', [0.8 0.2 0.6], ...
    'Marker', 'd', 'MarkerSize', 7, 'MarkerFaceColor', [0.8 0.2 0.6], ...
    'DisplayName', 'Turbo Code (R=1/3, 8 iter)');

% 标注区域
text(0.5, 1e-2, 'Waterfall', 'FontSize', 12, 'FontWeight', 'bold', 'Color', [0.2 0.4 0.8], 'Rotation', -55);
text(4.5, 3e-8, 'Error Floor', 'FontSize', 12, 'FontWeight', 'bold', 'Color', [0.8 0.2 0.2]);
yline(1e-7, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);

%% === 编码增益 ===
target_ber = 1e-5;
EbN0_uncoded = fzero(@(x) 0.5*erfc(sqrt(10^(x/10))) - target_ber, 9);
if min(ber) < target_ber && max(ber) > target_ber
    EbN0_turbo = interp1(log10(ber), snr_db, log10(target_ber));
    coding_gain = EbN0_uncoded - EbN0_turbo;
    text(-1, 1e-8, sprintf('Coding Gain: %.1f dB @ BER=10^{-5}', coding_gain), ...
        'FontSize', 11, 'FontWeight', 'bold', 'Color', [0.1 0.5 0.1], ...
        'BackgroundColor', 'w', 'EdgeColor', [0.1 0.5 0.1]);
end

xlabel('$E_b/N_0$ (dB)', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Bit Error Rate (BER)', 'FontSize', 14, 'Interpreter', 'latex');
title({'Turbo Code Complete BER Performance', ...
    'PCCC $R=1/3$, $k=1024$ bits, 8 iterations, AWGN'}, 'FontSize', 14, 'Interpreter', 'latex');
legend('Location', 'southwest', 'FontSize', 11);
grid on;
xlim([-2 7]);
ylim([1e-9 1]);
set(gca, 'FontSize', 12);

saveas(gcf, 'plot_complete.png');
print(gcf, 'plot_complete_hires.png', '-dpng', '-r300');
fprintf('已保存: plot_complete.png, plot_complete_hires.png\n');

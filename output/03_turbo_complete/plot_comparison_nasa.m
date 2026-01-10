%% Turbo Code: Our Implementation vs ALL NASA Baselines
% 对比你的 1024-bit 实现与所有 NASA 标准 block sizes

clear; close all; clc;

%% === 读取你的仿真数据 ===
your_data = readtable('ber_turbo_complete.csv', 'CommentStyle', '#');
your_snr = your_data.SNR_dB;
your_ber = your_data.BER;
your_ber(your_ber == 0) = NaN;

%% === NASA 基准数据 (硬编码避免读取问题) ===
% 1784 bits
nasa_1784_snr = [-0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0];
nasa_1784_ber = [9.9020e-01, 9.0090e-01, 6.8493e-01, 2.9762e-01, 4.7174e-02, 4.4583e-03, 9.2350e-05, 1.9100e-06];

% 3568 bits
nasa_3568_snr = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7];
nasa_3568_ber = [1.8065e-01, 5.1557e-02, 9.0463e-03, 9.5734e-04, 4.1120e-05, 4.5100e-06];

% 7136 bits
nasa_7136_snr = [0.1, 0.2, 0.3, 0.4, 0.5];
nasa_7136_ber = [4.3328e-01, 1.0761e-01, 1.0989e-02, 4.4099e-04, 1.0050e-05];

% 8920 bits
nasa_8920_snr = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7];
nasa_8920_ber = [8.3333e-01, 4.9505e-01, 9.7752e-02, 8.9847e-03, 2.0755e-04, 2.8730e-05, 1.4360e-05, 1.1490e-05];

%% === 理论曲线 ===
EbN0_theory = -1:0.1:3;
BER_uncoded = 0.5 * erfc(sqrt(10.^(EbN0_theory/10)));

%% === 绘图 ===
figure('Position', [50 50 1100 750], 'Color', 'white');

% Uncoded BPSK
semilogy(EbN0_theory, BER_uncoded, '-.', 'Color', [0.75 0.75 0.75], ...
    'LineWidth', 1.5, 'DisplayName', 'Uncoded BPSK');
hold on;

% NASA 1784 bits - 蓝色
semilogy(nasa_1784_snr, nasa_1784_ber, '-s', 'Color', [0.0 0.4 0.8], ...
    'LineWidth', 1.8, 'MarkerSize', 8, 'MarkerFaceColor', [0.0 0.4 0.8], ...
    'DisplayName', 'NASA 1784 bits');

% NASA 3568 bits - 绿色
semilogy(nasa_3568_snr, nasa_3568_ber, '-^', 'Color', [0.0 0.7 0.3], ...
    'LineWidth', 1.8, 'MarkerSize', 8, 'MarkerFaceColor', [0.0 0.7 0.3], ...
    'DisplayName', 'NASA 3568 bits');

% NASA 7136 bits - 橙色
semilogy(nasa_7136_snr, nasa_7136_ber, '-v', 'Color', [0.9 0.5 0.0], ...
    'LineWidth', 1.8, 'MarkerSize', 8, 'MarkerFaceColor', [0.9 0.5 0.0], ...
    'DisplayName', 'NASA 7136 bits');

% NASA 8920 bits - 紫色
semilogy(nasa_8920_snr, nasa_8920_ber, '-o', 'Color', [0.6 0.0 0.6], ...
    'LineWidth', 1.8, 'MarkerSize', 8, 'MarkerFaceColor', [0.6 0.0 0.6], ...
    'DisplayName', 'NASA 8920 bits');

% 你的实现 - 红色突出显示
semilogy(your_snr, your_ber, '-d', 'Color', [0.9 0.1 0.3], ...
    'LineWidth', 2.5, 'MarkerSize', 6, 'MarkerFaceColor', [0.9 0.1 0.3], ...
    'DisplayName', 'Our Impl. (1024 bits)');

%% === 图表美化 ===
xlabel('$E_b/N_0$ (dB)', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Bit Error Rate (BER)', 'FontSize', 14, 'Interpreter', 'latex');
title({'Turbo Code Performance: Our Implementation vs NASA CCSDS Baselines', ...
    'Rate 1/3, AWGN Channel'}, 'FontSize', 13, 'Interpreter', 'latex');

legend('Location', 'southwest', 'FontSize', 10);
grid on;
xlim([-0.5 2]);
ylim([1e-7 1]);
set(gca, 'FontSize', 12);

% Block size 趋势注释
annotation('textbox', [0.55 0.70 0.33 0.18], ...
    'String', {'Block Size Effect:', ...
    '  Larger block = Steeper waterfall', ...
    '  Larger block = Lower error floor', ...
    '', ...
    'Our: 1024 bits (smallest)', ...
    'NASA: 1784 - 8920 bits'}, ...
    'FontSize', 10, 'BackgroundColor', [1 1 0.95], 'EdgeColor', [0.3 0.3 0.3]);

%% === 保存 ===
saveas(gcf, 'comparison_nasa_all.png');
print(gcf, 'comparison_nasa_all_hires.png', '-dpng', '-r300');
fprintf('已保存: comparison_nasa_all.png 和 comparison_nasa_all_hires.png\n');

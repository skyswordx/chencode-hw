%% CCSDS Turbo Code K=8920 Validation
% 对比仿真结果与 NASA Table A-4 (K=8920) 标准数据

clear; close all; clc;

%% === 读取仿真数据 ===
csv_file = '../ber_Turbo_LogMAP_20260110_180444_merged.csv';
if ~exist(csv_file, 'file')
    error('Could not find CSV file: %s', csv_file);
end

data = readtable(csv_file, 'CommentStyle', '#');
fprintf('Reading data from: %s\n', csv_file);

% 提取数据列
my_snr = data{:, 1};   % Eb/N0 (dB)
my_ber = data{:, 4};   % BER
my_fer = data{:, 7};   % FER

% 过滤掉 0 值以便对数绘图
my_ber(my_ber == 0) = NaN;
my_fer(my_fer == 0) = NaN;

%% === NASA 基准数据 (Table A-4, K=8920, R=1/3) ===
nasa_8920_snr = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7];
nasa_8920_fer = [8.3333e-01, 4.9505e-01, 9.7752e-02, 8.9847e-03, 2.0755e-04, 2.8730e-05, 1.4360e-05, 1.1490e-05];

%% === 理论参考曲线 ===
EbN0_theory = -0.5:0.05:1.5;
BER_uncoded = 0.5 * erfc(sqrt(10.^(EbN0_theory/10)));

%% === 绘图 ===
figure('Position', [100 100 1200 800], 'Color', 'white');

% 1. Uncoded BPSK (Reference)
semilogy(EbN0_theory, BER_uncoded, '-.', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5, ...
    'DisplayName', 'Uncoded BPSK Theory');
hold on;

% 2. NASA Baseline (K=8920 FER) - 蓝色
semilogy(nasa_8920_snr, nasa_8920_fer, '--s', 'Color', [0.0 0.3 0.7], 'LineWidth', 2, ...
    'MarkerSize', 10, 'MarkerFaceColor', [0.0 0.3 0.7], ...
    'DisplayName', 'NASA K=8920 (FER)');

% 3. My Implementation (FER) - 绿色
semilogy(my_snr, my_fer, '-^', 'Color', [0.1 0.7 0.2], 'LineWidth', 2.5, ...
    'MarkerSize', 9, 'MarkerFaceColor', [0.1 0.7 0.2], ...
    'DisplayName', 'Our Implementation (FER)');

% 4. My Implementation (BER) - 红色虚线
semilogy(my_snr, my_ber, ':o', 'Color', [0.8 0.1 0.1], 'LineWidth', 2, ...
    'MarkerSize', 7, 'MarkerFaceColor', 'none', ...
    'DisplayName', 'Our Implementation (BER)');

%% === 图表装饰 ===
grid on;
grid minor;
title({'CCSDS Turbo Code (K=8920) Validation', 'Comparison with NASA CCSDS 131.0-B-5 Standard'}, ...
    'FontSize', 14, 'FontWeight', 'bold');
xlabel('Eb/N0 (dB)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Error Rate (BER / FER)', 'FontSize', 12, 'FontWeight', 'bold');

legend('Location', 'southwest', 'FontSize', 11);
axis([-0.2 1.2 1e-6 1]);

% 添加文本标注
annotation('textbox', [0.6 0.7 0.25 0.15], 'String', ...
    {'Configuration:', ' - Code: CCSDS Turbo (16-state)', ' - Rate: 1/3', ' - Block: K=8920 bits', ' - Iterations: 10'}, ...
    'FontSize', 10, 'EdgeColor', 'black', 'BackgroundColor', [1 1 0.9]);

%% === 保存与输出 ===
saveas(gcf, 'validation_8920_plot.png');
fprintf('Plot saved to: %s\n', fullfile(pwd, 'validation_8920_plot.png'));

%% === 关键点对比 ===
fprintf('\n=== Critical SNR Points Comparison (K=8920) ===\n');
fprintf('Eb/N0 (dB) | NASA FER     | Our FER      | Match?\n');
fprintf('-----------|--------------|--------------|---------\n');
for i = 1:length(nasa_8920_snr)
    snr = nasa_8920_snr(i);
    nasa_val = nasa_8920_fer(i);
    
    % 找最接近的 My Data 点
    [min_diff, idx] = min(abs(my_snr - snr));
    if min_diff < 0.05
        my_val = my_fer(idx);
        ratio = my_val / nasa_val;
        if ratio > 0.5 && ratio < 2.0
            match = 'OK';
        else
            match = 'CHECK';
        end
        fprintf('%5.1f      | %10.2e   | %10.2e   | %s\n', snr, nasa_val, my_val, match);
    else
        fprintf('%5.1f      | %10.2e   |      N/A     | -\n', snr, nasa_val);
    end
end

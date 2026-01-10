%% CCSDS Turbo Code Validation (Refined)
% 对比 CCSDS 16-state 实现 (High Precision) 与 NASA Table A-4 标准数据

clear; close all; clc;

%% === 读取仿真数据 ===
csv_file = 'ber_ccsds_refined.csv';
if ~exist(csv_file, 'file')
    error('Could not find CSV file: %s', csv_file);
end

data = readtable(csv_file, 'CommentStyle', '#');
fprintf('Reading data from: %s\n', csv_file);

% 提取数据列
my_snr = data{:, 1};   % Eb/N0 (dB)
my_ber = data{:, 4};   % BER

% 提取 FER (如果是7列格式自动提取, 否则置NaN)
if width(data) >= 7
    my_fer = data{:, 7};   % FER
else
    my_fer = NaN(size(my_snr));
    warning('FER data column missing.');
end

% 过滤掉 0 值以便对数绘图
my_ber(my_ber == 0) = NaN;
my_fer(my_fer == 0) = NaN;

%% === NASA 基准数据 (Table A-4, K=1784, R=1/3) ===
% Eb/N0 | FER
nasa_snr = [-0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0];
nasa_fer = [9.9020e-01, 9.0090e-01, 6.8493e-01, 2.9762e-01, 4.7174e-02, 4.4583e-03, 9.2350e-05, 1.9100e-06];

%% === 理论参考曲线 ===
EbN0_theory = -1:0.05:2;
BER_uncoded = 0.5 * erfc(sqrt(10.^(EbN0_theory/10)));

%% === 绘图 ===
figure('Position', [100 100 1200 800], 'Color', 'white');

% 1. Uncoded BPSK (Reference)
semilogy(EbN0_theory, BER_uncoded, '-.', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5, ...
    'DisplayName', 'Uncoded BPSK Theory');
hold on;

% 2. NASA Baseline (FER) - 蓝色填充区域/粗线
semilogy(nasa_snr, nasa_fer, '--s', 'Color', [0.0 0.3 0.7], 'LineWidth', 2, ...
    'MarkerSize', 8, 'MarkerFaceColor', [0.0 0.3 0.7], ...
    'DisplayName', 'NASA Standard (FER)');

% 3. My Implementation (FER) - 绿色实线
semilogy(my_snr, my_fer, '-^', 'Color', [0.1 0.7 0.2], 'LineWidth', 2.5, ...
    'MarkerSize', 9, 'MarkerFaceColor', [0.1 0.7 0.2], ...
    'DisplayName', 'Our CCSDS Impl (FER)');

% 4. My Implementation (BER) - 红色虚线
semilogy(my_snr, my_ber, ':o', 'Color', [0.8 0.1 0.1], 'LineWidth', 2, ...
    'MarkerSize', 7, 'MarkerFaceColor', 'none', ...
    'DisplayName', 'Our CCSDS Impl (BER)');

%% === 图表装饰 ===
grid on;
grid minor;
title({'CCSDS Turbo Code (K=1784) Performance Check', 'Comparison with NASA CCSDS 131.0-B-5 Standard'}, ...
    'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'none');
xlabel('Eb/N0 (dB)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Error Rate (BER / FER)', 'FontSize', 12, 'FontWeight', 'bold');

legend('Location', 'southwest', 'FontSize', 11);
axis([-0.6 1.3 1e-6 1]);

% 添加文本标注
annotation('textbox', [0.6 0.7 0.25 0.15], 'String', ...
    {'Configuration:', ' - Code: CCSDS Turbo (16-state)', ' - Rate: 1/3', ' - Block: K=1784 bits', ' - Iterations: 10'}, ...
    'FontSize', 10, 'EdgeColor', 'black', 'BackgroundColor', [1 1 0.9]);

%% === 保存与输出 ===
saveas(gcf, 'refined_comparison_plot.png');
fprintf('Plot saved to: %s\\refined_comparison_plot.png\n', pwd);

% 打印关键点对比
fprintf('\n=== Critical SNR Points Comparison ===\n');
fprintf('Eb/N0 (dB) | NASA FER     | Our FER      | Our BER\n');
fprintf('-----------|--------------|--------------|--------------\n');
check_points = [0.0, 0.4, 0.6, 0.8, 1.0];
for snr = check_points
    % 找最接近的 NASA 点
    [min_diff_nasa, idx_nasa] = min(abs(nasa_snr - snr));
    nasa_val = NaN;
    if min_diff_nasa < 0.05
        nasa_val = nasa_fer(idx_nasa);
    end
    
    % 找最接近的 My Data 点
    [min_diff_my, idx_my] = min(abs(my_snr - snr));
    my_fer_val = NaN;
    my_ber_val = NaN;
    if min_diff_my < 0.05
        my_fer_val = my_fer(idx_my);
        my_ber_val = my_ber(idx_my);
    end
    
    fprintf('%5.1f      | %10.2e   | %10.2e   | %10.2e\n', snr, nasa_val, my_fer_val, my_ber_val);
end

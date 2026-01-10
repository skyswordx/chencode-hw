%% CCSDS Turbo Code K=1784 R=1/2 Validation
% 对比 R=1/2 仿真结果与 NASA Table A-3 (R=1/2) 基准

clear; close all; clc;

%% === 读取合并后的仿真数据 ===
csv_file = 'ber_Turbo_LogMAP_20260110_203824_merged.csv';
if ~exist(csv_file, 'file')
    error('Could not find CSV file: %s', csv_file);
end

data = readtable(csv_file, 'CommentStyle', '#');
fprintf('Reading merged R=1/2 data from: %s\n', csv_file);

% 提取数据列
my_snr = data{:, 1};   % Eb/N0 (dB)
my_ber = data{:, 4};   % BER
my_fer = data{:, 7};   % FER

% 过滤掉 0 值
my_ber(my_ber == 0) = NaN;
my_fer(my_fer == 0) = NaN;

%% === NASA Table A-3 基准数据 (K=1784, R=1/2) ===
nasa_r12_snr = [0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2];
nasa_r12_fer = [7.5000e-01, 3.8931e-01, 7.5529e-02, 7.9605e-03, 3.2503e-04, ...
    1.1620e-05, 3.7500e-06, 2.2500e-06, 7.5000e-07];

%% === NASA Table A-4 基准数据 (K=1784, R=1/3 无打孔) ===
nasa_r13_snr = [-0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
nasa_r13_fer = [8.1169e-01, 4.9505e-01, 1.6026e-01, 2.4938e-02, 2.0964e-03, ...
    1.1696e-04, 8.4103e-06, 5.0116e-06, 4.0185e-06, 3.3272e-06, 2.8830e-06];

%% === 理论参考曲线 ===
EbN0_theory = -0.5:0.05:3.5;
BER_uncoded = 0.5 * erfc(sqrt(10.^(EbN0_theory/10)));

%% === 绘图 ===
figure('Position', [100 100 1200 800], 'Color', 'white');

% 1. Uncoded BPSK (Reference)
semilogy(EbN0_theory, BER_uncoded, '-.', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5, ...
    'DisplayName', 'Uncoded BPSK');
hold on;

% 2. NASA K=1784 R=1/3 (FER) - 蓝色虚线
semilogy(nasa_r13_snr, nasa_r13_fer, '--s', 'Color', [0.0 0.3 0.7], 'LineWidth', 2, ...
    'MarkerSize', 9, 'MarkerFaceColor', [0.0 0.3 0.7], ...
    'DisplayName', 'NASA K=1784 R=1/3 (FER)');

% 3. NASA K=1784 R=1/2 (FER) - 橙色虚线
semilogy(nasa_r12_snr, nasa_r12_fer, '--d', 'Color', [0.9 0.5 0.0], 'LineWidth', 2, ...
    'MarkerSize', 9, 'MarkerFaceColor', [0.9 0.5 0.0], ...
    'DisplayName', 'NASA K=1784 R=1/2 (FER)');

% 4. Our R=1/2 Implementation (FER) - 绿色实线
semilogy(my_snr, my_fer, '-^', 'Color', [0.1 0.7 0.2], 'LineWidth', 2.5, ...
    'MarkerSize', 8, 'MarkerFaceColor', [0.1 0.7 0.2], ...
    'DisplayName', 'Our R=1/2 Implementation (FER)');

% 5. Our R=1/2 Implementation (BER) - 红色实线
semilogy(my_snr, my_ber, '-o', 'Color', [0.8 0.1 0.1], 'LineWidth', 2, ...
    'MarkerSize', 7, 'MarkerFaceColor', 'none', ...
    'DisplayName', 'Our R=1/2 Implementation (BER)');

%% === 图表装饰 ===
grid on;
grid minor;
title({'CCSDS Turbo Code K=1784 Rate 1/2 Validation', ...
    'Comparison with NASA Table A-3 Baseline'}, ...
    'FontSize', 14, 'FontWeight', 'bold');
xlabel('Eb/N0 (dB)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Error Rate (BER / FER)', 'FontSize', 12, 'FontWeight', 'bold');

legend('Location', 'southwest', 'FontSize', 10);
axis([-0.5 3.5 1e-7 1]);

% 添加文本标注
annotation('textbox', [0.55 0.60 0.32 0.25], 'String', ...
    {'Configuration:', ...
    ' - Codec: CCSDS 16-state Turbo', ...
    ' - Block: K=1784 bits', ...
    ' - Rate: R=1/2 (punctured)', ...
    '', ...
    'Reference:', ...
    ' - NASA Table A-3 (R=1/2)', ...
    ' - NASA Table A-4 (R=1/3)'}, ...
    'FontSize', 10, 'EdgeColor', 'black', 'BackgroundColor', [1 1 0.9]);

%% === 保存 ===
saveas(gcf, 'r12_validation_nasa.png');
fprintf('\nPlot saved to: %s\n', fullfile(pwd, 'r12_validation_nasa.png'));

%% === NASA 对比表 ===
fprintf('\n=== Comparison with NASA Table A-3 (K=1784 R=1/2) ===\n');
fprintf('Eb/N0 (dB) | NASA FER     | Our FER      | Ratio\n');
fprintf('-----------+--------------+--------------+-------\n');
for i = 1:length(nasa_r12_snr)
    snr = nasa_r12_snr(i);
    nasa_val = nasa_r12_fer(i);
    
    idx = find(abs(my_snr - snr) < 0.05);
    if ~isempty(idx) && ~isnan(my_fer(idx))
        my_val = my_fer(idx);
        ratio = my_val / nasa_val;
        fprintf('%5.1f      | %10.2e   | %10.2e   | %.2f\n', snr, nasa_val, my_val, ratio);
    else
        fprintf('%5.1f      | %10.2e   |      N/A     | -\n', snr, nasa_val);
    end
end

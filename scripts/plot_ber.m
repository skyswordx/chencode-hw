%% Turbo Code BER Performance Plot
% 绘制完整的 Turbo 码仿真结果曲线（包含瀑布区和错误平层）
% 数据来源：ber_turbo_complete.csv

clear; close all; clc;

%% === 读取 CSV 数据 ===
data_dir = '.';
target_file = 'ber_turbo_complete.csv';
filepath = fullfile(data_dir, target_file);

fprintf('正在读取: %s\n', target_file);
data = readtable(filepath, 'CommentStyle', '#');

% 提取数据
snr_db = data.SNR_dB;
ber = data.BER;

% 过滤掉 BER=0 的点（无法在对数坐标显示）
valid_idx = ber > 0;
snr_db = snr_db(valid_idx);
ber = ber(valid_idx);

fprintf('数据范围: SNR = [%.1f, %.1f] dB, 共 %d 个有效数据点\n', ...
    min(snr_db), max(snr_db), length(snr_db));

%% === 理论曲线计算 ===
% Uncoded BPSK 理论曲线
EbN0_dB_theory = -2:0.2:8;
EbN0_theory = 10.^(EbN0_dB_theory/10);
BER_uncoded = 0.5 * erfc(sqrt(EbN0_theory));

%% === 绘图 ===
figure('Position', [100 100 1000 700], 'Color', 'white');

% 绘制 Uncoded BPSK 理论曲线
semilogy(EbN0_dB_theory, BER_uncoded, '-.', 'Color', [0.6 0.6 0.6], ...
    'LineWidth', 1.5, 'DisplayName', 'Uncoded BPSK (Theory)');
hold on;

% 绘制 Turbo 码仿真曲线
semilogy(snr_db, ber, '-', 'LineWidth', 2.5, 'Color', [0.8 0.2 0.6], ...
    'Marker', 'd', 'MarkerSize', 8, 'MarkerFaceColor', [0.8 0.2 0.6], ...
    'DisplayName', 'Turbo Code (R=1/3, 8 iter)');

%% === 标注瀑布区和错误平层 ===
% 瀑布区标注
text(0.5, 1e-2, 'Waterfall Region', 'FontSize', 12, 'FontWeight', 'bold', ...
    'Color', [0.2 0.4 0.8], 'Rotation', -60);

% 错误平层标注
text(4, 5e-8, 'Error Floor Region', 'FontSize', 12, 'FontWeight', 'bold', ...
    'Color', [0.8 0.2 0.2]);

% 绘制错误平层参考线
yline(1e-7, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, ...
    'Label', 'BER = 10^{-7}', 'LabelHorizontalAlignment', 'left');

%% === 图表美化 ===
xlabel('$E_b/N_0$ (dB)', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Bit Error Rate (BER)', 'FontSize', 14, 'Interpreter', 'latex');
title({'Turbo Code BER Performance with Error Floor', ...
    'PCCC $R=1/3$, Block size $k=1024$ bits, 8 iterations, AWGN Channel'}, ...
    'FontSize', 14, 'Interpreter', 'latex');

legend('Location', 'southwest', 'FontSize', 11);
grid on;

% 设置坐标轴范围
xlim([-2 8]);
ylim([1e-9 1]);

% 设置刻度
set(gca, 'YScale', 'log');
set(gca, 'FontSize', 12);
ax = gca;
ax.YMinorGrid = 'on';
ax.YMinorTick = 'on';
ax.XMinorGrid = 'on';

%% === 添加编码增益标注 ===
% 在 BER = 1e-5 处计算编码增益
target_ber = 1e-5;

% Uncoded: BER = 0.5*erfc(sqrt(EbN0)) => 需要数值求解
EbN0_uncoded_at_target = fzero(@(x) 0.5*erfc(sqrt(10^(x/10))) - target_ber, 9);

% 从仿真数据插值获取 Turbo 码在目标 BER 处的 SNR
if min(ber) < target_ber && max(ber) > target_ber
    EbN0_turbo_at_target = interp1(log10(ber), snr_db, log10(target_ber), 'linear');
    coding_gain = EbN0_uncoded_at_target - EbN0_turbo_at_target;
    
    % 添加编码增益文本
    text(-1, 2e-8, sprintf('Coding Gain \\approx %.1f dB @ BER=10^{-5}', coding_gain), ...
        'FontSize', 12, 'FontWeight', 'bold', 'Color', [0.1 0.5 0.1], ...
        'BackgroundColor', 'white', 'EdgeColor', [0.1 0.5 0.1]);
    
    fprintf('\n编码增益 @ BER=1e-5: %.2f dB\n', coding_gain);
end

%% === 保存图片 ===
output_fig = fullfile(data_dir, 'ber_turbo_complete_plot.png');
saveas(gcf, output_fig);

% 保存高分辨率版本
print(gcf, fullfile(data_dir, 'ber_turbo_complete_plot_hires.png'), '-dpng', '-r300');

fprintf('\n绘图完成！已保存到:\n  - %s\n  - %s\n', output_fig, ...
    fullfile(data_dir, 'ber_turbo_complete_plot_hires.png'));

%% === 打印数据摘要 ===
fprintf('\n=== 仿真数据摘要 ===\n');
fprintf('%-10s %-15s %-15s %-15s\n', 'SNR (dB)', 'Bit Errors', 'Total Bits', 'BER');
fprintf('%s\n', repmat('-', 1, 55));
for i = 1:length(snr_db)
    fprintf('%-10.1f %-15d %-15d %-15.2e\n', ...
        snr_db(i), data.Bit_Errors(valid_idx), data.Total_Bits(valid_idx), ber(i));
end

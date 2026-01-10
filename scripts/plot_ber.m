%% BER vs Eb/N0 Comparison Plot
% 信道编码仿真结果绘图脚本
% 使用方法：修改 filenames 数组，指向你的 CSV 文件
%
% CSV 文件格式要求：
%   - 表头行: SNR_dB,Bit_Errors,Total_Bits,BER
%   - 支持 # 开头的注释行

clear; close all; clc;

%% === 配置区域 ===
% 修改这里的文件名以匹配你的 CSV 文件
data_dir = '../output/';

% 自动搜索 output 目录下的所有 CSV 文件
csv_files = dir(fullfile(data_dir, 'ber_*.csv'));

% 手动指定文件（如果不想自动搜索）
% filenames = {
%     'ber_uncoded_*.csv',
%     'ber_hard_viterbi_*.csv', 
%     'ber_soft_viterbi_*.csv',
%     'ber_bcjr_*.csv',
%     'ber_turbo_*.csv'
% };

% 标签与样式
style_map = containers.Map();
style_map('uncoded')      = struct('label', 'Uncoded BPSK (Theory)', 'style', '--k', 'marker', 'none');
style_map('hard_viterbi') = struct('label', 'Hard Viterbi', 'style', '-r', 'marker', 'o');
style_map('soft_viterbi') = struct('label', 'Soft Viterbi', 'style', '-b', 'marker', 's');
style_map('bcjr')         = struct('label', 'BCJR (Log-MAP)', 'style', '-g', 'marker', '^');
style_map('turbo')        = struct('label', 'Turbo (8 iter)', 'style', '-m', 'marker', 'd');

%% === 理论曲线计算 ===
EbN0_dB = 0:0.5:10;
EbN0 = 10.^(EbN0_dB/10);
BER_theory = 0.5 * erfc(sqrt(EbN0));

%% === 绘图 ===
figure('Position', [100 100 900 650], 'Color', 'white');
hold on; grid on;

% 绘制理论曲线
semilogy(EbN0_dB, BER_theory, '-.', 'Color', [0.5 0.5 0.5], ...
    'LineWidth', 1.5, 'DisplayName', 'Uncoded BPSK (Theory)');

% 读取并绘制每个 CSV 文件
legend_entries = {'Uncoded BPSK (Theory)'};

for i = 1:length(csv_files)
    filepath = fullfile(data_dir, csv_files(i).name);
    
    try
        % 使用 readtable 读取 (支持注释行)
        data = readtable(filepath, 'CommentStyle', '#');
        
        % 提取译码器类型 (从文件名)
        name_parts = strsplit(csv_files(i).name, '_');
        if length(name_parts) >= 2
            decoder_type = name_parts{2};
            % 处理带时间戳的文件名 (ber_soft_viterbi_20260109_...)
            if any(strcmp(decoder_type, {'soft', 'hard'})) && length(name_parts) >= 3
                decoder_type = [name_parts{2} '_' name_parts{3}];
            end
        else
            decoder_type = 'unknown';
        end
        
        % 获取样式
        if isKey(style_map, decoder_type)
            s = style_map(decoder_type);
            semilogy(data.SNR_dB, data.BER, s.style, ...
                'LineWidth', 1.5, 'MarkerSize', 7, ...
                'MarkerFaceColor', 'auto', 'DisplayName', s.label);
            legend_entries{end+1} = s.label;
        else
            % 默认样式
            semilogy(data.SNR_dB, data.BER, '-', ...
                'LineWidth', 1.5, 'DisplayName', decoder_type);
            legend_entries{end+1} = decoder_type;
        end
        
        fprintf('已加载: %s\n', csv_files(i).name);
    catch ME
        fprintf('无法读取 %s: %s\n', csv_files(i).name, ME.message);
    end
end

%% === 图表美化 ===
xlabel('$E_b/N_0$ (dB)', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('BER', 'FontSize', 14);
title({'BER Performance: $(7,5)_8$ Convolutional Code, $R=1/2$', ...
       'Block size $k=1024$ bits | AWGN Channel | BPSK Modulation'}, ...
       'FontSize', 14, 'Interpreter', 'latex');

legend('Location', 'southwest', 'FontSize', 11);

% 设置坐标轴范围
xlim([0 10]);
ylim([1e-7 1]);

% 设置刻度
set(gca, 'YScale', 'log');
set(gca, 'FontSize', 12);
set(gca, 'LineWidth', 1);

% 添加 minor grid
ax = gca;
ax.YMinorGrid = 'on';
ax.YMinorTick = 'on';

%% === 保存图片 ===
% 保存为 PNG
saveas(gcf, fullfile(data_dir, 'ber_comparison.png'));

% 保存为 PDF (矢量图，适合论文)
% exportgraphics(gcf, fullfile(data_dir, 'ber_comparison.pdf'), 'ContentType', 'vector');

fprintf('\n绘图完成！已保存到: %s\n', fullfile(data_dir, 'ber_comparison.png'));

%% === 可选：生成表格数据 ===
% 打印 LaTeX 格式表格
fprintf('\n=== LaTeX 表格数据 ===\n');
fprintf('$E_b/N_0$ (dB) & Uncoded & Hard Viterbi & Soft Viterbi & BCJR \\\\ \\hline\n');

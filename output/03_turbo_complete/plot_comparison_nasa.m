%% Turbo Code FER: Our Implementation vs NASA Baselines
% 使用 FER (Frame Error Rate) 与 NASA CCSDS 标准对比

clear; close all; clc;

%% === 读取仿真数据 ===
% 新 CSV 格式: SNR_dB,Bit_Errors,Total_Bits,BER,Frame_Errors,Total_Frames,FER
data = readtable('ber_turbo_complete.csv', 'CommentStyle', '#');

our_snr = data.SNR_dB;

% 检查是否有 FER 列
if ismember('FER', data.Properties.VariableNames)
    our_fer = data.FER;
    fprintf('使用 CSV 中的实际 FER 数据\n');
else
    % 如果没有 FER 列,从 BER 估算
    our_ber = data.BER;
    block_size = 1024;
    our_fer = 1 - (1 - our_ber).^block_size;
    fprintf('从 BER 估算 FER\n');
end

our_fer(our_fer == 0) = NaN;

%% === NASA 基准数据 (FER) ===
nasa_1784_snr = [-0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0];
nasa_1784_fer = [9.9020e-01, 9.0090e-01, 6.8493e-01, 2.9762e-01, 4.7174e-02, 4.4583e-03, 9.2350e-05, 1.9100e-06];

nasa_3568_snr = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7];
nasa_3568_fer = [1.8065e-01, 5.1557e-02, 9.0463e-03, 9.5734e-04, 4.1120e-05, 4.5100e-06];

nasa_7136_snr = [0.1, 0.2, 0.3, 0.4, 0.5];
nasa_7136_fer = [4.3328e-01, 1.0761e-01, 1.0989e-02, 4.4099e-04, 1.0050e-05];

nasa_8920_snr = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7];
nasa_8920_fer = [8.3333e-01, 4.9505e-01, 9.7752e-02, 8.9847e-03, 2.0755e-04, 2.8730e-05, 1.4360e-05, 1.1490e-05];

%% === 绘图 ===
figure('Position', [50 50 1100 750], 'Color', 'white');

% NASA curves
semilogy(nasa_1784_snr, nasa_1784_fer, '-s', 'Color', [0.0 0.4 0.8], ...
    'LineWidth', 1.8, 'MarkerSize', 8, 'MarkerFaceColor', [0.0 0.4 0.8], ...
    'DisplayName', 'NASA 1784 bits');
hold on;

semilogy(nasa_3568_snr, nasa_3568_fer, '-^', 'Color', [0.0 0.7 0.3], ...
    'LineWidth', 1.8, 'MarkerSize', 8, 'MarkerFaceColor', [0.0 0.7 0.3], ...
    'DisplayName', 'NASA 3568 bits');

semilogy(nasa_7136_snr, nasa_7136_fer, '-v', 'Color', [0.9 0.5 0.0], ...
    'LineWidth', 1.8, 'MarkerSize', 8, 'MarkerFaceColor', [0.9 0.5 0.0], ...
    'DisplayName', 'NASA 7136 bits');

semilogy(nasa_8920_snr, nasa_8920_fer, '-o', 'Color', [0.6 0.0 0.6], ...
    'LineWidth', 1.8, 'MarkerSize', 8, 'MarkerFaceColor', [0.6 0.0 0.6], ...
    'DisplayName', 'NASA 8920 bits');

% Our implementation - RED highlight
semilogy(our_snr, our_fer, '-d', 'Color', [0.9 0.1 0.3], ...
    'LineWidth', 2.5, 'MarkerSize', 6, 'MarkerFaceColor', [0.9 0.1 0.3], ...
    'DisplayName', 'Our Impl. (1024 bits)');

%% === 图表美化 ===
xlabel('$E_b/N_0$ (dB)', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Frame Error Rate (FER)', 'FontSize', 14, 'Interpreter', 'latex');
title({'Turbo Code FER: Our Implementation vs NASA CCSDS Baselines', ...
    'Rate 1/3, AWGN Channel'}, 'FontSize', 13, 'Interpreter', 'latex');

legend('Location', 'southwest', 'FontSize', 10);
grid on;
xlim([-0.5 2]);
ylim([1e-6 1]);
set(gca, 'FontSize', 12);

annotation('textbox', [0.55 0.70 0.33 0.15], ...
    'String', {'Block Size Comparison:', ...
    '  Our: 1024 bits', ...
    '  NASA: 1784 - 8920 bits', ...
    '', ...
    'Larger block = Better FER'}, ...
    'FontSize', 10, 'BackgroundColor', [1 1 0.95], 'EdgeColor', [0.3 0.3 0.3]);

saveas(gcf, 'comparison_fer_nasa.png');
print(gcf, 'comparison_fer_nasa_hires.png', '-dpng', '-r300');
fprintf('已保存: comparison_fer_nasa.png\n');

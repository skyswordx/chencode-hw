%% Plot Fused Results with NASA Comparison
% Plots both K=1784 and K=8920 fused results and compares with NASA baselines.

clear; close all; clc;

%% === Configuration ===
SCRIPT_DIR = fileparts(mfilename('fullpath'));
DATA_DIR = fullfile(SCRIPT_DIR, 'fused_results'); % Fused data in subfolder
FIGURE_DIR = fullfile(SCRIPT_DIR, 'figures');

if ~exist(FIGURE_DIR, 'dir')
    mkdir(FIGURE_DIR);
end

%% === NASA Reference Data (FER from CCSDS 131.0-B-5 Table A-4) ===
% Rate 1/3, K=1784
nasa_1784_r13_snr = [-0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0];
nasa_1784_r13_fer = [9.9020e-01, 9.0090e-01, 6.8493e-01, 2.9762e-01, 4.7174e-02, 4.4583e-03, 9.2350e-05, 1.9100e-06];

% Rate 1/3, K=8920
nasa_8920_r13_snr = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7];
nasa_8920_r13_fer = [8.3333e-01, 4.9505e-01, 9.7752e-02, 8.9847e-03, 2.0755e-04, 2.8730e-05, 1.4360e-05, 1.1490e-05];

%% === Process Both Fused Files ===
files_to_plot = {
    'ber_Turbo_LogMAP_K1784_R1-3_FUSED.csv', ...
    'ber_Turbo_LogMAP_K8920_R1-3_FUSED.csv'
    };

failed = {};

for i = 1:length(files_to_plot)
    filename = files_to_plot{i};
    filepath = fullfile(DATA_DIR, filename);
    
    fprintf('Processing: %s\n', filename);
    
    if exist(filepath, 'file')
        try
            process_file(filepath, filename, FIGURE_DIR, ...
                nasa_1784_r13_snr, nasa_1784_r13_fer, ...
                nasa_8920_r13_snr, nasa_8920_r13_fer);
        catch ME
            fprintf('  Error: %s\n', ME.message);
            failed{end+1} = filename;
        end
    else
        fprintf('  Error: File not found at %s\n', filepath);
        failed{end+1} = filename;
    end
end

% Also create combined comparison plot
create_combined_plot(DATA_DIR, FIGURE_DIR, ...
    nasa_1784_r13_snr, nasa_1784_r13_fer, ...
    nasa_8920_r13_snr, nasa_8920_r13_fer);

if ~isempty(failed)
    fprintf('\n=== Failed files ===\n');
    for k = 1:length(failed)
        fprintf('  - %s\n', failed{k});
    end
end

fprintf('\n=== Done! Figures saved to: %s ===\n', FIGURE_DIR);

%% === Individual File Plot Function ===
function process_file(filepath, filename, fig_dir, ...
    nasa1784_r13_x, nasa1784_r13_y, ...
    nasa8920_r13_x, nasa8920_r13_y)

% Read Data
opts = detectImportOptions(filepath);
opts.VariableNamingRule = 'preserve';
data = readtable(filepath, opts);

% Extract Columns
if ~all(ismember({'Eb_N0', 'BER', 'FER'}, data.Properties.VariableNames))
    error('Missing required columns (Eb_N0, BER, FER).');
end

snr = data.Eb_N0;
ber = data.BER;
fer = data.FER;

% Replace zeros with NaN for log plot
ber(ber == 0) = NaN;
fer(fer == 0) = NaN;

% Setup Figure
fig = figure('Visible', 'off', 'Position', [100 100 1000 750], 'Color', 'white');

% Plot BER & FER
semilogy(snr, ber, '-o', 'Color', [0.8 0.1 0.1], 'LineWidth', 1.8, ...
    'MarkerSize', 8, 'MarkerFaceColor', 'none', 'DisplayName', 'Simulation BER');
hold on;
semilogy(snr, fer, '-^', 'Color', [0.1 0.6 0.2], 'LineWidth', 1.8, ...
    'MarkerSize', 8, 'MarkerFaceColor', [0.1 0.6 0.2], 'DisplayName', 'Simulation FER');

% Determine type and add NASA baseline
if contains(filename, 'K1784_R1-3', 'IgnoreCase', true)
    title_str = 'CCSDS Turbo Code ($K=1784, R=1/3$)';
    semilogy(nasa1784_r13_x, nasa1784_r13_y, '--s', 'Color', [0.0 0.4 0.8], ...
        'LineWidth', 2, 'MarkerSize', 9, 'MarkerFaceColor', [0.0 0.4 0.8], ...
        'DisplayName', 'NASA Baseline (FER)');
    
elseif contains(filename, 'K8920_R1-3', 'IgnoreCase', true)
    title_str = 'CCSDS Turbo Code ($K=8920, R=1/3$)';
    semilogy(nasa8920_r13_x, nasa8920_r13_y, '--s', 'Color', [0.6 0.0 0.6], ...
        'LineWidth', 2, 'MarkerSize', 9, 'MarkerFaceColor', [0.6 0.0 0.6], ...
        'DisplayName', 'NASA Baseline (FER)');
else
    title_str = strrep(filename, '_FUSED.csv', '');
    title_str = strrep(title_str, '_', ' ');
end

% Styling with LaTeX
grid on; grid minor;
set(gca, 'FontSize', 12);
xlabel('$E_b/N_0$ (dB)', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Error Rate', 'FontSize', 14, 'Interpreter', 'latex');
title(title_str, 'FontSize', 16, 'Interpreter', 'latex');

% Legend
lgd = legend('Location', 'southwest');
set(lgd, 'Interpreter', 'latex', 'FontSize', 11);

% Axis limits
xlim([min(snr)-0.2, max(snr)+0.2]);
ylim([1e-7, 1.1]);

% Save PNG
out_name = strrep(filename, '.csv', '.png');
saveas(fig, fullfile(fig_dir, out_name));
fprintf('  -> Saved: %s\n', out_name);

close(fig);
end

%% === Combined Comparison Plot Function ===
function create_combined_plot(data_dir, fig_dir, ...
    nasa1784_x, nasa1784_y, nasa8920_x, nasa8920_y)

fprintf('\nCreating combined comparison plot...\n');

% Read both files
file1 = fullfile(data_dir, 'ber_Turbo_LogMAP_K1784_R1-3_FUSED.csv');
file2 = fullfile(data_dir, 'ber_Turbo_LogMAP_K8920_R1-3_FUSED.csv');

if ~exist(file1, 'file') || ~exist(file2, 'file')
    fprintf('  Skipping combined plot: one or both files missing.\n');
    return;
end

opts1 = detectImportOptions(file1); opts1.VariableNamingRule = 'preserve';
opts2 = detectImportOptions(file2); opts2.VariableNamingRule = 'preserve';
data1 = readtable(file1, opts1);
data2 = readtable(file2, opts2);

% Extract FER data
snr1 = data1.Eb_N0; fer1 = data1.FER; fer1(fer1 == 0) = NaN;
snr2 = data2.Eb_N0; fer2 = data2.FER; fer2(fer2 == 0) = NaN;

% Create figure
fig = figure('Visible', 'off', 'Position', [100 100 1100 800], 'Color', 'white');

% K=1784 Simulation
semilogy(snr1, fer1, '-o', 'Color', [0.2 0.6 0.3], 'LineWidth', 2, ...
    'MarkerSize', 8, 'MarkerFaceColor', [0.2 0.6 0.3], ...
    'DisplayName', '$K=1784$ Simulation');
hold on;

% K=1784 NASA
semilogy(nasa1784_x, nasa1784_y, '--s', 'Color', [0.1 0.4 0.2], 'LineWidth', 1.5, ...
    'MarkerSize', 8, 'MarkerFaceColor', 'none', ...
    'DisplayName', '$K=1784$ NASA');

% K=8920 Simulation
semilogy(snr2, fer2, '-^', 'Color', [0.8 0.2 0.2], 'LineWidth', 2, ...
    'MarkerSize', 8, 'MarkerFaceColor', [0.8 0.2 0.2], ...
    'DisplayName', '$K=8920$ Simulation');

% K=8920 NASA
semilogy(nasa8920_x, nasa8920_y, '--d', 'Color', [0.5 0.1 0.1], 'LineWidth', 1.5, ...
    'MarkerSize', 8, 'MarkerFaceColor', 'none', ...
    'DisplayName', '$K=8920$ NASA');

% Styling
grid on; grid minor;
set(gca, 'FontSize', 12, 'GridAlpha', 0.3, 'MinorGridAlpha', 0.15);
xlabel('$E_b/N_0$ (dB)', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Frame Error Rate (FER)', 'FontSize', 14, 'Interpreter', 'latex');
title('CCSDS Turbo Code FER Comparison ($R=1/3$, 10 iterations)', ...
    'FontSize', 16, 'Interpreter', 'latex');

% Legend
lgd = legend('Location', 'southwest', 'NumColumns', 2);
set(lgd, 'Interpreter', 'latex', 'FontSize', 11);

% Axis
xlim([-0.6, 1.2]);
ylim([1e-7, 1.1]);

% Save
saveas(fig, fullfile(fig_dir, 'combined_fer_comparison.png'));
fprintf('  -> Saved: combined_fer_comparison.png\n');

close(fig);
end

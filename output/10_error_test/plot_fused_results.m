%% Plot Fused Results with NASA Comparison
% Automatically analyzes fused CSVs and compares Turbo codes with NASA baselines.

clear; close all; clc;

%% === Configuration ===
SCRIPT_DIR = fileparts(mfilename('fullpath'));
DATA_DIR = SCRIPT_DIR; % Data is in the same directory as script
FIGURE_DIR = fullfile(DATA_DIR, 'figures');

if ~exist(FIGURE_DIR, 'dir')
    mkdir(FIGURE_DIR);
end

%% === NASA Reference Data (FER only usually) ===
% Table A-4: CCSDS Turbo, Rate 1/3
nasa_1784_r13_snr = [-0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0];
nasa_1784_r13_fer = [9.9020e-01, 9.0090e-01, 6.8493e-01, 2.9762e-01, 4.7174e-02, 4.4583e-03, 9.2350e-05, 1.9100e-06];

nasa_8920_r13_snr = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7];
nasa_8920_r13_fer = [8.3333e-01, 4.9505e-01, 9.7752e-02, 8.9847e-03, 2.0755e-04, 2.8730e-05, 1.4360e-05, 1.1490e-05];

% Table A-3: CCSDS Turbo, Rate 1/2
nasa_1784_r12_snr = [0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2];
nasa_1784_r12_fer = [7.5000e-01, 3.8931e-01, 7.5529e-02, 7.9605e-03, 3.2503e-04, 1.1620e-05, 3.7500e-06, 2.2500e-06, 7.5000e-07];

%% === File Processing ===
failed = {};
% Specific file to process
% Users requested: only use ber_Turbo_LogMAP_K1784_R1-3_FUSED.csv
filename = 'ber_Turbo_LogMAP_K1784_R1-3_FUSED.csv';
filepath = fullfile(DATA_DIR, filename);

fprintf('Processing target file: %s\n', filename);

if exist(filepath, 'file')
    try
        process_file(filepath, filename, FIGURE_DIR, ...
            nasa_1784_r13_snr, nasa_1784_r13_fer, ...
            nasa_8920_r13_snr, nasa_8920_r13_fer, ...
            nasa_1784_r12_snr, nasa_1784_r12_fer);
    catch ME
        fprintf('Error processing %s: %s\n', filename, ME.message);
        failed{end+1} = filename;
    end
else
    fprintf('Error: Target file %s not found in %s\n', filename, DATA_DIR);
    failed{end+1} = filename;
end

if ~isempty(failed)
    fprintf('\nFailed files:\n');
    disp(failed);
end

%% === Helper Function ===
function process_file(filepath, filename, fig_dir, ...
    nasa1784_r13_x, nasa1784_r13_y, ...
    nasa8920_r13_x, nasa8920_r13_y, ...
    nasa1784_r12_x, nasa1784_r12_y)

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

% Plot Our BER & FER (Semilogy = Log Y axis)
h_ber = semilogy(snr, ber, '-o', 'Color', [0.8 0.1 0.1], 'LineWidth', 1.5, ...
    'MarkerFaceColor', 'none', 'DisplayName', 'Simulation BER');
hold on;
h_fer = semilogy(snr, fer, '-^', 'Color', [0.1 0.6 0.2], 'LineWidth', 1.5, ...
    'MarkerFaceColor', [0.1 0.6 0.2], 'DisplayName', 'Simulation FER');

% --- Determine Type & Plot Logic ---
is_turbo = false;
title_str = strrep(filename, '_FUSED.csv', '');
title_str = strrep(title_str, '_', ' ');

if contains(filename, 'BCJR_MAP', 'IgnoreCase', true)
    % BCJR Only - No NASA Comparison
    title_str = 'BCJR MAP Decoder Performance';
    % User specified: No NASA comparison
    
elseif contains(filename, 'K1784_R1-3', 'IgnoreCase', true)
    is_turbo = true;
    title_str = 'CCSDS Turbo Code ($K=1784, R=1/3$)';
    semilogy(nasa1784_r13_x, nasa1784_r13_y, '--s', 'Color', [0.0 0.4 0.8], ...
        'LineWidth', 1.5, 'MarkerFaceColor', [0.0 0.4 0.8], ...
        'DisplayName', 'NASA Baseline (FER)');
    
elseif contains(filename, 'K8920_R1-3', 'IgnoreCase', true)
    is_turbo = true;
    title_str = 'CCSDS Turbo Code ($K=8920, R=1/3$)';
    semilogy(nasa8920_r13_x, nasa8920_r13_y, '--s', 'Color', [0.6 0.0 0.6], ...
        'LineWidth', 1.5, 'MarkerFaceColor', [0.6 0.0 0.6], ...
        'DisplayName', 'NASA Baseline (FER)');
    
elseif contains(filename, 'Viterbi', 'IgnoreCase', true)
    % Viterbi
    title_str = ['Convolutional Code: ', title_str];
end

% Styling with LaTeX
grid on; grid minor;
xlabel('$E_b/N_0$ (dB)', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Error Rate (BER / FER)', 'FontSize', 14, 'Interpreter', 'latex');
title(title_str, 'FontSize', 16, 'Interpreter', 'latex');

% Legend Style
lgd = legend('Location', 'southwest');
set(lgd, 'Interpreter', 'latex', 'FontSize', 12);

% Adjust limits
if ~isempty(snr)
    xlim([min(snr)-0.5, max(snr)+0.5]);
end
ylim([1e-7, 1.1]);

% Save
out_name = strrep(filename, '.csv', '.png');
saveas(fig, fullfile(fig_dir, out_name));
fprintf('  -> Saved plot: %s\n', out_name);

close(fig);
end

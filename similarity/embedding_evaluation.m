%% MAIN ANALYSIS: Load Data and Create Plots
%  ========================================================================
%  This script:
%  1. Loads all experiment CSV files
%  2. Creates 3 main analysis plots
%  3. Prints key findings
%  ========================================================================

clear; clc;

fprintf('\n');
fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║                                                                ║\n');
fprintf('║  EXPERIMENT ANALYSIS PIPELINE                                  ║\n');
fprintf('║  Embedding Validation for DTW Approximation                    ║\n');
fprintf('║                                                                ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n');
fprintf('\n');


%% STEP 1: Load Data
% ========================================================================

fprintf('STEP 1: Loading experiment data...\n');
fprintf('────────────────────────────────────────\n\n');

% Load all CSV files from current directory
data = loadExperimentData('results/');

%% STEP 1.5: Smart Excel Export - Append Only New Data
% ========================================================================

fprintf('\nSTEP 1.5: Smart Excel Export (Append New Data Only)...\n');
fprintf('────────────────────────────────────────\n\n');

% Fixed filename
export_filename = 'results/experiment_data.xlsx';

% Get the new data to potentially add
new_data = data.combined;

fprintf('New data loaded: %d rows\n', height(new_data));

%% Check if Excel file exists
% ========================================================================

if exist(export_filename, 'file')
    fprintf('✓ Found existing file: %s\n', export_filename);
    
    % Read existing data from "All Data" sheet
    try
        fprintf('  Reading existing data from sheet "All Data"...\n');
        existing_data = readtable(export_filename, 'Sheet', 'All Data', 'VariableNamingRule', 'preserve');
        fprintf('  ✓ Existing data: %d rows\n', height(existing_data));
        
        % Check which timestamps are already in the file
        if ismember('Timestamp', existing_data.Properties.VariableNames) && ...
           ismember('Timestamp', new_data.Properties.VariableNames)
            
            % Get unique timestamps from both datasets
            existing_timestamps = unique(existing_data.Timestamp);
            new_timestamps = unique(new_data.Timestamp);
            
            fprintf('\n  Existing timestamps in file: %d unique\n', length(existing_timestamps));
            fprintf('  New timestamps in loaded data: %d unique\n', length(new_timestamps));
            
            % Find which timestamps are NOT in the existing file
            timestamps_to_add = setdiff(new_timestamps, existing_timestamps);
            
            if isempty(timestamps_to_add)
                fprintf('\n  ℹ All data already exists in Excel file - nothing to add!\n\n');
            else
                fprintf('\n  → Found %d NEW timestamp(s) to add:\n', length(timestamps_to_add));
                for i = 1:length(timestamps_to_add)
                    fprintf('      %d. %s\n', i, char(timestamps_to_add(i)));
                end
                
                % Filter new_data to only include rows with new timestamps
                if iscellstr(timestamps_to_add) || isstring(timestamps_to_add)
                    mask_to_add = ismember(new_data.Timestamp, timestamps_to_add);
                else
                    mask_to_add = false(height(new_data), 1);
                    for i = 1:length(timestamps_to_add)
                        mask_to_add = mask_to_add | strcmp(new_data.Timestamp, timestamps_to_add(i));
                    end
                end
                
                data_to_add = new_data(mask_to_add, :);
                
                fprintf('  → Rows to add: %d\n\n', height(data_to_add));
                
                % Append new data to existing data
                combined_data = [existing_data; data_to_add];
                
                fprintf('  Combined data: %d rows (was %d, added %d)\n\n', ...
                    height(combined_data), height(existing_data), height(data_to_add));
                
                % Write back to Excel
                fprintf('  Writing updated data back to Excel...\n');
                writetable(combined_data, export_filename, 'Sheet', 'All Data', 'WriteRowNames', false);
                fprintf('  ✓ Successfully appended new data!\n\n');
                
                fprintf('╔════════════════════════════════════════════════════════════════╗\n');
                fprintf('║  DATA APPENDED SUCCESSFULLY                                    ║\n');
                fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');
                fprintf('Excel file updated: %s\n', export_filename);
                fprintf('  Previous rows: %d\n', height(existing_data));
                fprintf('  Added rows: %d\n', height(data_to_add));
                fprintf('  Total rows now: %d\n\n', height(combined_data));
            end
            
        else
            fprintf('  ⚠ Warning: Timestamp column not found in existing or new data\n');
            fprintf('  Cannot determine which data is new - skipping append\n\n');
        end
        
    catch ME
        fprintf('  ✗ Error reading existing file: %s\n', ME.message);
        fprintf('  Will create new file instead\n\n');
        
        % Create new file (fallback)
        writetable(new_data, export_filename, 'Sheet', 'All Data', 'WriteRowNames', false);
        fprintf('  ✓ Created new Excel file with %d rows\n\n', height(new_data));
    end
    
else
    % File doesn't exist - create new
    fprintf('ℹ File does not exist - creating new Excel file...\n');
    fprintf('  Location: %s\n\n', export_filename);
    
    try
        writetable(new_data, export_filename, 'Sheet', 'All Data', 'WriteRowNames', false);
        fprintf('  ✓ Successfully created new file with %d rows\n\n', height(new_data));
        
        fprintf('╔════════════════════════════════════════════════════════════════╗\n');
        fprintf('║  NEW FILE CREATED                                              ║\n');
        fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');
        fprintf('Excel file created: %s\n', export_filename);
        fprintf('  Total rows: %d\n', height(new_data));
        fprintf('  Sheet: All Data\n\n');
        
    catch ME
        fprintf('  ✗ Error creating file: %s\n\n', ME.message);
    end
end

fprintf('Note: File contains ONLY "All Data" sheet.\n');
fprintf('      Trajectory/Segment filtering can be done in Excel.\n\n');

% ========================================================================
%% FINAL CORRECTED PAPER ANALYSIS
% ========================================================================
% Figure 1: Dimensionality Invariance - ONLY BASELINES (Joint only, Pos only)
% Figure 2: Weight Mode Contribution - ALL modes sorted by composite
% ========================================================================

clear; clc; close all;

fprintf('\n');
fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  FINAL CORRECTED PAPER ANALYSIS                                ║\n');
fprintf('║  • Figure 1: ONLY Baselines (pure dimensionality)              ║\n');
fprintf('║  • Figure 2: ALL modes (incremental features)                  ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n');
fprintf('\n');

% Configuration
excel_file = 'results/experiment_data.xlsx';
output_folder = 'results';
use_trajectory_only = true;
use_non_normalized_dtw = true;
save_figures = false;
figure_format = 'both';

% ========================================================================
%% STEP 1: LOAD AND FILTER
% ========================================================================

fprintf('Loading data...\n');
data_table = readtable(excel_file, 'Sheet', 'All Data', 'VariableNamingRule', 'preserve');

if use_trajectory_only
    data_table = data_table(strcmp(data_table.Level, 'Trajectory'), :);
end
if use_non_normalized_dtw
    data_table = data_table(data_table.DTW_Normalization == 0, :);
end

filtered = data_table;
fprintf('✓ Dataset: %d experiments\n\n', height(filtered));

% ========================================================================
%% STEP 2: SELECT CONFIGS
% ========================================================================

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  CONFIG SELECTION                                              ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

fprintf('Composite Score: 0.4×NDCG + 0.3×(1/Rank) + 0.2×R50 + 0.1×Spearman\n\n');

dtw_modes = {'joint_states', 'position'};
mode_labels = {'MOTION (Joint)', 'SHAPE (Cartesian)'};

baseline_weights = {'Joint only', 'Position only'};

fprintf('Figure 1 (Dimensionality Invariance):\n');
fprintf('  • Motion: Joint only (pure joint dimensions)\n');
fprintf('  • Space:  Position only (pure position dimensions)\n\n');

% Filter baselines for Figure 1
filtered_baselines = [];
for m = 1:length(dtw_modes)
    mode_name = dtw_modes{m};
    baseline = baseline_weights{m};
    
    mask = strcmp(filtered.DTW_Mode, mode_name) & ...
           strcmp(filtered.Weight_Mode, baseline);
    filtered_baselines = [filtered_baselines; filtered(mask, :)];
end

fprintf('✓ Figure 1: %d experiments (baselines only)\n\n', height(filtered_baselines));

% ========================================================================
%% STEP 3: FIGURE 1 - DIMENSIONALITY (BASELINES ONLY!)
% ========================================================================

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  FIGURE 1: DIMENSIONALITY INVARIANCE (BASELINES)              ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

fig1 = figure('Position', [100, 100, 900, 550], 'Color', 'w');
set(fig1, 'Name', 'Figure 1: Dimensionality Invariance (Baselines)');

color_motion = [0.86, 0.13, 0.15];
color_shape  = [0.15, 0.39, 0.91];
color_zone   = [0.09, 0.64, 0.29];

h_spearman = [];
legend_entries = {};
all_plot_data = struct([]);

% Prepare data with ACTUAL dimensions
for m = 1:length(dtw_modes)
    curr_dtw_mode = dtw_modes{m};
    baseline = baseline_weights{m};
    
    mode_data = filtered_baselines(strcmp(filtered_baselines.DTW_Mode, curr_dtw_mode), :);
    
    all_dims_per_param = unique(mode_data.Total_Dims);
    all_dims_per_param = sort(all_dims_per_param);
    
    % ACTUAL dimensions (pure!)
    if strcmp(curr_dtw_mode, 'joint_states')
        actual_dims = all_dims_per_param * 6;
        dim_label = '6×n';
    else
        actual_dims = all_dims_per_param * 3;
        dim_label = '3×n';
    end
    
    fprintf('%s (%s):\n', mode_labels{m}, baseline);
    fprintf('  Dims per param: %s\n', mat2str(all_dims_per_param'));
    fprintf('  Actual dims:    %s (%s)\n', mat2str(actual_dims'), dim_label);
    
    spearman_vals = zeros(length(all_dims_per_param), 1);
    r50_vals = zeros(length(all_dims_per_param), 1);
    
    for i = 1:length(all_dims_per_param)
        dim_data = mode_data(mode_data.Total_Dims == all_dims_per_param(i), :);
        spearman_vals(i) = mean(dim_data.Spearman_DTWvsEB);
        r50_vals(i) = mean(dim_data.('R@50_DTWvsEB'));  % DTWvsEB!
    end
    
    all_plot_data(m).actual_dims = actual_dims;
    all_plot_data(m).spearman = spearman_vals;
    all_plot_data(m).r50 = r50_vals;
    all_plot_data(m).mode = curr_dtw_mode;
    all_plot_data(m).baseline = baseline;
    
    fprintf('\n');
end

% Plot setup
yyaxis left
ylabel('Spearman \rho (Approximation Quality)', 'FontWeight', 'bold', 'FontSize', 11);
ylim([0.1, 0.3]);
ax = gca;
ax.YColor = 'k';
hold on;

% Sweet spot zones
sweet_start = 30;
sweet_end = 150;
plateau_start = 300;
plateau_end = max([all_plot_data(1).actual_dims; all_plot_data(2).actual_dims]);

fill([sweet_start, sweet_end, sweet_end, sweet_start], [0, 0, 1.2, 1.2], ...
     color_zone, 'FaceAlpha', 0.10, 'EdgeColor', 'none', 'HandleVisibility', 'off');
fill([plateau_start, plateau_end, plateau_end, plateau_start], [0, 0, 1.2, 1.2], ...
     [0.5, 0.5, 0.5], 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');

% Plot data
for m = 1:length(all_plot_data)
    if strcmp(all_plot_data(m).mode, 'joint_states')
        line_color = color_motion;
        marker_style = 'o';
        legend_label = 'Motion (6×n, Joint only)';
    else
        line_color = color_shape;
        marker_style = 's';
        legend_label = 'Space (3×n, Position only)';
    end
    
    yyaxis left
    h = plot(all_plot_data(m).actual_dims, all_plot_data(m).spearman, ...
        [marker_style '-'], 'LineWidth', 2.5, 'MarkerSize', 8, ...
        'Color', line_color, 'MarkerFaceColor', line_color);
    h_spearman = [h_spearman, h];
    legend_entries{end+1} = legend_label;
    
    yyaxis right
    plot(all_plot_data(m).actual_dims, all_plot_data(m).r50, ...
        [marker_style '--'], 'LineWidth', 2.0, 'MarkerSize', 8, ...
        'Color', line_color, 'MarkerFaceColor', 'w');
end

% Finalize
yyaxis left
set(gca, 'XScale', 'log');
xlabel('Total Embedding Dimensions (Position: 3×n, Joint: 6×n)', 'FontWeight', 'bold', 'FontSize', 11);
set(gca, 'XTick', [6, 12, 15, 30, 60, 90, 120, 180, 300, 600, 1200, 2400]);

yyaxis right
ylabel('R@50 (DTW Coverage)', 'FontWeight', 'bold', 'FontSize', 11);
ylim('auto');
ax.YColor = 'k';

grid on;
set(gca, 'FontSize', 10);

h_solid = plot(nan, nan, 'k-', 'LineWidth', 2);
h_dash  = plot(nan, nan, 'k--', 'LineWidth', 2);

legend([h_spearman(1), h_spearman(2), h_solid, h_dash], ...
       {legend_entries{1}, legend_entries{2}, 'Spearman \rho', 'R@50'}, ...
       'Location', 'southeast', 'FontSize', 10);
hold off;

% Save
if save_figures
    fig1_file = fullfile(output_folder, 'figure1_dimensionality_BASELINES');
    if strcmp(figure_format, 'png') || strcmp(figure_format, 'both')
        saveas(fig1, [fig1_file '.png']);
        fprintf('✓ Saved: %s.png\n', fig1_file);
    end
    if strcmp(figure_format, 'pdf') || strcmp(figure_format, 'both')
        saveas(fig1, [fig1_file '.pdf']);
        fprintf('✓ Saved: %s.pdf\n', fig1_file);
    end
end

fprintf('\n');

% ========================================================================
%% TABLE 1: DIMENSIONALITY ANALYSIS (BASELINES)
% ========================================================================

fprintf('Creating Table 1: Dimensionality Analysis...\n');

table1_data = cell(0, 10);

for m = 1:length(dtw_modes)
    mode_name = dtw_modes{m};
    mode_label = mode_labels{m};
    baseline = baseline_weights{m};
    
    mode_data = filtered_baselines(strcmp(filtered_baselines.DTW_Mode, mode_name), :);
    
    all_dims_per_param = unique(mode_data.Total_Dims);
    all_dims_per_param = sort(all_dims_per_param);
    
    % Calculate actual dimensions
    if strcmp(mode_name, 'joint_states')
        actual_dims = all_dims_per_param * 6;
    else
        actual_dims = all_dims_per_param * 3;
    end
    
    for i = 1:length(all_dims_per_param)
        dim = all_dims_per_param(i);
        actual_dim = actual_dims(i);
        dim_data = mode_data(mode_data.Total_Dims == dim, :);
        
        if height(dim_data) < 1
            continue;
        end
        
        % Calculate metrics
        n = height(dim_data);
        spearman_mean = mean(dim_data.Spearman_DTWvsEB);
        spearman_std = std(dim_data.Spearman_DTWvsEB);
        r50_dtw = mean(dim_data.('R@50_DTWvsEB'));
        r10_dtw = mean(dim_data.('R@10_DTWvsEB'));
        r50_gt = mean(dim_data.('R@50_GTvsEB'));
        r10_gt = mean(dim_data.('R@10_GTvsEB'));
        mean_rank_gt = mean(dim_data.Mean_GTvsEB_Rank);
        ndcg10_gt = mean(dim_data.('NDCG@10_GTvsEB'), 'omitnan');
        
        % Store
        table1_data(end+1, :) = {
            mode_label, baseline, dim, actual_dim, n, ...
            spearman_mean, spearman_std, ...
            r10_dtw, r50_dtw, mean_rank_gt
        };
    end
end

% Create table
T1 = cell2table(table1_data, 'VariableNames', ...
    {'Mode', 'Config', 'Dims_Per_Param', 'Actual_Dims', 'N', ...
     'Spearman_Mean', 'Spearman_Std', ...
     'R10_DTW', 'R50_DTW', 'MeanRank_GT'});

% Create output folder
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Save
table1_file = fullfile(output_folder, 'table1_dimensionality.csv');
writetable(T1, table1_file);
fprintf('✓ Saved: %s\n\n', table1_file);

% ========================================================================
%% STEP 4: FIGURE 2 - WEIGHT MODES (ALL MODES)
% ========================================================================

fprintf('\n╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  FIGURE 2: WEIGHT MODE CONTRIBUTION (ALL)                     ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

fig2 = figure('Position', [150, 150, 1000, 600], 'Color', 'w');
set(fig2, 'Name', 'Figure 2: Weight Mode Contribution');

for m = 1:length(dtw_modes)
    subplot(2, 1, m);
    
    mode_name = dtw_modes{m};
    mode_label = mode_labels{m};
    mode_data = filtered(strcmp(filtered.DTW_Mode, mode_name), :);
    
    % Get all weight modes
    weight_modes_unique = unique(mode_data.Weight_Mode);
    
    wm_results = struct([]);
    for w = 1:length(weight_modes_unique)
        wm = weight_modes_unique{w};
        wm_data = mode_data(strcmp(mode_data.Weight_Mode, wm), :);
        
        if height(wm_data) < 3
            continue;
        end
        
        wm_results(end+1).name = wm;
        wm_results(end).spearman = mean(wm_data.Spearman_DTWvsEB);
        wm_results(end).r50 = mean(wm_data.('R@50_DTWvsEB'));  % DTWvsEB!
        wm_results(end).rank = mean(wm_data.Mean_GTvsEB_Rank);
        
        if ismember('NDCG@10_GTvsEB', wm_data.Properties.VariableNames)
            wm_results(end).ndcg10 = mean(wm_data.('NDCG@10_GTvsEB'), 'omitnan');
        else
            wm_results(end).ndcg10 = 0;
        end
        
        % Get R@50_DTWvsEB for composite
        r50_dtw = mean(wm_data.('R@50_DTWvsEB'));
        
        wm_results(end).composite = 0.4*wm_results(end).ndcg10 + ...
                                   0.3*(1/wm_results(end).rank) + ...
                                   0.2*r50_dtw + ...
                                   0.1*wm_results(end).spearman;
    end
    
    % Manual logical sorting for incremental features
    if strcmp(mode_name, 'joint_states')
        % Motion mode order: Baseline → +Meta → +Velocity → +Orient → +Position → +All
        feature_order = {
            'Joint only',           % Baseline
            'Joint + Meta',         % +Meta
            'Joint + Velocity',     % +Velocity
            'Joint + Orient',       % +Orient
            'Joint + Position',          % +Position (cross-domain)
            'Joint + All'           % +All
        };
    else
        % Space mode order: Baseline → +Meta → +Velocity → +Orient → +Joint → +All
        feature_order = {
            'Position only',        % Baseline
            'Pos + Meta',           % +Meta
            'Pos + Velocity',       % +Velocity
            'Pos + Orient',         % +Orient
            'Pos + Joint',          % +Joint (cross-domain)
            'Pos + All'             % +All
        };
    end
    
    % Reorder wm_results according to feature_order
    sorted_indices = [];
    for f = 1:length(feature_order)
        feature_name = feature_order{f};
        idx = find(strcmp({wm_results.name}, feature_name), 1);
        if ~isempty(idx)
            sorted_indices(end+1) = idx;
        end
    end
    
    % Add any remaining modes not in the predefined order
    for w = 1:length(wm_results)
        if ~ismember(wm_results(w).name, feature_order)
            sorted_indices(end+1) = w;
        end
    end
    
    % Apply the sorting
    wm_results = wm_results(sorted_indices);
    
    % Shorten names
    x_labels = cell(length(wm_results), 1);
    for i = 1:length(wm_results)
        name = wm_results(i).name;
        name = strrep(name, 'Joint only', 'Baseline');
        name = strrep(name, 'Position only', 'Baseline');
        name = strrep(name, 'Joint + ', '+');
        name = strrep(name, 'Pos + ', '+');
        x_labels{i} = name;
    end
    
    x_pos = 1:length(wm_results);
    
    % Color
    if strcmp(mode_name, 'joint_states')
        line_color = color_motion;
        marker_style = 'o';
    else
        line_color = color_shape;
        marker_style = 's';
    end
    
    % Left: Spearman
    yyaxis left
    h_spear = plot(x_pos, [wm_results.spearman], ...
        [marker_style '-'], 'LineWidth', 2.5, 'MarkerSize', 8, ...
        'Color', line_color, 'MarkerFaceColor', line_color);
    ylabel('Spearman \rho', 'FontWeight', 'bold', 'FontSize', 11);
    ylim([0.15, 0.3]);
    ax = gca;
    ax.YColor = 'k';
    hold on;
    
    % Right: R@50 (DTW vs EB)
    yyaxis right
    h_r50 = plot(x_pos, [wm_results.r50], ...
        [marker_style '--'], 'LineWidth', 2.0, 'MarkerSize', 8, ...
        'Color', line_color, 'MarkerFaceColor', 'w');
    ylabel('R@50 (DTW Coverage)', 'FontWeight', 'bold', 'FontSize', 11);
    ylim([0.4, 0.8]);
    ax.YColor = 'k';
    
    % Styling
    yyaxis left
    set(gca, 'XTick', x_pos);
    set(gca, 'XTickLabel', x_labels);
    xtickangle(45);
    xlabel('Weight Mode', 'FontWeight', 'bold', 'FontSize', 10);
    title(mode_label, 'FontSize', 10, 'FontWeight', 'bold');
    grid on;
    set(gca, 'FontSize', 10);
    
    if m == 1
        legend([h_spear, h_r50], {'Spearman \rho (L)', 'R@50 (R)'}, ...
            'Location', 'best', 'FontSize', 10);
    end
    
    hold off;
end

% Save
if save_figures
    fig2_file = fullfile(output_folder, 'figure2_weight_modes');
    if strcmp(figure_format, 'png') || strcmp(figure_format, 'both')
        saveas(fig2, [fig2_file '.png']);
        fprintf('✓ Saved: %s.png\n', fig2_file);
    end
    if strcmp(figure_format, 'pdf') || strcmp(figure_format, 'both')
        saveas(fig2, [fig2_file '.pdf']);
        fprintf('✓ Saved: %s.pdf\n', fig2_file);
    end
end

fprintf('\n');

% ========================================================================
%% TABLE 2: WEIGHT MODE ANALYSIS (ALL MODES)
% ========================================================================

fprintf('Creating Table 2: Weight Mode Analysis...\n');

table2_data = cell(0, 12);

for m = 1:length(dtw_modes)
    mode_name = dtw_modes{m};
    mode_label = mode_labels{m};
    mode_data = filtered(strcmp(filtered.DTW_Mode, mode_name), :);
    
    weight_modes_unique = unique(mode_data.Weight_Mode);
    
    for w = 1:length(weight_modes_unique)
        wm = weight_modes_unique{w};
        wm_data = mode_data(strcmp(mode_data.Weight_Mode, wm), :);
        
        if height(wm_data) < 3
            continue;
        end
        
        % Calculate all metrics
        n = height(wm_data);
        spearman_mean = mean(wm_data.Spearman_DTWvsEB);
        spearman_std = std(wm_data.Spearman_DTWvsEB);
        r50_dtw = mean(wm_data.('R@50_DTWvsEB'));
        r10_dtw = mean(wm_data.('R@10_DTWvsEB'));
        r50_gt = mean(wm_data.('R@50_GTvsEB'));
        r10_gt = mean(wm_data.('R@10_GTvsEB'));
        mean_rank_gt = mean(wm_data.Mean_GTvsEB_Rank);
        ndcg10_gt = mean(wm_data.('NDCG@10_GTvsEB'), 'omitnan');
        
        % Composite score
        composite = 0.4*ndcg10_gt + 0.3*(1/mean_rank_gt) + 0.2*r50_gt + 0.1*spearman_mean;
        
        % Store
        table2_data(end+1, :) = {
            mode_label, wm, n, ...
            spearman_mean, spearman_std, ...
            r10_dtw, r50_dtw, ...
            r10_gt, r50_gt, mean_rank_gt, ndcg10_gt, ...
            composite
        };
    end
end

% Create table
T2 = cell2table(table2_data, 'VariableNames', ...
    {'Mode', 'Weight_Mode', 'N', ...
     'Spearman_Mean', 'Spearman_Std', ...
     'R10_DTW', 'R50_DTW', ...
     'R10_GT', 'R50_GT', 'MeanRank_GT', 'NDCG10_GT', ...
     'Composite_Score'});

% Sort by composite score within each mode
motion_rows = strcmp(T2.Mode, mode_labels{1});
shape_rows = strcmp(T2.Mode, mode_labels{2});

T2_motion = T2(motion_rows, :);
T2_shape = T2(shape_rows, :);

[~, sort_idx_motion] = sort(T2_motion.Composite_Score, 'descend');
[~, sort_idx_shape] = sort(T2_shape.Composite_Score, 'descend');

T2 = [T2_motion(sort_idx_motion, :); T2_shape(sort_idx_shape, :)];

% Save
table2_file = fullfile(output_folder, 'table2_weight_modes.csv');
writetable(T2, table2_file);
fprintf('✓ Saved: %s\n\n', table2_file);

% ========================================================================
%% COMPLETION
% ========================================================================

fprintf('\n╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  ✓ FINAL ANALYSIS COMPLETE!                                   ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

fprintf('Files Created:\n');
fprintf('───────────────────────────────────────────────────────────────\n');
fprintf('Tables:\n');
fprintf('  • table1_dimensionality.csv (Baselines, actual dimensions)\n');
fprintf('  • table2_weight_modes.csv (All modes, sorted by composite)\n\n');

fprintf('Figures:\n');
fprintf('  • figure1_dimensionality_BASELINES.%s\n', figure_format);
fprintf('  • figure2_weight_modes.%s\n\n', figure_format);

fprintf('Figure 1 (Dimensionality Invariance):\n');
fprintf('  • ONLY baselines (pure dimensions)\n');
fprintf('  • Motion: Joint only (6×n)\n');
fprintf('  • Space:  Position only (3×n)\n');
fprintf('  • X-axis: Actual dimensions\n');
fprintf('  • Y: Spearman + R@50 (DTW approximation)\n\n');

fprintf('Figure 2 (Weight Mode Contribution):\n');
fprintf('  • ALL weight modes\n');
fprintf('  • Sorted by Composite Score (NDCG=40%%, Rank=30%%)\n');
fprintf('  • Y: Spearman + R@50 (DTW approximation)\n');
fprintf('  • Shows incremental feature improvement\n\n');

fprintf('Table 1: Dimensionality sweep for baselines\n');
fprintf('Table 2: All weight modes with composite scores\n\n');
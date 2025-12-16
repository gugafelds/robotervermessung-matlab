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
export_filename = 'experiment_data.xlsx';

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
%% PAPER ANALYSIS - SECTIONS 1-3
% ========================================================================
% Complete standalone script for creating:
% - Table 1: Dimensionality Analysis
% - Table 2: Weight Mode Analysis  
% - Figure 1: Dimensionality Invariance
% - Figure 2: Weight Mode Contribution
% - Figure 3: Orientation Failure
% ========================================================================
% Author: Your Name
% Date: 2025-12-14
% ========================================================================

clear; clc; close all;

fprintf('\n');
fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║                                                                ║\n');
fprintf('║  PAPER ANALYSIS - EMBEDDING EVALUATION                         ║\n');
fprintf('║  Sections 1-3: Architecture & Feature Engineering             ║\n');
fprintf('║                                                                ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n');
fprintf('\n');

% ========================================================================
%% CONFIGURATION
% ========================================================================

% File paths
excel_file = 'results/experiment_data.xlsx';
output_folder = 'results';

% Queries to exclude (trivial 1:1 duplicates)
exclude_queries = [];

% Filter settings
use_trajectory_only = false;      % true = Trajectory level, false = include Segment
use_non_normalized_dtw = true;   % true = Non-normalized DTW only

% Figure settings
save_figures = false;
figure_format = 'both';  % 'png', 'pdf', or 'both'
figure_dpi = 300;

fprintf('Configuration:\n');
fprintf('  Excel file: %s\n', excel_file);
fprintf('  Output folder: %s\n', output_folder);
fprintf('  Excluded queries: %s\n', mat2str(exclude_queries));
fprintf('  Trajectory-only: %d\n', use_trajectory_only);
fprintf('  Non-normalized DTW: %d\n', use_non_normalized_dtw);
fprintf('\n');

% ========================================================================
%% STEP 1: LOAD AND FILTER DATA
% ========================================================================

fprintf('\n╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  STEP 1: LOADING AND FILTERING DATA                           ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

% Check if file exists
if ~isfile(excel_file)
    error('Excel file not found: %s\nPlease run the main analysis script first.', excel_file);
end

% Load data
fprintf('Loading data from Excel...\n');
data_table = readtable(excel_file, 'Sheet', 'All Data', 'VariableNamingRule', 'preserve');
fprintf('✓ Loaded %d experiments\n', height(data_table));
fprintf('  Columns: %d\n', width(data_table));
fprintf('  Unique queries: %d\n\n', length(unique(data_table.Query_Bahn_ID)));

% Filter out excluded queries
if ~isempty(exclude_queries)
    fprintf('Filtering out excluded queries...\n');
    initial_count = height(data_table);
    
    for i = 1:length(exclude_queries)
        n_excluded = sum(data_table.Query_Bahn_ID == exclude_queries(i));
        fprintf('  Query %d: %d experiments\n', exclude_queries(i), n_excluded);
    end
    
    data_table = data_table(~ismember(data_table.Query_Bahn_ID, exclude_queries), :);
    fprintf('✓ Removed %d experiments\n', initial_count - height(data_table));
    fprintf('  Remaining: %d experiments\n\n', height(data_table));
end

% Filter by level (Trajectory vs Segment)
if use_trajectory_only
    fprintf('Filtering to Trajectory-level only...\n');
    initial_count = height(data_table);
    data_table = data_table(strcmp(data_table.Level, 'Trajectory'), :);
    fprintf('✓ Removed %d Segment-level experiments\n', initial_count - height(data_table));
    fprintf('  Remaining: %d experiments\n\n', height(data_table));
end

% Filter by DTW normalization
if use_non_normalized_dtw
    fprintf('Filtering to Non-normalized DTW only...\n');
    initial_count = height(data_table);
    data_table = data_table(data_table.DTW_Normalization == 0, :);
    fprintf('✓ Removed %d Normalized DTW experiments\n', initial_count - height(data_table));
    fprintf('  Remaining: %d experiments\n\n', height(data_table));
end

% Final dataset summary
fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('FINAL FILTERED DATASET:\n');
fprintf('  Total experiments: %d\n', height(data_table));
fprintf('  Unique queries: %d\n', length(unique(data_table.Query_Bahn_ID)));
fprintf('  DTW modes: %s\n', strjoin(unique(data_table.DTW_Mode), ', '));
fprintf('  Dimensions tested: %s\n', mat2str(unique(data_table.Total_Dims)'));
fprintf('═══════════════════════════════════════════════════════════════\n\n');

filtered = data_table;  % Renamed for clarity

%% ========================================================================
%% STEP 2: TABLE 1 - DIMENSIONALITY ANALYSIS
%% ========================================================================

fprintf('\n╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  STEP 2: TABLE 1 - DIMENSIONALITY ANALYSIS                    ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

% Storage for table data
table1_data = cell(0, 9);

% Process each DTW mode
dtw_modes = {'joint_states', 'position'};
mode_labels = {'MOTION (Joint)', 'SHAPE (Cartesian)'};

for m = 1:length(dtw_modes)
    mode_name = dtw_modes{m};
    mode_label = mode_labels{m};
    
    % Filter data for this mode
    mode_data = filtered(strcmp(filtered.DTW_Mode, mode_name), :);
    
    fprintf('%s\n', mode_label);
    fprintf('──────────────────────────────────────────────────────────────────────────\n');
    fprintf('Dims │  N  │ Spearman ρ      │  R@10 (DTW)   │  R@100 (DTW)   │ Mean Rank │ DTW Rank │\n');
    fprintf('──────────────────────────────────────────────────────────────────────────\n');
    
    % Get unique dimensions
    all_dims = unique(mode_data.Total_Dims);
    all_dims = sort(all_dims);
    
    % Process each dimension
    for i = 1:length(all_dims)
        dim = all_dims(i);
        dim_data = mode_data(mode_data.Total_Dims == dim, :);
        
        if height(dim_data) < 3
            continue;  % Skip if too few samples
        end
        
        % Calculate metrics
        n = height(dim_data);
        spearman_mean = mean(dim_data.Spearman_DTWvsEB);
        spearman_std = std(dim_data.Spearman_DTWvsEB);
        r10_dtw = mean(dim_data.("R@10_DTWvsEB"));
        r100_dtw = mean(dim_data.("R@K_DTWvsEB"));
        mean_rank_eb = mean(dim_data.Mean_GTvsEB_Rank);
        mean_rank_dtw = mean(dim_data.Mean_GTvsDTW_Rank);
        
        % Mark sweet spot (10-25 dims)
        if dim >= 10 && dim <= 25
            marker = ' ★';
        else
            marker = '';
        end
        
        % Print row
        fprintf('%4d │%4d │ %5.3f ± %5.4f │ %6.3f │ %6.3f │  %6.2f   │  %5.2f   │%s\n', ...
            dim, n, spearman_mean, spearman_std, r10_dtw, r100_dtw, mean_rank_eb, mean_rank_dtw, marker);
        
        % Store for export
        table1_data(end+1, :) = {
            mode_label, dim, n, spearman_mean, spearman_std, ...
            r10_dtw, r100_dtw, mean_rank_eb, mean_rank_dtw
        };
    end
    
    % Calculate invariance statistics for 50-600 dims
    main_dims_data = mode_data(mode_data.Total_Dims >= 50, :);
    
    if height(main_dims_data) > 0
        % Group by dimension and calculate mean Spearman for each
        dims_unique = unique(main_dims_data.Total_Dims);
        spearman_by_dim = zeros(length(dims_unique), 1);
        
        for i = 1:length(dims_unique)
            dim_subset = main_dims_data(main_dims_data.Total_Dims == dims_unique(i), :);
            spearman_by_dim(i) = nanmean(dim_subset.Spearman_DTWvsEB);
        end
        
        % Statistics
        mean_spear = mean(spearman_by_dim);
        std_spear = std(spearman_by_dim);
        cv = 100 * std_spear / mean_spear;
        range_spear = max(spearman_by_dim) - min(spearman_by_dim);
        
        fprintf('\n');
        fprintf('Invariance Analysis (50-600 dims):\n');
        fprintf('  Mean Spearman:  %.4f\n', mean_spear);
        fprintf('  Std of means:   %.4f\n', std_spear);
        fprintf('  CV:             %.2f%%', cv);
        if cv < 1.0
            fprintf('  ✓ Excellent!\n');
        elseif cv < 5.0
            fprintf('  ✓ Good\n');
        else
            fprintf('\n');
        end
        fprintf('  Range:          %.4f\n', range_spear);
    end
    
    fprintf('\n');
end

% Convert to table and save
T1 = cell2table(table1_data, 'VariableNames', ...
    {'Mode', 'Dims', 'N', 'Spearman_Mean', 'Spearman_Std', ...
     'R10_GT', 'R50_GT', 'MeanRank_EB', 'MeanRank_DTW'});

% Create output folder if needed
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Save CSV
table1_file = fullfile(output_folder, 'table1_dimensionality.csv');
writetable(T1, table1_file);
fprintf('✓ Saved: %s\n\n', table1_file);

% ========================================================================
%% STEP 3: TABLE 2 - WEIGHT MODE ANALYSIS (ALL METRICS)
% ========================================================================

fprintf('\n╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  STEP 3: TABLE 2 - WEIGHT MODE ANALYSIS                       ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

% Storage for table data (expanded to include all metrics)
table2_data = cell(0, 21);  % Increased from 8 to 24 columns

% Process each DTW mode
for m = 1:length(dtw_modes)
    mode_name = dtw_modes{m};
    mode_label = mode_labels{m};
    
    % Filter data for this mode
    mode_data = filtered(strcmp(filtered.DTW_Mode, mode_name), :);
    
    fprintf('%s\n', mode_label);
    fprintf('────────────────────────────────────────────────────────────────────────────────\n');
    fprintf('Weight Mode              │  N  │ Spearman ρ      │ Δρ      │  R@10   │ Rank   │\n');
    fprintf('────────────────────────────────────────────────────────────────────────────────\n');
    
    % Get unique weight modes
    weight_modes = unique(mode_data.Weight_Mode);
    
    % Calculate metrics for each weight mode
    weight_results = struct([]);
    for i = 1:length(weight_modes)
        wm = weight_modes{i};
        wm_data = mode_data(strcmp(mode_data.Weight_Mode, wm), :);
        
        if height(wm_data) < 3
            continue;
        end
        
        % Store ALL metrics
        weight_results(end+1).name = wm;
        weight_results(end).n = height(wm_data);
        
        % DTW vs EB metrics
        weight_results(end).spearman_dtw = mean(wm_data.Spearman_DTWvsEB);
        weight_results(end).spearman_dtw_std = std(wm_data.Spearman_DTWvsEB);
        weight_results(end).rk_dtw = mean(wm_data.("R@K_DTWvsEB"));
        weight_results(end).r50_dtw = mean(wm_data.("R@50_DTWvsEB"));
        weight_results(end).r10_dtw = mean(wm_data.("R@10_DTWvsEB"));
        weight_results(end).r5_dtw = mean(wm_data.("R@5_DTWvsEB"));
        weight_results(end).r3_dtw = mean(wm_data.("R@3_DTWvsEB"));
        weight_results(end).r1_dtw = mean(wm_data.("R@1_DTWvsEB"));
        weight_results(end).p_dtw = mean(wm_data.P_DTWvsEB);
        
        % GT vs EB metrics
        weight_results(end).r50_gt = mean(wm_data.("R@50_GTvsEB"));
        weight_results(end).r10_gt = mean(wm_data.("R@10_GTvsEB"));
        weight_results(end).r5_gt = mean(wm_data.("R@5_GTvsEB"));
        weight_results(end).r3_gt = mean(wm_data.("R@3_GTvsEB"));
        weight_results(end).r1_gt = mean(wm_data.("R@1_GTvsEB"));
        weight_results(end).p_gt = mean(wm_data.P_GTvsEB);
        weight_results(end).mean_rank_gt = mean(wm_data.Mean_GTvsEB_Rank);
        weight_results(end).num_gt = mean(wm_data.Num_GT);
    end
    
    % Find baseline (contains "only")
    baseline_idx = find(contains({weight_results.name}, 'only'), 1);
    if ~isempty(baseline_idx)
        baseline_spear = weight_results(baseline_idx).spearman_dtw;
    else
        baseline_spear = 0;
        warning('No baseline (only) configuration found for %s', mode_label);
    end
    
    % Sort by Spearman descending
    [~, sort_idx] = sort([weight_results.spearman_dtw], 'descend');
    weight_results = weight_results(sort_idx);
    
    % Print table (console output remains the same)
    for i = 1:length(weight_results)
        w = weight_results(i);
        
        % Calculate delta
        delta = w.spearman_dtw - baseline_spear;
        
        % Determine delta string
        if abs(delta) < 0.001
            delta_str = '  --   ';
        else
            delta_str = sprintf('%+.3f', delta);
        end
        
        % Mark best
        if i == 1
            marker = '► ';
        else
            marker = '  ';
        end
        
        % Print row (using old variable names for display)
        fprintf('%s%-22s│%4d │ %5.3f ± %5.4f│ %s │ %6.3f │ %6.2f │\n', ...
            marker, w.name, w.n, w.spearman_dtw, w.spearman_dtw_std, ...
            delta_str, w.r10_gt, w.mean_rank_gt);
        
        % Store ALL metrics for export
        table2_data(end+1, :) = {
            mode_label, w.name, w.n, ...
            w.spearman_dtw, w.spearman_dtw_std, ...
            w.rk_dtw, w.r50_dtw, w.r10_dtw, w.r5_dtw, w.r3_dtw, w.r1_dtw, w.p_dtw, ...
            w.r50_gt, w.r10_gt, w.r5_gt, w.r3_gt, w.r1_gt, w.p_gt, ...
            w.mean_rank_gt, w.num_gt, ...
            delta
        };
    end
    
    fprintf('\n');
end

% Convert to table with ALL column names
T2 = cell2table(table2_data, 'VariableNames', ...
    {'Mode', 'Weight_Mode', 'N', ...
     'Spearman_Mean', 'Spearman_Std', ...
     'RK_DTW', 'R50_DTW', 'R10_DTW', 'R5_DTW', 'R3_DTW', 'R1_DTW', 'P_DTW', ...
     'R50_GT', 'R10_GT', 'R5_GT', 'R3_GT', 'R1_GT', 'P_GT', ...
     'MeanRank_EB', 'Num_GT', ...
     'Delta_Spearman'});

% Save CSV
table2_file = fullfile(output_folder, 'table2_weight_modes.csv');
writetable(T2, table2_file);
fprintf('✓ Saved: %s\n\n', table2_file);

% ========================================================================
%% STEP 4: FIGURE 1 - DIMENSIONALITY INVARIANCE (BEST WEIGHT MODE)
% ========================================================================
fprintf('\n╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  STEP 4: FIGURE 1 - DIMENSIONALITY INVARIANCE (OPTIMIZED)     ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

fig1 = figure('Position', [100, 100, 900, 550], 'Color', 'w');
set(fig1, 'Name', 'Figure 1: Dimensionality Invariance (Best Config)');

% Color scheme
color_motion = [0.86, 0.13, 0.15];  % Red (Joint States)
color_shape  = [0.15, 0.39, 0.91];  % Blue (Cartesian/Shape)
color_zone   = [0.09, 0.64, 0.29];  % Green (Sweet Spot)

% Arrays für Legende und Daten
h_spearman = [];
legend_entries = {};
best_modes_info = {};
all_plot_data = struct([]);

% =========================================================================
% PHASE 1: FIND BEST WEIGHT MODES (COMPOSITE SCORE)
% =========================================================================

fprintf('PHASE 1: Finding best weight modes using composite score...\n');
fprintf('═══════════════════════════════════════════════════════════════\n\n');

for m = 1:length(dtw_modes)
    curr_dtw_mode = dtw_modes{m};
    
    fprintf('Evaluating %s:\n', upper(curr_dtw_mode));
    fprintf('───────────────────────────────────────────────────────────\n');
    
    % Filter auf aktuellen DTW Mode
    base_data = filtered(strcmp(filtered.DTW_Mode, curr_dtw_mode), :);
    
    % ---------------------------------------------------------
    % COMPOSITE SCORE SELECTION (OPTION 3)
    % ---------------------------------------------------------
    available_weights = unique(base_data.Weight_Mode);
    best_weight_mode = '';
    best_score = -inf;
    
    for w = 1:length(available_weights)
        w_mode = available_weights{w};
        w_data = base_data(strcmp(base_data.Weight_Mode, w_mode), :);
        
        if height(w_data) < 3
            continue;
        end
        
        % Calculate metrics for composite score
        avg_spearman = mean(w_data.Spearman_DTWvsEB);
        avg_rank = mean(w_data.Mean_GTvsEB_Rank);
        avg_r100_dtw = mean(w_data.('R@K_DTWvsEB'));
        
        % COMPOSITE SCORE (normalized and combined)
        % Higher Spearman = better
        % Lower Rank = better → invert with 1/rank
        % Higher R@100 = better
        composite_score = avg_spearman + (1/avg_rank) + avg_r100_dtw;
        
        fprintf('  %-25s | ρ=%.3f | Rank=%.2f | R@100=%.3f | Score=%.3f\n', ...
            w_mode, avg_spearman, avg_rank, avg_r100_dtw, composite_score);
        
        if composite_score > best_score
            best_score = composite_score;
            best_weight_mode = w_mode;
            best_metrics.spearman = avg_spearman;
            best_metrics.rank = avg_rank;
            best_metrics.r100 = avg_r100_dtw;
            best_metrics.composite = composite_score;
        end
    end
    
    fprintf('───────────────────────────────────────────────────────────\n');
    fprintf('  ✓ SELECTED: %s\n', best_weight_mode);
    fprintf('      Spearman ρ:     %.3f\n', best_metrics.spearman);
    fprintf('      Mean Rank:      %.2f\n', best_metrics.rank);
    fprintf('      R@100 (DTW):    %.3f\n', best_metrics.r100);
    fprintf('      Composite:      %.3f\n\n', best_metrics.composite);
    
    best_modes_info{end+1} = best_weight_mode;
    
    % Filter zu bestem Weight Mode
    mode_data = base_data(strcmp(base_data.Weight_Mode, best_weight_mode), :);
    
    % ---------------------------------------------------------
    % BERECHNUNG DER KURVEN
    % ---------------------------------------------------------
    all_dims = unique(mode_data.Total_Dims);
    all_dims = sort(all_dims);
    
    spearman_vals = zeros(length(all_dims), 1);
    r100_dtw_vals = zeros(length(all_dims), 1);
    
    for i = 1:length(all_dims)
        dim_data = mode_data(mode_data.Total_Dims == all_dims(i), :);
        spearman_vals(i) = mean(dim_data.Spearman_DTWvsEB);
        r100_dtw_vals(i) = mean(dim_data.('R@K_DTWvsEB'));
    end
    
    % Store data for later processing
    all_plot_data(m).dims = all_dims;
    all_plot_data(m).spearman = spearman_vals;
    all_plot_data(m).r100 = r100_dtw_vals;
    all_plot_data(m).mode = curr_dtw_mode;
    all_plot_data(m).weight = best_weight_mode;
end

% =========================================================================
% PHASE 2: AUTOMATIC SWEET SPOT DETECTION
% =========================================================================

fprintf('\nPHASE 2: Automatic sweet spot detection...\n');
fprintf('═══════════════════════════════════════════════════════════════\n\n');

% Combine all Spearman values across both modes
all_spearman_combined = [];
all_dims_combined = [];

for m = 1:length(all_plot_data)
    all_spearman_combined = [all_spearman_combined; all_plot_data(m).spearman];
    all_dims_combined = [all_dims_combined; all_plot_data(m).dims];
end

% Get unique dimensions
unique_dims = unique(all_dims_combined);
unique_dims = sort(unique_dims);

% Calculate mean Spearman per dimension (across both modes)
mean_spearman_per_dim = zeros(length(unique_dims), 1);
for i = 1:length(unique_dims)
    dim = unique_dims(i);
    mask = all_dims_combined == dim;
    mean_spearman_per_dim(i) = mean(all_spearman_combined(mask));
end

% --- SWEET SPOT DETECTION LOGIC ---
% Sweet spot = dimensions where performance is HIGH and STABLE
% Strategy: Find where performance reaches 95% of maximum and stays there

max_spearman = max(mean_spearman_per_dim);
threshold = 0.95 * max_spearman;

% Find first dimension that exceeds threshold
sweet_start_idx = find(mean_spearman_per_dim >= threshold, 1, 'first');

% Sweet spot ends where plateau begins (consecutive similar values)
% Plateau = CV < 0.5% for dimensions >= 50
plateau_start = 50;
plateau_dims = unique_dims(unique_dims >= plateau_start);

if ~isempty(plateau_dims)
    % Find index of plateau start
    sweet_end_idx = find(unique_dims == plateau_start, 1, 'first') - 1;
    
    % If sweet_end is before sweet_start, just use sweet_start
    if sweet_end_idx < sweet_start_idx
        sweet_end_idx = sweet_start_idx + 2; % At least 2-3 points
    end
else
    % Fallback: sweet spot is middle range
    sweet_end_idx = floor(length(unique_dims) * 0.5);
end

% Ensure we have valid indices
if isempty(sweet_start_idx)
    sweet_start_idx = 2; % Default to second dimension
end

sweet_start_dim = unique_dims(sweet_start_idx);
sweet_end_dim = unique_dims(sweet_end_idx);

fprintf('Sweet Spot Detection Results:\n');
fprintf('  Max Spearman:       %.4f\n', max_spearman);
fprintf('  95%% Threshold:      %.4f\n', threshold);
fprintf('  Sweet Spot Range:   %d - %d dims\n', sweet_start_dim, sweet_end_dim);
fprintf('  Plateau Start:      %d dims\n\n', plateau_start);

% --- PLATEAU DETECTION ---
plateau_end_dim = max(unique_dims);

fprintf('Detected Zones:\n');
fprintf('  Sweet Spot: %d - %d dims\n', sweet_start_dim, sweet_end_dim);
fprintf('  Plateau:    %d - %d dims\n\n', plateau_start, plateau_end_dim);

% =========================================================================
% PHASE 3: PLOTTING
% =========================================================================

fprintf('PHASE 3: Creating plot...\n');
fprintf('═══════════════════════════════════════════════════════════════\n\n');

% --- SETUP AXES ---
yyaxis left
ylabel('Spearman \rho (Approximation Quality)', 'FontWeight', 'bold', 'FontSize', 11);
ylim([0.22, 0.42]); 
ax = gca;
ax.YColor = 'k'; 
hold on;

% --- BACKGROUND ZONES (AUTOMATIC) ---
y_max = 1.2;

% Sweet Spot Zone (GREEN)
fill([sweet_start_dim, sweet_end_dim, sweet_end_dim, sweet_start_dim], ...
     [0, 0, y_max, y_max], color_zone, ...
     'FaceAlpha', 0.10, 'EdgeColor', 'none', 'HandleVisibility', 'off');

% Plateau Zone (GRAY)
fill([plateau_start, plateau_end_dim, plateau_end_dim, plateau_start], ...
     [0, 0, y_max, y_max], [0.5, 0.5, 0.5], ...
     'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');

% --- PLOT DATA ---
for m = 1:length(all_plot_data)
    % Style Definition
    if strcmp(all_plot_data(m).mode, 'joint_states')
        line_color = color_motion;
        marker_style = 'o';
        legend_label = 'Motion (Joints)';
    else
        line_color = color_shape;
        marker_style = 's';
        legend_label = 'Space (Cartesian)';
    end
    
    % Left Axis: Spearman
    yyaxis left
    h = plot(all_plot_data(m).dims, all_plot_data(m).spearman, ...
        [marker_style '-'], 'LineWidth', 2.5, 'MarkerSize', 8, ...
        'Color', line_color, 'MarkerFaceColor', line_color);
    h_spearman = [h_spearman, h];
    legend_entries{end+1} = legend_label;
    
    % Right Axis: R@100
    yyaxis right
    plot(all_plot_data(m).dims, all_plot_data(m).r100, ...
        [marker_style '--'], 'LineWidth', 2.0, 'MarkerSize', 8, ...
        'Color', line_color, 'MarkerFaceColor', 'w');
end

% --- FINALIZING PLOT ---

% X-Axis
yyaxis left
set(gca, 'XScale', 'log');
xlabel('Dimensions per parameter', 'FontWeight', 'bold', 'FontSize', 11);

% Right Axis Label
yyaxis right
ylabel('R@100', 'FontWeight', 'bold', 'FontSize', 11);
ylim([0.27, 0.52]);
ax.YColor = 'k';

% Grid & Ticks
set(gca, 'XTick', [2, 5, 10, 15, 25, 50, 100, 200, 400, 600]);
set(gca, 'XTickLabel', {'2', '5', '10', '15', '25', '50', '100', '200', '400', '600'});
grid on;
set(gca, 'FontSize', 10);

% --- SMART LEGEND ---
h_solid = plot(nan, nan, 'k-', 'LineWidth', 2, 'DisplayName', 'Spearman \rho (L)');
h_dash  = plot(nan, nan, 'k--', 'LineWidth', 1.5, 'DisplayName', 'R@100 (R)');

legend([h_spearman(1), h_spearman(2), h_solid, h_dash], ...
       {'Motion (Joints)', 'Space (Cartesian)', ...
        'Spearman \rho', 'R@100'}, ...
       'Location', 'southeast', 'FontSize', 10);
hold off;

% =========================================================================
% SAVE FIGURE
% =========================================================================

if save_figures
    fig1_basename = fullfile(output_folder, 'figure1_dimensionality_invariance');
    
    if strcmp(figure_format, 'png') || strcmp(figure_format, 'both')
        saveas(fig1, [fig1_basename '.png']);
        fprintf('✓ Saved: %s.png\n', fig1_basename);
    end
    
    if strcmp(figure_format, 'pdf') || strcmp(figure_format, 'both')
        saveas(fig1, [fig1_basename '.pdf']);
        fprintf('✓ Saved: %s.pdf\n', fig1_basename);
    end
end

fprintf('\n');

% =========================================================================
% SUMMARY OUTPUT
% =========================================================================

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  SUMMARY: FIGURE 1 CONFIGURATION                              ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

fprintf('Selected Configurations (Composite Score):\n');
fprintf('─────────────────────────────────────────────────────────────────\n');
for m = 1:length(all_plot_data)
    fprintf('  %s: %s\n', upper(all_plot_data(m).mode), all_plot_data(m).weight);
end

fprintf('\nAutomatic Zone Detection:\n');
fprintf('─────────────────────────────────────────────────────────────────\n');
fprintf('  Sweet Spot: %d-%d dims (95%% of max performance)\n', ...
    sweet_start_dim, sweet_end_dim);
fprintf('  Plateau:    %d-%d dims (performance stabilized)\n', ...
    plateau_start, plateau_end_dim);
fprintf('\n');

% ========================================================================
%% STEP 5: FIGURE 2 - WEIGHT MODE CONTRIBUTION (ABSOLUTE METRICS)
% ========================================================================

fprintf('\n╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  STEP 5: FIGURE 2 - WEIGHT MODE CONTRIBUTION                  ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

fig2 = figure('Position', [150, 150, 900, 700], 'Color', 'w');
set(fig2, 'Name', 'Figure 2: Weight Mode Contribution');

% Colors (same as Figure 1)
color_motion = [0.86, 0.13, 0.15];  % Red (Joint States)
color_shape  = [0.15, 0.39, 0.91];  % Blue (Cartesian/Shape)

% --- FIND GLOBAL MIN/MAX FOR Y-LIMITS ---
all_spearman = T2.Spearman_Mean;
all_rk = T2.RK_DTW;

spearman_min = min(all_spearman);
spearman_max = max(all_spearman);
rk_min = min(all_rk);
rk_max = max(all_rk);

% Add small margin (5%)
spearman_margin = (spearman_max - spearman_min) * 0.10;
rk_margin = (rk_max - rk_min) * 0.10;

spearman_ylim = [spearman_min - spearman_margin, spearman_max + spearman_margin];
rk_ylim = [rk_min - rk_margin, rk_max + rk_margin];

fprintf('Y-axis limits determined:\n');
fprintf('  Spearman: [%.3f, %.3f]\n', spearman_ylim(1), spearman_ylim(2));
fprintf('  R@100:    [%.3f, %.3f]\n\n', rk_ylim(1), rk_ylim(2));

% --- BUILD FEATURE ORDER FROM FIRST MODE ---
mode1_data = T2(strcmp(T2.Mode, mode_labels{1}), :);

% Sort: Baseline first, then best to worst
baseline_idx = find(contains(mode1_data.Weight_Mode, 'only'), 1);
non_baseline = mode1_data(~contains(mode1_data.Weight_Mode, 'only'), :);
[~, sort_idx] = sort(non_baseline.RK_DTW, 'descend');
non_baseline = non_baseline(sort_idx, :);

if ~isempty(baseline_idx)
    sorted_mode1 = [mode1_data(baseline_idx, :); non_baseline];
else
    sorted_mode1 = non_baseline;
end

% Build feature order and labels
feature_keywords = {};
unified_labels = {};

for i = 1:height(sorted_mode1)
    wm = sorted_mode1.Weight_Mode{i};
    
    if contains(wm, 'only')
        feature_keywords{end+1} = 'only';
        unified_labels{end+1} = 'Baseline';
    elseif contains(wm, 'Position') || contains(wm, 'Pos + Joint')
        % Cross-domain: Joint + Position (for Motion) or Position + Joint (for Shape)
        if ~ismember('+Cross-Domain', unified_labels)
            feature_keywords{end+1} = 'cross';
            unified_labels{end+1} = '+Pos./Joint';
        end
    elseif contains(wm, 'Velocity')
        if ~ismember('+Velocity', unified_labels)
            feature_keywords{end+1} = 'Velocity';
            unified_labels{end+1} = '+Vel.';
        end
    elseif contains(wm, 'All')
        if ~ismember('+All', unified_labels)
            feature_keywords{end+1} = 'All';
            unified_labels{end+1} = '+All';
        end
    elseif contains(wm, 'Meta')
        if ~ismember('+Metadata', unified_labels)
            feature_keywords{end+1} = 'Meta';
            unified_labels{end+1} = '+Meta';
        end
    elseif contains(wm, 'Orient')
        if ~ismember('+Orient', unified_labels)
            feature_keywords{end+1} = 'Orient';
            unified_labels{end+1} = '+Orient.';
        end
    end
end

fprintf('Unified feature order (%d features):\n', length(unified_labels));
for i = 1:length(unified_labels)
    fprintf('  %d. %s (keyword: %s)\n', i, unified_labels{i}, feature_keywords{i});
end
fprintf('\n');

% --- PROCESS BOTH MODES WITH UNIFIED ORDER ---
for m = 1:length(dtw_modes)
    subplot(2, 1, m);
    
    % Get data for this mode
    mode_subset = T2(strcmp(T2.Mode, mode_labels{m}), :);
    
    % Map to unified features
    spearman_vals = nan(length(feature_keywords), 1);
    rk_vals = nan(length(feature_keywords), 1);
    
    for i = 1:length(feature_keywords)
        keyword = feature_keywords{i};
        
        match_idx = [];
        
        if strcmp(keyword, 'only')
            % Baseline
            match_idx = find(contains(mode_subset.Weight_Mode, 'only'), 1);
        elseif strcmp(keyword, 'cross')
            % Cross-domain: different for each mode
            if strcmp(dtw_modes{m}, 'joint_states')
                % Motion mode: look for Joint + Position
                match_idx = find(contains(mode_subset.Weight_Mode, 'Position'), 1);
            else
                % Shape mode: look for Position + Joint
                match_idx = find(contains(mode_subset.Weight_Mode, 'Joint'), 1);
            end
        else
            % Other features: direct keyword match
            match_idx = find(contains(mode_subset.Weight_Mode, keyword), 1);
        end
        
        if ~isempty(match_idx)
            spearman_vals(i) = mode_subset.Spearman_Mean(match_idx);
            rk_vals(i) = mode_subset.RK_DTW(match_idx);
            fprintf('  Mode %d, Feature %s: matched "%s" (Spearman=%.3f)\n', ...
                m, unified_labels{i}, mode_subset.Weight_Mode{match_idx}, spearman_vals(i));
        else
            fprintf('  Mode %d, Feature %s: NO MATCH\n', m, unified_labels{i});
        end
    end
    
    % X positions
    x_pos = 1:length(unified_labels);
    
    % Select color and marker
    if strcmp(dtw_modes{m}, 'joint_states')
        line_color = color_motion;
        marker_style = 'o';
    else
        line_color = color_shape;
        marker_style = 's';
    end
    
    % --- LEFT AXIS: Spearman ---
    yyaxis left
    h_spear = plot(x_pos, spearman_vals, ...
        [marker_style '-'], 'LineWidth', 2.5, 'MarkerSize', 8, ...
        'Color', line_color, 'MarkerFaceColor', line_color);
    ylabel('Spearman \rho', 'FontWeight', 'bold', 'FontSize', 11);
    ylim(spearman_ylim);
    ax = gca;
    ax.YColor = 'k';
    hold on;
    
    % --- RIGHT AXIS: R@100 ---
    yyaxis right
    h_rk = plot(x_pos, rk_vals, ...
        [marker_style '--'], 'LineWidth', 2.0, 'MarkerSize', 8, ...
        'Color', line_color, 'MarkerFaceColor', 'w');
    ylabel('R@100', 'FontWeight', 'bold', 'FontSize', 11);
    ylim(rk_ylim);
    ax.YColor = 'k';
    
    % --- STYLING ---
    yyaxis left
    set(gca, 'XTick', x_pos);
    set(gca, 'XTickLabel', unified_labels);
    xlabel('Weight Mode', 'FontWeight', 'bold', 'FontSize', 10);
    title(mode_labels{m}, 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
    set(gca, 'FontSize', 10);
    
    % Legend (only for first subplot)
    if m == 1
        legend([h_spear, h_rk], {'Spearman \rho (L)', 'R@100 (R)'}, ...
            'Location', 'best', 'FontSize', 10);
    end
    
    hold off;
end

% Save figure
if save_figures
    fig2_basename = fullfile(output_folder, 'figure2_weight_mode_contribution');
    
    if strcmp(figure_format, 'png') || strcmp(figure_format, 'both')
        saveas(fig2, [fig2_basename '.png']);
        fprintf('✓ Saved: %s.png\n', fig2_basename);
    end
    
    if strcmp(figure_format, 'pdf') || strcmp(figure_format, 'both')
        saveas(fig2, [fig2_basename '.pdf']);
        fprintf('✓ Saved: %s.pdf\n', fig2_basename);
    end
end

fprintf('\n');


% ========================================================================
%% STEP 7: SUMMARY STATISTICS
% ========================================================================

fprintf('\n╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  STEP 7: SUMMARY STATISTICS                                    ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('BEST CONFIGURATIONS (COMPOSITE SCORE)\n');
fprintf('═══════════════════════════════════════════════════════════════\n\n');

for m = 1:length(mode_labels)
    mode_weights = T2(strcmp(T2.Mode, mode_labels{m}), :);
    
    % --- COMPOSITE SCORE CALCULATION ---
    % Same as Figure 1: Spearman + (1/Rank) + R@100
    composite_scores = mode_weights.Spearman_Mean + ...
                      (1 ./ mode_weights.MeanRank_EB) + ...
                      mode_weights.RK_DTW;
    
    [best_score, best_idx] = max(composite_scores);
    best = mode_weights(best_idx, :);
    
    baseline = mode_weights(contains(mode_weights.Weight_Mode, 'only'), :);
    
    fprintf('%s:\n', mode_labels{m});
    fprintf('──────────────────────────────────────────────────────────────\n');
    fprintf('  Best Configuration: %s\n', best.Weight_Mode{1});
    fprintf('    Spearman ρ:       %.3f\n', best.Spearman_Mean);
    fprintf('    R@100 (DTW):      %.3f\n', best.RK_DTW);
    fprintf('    R@10 (GT):        %.3f\n', best.R10_GT);
    fprintf('    Mean Rank:        %.2f\n', best.MeanRank_EB);
    fprintf('    Composite Score:  %.3f\n', composite_scores(best_idx));
    
    if ~isempty(baseline)
        improvement_spearman = 100 * (best.Spearman_Mean / baseline.Spearman_Mean(1) - 1);
        improvement_r100 = 100 * (best.RK_DTW / baseline.RK_DTW(1) - 1);
        
        fprintf('\n  Baseline: %s\n', baseline.Weight_Mode{1});
        fprintf('    Spearman ρ:       %.3f\n', baseline.Spearman_Mean(1));
        fprintf('    R@100 (DTW):      %.3f\n', baseline.RK_DTW(1));
        fprintf('    Composite Score:  %.3f\n', ...
            baseline.Spearman_Mean(1) + (1/baseline.MeanRank_EB(1)) + baseline.RK_DTW(1));
        
        fprintf('\n  Improvement:\n');
        fprintf('    Spearman:         %+.1f%%\n', improvement_spearman);
        fprintf('    R@100:            %+.1f%%\n', improvement_r100);
        fprintf('    Composite:        %+.1f%%\n', ...
            100 * (composite_scores(best_idx) / ...
            (baseline.Spearman_Mean(1) + (1/baseline.MeanRank_EB(1)) + baseline.RK_DTW(1)) - 1));
    end
    fprintf('\n');
end

fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('BEST EMBEDDING CONFIGURATIONS\n');
fprintf('═══════════════════════════════════════════════════════════════\n\n');

% Analyze filtered data (which should have Embedding_Config column)
if ismember('Embedding_Config', filtered.Properties.VariableNames)
    
    for m = 1:length(dtw_modes)
        mode_name = dtw_modes{m};
        mode_label = mode_labels{m};
        
        % Get data for this mode
        mode_data = filtered(strcmp(filtered.DTW_Mode, mode_name), :);
        
        % Get unique embedding configs
        emb_configs = unique(mode_data.Embedding_Config);
        
        fprintf('%s:\n', mode_label);
        fprintf('──────────────────────────────────────────────────────────────\n');
        
        % Calculate composite score for each config
        config_results = struct([]);
        for i = 1:length(emb_configs)
            config = emb_configs{i};
            config_data = mode_data(strcmp(mode_data.Embedding_Config, config), :);
            
            if height(config_data) < 3
                continue;
            end
            
            avg_spearman = mean(config_data.Spearman_DTWvsEB);
            avg_rank = mean(config_data.Mean_GTvsEB_Rank);
            avg_r100 = mean(config_data.('R@K_DTWvsEB'));
            
            composite = avg_spearman + (1/avg_rank) + avg_r100;
            
            config_results(end+1).name = config;
            config_results(end).n = height(config_data);
            config_results(end).spearman = avg_spearman;
            config_results(end).r100 = avg_r100;
            config_results(end).rank = avg_rank;
            config_results(end).composite = composite;
        end
        
        % Sort by composite score
        [~, sort_idx] = sort([config_results.composite], 'descend');
        config_results = config_results(sort_idx);
        
        % Print top 5
        fprintf('  Top 5 Embedding Configurations:\n\n');
        for i = 1:min(5, length(config_results))
            c = config_results(i);
            fprintf('  %d. %s (N=%d)\n', i, c.name, c.n);
            fprintf('     Spearman ρ:      %.3f\n', c.spearman);
            fprintf('     R@100:           %.3f\n', c.r100);
            fprintf('     Mean Rank:       %.2f\n', c.rank);
            fprintf('     Composite Score: %.3f\n', c.composite);
            
            if i == 1
                fprintf('     ★ BEST OVERALL\n');
            end
            fprintf('\n');
        end
    end
    
else
    fprintf('⚠ Warning: Embedding_Config column not found in filtered data\n');
    fprintf('  Skipping embedding configuration analysis\n\n');
end

fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('DIMENSION INVARIANCE (50-600 dims)\n');
fprintf('═══════════════════════════════════════════════════════════════\n\n');

for m = 1:length(mode_labels)
    mode_dims = T1(strcmp(T1.Mode, mode_labels{m}), :);
    main = mode_dims(mode_dims.Dims >= 50, :);
    
    if height(main) > 0
        cv = 100 * std(main.Spearman_Mean) / mean(main.Spearman_Mean);
        
        fprintf('%s:\n', mode_labels{m});
        fprintf('  Mean Spearman:  %.4f\n', mean(main.Spearman_Mean));
        fprintf('  Std of means:   %.4f\n', std(main.Spearman_Mean));
        fprintf('  CV:             %.2f%%', cv);
        if cv < 1.0
            fprintf('  ✓ Excellent!\n');
        else
            fprintf('\n');
        end
        fprintf('  Range:          %.4f\n', max(main.Spearman_Mean) - min(main.Spearman_Mean));
        fprintf('\n');
    end
end

fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('ORIENTATION FAILURE ANALYSIS\n');
fprintf('═══════════════════════════════════════════════════════════════\n\n');

orient_rows = T2(contains(T2.Weight_Mode, 'Orient'), :);

for i = 1:height(orient_rows)
    row = orient_rows(i, :);
    baseline = T2(strcmp(T2.Mode, row.Mode{1}) & contains(T2.Weight_Mode, 'only'), :);
    
    fprintf('%s - %s:\n', row.Mode{1}, row.Weight_Mode{1});
    fprintf('  Spearman:  %.3f  (baseline: %.3f, Δ: %.3f)\n', ...
        row.Spearman_Mean, baseline.Spearman_Mean(1), row.Delta_Spearman);
    fprintf('  R@100:     %.3f  (baseline: %.3f, Δ: %.3f)\n', ...
        row.RK_DTW, baseline.RK_DTW(1), row.RK_DTW - baseline.RK_DTW(1));
    fprintf('  R@10:      %.3f  (baseline: %.3f)\n', ...
        row.R10_GT, baseline.R10_GT(1));
    fprintf('  Mean Rank: %.2f  (baseline: %.2f)\n', ...
        row.MeanRank_EB, baseline.MeanRank_EB(1));
    
    r10_drop = 100 * (1 - row.R10_GT / baseline.R10_GT(1));
    spearman_drop = 100 * (1 - row.Spearman_Mean / baseline.Spearman_Mean(1));
    fprintf('  💥 R@10 drops by %.0f%%!\n', r10_drop);
    fprintf('  💥 Spearman drops by %.0f%%!\n\n', spearman_drop);
end
% ========================================================================
%% COMPLETION
% ========================================================================

fprintf('\n╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║                                                                ║\n');
fprintf('║  ✓ ANALYSIS COMPLETE!                                          ║\n');
fprintf('║                                                                ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

fprintf('Files created in %s/:\n', output_folder);
fprintf('  Tables:\n');
fprintf('    - table1_dimensionality.csv\n');
fprintf('    - table2_weight_modes.csv\n');
fprintf('  Figures:\n');
fprintf('    - figure1_dimensionality_invariance.%s\n', figure_format);
fprintf('    - figure2_weight_mode_contribution.%s\n', figure_format);
fprintf('    - figure3_orientation_failure.%s\n', figure_format);
fprintf('\n');

fprintf('You can now customize the figures and tables as needed!\n\n');
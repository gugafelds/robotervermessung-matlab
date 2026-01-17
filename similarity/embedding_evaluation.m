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
        fprintf('  Reading existing data from sheet "All Data - Phase 1"...\n');
        existing_data = readtable(export_filename, 'Sheet', 'All Data - Phase 1', 'VariableNamingRule', 'preserve');
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
                writetable(combined_data, export_filename, 'Sheet', 'All Data - Phase 1', 'WriteRowNames', false);
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
        writetable(new_data, export_filename, 'Sheet', 'All Data - Phase 1', 'WriteRowNames', false);
        fprintf('  ✓ Created new Excel file with %d rows\n\n', height(new_data));
    end
    
else
    % File doesn't exist - create new
    fprintf('ℹ File does not exist - creating new Excel file...\n');
    fprintf('  Location: %s\n\n', export_filename);
    
    try
        writetable(new_data, export_filename, 'Sheet', 'All Data - Phase 1', 'WriteRowNames', false);
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

fprintf('Note: File contains ONLY "All Data - Phase 1" sheet.\n');
fprintf('      Trajectory/Segment filtering can be done in Excel.\n\n');

%% STEP 1.6: Smart Excel Export - Phase 2 (Two-Stage Retrieval)
% ========================================================================

fprintf('\nSTEP 1.6: Smart Excel Export - Phase 2 (Two-Stage)...\n');
fprintf('────────────────────────────────────────\n\n');

% Load Phase 2 CSV files
phase2_files = dir('results/two_stage_results_*.csv');

if isempty(phase2_files)
    fprintf('ℹ No Phase 2 CSV files found - skipping\n\n');
else
    fprintf('Found %d Phase 2 CSV file(s)\n', length(phase2_files));
    
    % Load and combine all Phase 2 CSVs
    phase2_data_all = [];
    for i = 1:length(phase2_files)
        filepath = fullfile(phase2_files(i).folder, phase2_files(i).name);
        temp_data = readtable(filepath, 'VariableNamingRule', 'preserve');
        
        % Extract timestamp from filename (two_stage_results_TIMESTAMP.csv)
        filename = phase2_files(i).name;
        timestamp_match = regexp(filename, 'two_stage_results_(.+)\.csv', 'tokens');
        if ~isempty(timestamp_match)
            timestamp_str = timestamp_match{1}{1};
            temp_data.Timestamp = repmat({timestamp_str}, height(temp_data), 1);
        end
        
        phase2_data_all = [phase2_data_all; temp_data];
        fprintf('  Loaded: %s (%d rows)\n', phase2_files(i).name, height(temp_data));
    end
    
    fprintf('\nTotal Phase 2 data loaded: %d rows\n', height(phase2_data_all));
    
    % Sheet name for Phase 2
    sheet_name_phase2 = 'All Data - Phase 2';
    
    %% Check if Phase 2 sheet exists
    if exist(export_filename, 'file')
        try
            fprintf('\n  Reading existing Phase 2 data from sheet "%s"...\n', sheet_name_phase2);
            existing_data_p2 = readtable(export_filename, 'Sheet', sheet_name_phase2, 'VariableNamingRule', 'preserve');
            fprintf('  ✓ Existing Phase 2 data: %d rows\n', height(existing_data_p2));
            
            % Check timestamps
            if ismember('Timestamp', existing_data_p2.Properties.VariableNames) && ...
               ismember('Timestamp', phase2_data_all.Properties.VariableNames)
                
                existing_timestamps_p2 = unique(existing_data_p2.Timestamp);
                new_timestamps_p2 = unique(phase2_data_all.Timestamp);
                
                fprintf('\n  Existing Phase 2 timestamps: %d unique\n', length(existing_timestamps_p2));
                fprintf('  New Phase 2 timestamps: %d unique\n', length(new_timestamps_p2));
                
                % Find new timestamps
                timestamps_to_add_p2 = setdiff(new_timestamps_p2, existing_timestamps_p2);
                
                if isempty(timestamps_to_add_p2)
                    fprintf('\n  ℹ All Phase 2 data already exists - nothing to add!\n\n');
                else
                    fprintf('\n  → Found %d NEW Phase 2 timestamp(s) to add:\n', length(timestamps_to_add_p2));
                    for i = 1:length(timestamps_to_add_p2)
                        fprintf('      %d. %s\n', i, char(timestamps_to_add_p2(i)));
                    end
                    
                    % Filter new data
                    if iscellstr(timestamps_to_add_p2) || isstring(timestamps_to_add_p2)
                        mask_to_add_p2 = ismember(phase2_data_all.Timestamp, timestamps_to_add_p2);
                    else
                        mask_to_add_p2 = false(height(phase2_data_all), 1);
                        for i = 1:length(timestamps_to_add_p2)
                            mask_to_add_p2 = mask_to_add_p2 | strcmp(phase2_data_all.Timestamp, timestamps_to_add_p2(i));
                        end
                    end
                    
                    data_to_add_p2 = phase2_data_all(mask_to_add_p2, :);
                    
                    fprintf('  → Rows to add: %d\n\n', height(data_to_add_p2));
                    
                    % Append
                    combined_data_p2 = [existing_data_p2; data_to_add_p2];
                    
                    fprintf('  Combined Phase 2 data: %d rows (was %d, added %d)\n\n', ...
                        height(combined_data_p2), height(existing_data_p2), height(data_to_add_p2));
                    
                    % Write
                    fprintf('  Writing updated Phase 2 data to Excel...\n');
                    writetable(combined_data_p2, export_filename, 'Sheet', sheet_name_phase2, 'WriteRowNames', false);
                    fprintf('  ✓ Successfully appended Phase 2 data!\n\n');
                    
                    fprintf('╔════════════════════════════════════════════════════════════════╗\n');
                    fprintf('║  PHASE 2 DATA APPENDED                                         ║\n');
                    fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');
                    fprintf('Sheet: %s\n', sheet_name_phase2);
                    fprintf('  Previous rows: %d\n', height(existing_data_p2));
                    fprintf('  Added rows: %d\n', height(data_to_add_p2));
                    fprintf('  Total rows now: %d\n\n', height(combined_data_p2));
                end
                
            else
                fprintf('  ⚠ Timestamp column not found - cannot determine new data\n\n');
            end
            
        catch ME
            fprintf('  ℹ Sheet "%s" does not exist yet - creating new\n', sheet_name_phase2);
            
            % Create new sheet
            writetable(phase2_data_all, export_filename, 'Sheet', sheet_name_phase2, 'WriteRowNames', false);
            fprintf('  ✓ Created new Phase 2 sheet with %d rows\n\n', height(phase2_data_all));
            
            fprintf('╔════════════════════════════════════════════════════════════════╗\n');
            fprintf('║  PHASE 2 SHEET CREATED                                         ║\n');
            fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');
            fprintf('Sheet: %s\n', sheet_name_phase2);
            fprintf('  Total rows: %d\n\n', height(phase2_data_all));
        end
        
    else
        fprintf('  ⚠ Excel file does not exist - Phase 1 must run first\n\n');
    end
end

%% STEP 1.7: Smart Excel Export - Phase 3 (pgvector Two-Stage)
% ========================================================================

fprintf('\nSTEP 1.7: Smart Excel Export - Phase 3 (pgvector)...\n');
fprintf('────────────────────────────────────────\n\n');

% Load Phase 3 CSV files
phase3_files = dir('results/two_stage_pgvector_*.csv');

if isempty(phase3_files)
    fprintf('ℹ No Phase 3 CSV files found - skipping\n\n');
else
    fprintf('Found %d Phase 3 CSV file(s)\n', length(phase3_files));
    
    % Load and combine all Phase 3 CSVs
    phase3_data_all = [];
    
    for i = 1:length(phase3_files)
        filepath = fullfile(phase3_files(i).folder, phase3_files(i).name);
        temp_data = readtable(filepath, 'VariableNamingRule', 'preserve');
        
        % Force segment_id to be cell array (string)
        if ismember('segment_id', temp_data.Properties.VariableNames)
            if ~iscell(temp_data.segment_id)
                temp_data.segment_id = cellstr(string(temp_data.segment_id));
            end
        end
        
        % Force bahn_id to be cell array (string) - falls gleiche Problem
        if ismember('bahn_id', temp_data.Properties.VariableNames)
            if ~iscell(temp_data.bahn_id)
                temp_data.bahn_id = cellstr(string(temp_data.bahn_id));
            end
        end
        
        % Extract timestamp from filename
        filename = phase3_files(i).name;
        timestamp_match = regexp(filename, 'two_stage_pgvector_(.+)\.csv', 'tokens');
        if ~isempty(timestamp_match)
            timestamp_str = timestamp_match{1}{1};
            temp_data.Timestamp = repmat({timestamp_str}, height(temp_data), 1);
        end
        
        if i == 1
            phase3_data_all = temp_data;
        else
            phase3_data_all = [phase3_data_all; temp_data];
        end
        
        fprintf('  Loaded: %s (%d rows)\n', phase3_files(i).name, height(temp_data));
    end
    
    fprintf('\nTotal Phase 3 data loaded: %d rows\n', height(phase3_data_all));
    
    % Sheet name for Phase 3
    sheet_name_phase3 = 'All Data - Phase 3';
    
    %% Check if Phase 3 sheet exists
    if exist(export_filename, 'file')
        try
            fprintf('\n  Reading existing Phase 3 data from sheet "%s"...\n', sheet_name_phase3);
            existing_data_p3 = readtable(export_filename, 'Sheet', sheet_name_phase3, 'VariableNamingRule', 'preserve');
            fprintf('  ✓ Existing Phase 3 data: %d rows\n', height(existing_data_p3));
            
            % Check timestamps
            if ismember('Timestamp', existing_data_p3.Properties.VariableNames) && ...
               ismember('Timestamp', phase3_data_all.Properties.VariableNames)
                
                existing_timestamps_p3 = unique(existing_data_p3.Timestamp);
                new_timestamps_p3 = unique(phase3_data_all.Timestamp);
                
                fprintf('\n  Existing Phase 3 timestamps: %d unique\n', length(existing_timestamps_p3));
                fprintf('  New Phase 3 timestamps: %d unique\n', length(new_timestamps_p3));
                
                % Find new timestamps
                timestamps_to_add_p3 = setdiff(new_timestamps_p3, existing_timestamps_p3);
                
                if isempty(timestamps_to_add_p3)
                    fprintf('\n  ℹ All Phase 3 data already exists - nothing to add!\n\n');
                else
                    fprintf('\n  → Found %d NEW Phase 3 timestamp(s) to add:\n', length(timestamps_to_add_p3));
                    for i = 1:length(timestamps_to_add_p3)
                        fprintf('      %d. %s\n', i, char(timestamps_to_add_p3(i)));
                    end
                    
                    % Filter new data
                    if iscellstr(timestamps_to_add_p3) || isstring(timestamps_to_add_p3)
                        mask_to_add_p3 = ismember(phase3_data_all.Timestamp, timestamps_to_add_p3);
                    else
                        mask_to_add_p3 = false(height(phase3_data_all), 1);
                        for i = 1:length(timestamps_to_add_p3)
                            mask_to_add_p3 = mask_to_add_p3 | strcmp(phase3_data_all.Timestamp, timestamps_to_add_p3(i));
                        end
                    end
                    
                    data_to_add_p3 = phase3_data_all(mask_to_add_p3, :);
                    
                    fprintf('  → Rows to add: %d\n\n', height(data_to_add_p3));
                    
                    % Append
                    combined_data_p3 = [existing_data_p3; data_to_add_p3];
                    
                    fprintf('  Combined Phase 3 data: %d rows (was %d, added %d)\n\n', ...
                        height(combined_data_p3), height(existing_data_p3), height(data_to_add_p3));
                    
                    % Write
                    fprintf('  Writing updated Phase 3 data to Excel...\n');
                    writetable(combined_data_p3, export_filename, 'Sheet', sheet_name_phase3, 'WriteRowNames', false);
                    fprintf('  ✓ Successfully appended Phase 3 data!\n\n');
                    
                    fprintf('╔════════════════════════════════════════════════════════════════╗\n');
                    fprintf('║  PHASE 3 DATA APPENDED                                         ║\n');
                    fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');
                    fprintf('Sheet: %s\n', sheet_name_phase3);
                    fprintf('  Previous rows: %d\n', height(existing_data_p3));
                    fprintf('  Added rows: %d\n', height(data_to_add_p3));
                    fprintf('  Total rows now: %d\n\n', height(combined_data_p3));
                end
                
            else
                fprintf('  ⚠ Timestamp column not found - cannot determine new data\n\n');
            end
            
        catch ME
            fprintf('  ℹ Sheet "%s" does not exist yet - creating new\n', sheet_name_phase3);
            
            % Create new sheet
            writetable(phase3_data_all, export_filename, 'Sheet', sheet_name_phase3, 'WriteRowNames', false);
            fprintf('  ✓ Created new Phase 3 sheet with %d rows\n\n', height(phase3_data_all));
            
            fprintf('╔════════════════════════════════════════════════════════════════╗\n');
            fprintf('║  PHASE 3 SHEET CREATED                                         ║\n');
            fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');
            fprintf('Sheet: %s\n', sheet_name_phase3);
            fprintf('  Total rows: %d\n\n', height(phase3_data_all));
        end
        
    else
        fprintf('  ⚠ Excel file does not exist - Phase 1 must run first\n\n');
    end
end

fprintf('Note: Excel file now contains:\n');
fprintf('  • All Data - Phase 1: Embedding Validation\n');
fprintf('  • All Data - Phase 2: Two-Stage Retrieval (MATLAB)\n');
fprintf('  • All Data - Phase 3: Two-Stage Retrieval (pgvector)\n\n');

fprintf('Note: Excel file now contains:\n');
fprintf('  • All Data - Phase 1: Embedding Validation)\n');
fprintf('  • All Data - Phase 2: (Phase 2: Two-Stage Retrieval)\n\n');

% ========================================================================
%% FINAL CORRECTED PAPER ANALYSIS - PHASE 1
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
figure_folder = 'figs';
use_trajectory_only = false;
use_segment_only = false;
use_non_normalized_dtw = true;
save_figures = true;
figure_format = 'both';

% ========================================================================
%% STEP 1: LOAD AND FILTER
% ========================================================================

fprintf('Loading data...\n');
data_table = readtable(excel_file, 'Sheet', 'All Data - Phase 1', 'VariableNamingRule', 'preserve');

if use_trajectory_only
    data_table = data_table(strcmp(data_table.Level, 'Trajectory'), :);
end
if use_segment_only
    data_table = data_table(strcmp(data_table.Level, 'Segment'), :);
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
%% STEP 3: FIGURE 1 - DIMENSIONALITY (BASELINES ONLY!) - COMPOSITE SCORE
% ========================================================================

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  FIGURE 1: DIMENSIONALITY INVARIANCE (BASELINES)              ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

% Filter für Database_Size = 1000
filtered_baselines_1000 = filtered_baselines(filtered_baselines.Database_Size == 1000, :);

fig1 = figure('Position', [100, 100, 750, 500], 'Color', 'w');
set(fig1, 'Name', 'Figure 1: Dimensionality Invariance (Baselines) - Composite');

color_motion = [0.86, 0.13, 0.15];
color_shape  = [0.15, 0.39, 0.91];
color_zone   = [0.09, 0.64, 0.29];

h_lines = [];
legend_entries = {};
all_plot_data = struct([]);

% Prepare data with ACTUAL dimensions
for m = 1:length(dtw_modes)
    curr_dtw_mode = dtw_modes{m};
    baseline = baseline_weights{m};
    
    mode_data = filtered_baselines_1000(strcmp(filtered_baselines_1000.DTW_Mode, curr_dtw_mode), :);
    
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
    
    composite_vals = zeros(length(all_dims_per_param), 1);
    
    for i = 1:length(all_dims_per_param)
        dim_data = mode_data(mode_data.Total_Dims == all_dims_per_param(i), :);
        
        % Berechne alle Metriken für Composite Score
        spearman_mean = mean(dim_data.Spearman_DTWvsEB);
        r50_dtw = mean(dim_data.('R@50_DTWvsEB'));
        r50_gt = mean(dim_data.('R@50_GTvsEB'));
        mean_rank_gt = mean(dim_data.Mean_GTvsEB_Rank);
        ndcg50_gt = mean(dim_data.('NDCG@50_GTvsEB'), 'omitnan');
        
        % Composite Score (gleiche Formel wie Tabelle 1)
        composite_vals(i) = calculate_composite_score(ndcg50_gt, r50_dtw, mean_rank_gt, r50_gt, spearman_mean);

    end
    
    all_plot_data(m).actual_dims = actual_dims;
    all_plot_data(m).composite = composite_vals;
    all_plot_data(m).mode = curr_dtw_mode;
    all_plot_data(m).baseline = baseline;
    
    fprintf('  Composite scores: %s\n', mat2str(composite_vals', 3));
    fprintf('\n');
end

% ========================================================================
% MANUAL Sweet Spot Configuration
% ========================================================================
sweet_start = 30;
sweet_end = 150;
plateau_start = 150;
plateau_end = max([all_plot_data(1).actual_dims; all_plot_data(2).actual_dims]);

fprintf('=== Manual Sweet Spot Configuration ===\n');
fprintf('Sweet Spot Zone: [%d, %d] dims\n', sweet_start, sweet_end);
fprintf('Plateau Zone: [%d, %d] dims\n\n', plateau_start, plateau_end);

% ========================================================================
% Plot
% ========================================================================
hold on;

% Get Y-axis range für Zonen (Auto-Skalierung)
all_scores = [all_plot_data(1).composite; all_plot_data(2).composite];
y_min = min(all_scores) * 0.99;
y_max = max(all_scores) * 1.01;

% Sweet spot zones (manuell konfiguriert!)
fill([sweet_start, sweet_end, sweet_end, sweet_start], [y_min, y_min, y_max, y_max], ...
     color_zone, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
fill([plateau_start, plateau_end, plateau_end, plateau_start], [y_min, y_min, y_max, y_max], ...
     [0.5, 0.5, 0.5], 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');

% Plot data
for m = 1:length(all_plot_data)
    if strcmp(all_plot_data(m).mode, 'joint_states')
        line_color = color_motion;
        marker_style = 'o';
        legend_label = 'Motion (6×n - Joint only)';
    else
        line_color = color_shape;
        marker_style = 's';
        legend_label = 'Space (3×n - Position only)';
    end
    
    h = plot(all_plot_data(m).actual_dims, all_plot_data(m).composite, ...
        [marker_style '-'], 'LineWidth', 3, 'MarkerSize', 11, ...  % Dickere Linie, größere Marker
        'Color', line_color, 'MarkerFaceColor', line_color);
    h_lines = [h_lines, h];
    legend_entries{end+1} = legend_label;
end

% Styling mit größeren Schriften
set(gca, 'XScale', 'log');
xlabel('Total embedding dimensions', 'FontWeight', 'bold', 'FontSize', 20);
ylabel('Composite score', 'FontWeight', 'bold', 'FontSize', 20);
set(gca, 'XTick', [6, 15, 30, 75, 150, 300, 600, 1200, 3600]);
set(gca, 'FontName', 'Courier New');
yticks([0.67,0.68,0.69,0.70,0.71,0.72,0.73])
yticklabels({'0.67','0.68','0.69','0.70','0.71','0.72','0.73'});
xtickangle(0);
set(gca, 'FontName', 'courier new', 'FontWeight', 'bold', 'FontSize', 18);

% Auto-Skalierung der Y-Achse
ylim([y_min, y_max]);


% Kompaktere Legend
legend(h_lines, legend_entries, 'Location', 'southeast', 'FontSize', 16);

% Engere Margins
ax = gca;
ax.Position = [0.14 0.14 0.84 0.84];
ax.LooseInset = [0, 0, 0, 0];  % Minimiert weiße Ränder
ax.YGrid = 'on';
ax.XGrid = 'off';

hold off;

% Save
if save_figures
    fig1_file = fullfile(figure_folder, 'dimensionality');
    if strcmp(figure_format, 'pdf') || strcmp(figure_format, 'both')
        exportgraphics(gca, [fig1_file '.pdf']);
        fprintf('✓ Saved: %s.pdf\n', fig1_file);
    end
end

fprintf('\n');
% ========================================================================
%% TABLE 1: DIMENSIONALITY ANALYSIS (BASELINES)
% ========================================================================
fprintf('Creating Table 1: Dimensionality Analysis...\n');
table1_data = cell(0, 10);  % Eine zusätzliche Spalte für Composite Score

filtered_baselines_1000 = filtered_baselines(filtered_baselines.Database_Size == 1000, :);

for m = 1:length(dtw_modes)
    mode_name = dtw_modes{m};
    mode_label = mode_labels{m};
    baseline = baseline_weights{m};
    mode_data = filtered_baselines_1000(strcmp(filtered_baselines_1000.DTW_Mode, mode_name), :);
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
        ndcg50_gt = mean(dim_data.('NDCG@50_GTvsEB'), 'omitnan');

        % Composite score (gleiche Gewichtung wie Tabelle 2)
        composite = calculate_composite_score(ndcg50_gt, r50_dtw, mean_rank_gt, r50_gt, spearman_mean);

        % Store
        table1_data(end+1, :) = {
            mode_label, baseline, dim, actual_dim, n, ...
            spearman_mean, ...
            r50_dtw, r50_gt, ndcg50_gt, ...
            composite
        };
    end
end

% Create table
T1 = cell2table(table1_data, 'VariableNames', ...
    {'Mode', 'Config', 'Dims_Per_Param', 'Actual_Dims', 'N', ...
    'Spearman', ...
    'R@50 (DTW)', 'R@50 (GT)', 'NDCG@50 (GT)', ...
    'Composite_Score'});

% Sort by composite score within each mode
motion_rows = strcmp(T1.Mode, mode_labels{1});
shape_rows = strcmp(T1.Mode, mode_labels{2});
T1_motion = T1(motion_rows, :);
T1_shape = T1(shape_rows, :);
[~, sort_idx_motion] = sort(T1_motion.Composite_Score, 'descend');
[~, sort_idx_shape] = sort(T1_shape.Composite_Score, 'descend');
T1 = [T1_motion(sort_idx_motion, :); T1_shape(sort_idx_shape, :)];

% Create output folder
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Save
table1_file = fullfile(output_folder, 'table1_dimensionality.csv');
writetable(T1, table1_file);
fprintf('✓ Saved: %s\n\n', table1_file);

% ========================================================================
%% STEP 4: FIGURE 2 - WEIGHT MODES (COMBINED) - COMPOSITE SCORE
% ========================================================================

fprintf('\n╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  FIGURE 2: WEIGHT MODE CONTRIBUTION (COMBINED)                ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

% Filter für Database_Size = 1000
filtered_1000 = filtered(filtered.Database_Size == 1000, :);

fig2 = figure('Position', [100, 100, 750, 500], 'Color', 'w');
set(fig2, 'Name', 'Figure 2: Weight Mode Contribution - Composite');

hold on;

all_composite_vals = [];
h_lines = [];
legend_entries = {};

for m = 1:length(dtw_modes)
    mode_name = dtw_modes{m};
    mode_label = mode_labels{m};
    mode_data = filtered_1000(strcmp(filtered_1000.DTW_Mode, mode_name), :);
    
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
        
        % Berechne alle Metriken für Composite Score
        spearman_mean = mean(wm_data.Spearman_DTWvsEB);
        r50_dtw = mean(wm_data.('R@50_DTWvsEB'));
        r50_gt = mean(wm_data.('R@50_GTvsEB'));
        mean_rank_gt = mean(wm_data.Mean_GTvsEB_Rank);
        
        if ismember('NDCG@10_GTvsEB', wm_data.Properties.VariableNames)
            ndcg10_gt = mean(wm_data.('NDCG@10_GTvsEB'), 'omitnan');
        else
            ndcg10_gt = 0;
        end
        
        % Composite Score (gleiche Formel wie Tabelle 2)
        wm_results(end).composite = calculate_composite_score(ndcg50_gt, r50_dtw, mean_rank_gt, r50_gt, spearman_mean);

    end
    
    % Manual logical sorting for incremental features
    if strcmp(mode_name, 'joint_states')
        % Motion mode order: Baseline → +Meta → +Velocity → +Orient → +Position → +All
        feature_order = {
            'Joint only',           % Baseline
            'Joint + Meta',         % +Meta
            'Joint + Velocity',     % +Velocity
            'Joint + Orient',       % +Orient
            'Joint + Position',     % +Position (cross-domain)
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
    
    % Shorten names (only for this mode)
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
    
    % Color and style
    if strcmp(mode_name, 'joint_states')
        line_color = color_motion;
        marker_style = 'o';
        legend_label = 'Motion';
    else
        line_color = color_shape;
        marker_style = 's';
        legend_label = 'Space';
    end
    
    % Plot Composite Score
    composite_vals = [wm_results.composite];
    all_composite_vals = [all_composite_vals; composite_vals'];
    
    h = plot(x_pos, composite_vals, ...
        [marker_style '-'], 'LineWidth', 3, 'MarkerSize', 11, ...
        'Color', line_color, 'MarkerFaceColor', line_color);
    
    h_lines = [h_lines, h];
    legend_entries{end+1} = legend_label;
end

% Auto-Skalierung
y_min = min(all_composite_vals) * 0.98;
y_max = max(all_composite_vals) * 1.02;

ylabel('Composite score', 'FontWeight', 'bold', 'FontSize', 20);
ylim([y_min, y_max]);

% X-axis: Generische Labels (1-6)
set(gca, 'XTick', 1:6);
set(gca, 'XTickLabel', {'Baseline', '+Meta.', '+Vel.', '+Orient.', '+Cross', '+All'});
xlabel('Incremental feature addition', 'FontWeight', 'bold', 'FontSize', 20);
yticks([0.60,0.64,0.68,0.72,0.76]);
yticklabels({'0.60','0.64','0.68','0.72','0.76'});
xtickangle(0);
set(gca, 'FontName', 'courier new');


grid on;
set(gca, 'FontSize', 18, 'FontWeight', 'bold');
ax = gca;
ax.Position = [0.14 0.14 0.82 0.84];
ax.LooseInset = [0, 0, 0, 0];  % Minimiert weiße Ränder


% Legend
legend(h_lines, legend_entries, 'Location', 'southeast', 'FontSize', 16);

hold off;

% Save
if save_figures
    fig2_file = fullfile(figure_folder, 'incremental_feature');
    if strcmp(figure_format, 'pdf') || strcmp(figure_format, 'both')
        exportgraphics(gca, 'figs/incremental_add.pdf');
        fprintf('✓ Saved: %s.pdf\n', fig2_file);
    end
end

fprintf('\n');

% ========================================================================
%% TABLE 2: WEIGHT MODE ANALYSIS (ALL MODES)
% ========================================================================

fprintf('Creating Table 2: Weight Mode Analysis...\n');

table2_data = cell(0, 12);

filtered_1000 = filtered(filtered.Database_Size == 1000, :);

for m = 1:length(dtw_modes)
    mode_name = dtw_modes{m};
    mode_label = mode_labels{m};
    mode_data = filtered_1000(strcmp(filtered_1000.DTW_Mode, mode_name), :);
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
        ndcg50_gt = mean(wm_data.('NDCG@50_GTvsEB'), 'omitnan');
        
        % Composite score
        composite = calculate_composite_score(ndcg50_gt, r50_dtw, mean_rank_gt, r50_gt, spearman_mean);

        
        % Store
        table2_data(end+1, :) = {
            mode_label, wm, n, ...
            spearman_mean, spearman_std, ...
            r10_dtw, r50_dtw, ...
            r10_gt, r50_gt, mean_rank_gt, ndcg50_gt, ...
            composite
        };
    end
end

% Create table
T2 = cell2table(table2_data, 'VariableNames', ...
    {'Mode', 'Weight_Mode', 'N', ...
     'Spearman_Mean', 'Spearman_Std', ...
     'R10_DTW', 'R50_DTW', ...
     'R10_GT', 'R50_GT', 'MeanRank_GT', 'NDCG50_GT', ...
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

% ========================================================================
%% COMPOSITE SCORE CONFIGURATION
% ========================================================================

function score = calculate_composite_score(ndcg, r50_dtw, mean_rank_gt, r50_gt, spearman)
    % Zentrale Composite Score Berechnung
    % Alle Gewichte hier ändern, um sie global anzupassen
    
    % Gewichtungen (müssen zu 1.0 summieren)
    w_ndcg = 0.25;
    w_r50_dtw = 0.25;
    w_rank = 0.0;       % wird als 1/mean_rank verwendet
    w_r50_gt = 0.25;
    w_spearman = 0.25;
    
    % Berechnung
    score = w_ndcg * ndcg + ...
            w_r50_dtw * r50_dtw + ...
            w_rank * (1/mean_rank_gt) + ...
            w_r50_gt * r50_gt + ...
            w_spearman * spearman;
end
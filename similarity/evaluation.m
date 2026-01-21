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

%% STEP 1.7: Smart Excel Export - Similarity Search Multi-Query
% ========================================================================

fprintf('\nSTEP 1.7: Smart Excel Export - Similarity Search...\n');
fprintf('────────────────────────────────────────\n\n');

% Load Similarity Search CSV files
similarity_files = dir('results/similarity_search_*.csv');

if isempty(similarity_files)
    fprintf('ℹ No Similarity Search CSV files found - skipping\n\n');
else
    fprintf('Found %d Similarity Search CSV file(s)\n', length(similarity_files));
    
    % Load and combine all CSVs
    similarity_data_all = [];
    
    for i = 1:length(similarity_files)
        filepath = fullfile(similarity_files(i).folder, similarity_files(i).name);
        temp_data = readtable(filepath, 'VariableNamingRule', 'preserve');
        
        % Force query_id to be cell array (string)
        if ismember('query_id', temp_data.Properties.VariableNames)
            if ~iscell(temp_data.query_id)
                temp_data.query_id = cellstr(string(temp_data.query_id));
            end
        end
        
        % Force segment_id to be cell array (string)
        if ismember('segment_id', temp_data.Properties.VariableNames)
            if ~iscell(temp_data.segment_id)
                temp_data.segment_id = cellstr(string(temp_data.segment_id));
            end
        end
        
        % Force level to be cell array (string)
        if ismember('level', temp_data.Properties.VariableNames)
            if ~iscell(temp_data.level)
                temp_data.level = cellstr(string(temp_data.level));
            end
        end
        
        % Extract timestamp from filename
        filename = similarity_files(i).name;
        timestamp_match = regexp(filename, 'similarity_search_(.+)\.csv', 'tokens');
        if ~isempty(timestamp_match)
            timestamp_str = timestamp_match{1}{1};
            temp_data.Timestamp = repmat({timestamp_str}, height(temp_data), 1);
        end
        
        if i == 1
            similarity_data_all = temp_data;
        else
            similarity_data_all = [similarity_data_all; temp_data];
        end
        
        fprintf('  Loaded: %s (%d rows)\n', similarity_files(i).name, height(temp_data));
    end
    
    fprintf('\nTotal Similarity Search data loaded: %d rows\n', height(similarity_data_all));
    
    % Sheet name
    sheet_name_similarity = 'All Data - Similarity Search';
    
    %% Check if sheet exists
    if exist(export_filename, 'file')
        try
            fprintf('\n  Reading existing data from sheet "%s"...\n', sheet_name_similarity);
            existing_data_sim = readtable(export_filename, 'Sheet', sheet_name_similarity, 'VariableNamingRule', 'preserve');
            fprintf('  ✓ Existing data: %d rows\n', height(existing_data_sim));
            
            % Check timestamps
            if ismember('Timestamp', existing_data_sim.Properties.VariableNames) && ...
               ismember('Timestamp', similarity_data_all.Properties.VariableNames)
                
                existing_timestamps_sim = unique(existing_data_sim.Timestamp);
                new_timestamps_sim = unique(similarity_data_all.Timestamp);
                
                fprintf('\n  Existing timestamps: %d unique\n', length(existing_timestamps_sim));
                fprintf('  New timestamps: %d unique\n', length(new_timestamps_sim));
                
                % Find new timestamps
                timestamps_to_add_sim = setdiff(new_timestamps_sim, existing_timestamps_sim);
                
                if isempty(timestamps_to_add_sim)
                    fprintf('\n  ℹ All data already exists - nothing to add!\n\n');
                else
                    fprintf('\n  → Found %d NEW timestamp(s) to add:\n', length(timestamps_to_add_sim));
                    for i = 1:length(timestamps_to_add_sim)
                        fprintf('      %d. %s\n', i, char(timestamps_to_add_sim(i)));
                    end
                    
                    % Filter new data
                    if iscellstr(timestamps_to_add_sim) || isstring(timestamps_to_add_sim)
                        mask_to_add_sim = ismember(similarity_data_all.Timestamp, timestamps_to_add_sim);
                    else
                        mask_to_add_sim = false(height(similarity_data_all), 1);
                        for i = 1:length(timestamps_to_add_sim)
                            mask_to_add_sim = mask_to_add_sim | strcmp(similarity_data_all.Timestamp, timestamps_to_add_sim(i));
                        end
                    end
                    
                    data_to_add_sim = similarity_data_all(mask_to_add_sim, :);
                    
                    fprintf('  → Rows to add: %d\n\n', height(data_to_add_sim));
                    
                    % Append
                    combined_data_sim = [existing_data_sim; data_to_add_sim];
                    
                    fprintf('  Combined data: %d rows (was %d, added %d)\n\n', ...
                        height(combined_data_sim), height(existing_data_sim), height(data_to_add_sim));
                    
                    % Write
                    fprintf('  Writing updated data to Excel...\n');
                    writetable(combined_data_sim, export_filename, 'Sheet', sheet_name_similarity, 'WriteRowNames', false);
                    fprintf('  ✓ Successfully appended data!\n\n');
                    
                    fprintf('╔════════════════════════════════════════════════════════════════╗\n');
                    fprintf('║  SIMILARITY SEARCH DATA APPENDED                               ║\n');
                    fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');
                    fprintf('Sheet: %s\n', sheet_name_similarity);
                    fprintf('  Previous rows: %d\n', height(existing_data_sim));
                    fprintf('  Added rows: %d\n', height(data_to_add_sim));
                    fprintf('  Total rows now: %d\n\n', height(combined_data_sim));
                end
                
            else
                fprintf('  ⚠ Timestamp column not found - cannot determine new data\n\n');
            end
            
        catch ME
            fprintf('  ℹ Sheet "%s" does not exist yet - creating new\n', sheet_name_similarity);
            
            % Create new sheet
            writetable(similarity_data_all, export_filename, 'Sheet', sheet_name_similarity, 'WriteRowNames', false);
            fprintf('  ✓ Created new sheet with %d rows\n\n', height(similarity_data_all));
            
            fprintf('╔════════════════════════════════════════════════════════════════╗\n');
            fprintf('║  SIMILARITY SEARCH SHEET CREATED                               ║\n');
            fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');
            fprintf('Sheet: %s\n', sheet_name_similarity);
            fprintf('  Total rows: %d\n\n', height(similarity_data_all));
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
use_normalized_dtw = false;
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
if use_normalized_dtw
    data_table = data_table(data_table.DTW_Normalization == 1, :);
end

data_table = data_table(contains(data_table.Timestamp, '2026'), :);

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


matlab
% ========================================================================
%% STEP 3: FIGURE 1 - DIMENSIONALITY (BASELINES ONLY!) - COMPOSITE SCORE
% ========================================================================

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  FIGURE 1: DIMENSIONALITY INVARIANCE (BASELINES)              ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

% Filter für Database_Size = 1000
filtered_baselines_1000 = filtered_baselines(filtered_baselines.Database_Size == 1000, :);
filtered_baselines_1000 = filtered_baselines_1000(contains(filtered_baselines_1000.Timestamp, '2026'), :);

fig1 = figure('Position', [100, 100, 750, 500], 'Color', 'w');
set(fig1, 'Name', 'Figure 1: Dimensionality Invariance (Baselines) - Composite');

% Farben
color_motion = [0.86, 0.13, 0.15];
color_shape  = [0.15, 0.39, 0.91];
color_zone   = [0.09, 0.64, 0.29];
color_motion_dark = color_motion * 0.6;
color_shape_dark = color_shape * 0.6;

h_lines = [];
legend_entries = {};
all_plot_data = struct([]);
plot_idx = 0;

levels = {'Trajectory', 'Segment'};

% ========================================================================
% Daten vorbereiten für beide Levels
% ========================================================================

for lvl = 1:length(levels)
    level_name = levels{lvl};
    filtered_level = filtered_baselines_1000(strcmp(filtered_baselines_1000.Level, level_name), :);
    
    for m = 1:length(dtw_modes)
        plot_idx = plot_idx + 1;
        curr_dtw_mode = dtw_modes{m};
        baseline = baseline_weights{m};
        
        mode_data = filtered_level(strcmp(filtered_level.DTW_Mode, curr_dtw_mode), :);
        
        all_dims_per_param = unique(mode_data.Total_Dims);
        all_dims_per_param = sort(all_dims_per_param);
        
        % ACTUAL dimensions
        if strcmp(curr_dtw_mode, 'joint_states')
            actual_dims = all_dims_per_param * 6;
            dim_label = '6×n';
        else
            actual_dims = all_dims_per_param * 3;
            dim_label = '3×n';
        end
        
        fprintf('%s - %s (%s):\n', level_name, mode_labels{m}, baseline);
        fprintf('  Dims per param: %s\n', mat2str(all_dims_per_param'));
        fprintf('  Actual dims:    %s (%s)\n', mat2str(actual_dims'), dim_label);
        
        composite_vals = zeros(length(all_dims_per_param), 1);
        
        for i = 1:length(all_dims_per_param)
            dim_data = mode_data(mode_data.Total_Dims == all_dims_per_param(i), :);
            
            spearman_mean = mean(dim_data.Spearman_DTWvsEB);
            r50_dtw = mean(dim_data.('R@50_DTWvsEB'));
            r50_gt = mean(dim_data.('R@50_GTvsEB'));
            mean_rank_gt = mean(dim_data.Mean_GTvsEB_Rank);
            ndcg50_gt = mean(dim_data.('NDCG@50_GTvsEB'), 'omitnan');
            
            composite_vals(i) = calculate_composite_score(ndcg50_gt, r50_dtw, mean_rank_gt, r50_gt, spearman_mean);
        end
        
        all_plot_data(plot_idx).actual_dims = actual_dims;
        all_plot_data(plot_idx).composite = composite_vals;
        all_plot_data(plot_idx).mode = curr_dtw_mode;
        all_plot_data(plot_idx).baseline = baseline;
        all_plot_data(plot_idx).level = level_name;
        
        fprintf('  Composite scores: %s\n\n', mat2str(composite_vals', 3));
    end
end

% ========================================================================
% Sweet Spot Configuration
% ========================================================================
sweet_start = 30;
sweet_end = 150;
plateau_start = 150;
all_dims_combined = [];
for p = 1:length(all_plot_data)
    all_dims_combined = [all_dims_combined; all_plot_data(p).actual_dims];
end
plateau_end = max(all_dims_combined);

fprintf('=== Manual Sweet Spot Configuration ===\n');
fprintf('Sweet Spot Zone: [%d, %d] dims\n', sweet_start, sweet_end);
fprintf('Plateau Zone: [%d, %d] dims\n\n', plateau_start, plateau_end);

% ========================================================================
% Plot
% ========================================================================
hold on;

% Y-Achsen Range
all_scores = [];
for p = 1:length(all_plot_data)
    all_scores = [all_scores; all_plot_data(p).composite];
end
y_min = min(all_scores) * 0.991;
y_max = max(all_scores) * 1.015;

% Sweet spot zones
fill([sweet_start, sweet_end, sweet_end, sweet_start], [y_min, y_min, y_max, y_max], ...
     color_zone, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
fill([plateau_start, plateau_end, plateau_end, plateau_start], [y_min, y_min, y_max, y_max], ...
     [0.5, 0.5, 0.5], 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');

% Plot alle Daten
for p = 1:length(all_plot_data)
    is_trajectory = strcmp(all_plot_data(p).level, 'Trajectory');
    is_motion = strcmp(all_plot_data(p).mode, 'joint_states');
    
    % Farbe und Stil bestimmen
    if is_motion
        if is_trajectory
            line_color = color_motion;
            line_style = '-';
            marker_style = 'o';
            line_width = 3;
            marker_size = 11;
            legend_label = 'Motion (Traj.)';
        else
            line_color = color_motion_dark;
            line_style = '--';
            marker_style = 'o';
            line_width = 2;
            marker_size = 8;
            legend_label = 'Motion (Seg.)';
        end
    else
        if is_trajectory
            line_color = color_shape;
            line_style = '-';
            marker_style = 's';
            line_width = 3;
            marker_size = 11;
            legend_label = 'Shape (Traj.)';
        else
            line_color = color_shape_dark;
            line_style = '--';
            marker_style = 's';
            line_width = 2;
            marker_size = 8;
            legend_label = 'Shape (Seg.)';
        end
    end
    
    h = plot(all_plot_data(p).actual_dims, all_plot_data(p).composite, ...
        [marker_style line_style], 'LineWidth', line_width, 'MarkerSize', marker_size, ...
        'Color', line_color, 'MarkerFaceColor', line_color);
    h_lines = [h_lines, h];
    legend_entries{end+1} = legend_label;
end

% Styling
set(gca, 'XScale', 'log');
xlabel('Total embedding dimensions', 'FontWeight', 'bold', 'FontSize', 20);
ylabel('Composite score', 'FontWeight', 'bold', 'FontSize', 20);
set(gca, 'XTick', [6, 15, 30, 75, 150, 300, 600, 1200, 3600]);
set(gca, 'FontName', 'Courier New', 'FontWeight', 'bold', 'FontSize', 18);
xtickangle(0);
ylim([y_min, y_max]);

legend(h_lines, legend_entries, 'Location', 'best', 'FontSize', 14);

ax = gca;
ax.Position = [0.14 0.14 0.84 0.84];
ax.LooseInset = [0, 0, 0, 0];
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
%% TABLE 1: DIMENSIONALITY ANALYSIS (BASELINES) - COMBINED
% ========================================================================
fprintf('Creating Table 1: Dimensionality Analysis (Combined)...\n');
table1_data = cell(0, 11);

for lvl = 1:length(levels)
    level_name = levels{lvl};
    filtered_level = filtered_baselines_1000(strcmp(filtered_baselines_1000.Level, level_name), :);
    
    for m = 1:length(dtw_modes)
        mode_name = dtw_modes{m};
        mode_label = mode_labels{m};
        baseline = baseline_weights{m};
        mode_data = filtered_level(strcmp(filtered_level.DTW_Mode, mode_name), :);
        
        all_dims_per_param = unique(mode_data.Total_Dims);
        all_dims_per_param = sort(all_dims_per_param);
        
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
            
            n = height(dim_data);
            spearman_mean = mean(dim_data.Spearman_DTWvsEB);
            r50_dtw = mean(dim_data.('R@50_DTWvsEB'));
            r50_gt = mean(dim_data.('R@50_GTvsEB'));
            mean_rank_gt = mean(dim_data.Mean_GTvsEB_Rank);
            ndcg50_gt = mean(dim_data.('NDCG@50_GTvsEB'), 'omitnan');
            
            composite = calculate_composite_score(ndcg50_gt, r50_dtw, mean_rank_gt, r50_gt, spearman_mean);
            
            table1_data(end+1, :) = {
                level_name, ...
                mode_label, baseline, dim, actual_dim, n, ...
                spearman_mean, ...
                r50_dtw, r50_gt, ndcg50_gt, ...
                composite
            };
        end
    end
end

T1 = cell2table(table1_data, 'VariableNames', ...
    {'Level', 'Mode', 'Config', 'Dims_Per_Param', 'Actual_Dims', 'N', ...
    'Spearman', ...
    'R@50 (DTW)', 'R@50 (GT)', 'NDCG@50 (GT)', ...
    'Composite_Score'});

% Sort: Level -> Mode -> Composite (descending)
T1 = sortrows(T1, {'Level', 'Mode', 'Composite_Score'}, {'ascend', 'ascend', 'descend'});

% Save
table1_file = fullfile(output_folder, 'table1_dimensionality.csv');
writetable(T1, table1_file);
fprintf('✓ Saved: %s\n\n', table1_file);

% ========================================================================
%% TABLE 1: DIMENSIONALITY ANALYSIS (BASELINES) - COMBINED
% ========================================================================
fprintf('Creating Table 1: Dimensionality Analysis (Combined)...\n');
table1_data = cell(0, 11);

for lvl = 1:length(levels)
    level_name = levels{lvl};
    filtered_level = filtered_baselines_1000(strcmp(filtered_baselines_1000.Level, level_name), :);
    
    for m = 1:length(dtw_modes)
        mode_name = dtw_modes{m};
        mode_label = mode_labels{m};
        baseline = baseline_weights{m};
        mode_data = filtered_level(strcmp(filtered_level.DTW_Mode, mode_name), :);
        
        all_dims_per_param = unique(mode_data.Total_Dims);
        all_dims_per_param = sort(all_dims_per_param);
        
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
            
            n = height(dim_data);
            spearman_mean = mean(dim_data.Spearman_DTWvsEB);
            r50_dtw = mean(dim_data.('R@50_DTWvsEB'));
            r50_gt = mean(dim_data.('R@50_GTvsEB'));
            mean_rank_gt = mean(dim_data.Mean_GTvsEB_Rank);
            ndcg50_gt = mean(dim_data.('NDCG@50_GTvsEB'), 'omitnan');
            
            composite = calculate_composite_score(ndcg50_gt, r50_dtw, mean_rank_gt, r50_gt, spearman_mean);
            
            table1_data(end+1, :) = {
                level_name, ...
                mode_label, baseline, dim, actual_dim, n, ...
                spearman_mean, ...
                r50_dtw, r50_gt, ndcg50_gt, ...
                composite
            };
        end
    end
end

T1 = cell2table(table1_data, 'VariableNames', ...
    {'Level', 'Mode', 'Config', 'Dims_Per_Param', 'Actual_Dims', 'N', ...
    'Spearman', ...
    'R@50 (DTW)', 'R@50 (GT)', 'NDCG@50 (GT)', ...
    'Composite_Score'});

% Sort: Level -> Mode -> Composite (descending)
T1 = sortrows(T1, {'Level', 'Mode', 'Composite_Score'}, {'ascend', 'ascend', 'descend'});

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
filtered_1000 = filtered_1000(contains(filtered_1000.Timestamp, '2026'), :);
filtered_1000 = filtered_1000(contains(filtered_1000.Embedding_Config, 'Multi-Balanced-25'), :);

% Farben
color_motion = [0.86, 0.13, 0.15];
color_shape  = [0.15, 0.39, 0.91];
color_motion_dark = color_motion * 0.6;
color_shape_dark = color_shape * 0.6;

fig2 = figure('Position', [100, 100, 750, 500], 'Color', 'w');
set(fig2, 'Name', 'Figure 2: Weight Mode Contribution - Composite');

hold on;

all_composite_vals = [];
h_lines = [];
legend_entries = {};

levels = {'Trajectory', 'Segment'};

for lvl = 1:length(levels)
    level_name = levels{lvl};
    filtered_level = filtered_1000(strcmp(filtered_1000.Level, level_name), :);
    
    for m = 1:length(dtw_modes)
        mode_name = dtw_modes{m};
        mode_label = mode_labels{m};
        mode_data = filtered_level(strcmp(filtered_level.DTW_Mode, mode_name), :);
        
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
            ndcg50_gt = mean(wm_data.('NDCG@50_GTvsEB'), 'omitnan');
            
            % Composite Score
            wm_results(end).composite = calculate_composite_score(ndcg50_gt, r50_dtw, mean_rank_gt, r50_gt, spearman_mean);
        end
        
        % Manual logical sorting for incremental features
        if strcmp(mode_name, 'joint_states')
            feature_order = {
                'Joint only',
                'Joint + Meta',
                'Joint + Velocity',
                'Joint + Orient',
                'Joint + Position',
                'Joint + All'
            };
        else
            feature_order = {
                'Position only',
                'Pos + Meta',
                'Pos + Velocity',
                'Pos + Orient',
                'Pos + Joint',
                'Pos + All'
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
        
        for w = 1:length(wm_results)
            if ~ismember(wm_results(w).name, feature_order)
                sorted_indices(end+1) = w;
            end
        end
        
        wm_results = wm_results(sorted_indices);
        
        x_pos = 1:length(wm_results);
        
        % Farbe und Stil bestimmen
        is_trajectory = strcmp(level_name, 'Trajectory');
        is_motion = strcmp(mode_name, 'joint_states');
        
        if is_motion
            if is_trajectory
                line_color = color_motion;
                line_style = '-';
                marker_style = 'o';
                line_width = 3;
                marker_size = 11;
                legend_label = 'Motion (Traj.)';
            else
                line_color = color_motion_dark;
                line_style = '--';
                marker_style = 'o';
                line_width = 2;
                marker_size = 8;
                legend_label = 'Motion (Seg.)';
            end
        else
            if is_trajectory
                line_color = color_shape;
                line_style = '-';
                marker_style = 's';
                line_width = 3;
                marker_size = 11;
                legend_label = 'Shape (Traj.)';
            else
                line_color = color_shape_dark;
                line_style = '--';
                marker_style = 's';
                line_width = 2;
                marker_size = 8;
                legend_label = 'Shape (Seg.)';
            end
        end
        
        % Plot Composite Score
        composite_vals = [wm_results.composite];
        all_composite_vals = [all_composite_vals; composite_vals'];
        
        h = plot(x_pos, composite_vals, ...
            [marker_style line_style], 'LineWidth', line_width, 'MarkerSize', marker_size, ...
            'Color', line_color, 'MarkerFaceColor', line_color);
        
        h_lines = [h_lines, h];
        legend_entries{end+1} = legend_label;
    end
end

% Auto-Skalierung
y_min = min(all_composite_vals) * 0.98;
y_max = max(all_composite_vals) * 1.02;

ylabel('Composite score', 'FontWeight', 'bold', 'FontSize', 20);
ylim([y_min, y_max]);

% X-axis
set(gca, 'XTick', 1:6);
set(gca, 'XTickLabel', {'Baseline', '+Meta.', '+Vel.', '+Orient.', '+Cross', '+All'});
xlabel('Incremental feature addition', 'FontWeight', 'bold', 'FontSize', 20);
xtickangle(0);
set(gca, 'FontName', 'courier new');
%yticks([0.60,0.64,0.68,0.72,0.76]);
%yticklabels({'0.60','0.64','0.68','0.72','0.76'});


grid on;
set(gca, 'FontSize', 18, 'FontWeight', 'bold');
ax = gca;
ax.Position = [0.14 0.14 0.82 0.84];
ax.LooseInset = [0, 0, 0, 0];

legend(h_lines, legend_entries, 'Location', 'southwest', 'FontSize', 14);

hold off;

% Save
if save_figures
    fig2_file = fullfile(figure_folder, 'incremental_feature');
    if strcmp(figure_format, 'pdf') || strcmp(figure_format, 'both')
        exportgraphics(gca, [fig2_file '.pdf']);
        fprintf('✓ Saved: %s.pdf\n', fig2_file);
    end
end

fprintf('\n');

% ========================================================================
%% TABLE 2: WEIGHT MODE ANALYSIS (ALL MODES) - COMBINED
% ========================================================================

fprintf('Creating Table 2: Weight Mode Analysis (Combined)...\n');

table2_data = cell(0, 13);

for lvl = 1:length(levels)
    level_name = levels{lvl};
    filtered_level = filtered_1000(strcmp(filtered_1000.Level, level_name), :);
    
    for m = 1:length(dtw_modes)
        mode_name = dtw_modes{m};
        mode_label = mode_labels{m};
        mode_data = filtered_level(strcmp(filtered_level.DTW_Mode, mode_name), :);
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
                level_name, ...
                mode_label, wm, n, ...
                spearman_mean, spearman_std, ...
                r10_dtw, r50_dtw, ...
                r10_gt, r50_gt, mean_rank_gt, ndcg50_gt, ...
                composite
            };
        end
    end
end

% Create table
T2 = cell2table(table2_data, 'VariableNames', ...
    {'Level', 'Mode', 'Weight_Mode', 'N', ...
     'Spearman_Mean', 'Spearman_Std', ...
     'R10_DTW', 'R50_DTW', ...
     'R10_GT', 'R50_GT', 'MeanRank_GT', 'NDCG50_GT', ...
     'Composite_Score'});

% Sort: Level -> Mode -> Composite (descending)
T2 = sortrows(T2, {'Level', 'Mode', 'Composite_Score'}, {'ascend', 'ascend', 'descend'});

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
    w_ndcg = 0.50;
    w_r50_dtw = 0.25;
    w_rank = 0.0;       % wird als 1/mean_rank verwendet
    w_r50_gt = 0.0;
    w_spearman = 0.25;
    
    % Berechnung
    score = w_ndcg * ndcg + ...
            w_r50_dtw * r50_dtw + ...
            w_rank * (1/mean_rank_gt) + ...
            w_r50_gt * r50_gt + ...
            w_spearman * spearman;
end
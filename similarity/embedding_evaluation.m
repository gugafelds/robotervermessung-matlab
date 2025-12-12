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
data = loadExperimentData('.');

%% STEP 1.5: Smart Excel Export - Append Only New Data
% ========================================================================

fprintf('\nSTEP 1.5: Smart Excel Export (Append New Data Only)...\n');
fprintf('────────────────────────────────────────\n\n');

% Create export folder if it doesn't exist
export_folder = 'exports';
if ~exist(export_folder, 'dir')
    mkdir(export_folder);
    fprintf('Created folder: %s\n', export_folder);
end

% Fixed filename
export_filename = fullfile(export_folder, 'experiment_data.xlsx');

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
                fprintf('Press Enter to continue to analysis...\n');
                pause;
                % Skip rest of export
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

fprintf('Press Enter to continue to analysis...\n');
pause;
%% STEP 2: Find Best Embedding Config + Weight Mode
% ========================================================================

fprintf('\nSTEP 2: Finding Best Embedding Config + Weight Mode...\n');
fprintf('────────────────────────────────────────────────────────────────\n\n');

% Focus on Trajectory level for main analysis
traj_data = data.combined;

fprintf('Analyzing %d trajectory-level experiments...\n\n', height(traj_data));

%% ========================================================================
%  DTW CONFIGURATION FILTER
%  ========================================================================

fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('DTW CONFIGURATION FILTER\n');
fprintf('═══════════════════════════════════════════════════════════════\n\n');

% Check what configurations are available
unique_dtw_norm = unique(traj_data.DTW_Normalization);
unique_dtw_rot = unique(traj_data.DTW_RotationAlign);
unique_lb_keogh = unique(traj_data.LB_Keogh_Candidates);
unique_db_size = unique(traj_data.Database_Size);

fprintf('Available configurations:\n');
fprintf('  DTW_Normalization:   %s\n', mat2str(unique_dtw_norm));
fprintf('  DTW_RotationAlign:   %s\n', mat2str(unique_dtw_rot));
fprintf('  LB_Keogh_Candidates: %s\n', mat2str(unique_lb_keogh));
fprintf('  Database_Size:       %s\n\n', mat2str(unique_db_size));

% ========================================================================
% USER FILTER: Set your desired configuration here
% ========================================================================

USE_FILTER = true;  % Set to false to analyze all data

if USE_FILTER
    % Filter settings
    FILTER_DTW_NORM = 0;        % 0 = no normalization, 1 = normalized
    FILTER_DTW_ROT = 0;         % 0 = no rotation, 1 = with rotation alignment
    FILTER_LB_KEOGH = 200;       % Leave empty [] to include all, or set to 100 or 200
    FILTER_DB_SIZE = [];        % Leave empty [] to include all, or set to specific size (e.g., 750, 1000)
    
    fprintf('FILTER ACTIVE:\n');
    fprintf('  DTW_Normalization = %d\n', FILTER_DTW_NORM);
    fprintf('  DTW_RotationAlign = %d\n', FILTER_DTW_ROT);
    
    if ~isempty(FILTER_LB_KEOGH)
        fprintf('  LB_Keogh_Candidates = %d\n', FILTER_LB_KEOGH);
    else
        fprintf('  LB_Keogh_Candidates = ALL\n');
    end
    
    if ~isempty(FILTER_DB_SIZE)
        fprintf('  Database_Size = %d\n', FILTER_DB_SIZE);
    else
        fprintf('  Database_Size = ALL\n');
    end
    
    % Apply filter
    dtw_filter_mask = (traj_data.DTW_Normalization == FILTER_DTW_NORM) & ...
                      (traj_data.DTW_RotationAlign == FILTER_DTW_ROT);
    
    % Add LB_Keogh filter if specified
    if ~isempty(FILTER_LB_KEOGH)
        dtw_filter_mask = dtw_filter_mask & (traj_data.LB_Keogh_Candidates == FILTER_LB_KEOGH);
    end
    
    % Add Database_Size filter if specified
    if ~isempty(FILTER_DB_SIZE)
        dtw_filter_mask = dtw_filter_mask & (traj_data.Database_Size == FILTER_DB_SIZE);
    end
    
    traj_data_filtered = traj_data(dtw_filter_mask, :);
    
    fprintf('\n  Filtered: %d → %d experiments (%.1f%%)\n\n', ...
        height(traj_data), height(traj_data_filtered), ...
        100 * height(traj_data_filtered) / height(traj_data));
    
    if height(traj_data_filtered) == 0
        error('No data matches the filter criteria! Adjust filter settings.');
    end
    
    % Use filtered data
    traj_data = traj_data_filtered;
    
else
    fprintf('FILTER DISABLED: Analyzing ALL configurations\n\n');
end

fprintf('═══════════════════════════════════════════════════════════════\n\n');

%% Aggregate Performance by Config + Weight Mode
% ========================================================================

configs = data.metadata.unique_embedding_configs;
weight_modes = data.metadata.unique_weight_modes;
dtw_modes = {'joint_states', 'position'};
levels = {'Trajectory', 'Segment'};

% ========================================================================
% USER CONFIG: Composite Score Weights
% ========================================================================

WEIGHT_GT_RETRIEVAL = 0.5;      % GT Coverage (R@50_GTvsEB)
WEIGHT_CONSISTENCY = 0.3;       % Ranking Consistency (|Rank_EB - Rank_DTW|)
WEIGHT_DTW_APPROX = 0.2;        % DTW Approximation (Spearman_DTWvsEB)

fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('COMPOSITE SCORE CONFIGURATION\n');
fprintf('═══════════════════════════════════════════════════════════════\n\n');
fprintf('Score = %.1f × R@50_GT + %.1f × Consistency + %.1f × Spearman_DTW\n', ...
    WEIGHT_GT_RETRIEVAL, WEIGHT_CONSISTENCY, WEIGHT_DTW_APPROX);
fprintf('where Consistency = 1 / (1 + |Rank_EB - Rank_DTW|)\n\n');

fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('BEST CONFIGURATIONS - MULTI-OBJECTIVE RANKING\n');
fprintf('═══════════════════════════════════════════════════════════════\n\n');

% Use the filtered data (contains BOTH Trajectory and Segment levels!)
all_filtered_data = traj_data;  % ← ÄNDERUNG: Umbenennen für Klarheit

% For each level (Trajectory / Segment)
for level_idx = 1:length(levels)
    level = levels{level_idx};
    
    fprintf('\n╔═══════════════════════════════════════════════════════════════╗\n');
    fprintf('║  %s LEVEL                                                  \n', upper(level));
    fprintf('╚═══════════════════════════════════════════════════════════════╝\n\n');
    
    % Filter for this level
    level_mask = strcmp(all_filtered_data.Level, level);  % ← ÄNDERUNG
    level_data = all_filtered_data(level_mask, :);        % ← ÄNDERUNG
    
    if isempty(level_data)
        fprintf('  No data for this level.\n\n');
        continue;
    end
    
    fprintf('Analyzing %d experiments at %s level...\n\n', height(level_data), level);
    
    % For each DTW mode
    for dtw_idx = 1:length(dtw_modes)
        dtw_mode = dtw_modes{dtw_idx};
        
        fprintf('───────────────────────────────────────────────────────────────\n');
        if strcmp(dtw_mode, 'joint_states')
            fprintf('JOINT STATES SPACE\n');
        else
            fprintf('CARTESIAN SPACE\n');
        end
        fprintf('───────────────────────────────────────────────────────────────\n\n');
        
        % Filter for this DTW mode
        mode_mask = strcmp(level_data.DTW_Mode, dtw_mode);
        mode_data = level_data(mode_mask, :);
        
        if isempty(mode_data)
            fprintf('  No data for this DTW mode.\n\n');
            continue;
        end
        
        % Create results array
        results = [];
        
        for c = 1:length(configs)
            for w = 1:length(weight_modes)
                config = configs{c};
                weight_mode = weight_modes{w};
                
                % Filter for this config + weight combination
                combo_mask = strcmp(mode_data.Embedding_Config, config) & ...
                             strcmp(mode_data.Weight_Mode, weight_mode);
                combo_data = mode_data(combo_mask, :);
                
                if ~isempty(combo_data)
                    % Metric 1: DTW Approximation
                    spearman_dtw = nanmean(combo_data.Spearman_DTWvsEB);
                    
                    % Metric 2: GT Retrieval Quality
                    r50_gt = nanmean(combo_data.("R@50_GTvsEB"));
                    r10_gt = nanmean(combo_data.("R@10_GTvsEB"));
                    mean_rank_eb = nanmean(combo_data.Mean_GTvsEB_Rank);
                    mean_rank_dtw = nanmean(combo_data.Mean_GTvsDTW_Rank);
                    
                    % Metric 3: Ranking Consistency
                    rank_diff = abs(mean_rank_eb - mean_rank_dtw);
                    consistency_score = 1 / (1 + rank_diff);
                    
                    % Composite Score (user-configurable weights)
                    composite_score = WEIGHT_GT_RETRIEVAL * r50_gt + ...
                                      WEIGHT_CONSISTENCY * consistency_score + ...
                                      WEIGHT_DTW_APPROX * spearman_dtw;
                    
                    % Store result
                    results(end+1).config = config;
                    results(end).weight_mode = weight_mode;
                    results(end).spearman_dtw = spearman_dtw;
                    results(end).r50_gt = r50_gt;
                    results(end).r10_gt = r10_gt;
                    results(end).mean_rank_eb = mean_rank_eb;
                    results(end).mean_rank_dtw = mean_rank_dtw;
                    results(end).rank_diff = rank_diff;
                    results(end).consistency_score = consistency_score;
                    results(end).composite_score = composite_score;
                    results(end).n = height(combo_data);
                end
            end
        end
        
        if isempty(results)
            fprintf('  No results found.\n\n');
            continue;
        end
        
        % Sort by composite score (descending)
        [~, sort_idx] = sort([results.composite_score], 'descend');
        results = results(sort_idx);
        
        % Print top 10
        fprintf('Rank | Config       | Weight Mode          | Spear | R@50  | Rank | Rank  | Cons | Score | n\n');
        fprintf('     |              |                      | (DTW) | (GT)  | (EB) | Diff  |      |       |  \n');
        fprintf('-----|--------------|----------------------|-------|-------|------|-------|------|-------|----\n');
        
        for i = 1:min(10, length(results))
            r = results(i);
            fprintf(' %2d  | %-12s | %-20s | %.3f | %.3f | %4.1f | %4.1f | %.2f | %.3f | %d\n', ...
                i, r.config, r.weight_mode, r.spearman_dtw, r.r50_gt, ...
                r.mean_rank_eb, r.rank_diff, r.consistency_score, r.composite_score, r.n);
        end
        
        fprintf('\n');
        
        % Highlight best
        best = results(1);
        fprintf('★ BEST OVERALL for %s - %s:\n', level, upper(dtw_mode));
        fprintf('  Config: %s\n', best.config);
        fprintf('  Weight Mode: %s\n', best.weight_mode);
        fprintf('  Composite Score: %.3f\n', best.composite_score);
        fprintf('  ├─ Spearman (DTW): %.3f (weight: %.1f)\n', best.spearman_dtw, WEIGHT_DTW_APPROX);
        fprintf('  ├─ R@50 (GT): %.3f (weight: %.1f)\n', best.r50_gt, WEIGHT_GT_RETRIEVAL);
        fprintf('  ├─ R@10 (GT): %.3f\n', best.r10_gt);
        fprintf('  ├─ Mean Rank (EB): %.1f\n', best.mean_rank_eb);
        fprintf('  ├─ Mean Rank (DTW): %.1f\n', best.mean_rank_dtw);
        fprintf('  ├─ Rank Difference: %.1f\n', best.rank_diff);
        fprintf('  └─ Consistency Score: %.2f (weight: %.1f)\n', best.consistency_score, WEIGHT_CONSISTENCY);
        fprintf('  Based on %d experiments\n\n', best.n);
        
        % Show Pareto-optimal configs
        fprintf('PARETO-OPTIMAL CONFIGS:\n\n');
        
        % Best DTW approximation
        [~, best_dtw_idx] = max([results.spearman_dtw]);
        if best_dtw_idx ~= 1
            fprintf('  • Best DTW Approximation: %s + %s (Spearman=%.3f)\n', ...
                results(best_dtw_idx).config, results(best_dtw_idx).weight_mode, ...
                results(best_dtw_idx).spearman_dtw);
        end
        
        % Best GT retrieval
        [~, best_gt_idx] = max([results.r50_gt]);
        if best_gt_idx ~= 1
            fprintf('  • Best GT Retrieval: %s + %s (R@50=%.3f)\n', ...
                results(best_gt_idx).config, results(best_gt_idx).weight_mode, ...
                results(best_gt_idx).r50_gt);
        end
        
        % Best ranking consistency
        [~, best_cons_idx] = min([results.rank_diff]);
        if best_cons_idx ~= 1
            fprintf('  • Best Ranking Consistency: %s + %s (Diff=%.1f)\n', ...
                results(best_cons_idx).config, results(best_cons_idx).weight_mode, ...
                results(best_cons_idx).rank_diff);
        end
        
        fprintf('\n');
    end
end

fprintf('\n═══════════════════════════════════════════════════════════════\n');
fprintf('ANALYSIS COMPLETE\n');
fprintf('═══════════════════════════════════════════════════════════════\n\n');

fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('SCORING FORMULA\n');
fprintf('═══════════════════════════════════════════════════════════════\n\n');
fprintf('Composite Score = %.1f × R@50_GT + %.1f × Consistency + %.1f × Spearman_DTW\n', ...
    WEIGHT_GT_RETRIEVAL, WEIGHT_CONSISTENCY, WEIGHT_DTW_APPROX);
fprintf('where Consistency = 1 / (1 + |Rank_EB - Rank_DTW|)\n\n');
fprintf('Rationale:\n');
fprintf('  • GT Retrieval (%.0f%%): Most important - find the right trajectories\n', WEIGHT_GT_RETRIEVAL * 100);
fprintf('  • Ranking Consistency (%.0f%%): Embeddings rank similar to DTW\n', WEIGHT_CONSISTENCY * 100);
fprintf('  • DTW Approximation (%.0f%%): Nice-to-have bonus\n\n', WEIGHT_DTW_APPROX * 100);
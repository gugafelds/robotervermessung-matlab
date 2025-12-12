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
%% STEP 2: Incremental Decision Making - Finding Best Configuration
% ========================================================================
% Goal: Show step-by-step how we arrive at the best embedding configuration
% through systematic analysis and justification of each decision
% ========================================================================

fprintf('\n╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  STEP 2: CONFIGURATION SELECTION ANALYSIS                      ║\n');
fprintf('║  Incremental decision-making with visual evidence              ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

% Use all data (no filtering yet - we'll decide what to filter based on analysis)
analysis_data = data.combined;

fprintf('Total experiments available: %d\n', height(analysis_data));
fprintf('  - Trajectory level: %d\n', sum(strcmp(analysis_data.Level, 'Trajectory')));
fprintf('  - Segment level: %d\n\n', sum(strcmp(analysis_data.Level, 'Segment')));

% ========================================================================
%% DECISION 1: Which Level? (Trajectory vs Segment)
% ========================================================================

fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('DECISION 1: Trajectory-Level vs. Segment-Level Analysis\n');
fprintf('═══════════════════════════════════════════════════════════════\n\n');

% Compare performance at both levels
traj_mask = strcmp(analysis_data.Level, 'Trajectory');
seg_mask = strcmp(analysis_data.Level, 'Segment');

fprintf('Performance Comparison:\n\n');
fprintf('                          | Trajectory | Segment  | Difference\n');
fprintf('--------------------------|------------|----------|------------\n');

% Key metrics
metrics_to_compare = {
    'Mean_GTvsEB_Rank',  'GT Mean Rank (EB)';
    'Mean_GTvsDTW_Rank', 'GT Mean Rank (DTW)';
    'R@10_GTvsEB',       'R@10 (GT vs EB)';
    'R@10_GTvsDTW',      'R@10 (GT vs DTW)';
    'Spearman_DTWvsEB',  'Spearman (DTW vs EB)';
};

for i = 1:size(metrics_to_compare, 1)
    metric = metrics_to_compare{i, 1};
    label = metrics_to_compare{i, 2};
    
    traj_val = nanmean(analysis_data.(metric)(traj_mask));
    seg_val = nanmean(analysis_data.(metric)(seg_mask));
    diff = ((traj_val - seg_val) / seg_val) * 100;
    
    fprintf('%-25s | %10.3f | %8.3f | %+7.1f%%\n', label, traj_val, seg_val, diff);
end

fprintf('\n📊 DECISION: Focus on TRAJECTORY level\n');
fprintf('Rationale: [Add your reasoning based on results]\n\n');

% Filter to trajectory level for remaining analysis
analysis_data = analysis_data(traj_mask, :);

% ========================================================================
%% DECISION 2: DTW Configuration (Normalized vs Non-Normalized)
% ========================================================================

fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('DECISION 2: DTW Normalization\n');
fprintf('═══════════════════════════════════════════════════════════════\n\n');

norm0_mask = analysis_data.DTW_Normalization == 0;
norm1_mask = analysis_data.DTW_Normalization == 1;

fprintf('Performance Comparison:\n\n');
fprintf('                          | Non-Norm   | Normalized | Difference\n');
fprintf('--------------------------|------------|------------|------------\n');

for i = 1:size(metrics_to_compare, 1)
    metric = metrics_to_compare{i, 1};
    label = metrics_to_compare{i, 2};
    
    norm0_val = nanmean(analysis_data.(metric)(norm0_mask));
    norm1_val = nanmean(analysis_data.(metric)(norm1_mask));
    diff = ((norm0_val - norm1_val) / norm1_val) * 100;
    
    fprintf('%-25s | %10.3f | %10.3f | %+7.1f%%\n', label, norm0_val, norm1_val, diff);
end

fprintf('\n📊 DECISION: Use NON-NORMALIZED DTW (normalize=0)\n');
fprintf('Rationale: Trajectory length is a discriminative feature\n\n');

% Filter to non-normalized
analysis_data = analysis_data(norm0_mask, :);

% ========================================================================
%% DECISION 3: DTW Mode (Joint States vs Position)
% ========================================================================

fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('DECISION 3: DTW Space (Joint States vs. Cartesian Position)\n');
fprintf('═══════════════════════════════════════════════════════════════\n\n');

joint_mask = strcmp(analysis_data.DTW_Mode, 'joint_states');
pos_mask = strcmp(analysis_data.DTW_Mode, 'position');

fprintf('Performance Comparison:\n\n');
fprintf('                          | Joint States | Position   | Difference\n');
fprintf('--------------------------|--------------|------------|------------\n');

for i = 1:size(metrics_to_compare, 1)
    metric = metrics_to_compare{i, 1};
    label = metrics_to_compare{i, 2};
    
    joint_val = nanmean(analysis_data.(metric)(joint_mask));
    pos_val = nanmean(analysis_data.(metric)(pos_mask));
    diff = ((joint_val - pos_val) / pos_val) * 100;
    
    fprintf('%-25s | %12.3f | %10.3f | %+7.1f%%\n', label, joint_val, pos_val, diff);
end

fprintf('\n📊 DECISION: Use JOINT STATES space\n');
fprintf('Rationale: [Add your reasoning based on results]\n\n');

% Filter to joint states
analysis_data = analysis_data(joint_mask, :);

% ========================================================================
%% ANALYSIS 4: Embedding Architecture Comparison
% ========================================================================

fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('ANALYSIS 4: Embedding Architecture Impact\n');
fprintf('═══════════════════════════════════════════════════════════════\n\n');

configs = data.metadata.unique_embedding_configs;

fprintf('Performance by Embedding Config:\n\n');
fprintf('Config            | Dims | Mean Rank | R@10  | Spearman | n\n');
fprintf('                  |      | (GTvsEB)  |(GTvsEB)| (DTWvsEB)|  \n');
fprintf('------------------|------|-----------|--------|----------|----\n');

config_results = [];
for c = 1:length(configs)
    config = configs{c};
    config_mask = strcmp(analysis_data.Embedding_Config, config);
    
    if any(config_mask)
        mean_rank = nanmean(analysis_data.Mean_GTvsEB_Rank(config_mask));
        r10 = nanmean(analysis_data.("R@10_GTvsEB")(config_mask));
        spearman = nanmean(analysis_data.Spearman_DTWvsEB(config_mask));
        dims = analysis_data.Total_Dims(find(config_mask, 1));
        n = sum(config_mask);
        
        fprintf('%-17s | %4d | %9.2f | %6.3f | %8.3f | %d\n', ...
            config, dims, mean_rank, r10, spearman, n);
        
        config_results(end+1).name = config;
        config_results(end).dims = dims;
        config_results(end).mean_rank = mean_rank;
        config_results(end).r10 = r10;
        config_results(end).spearman = spearman;
    end
end

fprintf('\n📊 OBSERVATION: All configs perform nearly identically!\n');
fprintf('→ Single-Fine-75 (75 dims) achieves same performance as Multi-Dense-200 (200 dims)\n');
fprintf('→ Dimensionality and architecture have minimal impact\n\n');

fprintf('📊 DECISION: Use Single-Fine-75 (simplest, lowest compute)\n\n');

% Filter to best config for next analysis
best_config = 'Single-Fine-75';
analysis_data = analysis_data(strcmp(analysis_data.Embedding_Config, best_config), :);

% ========================================================================
%% ANALYSIS 5: Weight Mode Impact
% ========================================================================

fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('ANALYSIS 5: Weight Mode Combinations\n');
fprintf('═══════════════════════════════════════════════════════════════\n\n');

weight_modes = data.metadata.unique_weight_modes;

fprintf('Performance by Weight Mode:\n\n');
fprintf('Weight Mode              | Mean Rank | R@10  | R@1   | Spearman | n\n');
fprintf('                         | (GTvsEB)  |(GTvsEB)|(GTvsEB)| (DTWvsEB)|  \n');
fprintf('-------------------------|-----------|-------|-------|----------|----\n');

weight_results = [];
for w = 1:length(weight_modes)
    weight_mode = weight_modes{w};
    weight_mask = strcmp(analysis_data.Weight_Mode, weight_mode);
    
    if any(weight_mask)
        mean_rank = nanmean(analysis_data.Mean_GTvsEB_Rank(weight_mask));
        r10 = nanmean(analysis_data.("R@10_GTvsEB")(weight_mask));
        r1 = nanmean(analysis_data.("R@1_GTvsEB")(weight_mask));
        spearman = nanmean(analysis_data.Spearman_DTWvsEB(weight_mask));
        n = sum(weight_mask);
        
        fprintf('%-24s | %9.2f | %5.3f | %5.3f | %8.3f | %d\n', ...
            weight_mode, mean_rank, r10, r1, spearman, n);
        
        weight_results(end+1).name = weight_mode;
        weight_results(end).mean_rank = mean_rank;
        weight_results(end).r10 = r10;
        weight_results(end).r1 = r1;
        weight_results(end).spearman = spearman;
    end
end

% Sort by mean rank
[~, sort_idx] = sort([weight_results.mean_rank]);
best_weight = weight_results(sort_idx(1));

fprintf('\n📊 BEST Weight Mode: %s\n', best_weight.name);
fprintf('  Mean Rank: %.2f\n', best_weight.mean_rank);
fprintf('  R@10: %.3f\n', best_weight.r10);
fprintf('  R@1: %.3f\n', best_weight.r1);
fprintf('  Spearman: %.3f\n\n', best_weight.spearman);

% ========================================================================
%% FINAL CONFIGURATION SUMMARY
% ========================================================================

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  FINAL CONFIGURATION                                           ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

fprintf('Selected Configuration:\n');
fprintf('  ✓ Level: Trajectory\n');
fprintf('  ✓ DTW Normalization: False (Non-normalized)\n');
fprintf('  ✓ DTW Mode: joint_states\n');
fprintf('  ✓ Embedding Config: %s (75 dims)\n', best_config);
fprintf('  ✓ Weight Mode: %s\n\n', best_weight.name);

fprintf('Performance (Average across %d queries):\n', length(data.metadata.unique_queries));
fprintf('  GT Mean Rank: %.2f (DTW: %.2f)\n', ...
    nanmean(analysis_data.Mean_GTvsEB_Rank), ...
    nanmean(analysis_data.Mean_GTvsDTW_Rank));
fprintf('  R@10: %.3f (DTW: %.3f)\n', ...
    nanmean(analysis_data.("R@10_GTvsEB")), ...
    nanmean(analysis_data.("R@10_GTvsDTW")));
fprintf('  R@1: %.3f (DTW: %.3f)\n', ...
    nanmean(analysis_data.("R@1_GTvsEB")), ...
    nanmean(analysis_data.("R@1_GTvsDTW")));
fprintf('  Spearman (DTW vs EB): %.3f\n\n', ...
    nanmean(analysis_data.Spearman_DTWvsEB));

fprintf(['═════════════════════════════════' ...
    '══════════════════════════════\n']);
fprintf('DECISION-MAKING COMPLETE\n');
fprintf('═══════════════════════════════════════════════════════════════\n\n');

%% ========================================================================
%  PLOT 5: Weight Mode Impact (BOTH DTW Spaces)
%  ========================================================================

fprintf('Creating Plot 5: Weight Mode Performance (Both DTW Spaces)...\n');

figure('Position', [100, 100, 1400, 600]);

% Build base filtered data (Trajectory + Non-Norm only)
filtered_data = data.combined;
filtered_data = filtered_data(strcmp(filtered_data.Level, 'Trajectory'), :);
filtered_data = filtered_data(filtered_data.DTW_Normalization == 0, :);

% Split by DTW Mode
dtw_modes = {'joint_states', 'position'};

for mode_idx = 1:2
    dtw_mode = dtw_modes{mode_idx};
    
    subplot(1, 2, mode_idx);
    
    % Filter for this DTW mode
    mode_data = filtered_data(strcmp(filtered_data.DTW_Mode, dtw_mode), :);
    
    % Use best config
    best_config_data = mode_data(strcmp(mode_data.Embedding_Config, 'Single-Fine-75'), :);
    
    if isempty(best_config_data)
        best_config_data = mode_data;
    end
    
    weight_modes = data.metadata.unique_weight_modes;
    
    weight_results = [];
    for w = 1:length(weight_modes)
        weight_mode = weight_modes{w};
        weight_mask = strcmp(best_config_data.Weight_Mode, weight_mode);
        
        if any(weight_mask)
            weight_results(end+1).name = weight_mode;
            weight_results(end).spearman = nanmean(best_config_data.Spearman_DTWvsEB(weight_mask));
            weight_results(end).r50 = nanmean(best_config_data.("R@50_DTWvsEB")(weight_mask));
            weight_results(end).r10 = nanmean(best_config_data.("R@10_DTWvsEB")(weight_mask));
        end
    end
    
    if isempty(weight_results)
        continue;
    end
    
    % Sort by Spearman
    [~, sort_idx] = sort([weight_results.spearman], 'descend');
    weight_results = weight_results(sort_idx);
    
    n_weights = length(weight_results);
    names = {weight_results.name};
    spearman_vals = [weight_results.spearman];
    r10_vals = [weight_results.r10];
    r50_vals = [weight_results.r50];
    
    x = 1:n_weights;
    
    % Dual axis plot
    yyaxis left
    bar(x, spearman_vals, 'FaceColor', [0.2, 0.5, 0.8], 'FaceAlpha', 0.8);
    ylabel('Spearman Correlation', 'FontWeight', 'bold');
    ylim([0, max(spearman_vals)*1.2]);
    
    for i = 1:n_weights
        text(x(i), spearman_vals(i)+0.01, sprintf('%.2f', spearman_vals(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 8);
    end
    
    yyaxis right
    hold on;
    plot(x, r10_vals, 'o-', 'LineWidth', 2, 'MarkerSize', 8, ...
        'Color', [0.8, 0.3, 0.2], 'MarkerFaceColor', [0.8, 0.3, 0.2]);
    plot(x, r50_vals, 's--', 'LineWidth', 2, 'MarkerSize', 7, ...
        'Color', [0.2, 0.7, 0.3], 'MarkerFaceColor', [0.2, 0.7, 0.3]);
    hold off;
    ylabel('Recall@K', 'FontWeight', 'bold');
    ylim([0, 1]);
    
    set(gca, 'XTick', x, 'XTickLabel', names, 'XTickLabelRotation', 45);
    
    if strcmp(dtw_mode, 'joint_states')
        title('Joint States Space', 'FontSize', 12, 'FontWeight', 'bold');
    else
        title('Cartesian Space', 'FontSize', 12, 'FontWeight', 'bold');
    end
    
    legend({'Spearman', 'R@10', 'R@50'}, 'Location', 'northwest');
    grid on;
end

sgtitle('Weight Mode Impact on DTW Approximation (Both Spaces)', ...
    'FontSize', 14, 'FontWeight', 'bold');


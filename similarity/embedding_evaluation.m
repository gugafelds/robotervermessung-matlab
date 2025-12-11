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

%% STEP 1.5: Export Filtered Data to Excel
% ========================================================================

fprintf('\nSTEP 1.5: Exporting data to Excel...\n');
fprintf('────────────────────────────────────────\n\n');

% Create export folder if it doesn't exist
export_folder = 'exports';
if ~exist(export_folder, 'dir')
    mkdir(export_folder);
    fprintf('Created folder: %s\n', export_folder);
end

% Generate filename with timestamp
timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd_HHmmss'));
export_filename = fullfile(export_folder, sprintf('experiment_data_%s.xlsx', timestamp));

% Get the combined data
export_data = data.combined;

% Convert numeric columns to use comma as decimal separator for European format
% Note: Excel will handle this based on system locale, but we prepare the data

fprintf('Preparing data for export...\n');
fprintf('  Total rows: %d\n', height(export_data));
fprintf('  Total columns: %d\n\n', width(export_data));

% Write to Excel with multiple sheets
fprintf('Writing to Excel file: %s\n\n', export_filename);

try
    % Sheet 1: All data
    writetable(export_data, export_filename, 'Sheet', 'All Data', 'WriteRowNames', false);
    fprintf('  ✓ Sheet 1: All Data (%d rows)\n', height(export_data));
    
    % Sheet 2: Trajectory level only
    traj_data_export = export_data(strcmp(export_data.Level, 'Trajectory'), :);
    writetable(traj_data_export, export_filename, 'Sheet', 'Trajectory Level', 'WriteRowNames', false);
    fprintf('  ✓ Sheet 2: Trajectory Level (%d rows)\n', height(traj_data_export));
    
    % Sheet 3: Segment level only
    seg_data_export = export_data(strcmp(export_data.Level, 'Segment'), :);
    writetable(seg_data_export, export_filename, 'Sheet', 'Segment Level', 'WriteRowNames', false);
    fprintf('  ✓ Sheet 3: Segment Level (%d rows)\n', height(seg_data_export));
    
    % Sheet 4: Summary statistics
    summary_table = table();
    summary_table.Metric = {
        'Total Experiments';
        'Queries';
        'Embedding Configs';
        'Weight Modes';
        'DTW Modes';
        'Trajectory Rows';
        'Segment Rows';
        'Mean Spearman (Traj)';
        'Mean Spearman (Seg)';
        'Mean R@50 (Traj)';
        'Mean R@50 (Seg)'
    };
    
    summary_table.Value = {
        height(export_data);
        data.metadata.num_queries;
        data.metadata.num_embedding_configs;
        data.metadata.num_weight_modes;
        data.metadata.num_dtw_modes;
        height(traj_data_export);
        height(seg_data_export);
        nanmean(traj_data_export.Spearman_DTWvsEB);
        nanmean(seg_data_export.Spearman_DTWvsEB);
        nanmean(traj_data_export.("R@50_DTWvsEB"));
        nanmean(seg_data_export.("R@50_DTWvsEB"))
    };
    
    writetable(summary_table, export_filename, 'Sheet', 'Summary', 'WriteRowNames', false);
    fprintf('  ✓ Sheet 4: Summary Statistics\n\n');
    
    fprintf('╔════════════════════════════════════════════════════════════════╗\n');
    fprintf('║  EXPORT SUCCESSFUL                                             ║\n');
    fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');
    fprintf('Excel file saved: %s\n', export_filename);
    fprintf('File size: %.2f MB\n\n', dir(export_filename).bytes / 1e6);
    
    fprintf('Excel file contains:\n');
    fprintf('  Sheet 1: All Data - Complete dataset\n');
    fprintf('  Sheet 2: Trajectory Level - Filtered view\n');
    fprintf('  Sheet 3: Segment Level - Filtered view\n');
    fprintf('  Sheet 4: Summary - Key statistics\n\n');
    
    fprintf('Note: Open in Excel and set your locale to use comma as decimal separator.\n');
    fprintf('      (File → Options → Advanced → Use system separators)\n\n');
    
catch ME
    fprintf('✗ ERROR: Could not write Excel file\n');
    fprintf('  Error message: %s\n\n', ME.message);
    
    % Fallback: Try CSV export with European format
    fprintf('Attempting CSV export with European format (comma as decimal)...\n');
    
    csv_filename = fullfile(export_folder, sprintf('experiment_data_%s.csv', timestamp));
    
    % Convert to European CSV format (semicolon separator, comma decimal)
    % This is tricky in MATLAB, so we use standard writetable
    writetable(export_data, csv_filename, 'Delimiter', ';');
    
    fprintf('  ✓ CSV exported: %s\n', csv_filename);
    fprintf('  Note: Use semicolon as separator when importing to Excel\n\n');
end

fprintf('Press Enter to continue to analysis...\n');
pause;


%% STEP 2: Find Best Embedding Config + Weight Mode
% ========================================================================

fprintf('\nSTEP 2: Finding Best Embedding Config + Weight Mode...\n');
fprintf('────────────────────────────────────────────────────────────────\n\n');

% Focus on Trajectory level for main analysis
traj_data = data.combined(strcmp(data.combined.Level, 'Trajectory'), :);

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
    FILTER_DTW_NORM = 1;        % 0 = no normalization, 1 = normalized
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

% Main metric for ranking
metric_name = 'R@50_GTvsDTW';

fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('BEST CONFIGURATIONS (Sorted by %s)\n', metric_name);
fprintf('═══════════════════════════════════════════════════════════════\n\n');

% For each DTW mode, find best combinations
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
    mode_mask = strcmp(traj_data.DTW_Mode, dtw_mode);
    mode_data = traj_data(mode_mask, :);
    
    if isempty(mode_data)
        fprintf('  No data for this DTW mode after filtering.\n\n');
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
                % Calculate mean over queries
                mean_r50 = nanmean(combo_data.(metric_name));
                std_r50 = nanstd(combo_data.(metric_name));
                mean_r10 = nanmean(combo_data.("R@10_DTWvsEB"));
                mean_spearman = nanmean(combo_data.Spearman_DTWvsEB);
                n_experiments = height(combo_data);
                
                % Store result
                results(end+1).config = config;
                results(end).weight_mode = weight_mode;
                results(end).mean_r50 = mean_r50;
                results(end).std_r50 = std_r50;
                results(end).mean_r10 = mean_r10;
                results(end).mean_spearman = mean_spearman;
                results(end).n = n_experiments;
            end
        end
    end
    
    if isempty(results)
        fprintf('  No results found.\n\n');
        continue;
    end
    
    % Sort by R@50 (descending)
    [~, sort_idx] = sort([results.mean_r50], 'descend');
    results = results(sort_idx);
    
    % Print top 10
    fprintf('Rank | Config            | Weight Mode              | R@50   | R@10   | Spearman | n\n');
    fprintf('-----|-------------------|--------------------------|--------|--------|----------|----\n');
    
    for i = 1:min(10, length(results))
        r = results(i);
        fprintf(' %2d  | %-17s | %-24s | %.3f  | %.3f  |  %.3f    | %d\n', ...
            i, r.config, r.weight_mode, r.mean_r50, r.mean_r10, r.mean_spearman, r.n);
    end
    
    fprintf('\n');
    
    % Highlight best
    best = results(1);
    fprintf('★ BEST for %s:\n', upper(dtw_mode));
    fprintf('  Config: %s\n', best.config);
    fprintf('  Weight Mode: %s\n', best.weight_mode);
    fprintf('  R@50: %.3f (±%.3f)\n', best.mean_r50, best.std_r50);
    fprintf('  R@10: %.3f\n', best.mean_r10);
    fprintf('  Spearman: %.3f\n', best.mean_spearman);
    fprintf('  Based on %d experiments\n\n', best.n);
end

%% Compare Embedding Configs (aggregated over weight modes)
% ========================================================================

fprintf('\n═══════════════════════════════════════════════════════════════\n');
fprintf('EMBEDDING CONFIGS COMPARISON (Aggregated over all Weight Modes)\n');
fprintf('═══════════════════════════════════════════════════════════════\n\n');

fprintf('Config            | Dims |     Joint States      |    Cartesian Space    \n');
fprintf('                  |      |  R@50  | R@10 | Spear |  R@50  | R@10 | Spear\n');
fprintf('------------------|------|-----------------------|----------------------\n');

for c = 1:length(configs)
    config = configs{c};
    
    % Get dims
    config_mask = strcmp(traj_data.Embedding_Config, config);
    if any(config_mask)
        dims = traj_data.Total_Dims(find(config_mask, 1));
    else
        dims = 0;
    end
    
    % Joint states
    joint_mask = strcmp(traj_data.Embedding_Config, config) & ...
                 strcmp(traj_data.DTW_Mode, 'joint_states');
    if any(joint_mask)
        joint_r50 = nanmean(traj_data.(metric_name)(joint_mask));
        joint_r10 = nanmean(traj_data.("R@10_DTWvsEB")(joint_mask));
        joint_spear = nanmean(traj_data.Spearman_DTWvsEB(joint_mask));
    else
        joint_r50 = NaN; joint_r10 = NaN; joint_spear = NaN;
    end
    
    % Position
    pos_mask = strcmp(traj_data.Embedding_Config, config) & ...
               strcmp(traj_data.DTW_Mode, 'position');
    if any(pos_mask)
        pos_r50 = nanmean(traj_data.(metric_name)(pos_mask));
        pos_r10 = nanmean(traj_data.("R@10_DTWvsEB")(pos_mask));
        pos_spear = nanmean(traj_data.Spearman_DTWvsEB(pos_mask));
    else
        pos_r50 = NaN; pos_r10 = NaN; pos_spear = NaN;
    end
    
    fprintf('%-17s | %4d | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f\n', ...
        config, dims, joint_r50, joint_r10, joint_spear, pos_r50, pos_r10, pos_spear);
end

fprintf('\n');

%% Key Insights
% ========================================================================

fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('KEY INSIGHTS\n');
fprintf('═══════════════════════════════════════════════════════════════\n\n');

% Best overall config (across both spaces)
all_r50 = [];
all_configs = {};
for c = 1:length(configs)
    config_mask = strcmp(traj_data.Embedding_Config, configs{c});
    if any(config_mask)
        all_r50(end+1) = nanmean(traj_data.(metric_name)(config_mask));
        all_configs{end+1} = configs{c};
    end
end

[best_overall_r50, best_idx] = max(all_r50);
fprintf('1. Best Overall Embedding Config: %s (R@50 = %.3f)\n', ...
    all_configs{best_idx}, best_overall_r50);

% Joint vs Position
joint_all = strcmp(traj_data.DTW_Mode, 'joint_states');
pos_all = strcmp(traj_data.DTW_Mode, 'position');

if any(joint_all) && any(pos_all)
    mean_joint = nanmean(traj_data.(metric_name)(joint_all));
    mean_pos = nanmean(traj_data.(metric_name)(pos_all));
    
    fprintf('2. Joint States vs Position: Joint is %.1f%% better (%.3f vs %.3f)\n', ...
        100*(mean_joint - mean_pos)/mean_pos, mean_joint, mean_pos);
end

% Multi-Scale vs Single-Scale
multi_mask = contains(traj_data.Embedding_Config, 'Multi');
single_mask = ~multi_mask;

if any(multi_mask) && any(single_mask)
    mean_multi = nanmean(traj_data.(metric_name)(multi_mask));
    mean_single = nanmean(traj_data.(metric_name)(single_mask));
    
    fprintf('3. Multi-Scale vs Single-Scale: Multi is %.1f%% better (%.3f vs %.3f)\n', ...
        100*(mean_multi - mean_single)/mean_single, mean_multi, mean_single);
end

% Dimensionality effect
dims = unique(traj_data.Total_Dims);
fprintf('4. Dimensionality Effect:\n');
for d = 1:length(dims)
    dim_mask = (traj_data.Total_Dims == dims(d));
    dim_perf = nanmean(traj_data.(metric_name)(dim_mask));
    fprintf('   %3d dims: R@50 = %.3f\n', dims(d), dim_perf);
end

fprintf('\n');

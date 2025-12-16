function data = loadExperimentData(csv_folder)
% LOADEXPERIMENTDATA - Load and organize all experiment CSV files
%
% Usage:
%   data = loadExperimentData('similarity/results/')
%   data = loadExperimentData()  % Uses current directory
%
% Input:
%   csv_folder - Path to folder containing CSV files (optional)
%
% Output:
%   data - Struct with organized experiment data:
%       .raw_tables     - Cell array of all loaded tables
%       .filenames      - Cell array of filenames
%       .combined       - Single combined table with all data
%       .summary        - Summary statistics
%       .metadata       - Experiment metadata structure
%
% Features:
%   - Automatically finds all CSV files in folder
%   - Validates column structure
%   - Combines data from multiple experiments
%   - Provides summary statistics
%   - Handles missing data gracefully

%% ========================================================================
%  INPUT VALIDATION & SETUP
%  ========================================================================

if nargin < 1 || isempty(csv_folder)
    csv_folder = pwd;
    fprintf('No folder specified, using current directory: %s\n', csv_folder);
end

if ~exist(csv_folder, 'dir')
    error('Folder does not exist: %s', csv_folder);
end

fprintf('\n╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  LOADING EXPERIMENT DATA                                       ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

fprintf('Scanning folder: %s\n', csv_folder);

%% ========================================================================
%  FIND ALL CSV FILES
%  ========================================================================

csv_files = dir(fullfile(csv_folder, 'embedding_validation_*.csv'));

if isempty(csv_files)
    error('No CSV files found in: %s', csv_folder);
end

num_files = length(csv_files);
fprintf('✓ Found %d CSV file(s)\n\n', num_files);

%% ========================================================================
%  LOAD ALL CSV FILES
%  ========================================================================

fprintf('Loading CSV files...\n');

raw_tables = cell(num_files, 1);
filenames = cell(num_files, 1);
file_info = struct();

for i = 1:num_files
    filename = csv_files(i).name;
    filepath = fullfile(csv_folder, filename);
    
    fprintf('  [%d/%d] %s ... ', i, num_files, filename);
    
    try
        % Read table
        tbl = readtable(filepath, 'VariableNamingRule', 'preserve');
        
        % Store
        raw_tables{i} = tbl;
        filenames{i} = filename;
        
        % Extract info
        file_info(i).filename = filename;
        file_info(i).rows = height(tbl);
        file_info(i).cols = width(tbl);
        
        % Try to extract timestamp from filename or table
        if ismember('Timestamp', tbl.Properties.VariableNames)
            file_info(i).timestamp = tbl.Timestamp{1};
        else
            file_info(i).timestamp = 'unknown';
        end
        
        fprintf('✓ (%d rows, %d cols)\n', height(tbl), width(tbl));
        
    catch ME
        fprintf('✗ ERROR: %s\n', ME.message);
        raw_tables{i} = [];
        filenames{i} = filename;
        file_info(i).filename = filename;
        file_info(i).rows = 0;
        file_info(i).cols = 0;
        file_info(i).timestamp = 'error';
    end
end

% Remove empty tables
valid_idx = ~cellfun(@isempty, raw_tables);
raw_tables = raw_tables(valid_idx);
filenames = filenames(valid_idx);
file_info = file_info(valid_idx);
num_valid_files = sum(valid_idx);

fprintf('\n✓ Successfully loaded %d/%d files\n\n', num_valid_files, num_files);

if num_valid_files == 0
    error('No valid CSV files could be loaded');
end

%% ========================================================================
%  VALIDATE COLUMN STRUCTURE
%  ========================================================================

fprintf('Validating column structure...\n');

% Expected columns (from your CSV structure)
expected_cols = {
    'Timestamp', 'DTW_Normalization', 'DTW_RotationAlign', ...
    'LB_Kim_Ratio', 'LB_Keogh_Candidates', 'Database_Size', ...
    'Num_Trajectories', 'Num_Segments', 'Top_K', ...
    'Query_Bahn_ID', 'Level', 'Embedding_Config', 'Weight_Mode', ...
    'N_Coarse', 'N_Fine', 'Total_Dims', 'DTW_Mode', ...
    'Spearman_DTWvsEB', 'R@K_DTWvsEB', 'R@50_DTWvsEB', ...
    'R@10_DTWvsEB', 'R@5_DTWvsEB', 'R@3_DTWvsEB', 'R@1_DTWvsEB', ...
    'P_DTWvsEB', ...
    'R@50_GTvsEB', 'R@10_GTvsEB', 'R@5_GTvsEB', ...
    'R@3_GTvsEB', 'R@1_GTvsEB', 'P_GTvsEB', 'Mean_GTvsEB_Rank', ...
    'R@50_GTvsDTW', 'R@10_GTvsDTW', 'R@5_GTvsDTW', ...
    'R@3_GTvsDTW', 'R@1_GTvsDTW', 'P_GTvsDTW', 'Mean_GTvsDTW_Rank', ...
    'Num_GT'
};

% Check each file
all_columns_match = true;
for i = 1:num_valid_files
    tbl = raw_tables{i};
    tbl_cols = tbl.Properties.VariableNames;
    
    % Check if all expected columns exist
    missing_cols = setdiff(expected_cols, tbl_cols);
    extra_cols = setdiff(tbl_cols, expected_cols);
    
    if ~isempty(missing_cols)
        fprintf('  ⚠ File %d missing columns: %s\n', i, strjoin(missing_cols, ', '));
        all_columns_match = false;
    end
    
    if ~isempty(extra_cols)
        fprintf('  ℹ File %d has extra columns: %s\n', i, strjoin(extra_cols, ', '));
    end
end

if all_columns_match
    fprintf('  ✓ All files have consistent column structure\n\n');
else
    fprintf('\n  ⚠ Warning: Some files have different columns\n');
    fprintf('  Continuing anyway - combined table will have all columns\n\n');
end

%% ========================================================================
%  COMBINE ALL TABLES
%  ========================================================================

fprintf('Combining all data into single table...\n');

if num_valid_files == 1
    combined_table = raw_tables{1};
else
    % Vertically concatenate all tables
    combined_table = vertcat(raw_tables{:});
end

fprintf('✓ Combined table: %d rows × %d columns\n\n', ...
    height(combined_table), width(combined_table));

%% ========================================================================
%  EXTRACT METADATA & STATISTICS
%  ========================================================================

fprintf('Extracting metadata and statistics...\n');

metadata = struct();

% Unique values for key dimensions
if ismember('Query_Bahn_ID', combined_table.Properties.VariableNames)
    unique_vals = unique(combined_table.Query_Bahn_ID);
    % Check if numeric or string
    if isnumeric(unique_vals)
        % Convert numeric to cell array of strings
        metadata.unique_queries = arrayfun(@(x) num2str(x), unique_vals, 'UniformOutput', false);
    else
        % Already string, just ensure cell array
        metadata.unique_queries = cellstr(unique_vals);
    end
    metadata.num_queries = length(metadata.unique_queries);
else
    metadata.unique_queries = {};
    metadata.num_queries = 0;
end

if ismember('Embedding_Config', combined_table.Properties.VariableNames)
    metadata.unique_embedding_configs = cellstr(unique(combined_table.Embedding_Config));
    metadata.num_embedding_configs = length(metadata.unique_embedding_configs);
else
    metadata.unique_embedding_configs = {};
    metadata.num_embedding_configs = 0;
end

if ismember('Weight_Mode', combined_table.Properties.VariableNames)
    metadata.unique_weight_modes = cellstr(unique(combined_table.Weight_Mode));
    metadata.num_weight_modes = length(metadata.unique_weight_modes);
else
    metadata.unique_weight_modes = {};
    metadata.num_weight_modes = 0;
end

if ismember('DTW_Mode', combined_table.Properties.VariableNames)
    metadata.unique_dtw_modes = cellstr(unique(combined_table.DTW_Mode));
    metadata.num_dtw_modes = length(metadata.unique_dtw_modes);
else
    metadata.unique_dtw_modes = {};
    metadata.num_dtw_modes = 0;
end

if ismember('Level', combined_table.Properties.VariableNames)
    metadata.unique_levels = cellstr(unique(combined_table.Level));
    metadata.num_levels = length(metadata.unique_levels);
else
    metadata.unique_levels = {};
    metadata.num_levels = 0;
end

% Total experiments
metadata.total_rows = height(combined_table);
metadata.total_experiments = metadata.num_queries * metadata.num_embedding_configs * ...
    metadata.num_weight_modes * metadata.num_dtw_modes * metadata.num_levels;

fprintf('  Queries: %d\n', metadata.num_queries);
fprintf('  Embedding Configs: %d\n', metadata.num_embedding_configs);
fprintf('  Weight Modes: %d\n', metadata.num_weight_modes);
fprintf('  DTW Modes: %d\n', metadata.num_dtw_modes);
fprintf('  Levels: %d\n', metadata.num_levels);
fprintf('  Total rows: %d\n', metadata.total_rows);
fprintf('  Expected experiments: %d\n\n', metadata.total_experiments);

%% ========================================================================
%  CREATE SUMMARY STATISTICS
%  ========================================================================

fprintf('Computing summary statistics...\n');

summary = struct();

% Separate by level
if ismember('Level', combined_table.Properties.VariableNames)
    traj_mask = strcmp(combined_table.Level, 'Trajectory');
    seg_mask = strcmp(combined_table.Level, 'Segment');
    
    summary.num_trajectory_rows = sum(traj_mask);
    summary.num_segment_rows = sum(seg_mask);
    
    % Mean performance metrics (Trajectory level)
    if any(traj_mask) && ismember('Spearman_DTWvsEB', combined_table.Properties.VariableNames)
        traj_data = combined_table(traj_mask, :);
        
        summary.trajectory.mean_spearman = nanmean(traj_data.Spearman_DTWvsEB);
        summary.trajectory.std_spearman = nanstd(traj_data.Spearman_DTWvsEB);
        
        if ismember('P_DTWvsEB', combined_table.Properties.VariableNames)
            summary.trajectory.mean_p_dtw = nanmean(traj_data.P_DTWvsEB);
        end
        
        if ismember('R@10_DTWvsEB', combined_table.Properties.VariableNames)
            summary.trajectory.mean_r10_dtw = nanmean(traj_data.("R@10_DTWvsEB"));
        end
    end
    
    % Mean performance metrics (Segment level)
    if any(seg_mask) && ismember('Spearman_DTWvsEB', combined_table.Properties.VariableNames)
        seg_data = combined_table(seg_mask, :);
        
        summary.segment.mean_spearman = nanmean(seg_data.Spearman_DTWvsEB);
        summary.segment.std_spearman = nanstd(seg_data.Spearman_DTWvsEB);
        
        if ismember('P_DTWvsEB', combined_table.Properties.VariableNames)
            summary.segment.mean_p_dtw = nanmean(seg_data.P_DTWvsEB);
        end
        
        if ismember('R@10_DTWvsEB', combined_table.Properties.VariableNames)
            summary.segment.mean_r10_dtw = nanmean(seg_data.("R@10_DTWvsEB"));
        end
    end
end

fprintf('  Trajectory rows: %d\n', summary.num_trajectory_rows);
fprintf('  Segment rows: %d\n\n', summary.num_segment_rows);

if isfield(summary, 'trajectory')
    fprintf('  Trajectory Level - Mean Spearman: %.4f (±%.4f)\n', ...
        summary.trajectory.mean_spearman, summary.trajectory.std_spearman);
end

if isfield(summary, 'segment')
    fprintf('  Segment Level - Mean Spearman: %.4f (±%.4f)\n', ...
        summary.segment.mean_spearman, summary.segment.std_spearman);
end

fprintf('\n');

%% ========================================================================
%  PACKAGE OUTPUT
%  ========================================================================

data = struct();
data.raw_tables = raw_tables;
data.filenames = filenames;
data.file_info = file_info;
data.combined = combined_table;
data.metadata = metadata;
data.summary = summary;

%% ========================================================================
%  FINAL SUMMARY
%  ========================================================================

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  DATA LOADING COMPLETE                                         ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

fprintf('Output structure:\n');
fprintf('  data.raw_tables     - Cell array of %d tables\n', length(data.raw_tables));
fprintf('  data.filenames      - Cell array of %d filenames\n', length(data.filenames));
fprintf('  data.combined       - Combined table (%d rows)\n', height(data.combined));
fprintf('  data.metadata       - Experiment metadata\n');
fprintf('  data.summary        - Summary statistics\n\n');

fprintf('Quick access examples:\n');
fprintf('  all_data = data.combined;\n');
fprintf('  traj_data = data.combined(strcmp(data.combined.Level, ''Trajectory''), :);\n');
fprintf('  queries = data.metadata.unique_queries;\n\n');

fprintf('✓ Ready for analysis!\n\n');

end
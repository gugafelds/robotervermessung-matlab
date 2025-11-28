%% ========================================================================
%  EXPERIMENT 1C: EMBEDDING ARCHITECTURE VALIDATION
%  ========================================================================
%  Goal: Validate optimal embedding configuration across multiple queries
%  
%  Tests different embedding architectures (single-scale vs multi-scale)
%  with all weight-mode combinations from multimodality ablation
%
%  Dimensions:
%  - 4 Embedding Configurations
%  - 4 Query Trajectories  
%  - 8 Weight-Mode Combinations (4 joint-based + 4 position-based)
%  = 128 Total Experiments
%
%  Runtime: ~4-5 hours (15-20 min preload + 128 experiments × 2 min)
%  Output: results/embedding_validation_YYYY-MM-DD.csv
%  ========================================================================

clear; clc;

% Add current directory and subdirectories
addpath(genpath(pwd));

% Or add specific folders:
addpath(genpath('../main'));
addpath(genpath('../lasertracker'));
addpath(genpath('../methods'));

fprintf('✓ Paths added successfully\n\n');

fprintf('========================================\n');
fprintf('EXPERIMENT: EMBEDDING ARCHITECTURE VALIDATION\n');
fprintf('========================================\n');

% === Base Configuration ===
base_config = struct();
base_config.database_sample_size = 1500;  % Fixed for fair comparison
base_config.random_seed = 42;
base_config.top_k_trajectories = 500;     % Fixed

% === DIMENSION 1: Embedding Architectures ===
embedding_configs = {
    % Name,                  n_coarse, n_fine, use_multi_scale
    'Single-Fine-100',        0,        100,     false;
    'Multi-Balanced-150',    50,       100,     true;
    'Single-Fine-300',       0,        300,    false;
    'Multi-Dense-350',       50,       300,    true;
};

% === DIMENSION 2: Query Trajectories ===
query_ids = {
    '1764336135';   % Baseline (currently used in other experiments)
    %'1764336030';   % First in original list
    %'1764169991';   % Middle of list
    %'1764334766';   % Last in original list
};

% === DIMENSION 3: DTW Mode + Weight Combinations ===
% (Already intelligently paired - joint weights with joint_states, position weights with position)
weight_mode_configs = {
    % Name,                      DTW_Mode,        [Pos, Joint, Orient, Vel, Meta]
    % Joint-based DTW experiments
    'Joint only',                'joint_states',  [0,   1,     0,      0,   0];
    'Joint + Orientation',       'joint_states',  [0,   1,     1,      0,   0];
    'Joint + Metadata',          'joint_states',  [0,   1,     0,      0,   1];
    'Joint + Orient + Meta',     'joint_states',  [0,   1,     1,      0,   1];
    
    % Position-based DTW experiments
    'Position only',             'position',      [1,   0,     0,      0,   0];
    'Pos + Velocity',            'position',      [1,   0,     0,      1,   0];
    'Pos + Metadata',            'position',      [1,   0,     0,      0,   1];
    'Pos + Vel + Meta',          'position',      [1,   0,     0,      1,   1];
};

% === Calculate Total Experiments ===
num_embeddings = size(embedding_configs, 1);      % 4
num_queries = size(query_ids, 1);                 % 4
num_weight_modes = size(weight_mode_configs, 1);  % 8

total_experiments = num_embeddings * num_queries * num_weight_modes;  % 128

fprintf('Total experiments: %d\n', total_experiments);
fprintf('  - Embedding configs: %d\n', num_embeddings);
fprintf('  - Query trajectories: %d\n', num_queries);
fprintf('  - Weight-mode combinations: %d\n', num_weight_modes);
fprintf('Estimated runtime WITH pre-loading: %.1f hours\n', (20/60) + (total_experiments * 2 / 60));
fprintf('  (15-20 min preload + 128 × ~2 min compute)\n\n');

%% ========================================================================
%  DATABASE CONNECTION & SAMPLING
%  ========================================================================

fprintf('=== Connecting to Database ===\n');

conn = connectingToPostgres;

if isopen(conn)
    fprintf('✓ Database connection successful\n\n');
else
    error('✗ Database connection failed');
end

% === Database Sampling ===
fprintf('=== Database Sampling ===\n');

schema = 'bewegungsdaten';

% Load full database metadata
full_db_query = sprintf(...
    ['SELECT bahn_id FROM robotervermessung.%s.bahn_metadata ' ...
     'WHERE bahn_id = segment_id'], schema);

full_db_metadata = fetch(conn, full_db_query);
num_total_trajectories = height(full_db_metadata);

fprintf('  Total trajectories in database: %d\n', num_total_trajectories);
fprintf('  Target sample size: %d (%.1f%%)\n', ...
    base_config.database_sample_size, ...
    100 * base_config.database_sample_size / num_total_trajectories);
fprintf('  Random seed: %d\n', base_config.random_seed);

% Set random seed
rng(base_config.random_seed);

% Random sample
sample_indices = randperm(num_total_trajectories, base_config.database_sample_size);
sampled_metadata = full_db_metadata(sample_indices, :);
candidate_ids = sampled_metadata.bahn_id;

fprintf('  ✓ Sampled %d trajectories\n\n', length(candidate_ids));

%% ========================================================================
%  PRE-LOAD ALL DATA (ONE-TIME COST)
%  ========================================================================

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  PRE-LOADING ALL EXPERIMENT DATA                               ║\n');
fprintf('║  This loads all data ONCE before running experiments          ║\n');
fprintf('║  Estimated time: 15-20 minutes                                 ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

preload_tic = tic;

chunk_size = 50;  % For batch loading

data_cache = loadDataExperiment(conn, schema, candidate_ids, query_ids, chunk_size);

preload_time = toc(preload_tic);

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  DATA PRE-LOADING COMPLETED                                    ║\n');
fprintf('║  Time: %.1f minutes                                            ║\n', preload_time/60);
fprintf('║  Memory: %.1f MB                                               ║\n', whos('data_cache').bytes / 1e6);
fprintf('║  Ready to run %d experiments!                                  ║\n', total_experiments);
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

%% ========================================================================
%  PRE-COMPUTE ALL DTW (ONE-TIME COST)
%  ========================================================================

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  PRE-COMPUTING DTW FOR ALL QUERIES & MODES                    ║\n');
fprintf('║  This computes DTW ONCE for reuse in all experiments          ║\n');
fprintf('║  Estimated time: 10-15 minutes                                 ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

dtw_tic = tic;

% Prepare config for DTW pre-computation
dtw_config = struct();
dtw_config.top_k_trajectories = base_config.top_k_trajectories;
dtw_config.lb_kim_keep_ratio = 0.8;
dtw_config.lb_keogh_candidates = 500;
dtw_config.cdtw_window = 0.15;
dtw_config.normalize_dtw = true;
dtw_config.use_rotation_alignment = false;

dtw_cache = precomputeDTW(data_cache, query_ids, dtw_config);

dtw_time = toc(dtw_tic);

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  DTW PRE-COMPUTATION COMPLETED                                 ║\n');
fprintf('║  Time: %.1f minutes                                            ║\n', dtw_time/60);
fprintf('║  Memory: %.1f MB                                               ║\n', whos('dtw_cache').bytes / 1e6);
fprintf('║  All experiments will skip DTW and only compute embeddings!    ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

% Close database connection (no longer needed!)
close(conn);
fprintf('✓ Database connection closed (no longer needed)\n\n');

%% ========================================================================
%  RUN EXPERIMENTS (USING PRE-LOADED DATA)
%  ========================================================================

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  STARTING EXPERIMENTS WITH PRE-LOADED DATA                     ║\n');
fprintf('║  No database queries will be made during experiments!          ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

% === Storage ===
all_results = cell(total_experiments, 1);
experiment_counter = 0;

% === Run Experiments ===
experiment_start = tic;

for emb_idx = 1:num_embeddings
    for query_idx = 1:num_queries
        for wm_idx = 1:num_weight_modes
            
            experiment_counter = experiment_counter + 1;
            
            % Build experiment name
            exp_name = sprintf('%s | Query_%s | %s', ...
                embedding_configs{emb_idx, 1}, ...
                query_ids{query_idx}, ...
                weight_mode_configs{wm_idx, 1});
            
            fprintf('\n╔════════════════════════════════════════════════════════════════╗\n');
            fprintf('║  EXPERIMENT %3d/%3d                                           ║\n', ...
                experiment_counter, total_experiments);
            fprintf('╠════════════════════════════════════════════════════════════════╣\n');
            fprintf('║  Embedding:   %-48s ║\n', embedding_configs{emb_idx, 1});
            fprintf('║  Query ID:    %-48s ║\n', query_ids{query_idx});
            fprintf('║  Weight Mode: %-48s ║\n', weight_mode_configs{wm_idx, 1});
            fprintf('╚════════════════════════════════════════════════════════════════╝\n');
            
            % Build config
            config = base_config;
            config.exp_name = exp_name;
            
            % Query
            config.query_bahn_id = query_ids{query_idx};
            
            % Embedding architecture
            config.n_coarse = embedding_configs{emb_idx, 2};
            config.n_fine = embedding_configs{emb_idx, 3};
            config.use_multi_scale = embedding_configs{emb_idx, 4};
            
            % DTW mode and weights
            config.dtw_mode = weight_mode_configs{wm_idx, 2};
            config.weights = weight_mode_configs{wm_idx, 3};
            
            % Metadata for results
            config.embedding_config_name = embedding_configs{emb_idx, 1};
            config.weight_mode_name = weight_mode_configs{wm_idx, 1};
            
            % ⭐ CRITICAL: Pass pre-loaded caches
            config.data_cache = data_cache;
            config.dtw_cache = dtw_cache;
            
            fprintf('Config: n_coarse=%d, n_fine=%d, multi_scale=%d, mode=%s\n', ...
                config.n_coarse, config.n_fine, config.use_multi_scale, config.dtw_mode);
            fprintf('Weights: [%.1f, %.1f, %.1f, %.1f, %.1f]\n', config.weights);
            fprintf('Using pre-loaded data cache: YES\n');
            fprintf('Using pre-computed DTW cache: YES\n');
            
            % Run experiment
            exp_tic = tic;
            results = runExperiment(config);
            exp_time = toc(exp_tic);
            
            % Store metadata in results
            results.embedding_config = embedding_configs{emb_idx, 1};
            results.query_id = query_ids{query_idx};
            results.weight_mode = weight_mode_configs{wm_idx, 1};
            results.exp_runtime = exp_time;
            
            all_results{experiment_counter} = results;
            
            % Progress update
            elapsed = toc(experiment_start);
            avg_time_per_exp = elapsed / experiment_counter;
            remaining_exp = total_experiments - experiment_counter;
            eta_seconds = remaining_exp * avg_time_per_exp;
            eta_hours = eta_seconds / 3600;
            
            fprintf('✓ Experiment %d/%d completed in %.1f min\n', ...
                experiment_counter, total_experiments, exp_time/60);
            fprintf('  Progress: %.1f%% | Elapsed: %.1fh | ETA: %.1fh\n\n', ...
                100*experiment_counter/total_experiments, elapsed/3600, eta_hours);
        end
    end
end

total_time = toc(experiment_start);

%% ========================================================================
%  CREATE RESULTS TABLES (TRAJECTORY + SEGMENT LEVEL)
%  ========================================================================

fprintf('\n=== Creating Results Tables ===\n');

% Preallocate matrices
traj_results_matrix = zeros(total_experiments, 6);  % 6 metrics
seg_results_matrix = zeros(total_experiments, 6);   % 6 metrics

% Metadata arrays
embedding_names = cell(total_experiments, 1);
query_bahn_ids = cell(total_experiments, 1);
weight_modes = cell(total_experiments, 1);
dtw_modes = cell(total_experiments, 1);
n_coarse_arr = zeros(total_experiments, 1);
n_fine_arr = zeros(total_experiments, 1);
total_dims_arr = zeros(total_experiments, 1);
multi_scale_arr = false(total_experiments, 1);
runtimes = zeros(total_experiments, 1);

% Weights arrays
w_pos = zeros(total_experiments, 1);
w_joint = zeros(total_experiments, 1);
w_orient = zeros(total_experiments, 1);
w_vel = zeros(total_experiments, 1);
w_meta = zeros(total_experiments, 1);

for i = 1:total_experiments
    r = all_results{i};
    
    % Metadata
    embedding_names{i} = r.embedding_config;
    query_bahn_ids{i} = r.query_id;
    weight_modes{i} = r.weight_mode;
    dtw_modes{i} = r.dtw_mode;
    n_coarse_arr(i) = r.n_coarse;
    n_fine_arr(i) = r.n_fine;
    total_dims_arr(i) = r.total_dims;
    
    % Multi-scale detection
    multi_scale_arr(i) = (r.n_coarse > 0);
    
    runtimes(i) = r.exp_runtime / 60;  % Convert to minutes
    
    % Weights
    w_pos(i) = r.weight_pos;
    w_joint(i) = r.weight_joint;
    w_orient(i) = r.weight_orient;
    w_vel(i) = r.weight_vel;
    w_meta(i) = r.weight_meta;
    
    % Trajectory-level metrics
    traj_results_matrix(i, :) = [
        r.spearman, r.p_at_k, r.p_at_10, r.p_at_5, r.p_at_3, r.p_at_1
    ];
    
    % Segment-level metrics
    seg_results_matrix(i, :) = [
        r.seg_spearman, r.seg_p_at_k, r.seg_p_at_10, ...
        r.seg_p_at_5, r.seg_p_at_3, r.seg_p_at_1
    ];
end

%% ========================================================================
%  BUILD TRAJECTORY-LEVEL TABLE
%  ========================================================================

traj_table = array2table(traj_results_matrix, ...
    'VariableNames', {'Spearman', sprintf('P@%d', base_config.top_k_trajectories), ...
                      'P@10', 'P@5', 'P@3', 'P@1'});

% Add configuration columns (from right to left to maintain order)
traj_table = addvars(traj_table, runtimes, 'Before', 'Spearman', ...
    'NewVariableNames', 'Runtime_min');
traj_table = addvars(traj_table, w_meta, 'Before', 'Runtime_min', ...
    'NewVariableNames', 'W_Meta');
traj_table = addvars(traj_table, w_vel, 'Before', 'W_Meta', ...
    'NewVariableNames', 'W_Vel');
traj_table = addvars(traj_table, w_orient, 'Before', 'W_Vel', ...
    'NewVariableNames', 'W_Orient');
traj_table = addvars(traj_table, w_joint, 'Before', 'W_Orient', ...
    'NewVariableNames', 'W_Joint');
traj_table = addvars(traj_table, w_pos, 'Before', 'W_Joint', ...
    'NewVariableNames', 'W_Pos');
traj_table = addvars(traj_table, multi_scale_arr, 'Before', 'W_Pos', ...
    'NewVariableNames', 'Multi_Scale');
traj_table = addvars(traj_table, total_dims_arr, 'Before', 'Multi_Scale', ...
    'NewVariableNames', 'Total_Dims');
traj_table = addvars(traj_table, n_fine_arr, 'Before', 'Total_Dims', ...
    'NewVariableNames', 'N_Fine');
traj_table = addvars(traj_table, n_coarse_arr, 'Before', 'N_Fine', ...
    'NewVariableNames', 'N_Coarse');
traj_table = addvars(traj_table, dtw_modes, 'Before', 'N_Coarse', ...
    'NewVariableNames', 'DTW_Mode');
traj_table = addvars(traj_table, weight_modes, 'Before', 'DTW_Mode', ...
    'NewVariableNames', 'Weight_Mode');
traj_table = addvars(traj_table, query_bahn_ids, 'Before', 'Weight_Mode', ...
    'NewVariableNames', 'Query_Bahn_ID');
traj_table = addvars(traj_table, embedding_names, 'Before', 'Query_Bahn_ID', ...
    'NewVariableNames', 'Embedding_Config');

% Add level identifier
level_traj = repmat({'Trajectory'}, total_experiments, 1);
traj_table = addvars(traj_table, level_traj, 'Before', 'Embedding_Config', ...
    'NewVariableNames', 'Level');

%% ========================================================================
%  BUILD SEGMENT-LEVEL TABLE
%  ========================================================================

seg_table = array2table(seg_results_matrix, ...
    'VariableNames', {'Spearman', sprintf('P@%d', base_config.top_k_trajectories), ...
                      'P@10', 'P@5', 'P@3', 'P@1'});

% Add configuration columns (same structure as trajectory table)
seg_table = addvars(seg_table, runtimes, 'Before', 'Spearman', ...
    'NewVariableNames', 'Runtime_min');
seg_table = addvars(seg_table, w_meta, 'Before', 'Runtime_min', ...
    'NewVariableNames', 'W_Meta');
seg_table = addvars(seg_table, w_vel, 'Before', 'W_Meta', ...
    'NewVariableNames', 'W_Vel');
seg_table = addvars(seg_table, w_orient, 'Before', 'W_Vel', ...
    'NewVariableNames', 'W_Orient');
seg_table = addvars(seg_table, w_joint, 'Before', 'W_Orient', ...
    'NewVariableNames', 'W_Joint');
seg_table = addvars(seg_table, w_pos, 'Before', 'W_Joint', ...
    'NewVariableNames', 'W_Pos');
seg_table = addvars(seg_table, multi_scale_arr, 'Before', 'W_Pos', ...
    'NewVariableNames', 'Multi_Scale');
seg_table = addvars(seg_table, total_dims_arr, 'Before', 'Multi_Scale', ...
    'NewVariableNames', 'Total_Dims');
seg_table = addvars(seg_table, n_fine_arr, 'Before', 'Total_Dims', ...
    'NewVariableNames', 'N_Fine');
seg_table = addvars(seg_table, n_coarse_arr, 'Before', 'N_Fine', ...
    'NewVariableNames', 'N_Coarse');
seg_table = addvars(seg_table, dtw_modes, 'Before', 'N_Coarse', ...
    'NewVariableNames', 'DTW_Mode');
seg_table = addvars(seg_table, weight_modes, 'Before', 'DTW_Mode', ...
    'NewVariableNames', 'Weight_Mode');
seg_table = addvars(seg_table, query_bahn_ids, 'Before', 'Weight_Mode', ...
    'NewVariableNames', 'Query_Bahn_ID');
seg_table = addvars(seg_table, embedding_names, 'Before', 'Query_Bahn_ID', ...
    'NewVariableNames', 'Embedding_Config');

% Add level identifier
level_seg = repmat({'Segment'}, total_experiments, 1);
seg_table = addvars(seg_table, level_seg, 'Before', 'Embedding_Config', ...
    'NewVariableNames', 'Level');

%% ========================================================================
%  COMBINE TABLES
%  ========================================================================

combined_table = [traj_table; seg_table];

%% ========================================================================
%  SAVE RESULTS
%  ========================================================================

output_dir = 'similarity/results';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd''T''HHmmss'));
output_file = fullfile(output_dir, sprintf('embedding_validation_%s.csv', timestamp));
writetable(combined_table, output_file, 'WriteRowNames', false);

fprintf('✓ Results saved to: %s\n\n', output_file);

%% ========================================================================
%  ANALYSIS & SUMMARY
%  ========================================================================

fprintf('========================================\n');
fprintf('EXPERIMENT COMPLETED\n');
fprintf('========================================\n');
fprintf('Pre-loading time: %.1f minutes\n', preload_time/60);
fprintf('DTW pre-computation time: %.1f minutes\n', dtw_time/60);
fprintf('Experiments runtime: %.2f hours (%.1f minutes)\n', total_time/3600, total_time/60);
fprintf('Total runtime: %.2f hours\n', (preload_time + dtw_time + total_time)/3600);
fprintf('Total Experiments: %d\n', total_experiments);
fprintf('Average time per experiment: %.2f minutes (embeddings only!)\n', total_time/60/total_experiments);
fprintf('Results file: %s\n\n', output_file);


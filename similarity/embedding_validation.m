%  EMBEDDING VALIDATION
%  ========================================================================
%  Goal: Validate optimal embedding configuration across multiple queries
%  
%  Tests different embedding architectures (single-scale vs multi-scale)
%  with all weight-mode combinations from multimodality ablation
%
%  - Pre-loads data ONCE (15-20 min)
%  - Pre-computes DTW ONCE (10-15 min)
%  - Pre-computes EMBEDDINGS ONCE (10-15 min)
%  - Experiments only apply different WEIGHTS via RRF fusion (~5 seconds each!)
%
%  Dimensions:
%  - 4 Embedding Configurations
%  - 7 Query Trajectories  
%  - 12 Weight-Mode Combinations
%  = 336 Total Experiments
%
%  Runtime: ~1 hour total! (40 min preload + 336 × 5 sec = 28 min)
%  Output: results/embedding_validation_YYYY-MM-DD.csv
%  ========================================================================

clear; clc;

% Add current directory and subdirectories
addpath(genpath(pwd));

% Or add specific folders:
addpath(genpath('../main'));
addpath(genpath('../lasertracker'));
addpath(genpath('../methods'));
addpath(genpath('results'));

fprintf('✓ Paths added successfully\n\n');

fprintf('========================================\n');
fprintf('EXPERIMENT: EMBEDDING ARCHITECTURE VALIDATION\n');
fprintf('WITH EMBEDDINGS PRE-COMPUTATION\n');
fprintf('========================================\n');

use_ground_truth = true;  % Set to false to disable

% === Base Configuration ===
base_config = struct();
base_config.database_sample_size = 1170;  % Fixed for fair comparison
base_config.random_seed = 42;
base_config.top_k_trajectories = 100;     % Fixed

% === DIMENSION 1: Embedding Architectures ===
embedding_configs = {
    % Name,                  n_coarse, n_fine, use_multi_scale
    'Single-Fine-75',        0,        75,     false;
    %'Multi-Balanced-100',    25,       75,     true;
    %'Single-Fine-150',       0,        150,    false;
    %'Multi-Dense-200',       50,       150,    true;
};

% === DIMENSION 2: Query Trajectories ===
query_ids = {
    '1763567277';   % gt = 9
    '1763567148';   % gt = 9
    '1763567026';   % gt = 9

};

% === DIMENSION 3: DTW Mode + Weight Combinations ===
weight_mode_configs = {
    % Name,                      DTW_Mode,        [Pos, Joint, Orient, Vel, Meta]
    % Joint-based DTW experiments
    'Joint only',                'joint_states',  [0,   1,     0,      0,   0];
    'Joint + Orientation',       'joint_states',  [0,   1,     1,      0,   0];
    'Joint + Metadata',          'joint_states',  [0,   1,     0,      0,   1];
    'Joint + Orient + Meta',     'joint_states',  [0,   1,     1,      0,   1];
    'Joint + Position',          'joint_states',  [1,   1,     0,      0,   0];
    'Meta',                      'joint_states',  [0,   0,     0,      0,   1];
    
    % Position-based DTW experiments
    'Position only',             'position',      [1,   0,     0,      0,   0];
    'Pos + Velocity',            'position',      [1,   0,     0,      1,   0];
    'Pos + Metadata',            'position',      [1,   0,     0,      0,   1];
    'Pos + Vel + Meta',          'position',      [1,   0,     0,      1,   1];
    'Position + Joint',          'position',      [1,   1,     0,      0,   0];
    'Meta',                      'position',      [0,   0,     0,      0,   1];
};

% === Calculate Total Experiments ===
num_embeddings = size(embedding_configs, 1);      % 4
num_queries = size(query_ids, 1);                 % 7
num_weight_modes = size(weight_mode_configs, 1);  % 12

total_experiments = num_embeddings * num_queries * num_weight_modes;  % 336

fprintf('Total experiments: %d\n', total_experiments);
fprintf('  - Embedding configs: %d\n', num_embeddings);
fprintf('  - Query trajectories: %d\n', num_queries);
fprintf('  - Weight-mode combinations: %d\n', num_weight_modes);
fprintf('\n⭐ OPTIMIZED RUNTIME ESTIMATE:\n');
fprintf('  Pre-loading: ~20 min\n');
fprintf('  DTW pre-compute: ~15 min\n');
fprintf('  Embeddings pre-compute: ~15 min\n');
fprintf('  Experiments: ~%.0f min (%d × 5 sec)\n', total_experiments * 5 / 60, total_experiments);
fprintf('  TOTAL: ~%.1f hours\n\n', (20 + 15 + 15 + (total_experiments * 5 / 60)) / 60);

% ========================================================================
%%  DATABASE CONNECTION & SAMPLING
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

% ========================================================================
%%  GROUND TRUTH ADDITION (OPTIONAL)
%  ========================================================================

if use_ground_truth
    % Get ground truth trajectories from same recordings as queries
    [ground_truth_ids, ground_truth_map] = getGTCandidates(conn, schema, query_ids);
    
    if ~isempty(ground_truth_ids)
        % Remove ground truth from random sample to avoid duplicates
        candidate_ids = setdiff(candidate_ids, ground_truth_ids);
        
        % ⭐ ADD GROUND TRUTH TO CANDIDATES - THIS LINE WAS MISSING!
        candidate_ids = [candidate_ids; ground_truth_ids];
        
        fprintf('=== Updated Candidate Pool ===\n');
        fprintf('  Random candidates: %d\n', length(candidate_ids) - length(ground_truth_ids));
        fprintf('  Ground truth: %d\n', length(ground_truth_ids));
        fprintf('  Total candidates: %d\n\n', length(candidate_ids));
        
        % Store ground truth info for later analysis
        base_config.ground_truth_map = ground_truth_map;
        base_config.has_ground_truth = true;
    else
        fprintf('⚠ No ground truth found - continuing with random sample only\n\n');
        base_config.has_ground_truth = false;
    end
else
    fprintf('=== Ground Truth Disabled ===\n');
    fprintf('  Using only random sample: %d trajectories\n\n', length(candidate_ids));
    base_config.has_ground_truth = false;
end

% ========================================================================
%%  PRE-LOAD ALL DATA (ONE-TIME COST)
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
fprintf('║  Ready to pre-compute DTW and embeddings!                      ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

% ========================================================================
%%  PRE-COMPUTE ALL DTW (ONE-TIME COST)
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
dtw_config.lb_kim_keep_ratio = 0.6;
dtw_config.lb_keogh_candidates = 200;
dtw_config.cdtw_window = 0.10;
dtw_config.normalize_dtw = true;
dtw_config.use_rotation_alignment = false;

dtw_cache = precomputeDTW(data_cache, query_ids, dtw_config);

dtw_time = toc(dtw_tic);

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  DTW PRE-COMPUTATION COMPLETED                                 ║\n');
fprintf('║  Time: %.1f minutes                                            ║\n', dtw_time/60);
fprintf('║  Memory: %.1f MB                                               ║\n', whos('dtw_cache').bytes / 1e6);
fprintf('║  All experiments will skip DTW computation!                    ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

% ========================================================================
%%  PRE-COMPUTE ALL EMBEDDINGS (ONE-TIME COST)
%  ========================================================================

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  PRE-COMPUTING EMBEDDINGS FOR ALL QUERIES & CONFIGS           ║\n');
fprintf('║  This computes embeddings ONCE for reuse in all experiments   ║\n');
fprintf('║  Estimated time: 10-15 minutes                                 ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

emb_tic = tic;

% Prepare config for embeddings pre-computation
emb_config = struct();
emb_config.norm_strategy = 'max_extent';

embeddings_cache = precomputeEmbeddings(data_cache, query_ids, embedding_configs, emb_config);

emb_time = toc(emb_tic);

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  EMBEDDINGS PRE-COMPUTATION COMPLETED                          ║\n');
fprintf('║  Time: %.1f minutes                                            ║\n', emb_time/60);
fprintf('║  Memory: %.1f MB                                               ║\n', whos('embeddings_cache').bytes / 1e6);
fprintf('║  All experiments will skip embedding computation!              ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

% Close database connection (no longer needed!)
close(conn);
fprintf('✓ Database connection closed (no longer needed)\n\n');

% ========================================================================
%%  RUN EXPERIMENTS (USING PRE-LOADED DATA + DTW + EMBEDDINGS)
%  ========================================================================

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  STARTING EXPERIMENTS WITH ALL PRE-COMPUTED CACHES             ║\n');
fprintf('║  No database queries, no DTW, no embedding computation!        ║\n');
fprintf('║  Only applying different weights via RRF fusion!               ║\n');
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
            
            % ⭐ CRITICAL: Pass all pre-loaded caches
            config.data_cache = data_cache;
            config.dtw_cache = dtw_cache;
            config.embeddings_cache = embeddings_cache;
            
            fprintf('Config: n_coarse=%d, n_fine=%d, multi_scale=%d, mode=%s\n', ...
                config.n_coarse, config.n_fine, config.use_multi_scale, config.dtw_mode);
            fprintf('Weights: [%.1f, %.1f, %.1f, %.1f, %.1f]\n', config.weights);
            fprintf('Using pre-loaded data cache: YES\n');
            fprintf('Using pre-computed DTW cache: YES\n');
            fprintf('Using pre-computed embeddings cache: YES\n');
            
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
            eta_minutes = eta_seconds / 60;
            
            fprintf('✓ Experiment %d/%d completed in %.1f seconds\n', ...
                experiment_counter, total_experiments, exp_time);
            fprintf('  Progress: %.1f%% | Elapsed: %.1f min | ETA: %.1f min\n\n', ...
                100*experiment_counter/total_experiments, elapsed/60, eta_minutes);
        end
    end
end

total_time = toc(experiment_start);

% ========================================================================
%%  CREATE RESULTS TABLES (TRAJECTORY + SEGMENT LEVEL)
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

% ========================================================================
% DTW vs Embedding Coverage metrics arrays (trajectory level) - EXISTING
% ========================================================================
coverage_10 = zeros(total_experiments, 1);
coverage_50 = zeros(total_experiments, 1);
expansion_10 = zeros(total_experiments, 1);
expansion_50 = zeros(total_experiments, 1);
recall_cov_10 = zeros(total_experiments, 1);
recall_cov_50 = zeros(total_experiments, 1);

% DTW vs Embedding Coverage metrics arrays (segment level) - EXISTING
seg_coverage_10 = zeros(total_experiments, 1);
seg_coverage_50 = zeros(total_experiments, 1);
seg_expansion_10 = zeros(total_experiments, 1);
seg_expansion_50 = zeros(total_experiments, 1);
seg_recall_cov_10 = zeros(total_experiments, 1);
seg_recall_cov_50 = zeros(total_experiments, 1);

% ========================================================================
% GT Coverage metrics arrays (trajectory level) - NEW!
% ========================================================================
gt_recall_10 = zeros(total_experiments, 1);
gt_recall_50 = zeros(total_experiments, 1);
gt_expansion_10 = zeros(total_experiments, 1);
gt_expansion_50 = zeros(total_experiments, 1);

% GT Coverage metrics arrays (segment level) - NEW!
seg_gt_recall_10 = zeros(total_experiments, 1);
seg_gt_recall_50 = zeros(total_experiments, 1);
seg_gt_expansion_10 = zeros(total_experiments, 1);
seg_gt_expansion_50 = zeros(total_experiments, 1);

% Additional GT statistics arrays - NEW!
num_gt_arr = zeros(total_experiments, 1);
mean_gt_rank_arr = zeros(total_experiments, 1);

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
    
    % ====================================================================
    % DTW vs Embedding Coverage metrics (trajectory level) - EXISTING
    % ====================================================================
    if isfield(r, 'coverage_at_10')
        coverage_10(i) = r.coverage_at_10;
        expansion_10(i) = r.expansion_ratio_10;
        recall_cov_10(i) = r.recall_at_10;
    else
        coverage_10(i) = NaN;
        expansion_10(i) = NaN;
        recall_cov_10(i) = NaN;
    end
    
    if isfield(r, 'coverage_at_50')
        coverage_50(i) = r.coverage_at_50;
        expansion_50(i) = r.expansion_ratio_50;
        recall_cov_50(i) = r.recall_at_50;
    else
        coverage_50(i) = NaN;
        expansion_50(i) = NaN;
        recall_cov_50(i) = NaN;
    end
    
    % ====================================================================
    % DTW vs Embedding Coverage metrics (segment level) - EXISTING
    % ====================================================================
    if isfield(r, 'seg_coverage_at_10')
        seg_coverage_10(i) = r.seg_coverage_at_10;
        seg_expansion_10(i) = r.seg_expansion_ratio_10;
        seg_recall_cov_10(i) = r.seg_recall_at_10;
    else
        seg_coverage_10(i) = NaN;
        seg_expansion_10(i) = NaN;
        seg_recall_cov_10(i) = NaN;
    end
    
    if isfield(r, 'seg_coverage_at_50')
        seg_coverage_50(i) = r.seg_coverage_at_50;
        seg_expansion_50(i) = r.seg_expansion_ratio_50;
        seg_recall_cov_50(i) = r.seg_recall_at_50;
    else
        seg_coverage_50(i) = NaN;
        seg_expansion_50(i) = NaN;
        seg_recall_cov_50(i) = NaN;
    end
    
    % ====================================================================
    % GT Coverage metrics (trajectory level) - NEW!
    % ====================================================================
    if isfield(r, 'gt_recall_at_10')
        gt_recall_10(i) = r.gt_recall_at_10;
        gt_expansion_10(i) = r.gt_expansion_ratio_10;
    else
        gt_recall_10(i) = NaN;
        gt_expansion_10(i) = NaN;
    end
    
    if isfield(r, 'gt_recall_at_50')
        gt_recall_50(i) = r.gt_recall_at_50;
        gt_expansion_50(i) = r.gt_expansion_ratio_50;
    else
        gt_recall_50(i) = NaN;
        gt_expansion_50(i) = NaN;
    end
    
    % ====================================================================
    % GT Coverage metrics (segment level) - NEW!
    % ====================================================================
    if isfield(r, 'seg_gt_recall_at_10')
        seg_gt_recall_10(i) = r.seg_gt_recall_at_10;
        seg_gt_expansion_10(i) = r.seg_gt_expansion_ratio_10;
    else
        seg_gt_recall_10(i) = NaN;
        seg_gt_expansion_10(i) = NaN;
    end
    
    if isfield(r, 'seg_gt_recall_at_50')
        seg_gt_recall_50(i) = r.seg_gt_recall_at_50;
        seg_gt_expansion_50(i) = r.seg_gt_expansion_ratio_50;
    else
        seg_gt_recall_50(i) = NaN;
        seg_gt_expansion_50(i) = NaN;
    end
    
    % ====================================================================
    % Additional GT statistics - NEW!
    % ====================================================================
    if isfield(r, 'num_gt')
        num_gt_arr(i) = r.num_gt;
    else
        num_gt_arr(i) = NaN;
    end
    
    if isfield(r, 'mean_gt_rank')
        mean_gt_rank_arr(i) = r.mean_gt_rank;
    else
        mean_gt_rank_arr(i) = NaN;
    end
end

% ========================================================================
%%  BUILD TRAJECTORY-LEVEL TABLE
%  ========================================================================

traj_table = array2table(traj_results_matrix, ...
    'VariableNames', {'Spearman', sprintf('P@%d', base_config.top_k_trajectories), ...
                      'P@10', 'P@5', 'P@3', 'P@1'});

% Runtime (rightmost column)
traj_table = addvars(traj_table, runtimes, 'After', 'P@1', ...
    'NewVariableNames', 'Runtime_min');

% GT Statistics
traj_table = addvars(traj_table, mean_gt_rank_arr, 'Before', 'Runtime_min', ...
    'NewVariableNames', 'Mean_GT_Rank');
traj_table = addvars(traj_table, num_gt_arr, 'Before', 'Mean_GT_Rank', ...
    'NewVariableNames', 'Num_GT');

% GT Coverage - K=50 metrics (before Num_GT)
traj_table = addvars(traj_table, gt_expansion_50, 'Before', 'Num_GT', ...
    'NewVariableNames', 'ER@50_GT');
traj_table = addvars(traj_table, gt_recall_50, 'Before', 'ER@50_GT', ...
    'NewVariableNames', 'R@50_GT');

% GT Coverage - K=10 metrics (before R@50_GT)
traj_table = addvars(traj_table, gt_expansion_10, 'Before', 'R@50_GT', ...
    'NewVariableNames', 'ER@10_GT');
traj_table = addvars(traj_table, gt_recall_10, 'Before', 'ER@10_GT', ...
    'NewVariableNames', 'R@10_GT');

% DTW vs Embedding Coverage - K=50 metrics (before R@10_GT)
traj_table = addvars(traj_table, expansion_50, 'Before', 'R@10_GT', ...
    'NewVariableNames', 'ER@50_all');
traj_table = addvars(traj_table, recall_cov_50, 'Before', 'ER@50_all', ...
    'NewVariableNames', 'R@50_all');

% DTW vs Embedding Coverage - K=10 metrics (before R@50_all)
traj_table = addvars(traj_table, expansion_10, 'Before', 'R@50_all', ...
    'NewVariableNames', 'ER@10_all');
traj_table = addvars(traj_table, recall_cov_10, 'Before', 'ER@10_all', ...
    'NewVariableNames', 'R@10_all');

% Add configuration columns (from right to left to maintain order)
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

% ========================================================================
%%  BUILD SEGMENT-LEVEL TABLE
%  ========================================================================

seg_table = array2table(seg_results_matrix, ...
    'VariableNames', {'Spearman', sprintf('P@%d', base_config.top_k_trajectories), ...
                      'P@10', 'P@5', 'P@3', 'P@1'});

% ========================================================================
% Add columns in SAME ORDER as trajectory table
% Use segment data but SAME column names (Level distinguishes!)
% ========================================================================

% Runtime
seg_table = addvars(seg_table, runtimes, 'After', 'P@1', ...
    'NewVariableNames', 'Runtime_min');

% GT Statistics (use segment GT stats where available, otherwise use traj stats)
seg_table = addvars(seg_table, mean_gt_rank_arr, 'Before', 'Runtime_min', ...
    'NewVariableNames', 'Mean_GT_Rank');
seg_table = addvars(seg_table, num_gt_arr, 'Before', 'Mean_GT_Rank', ...
    'NewVariableNames', 'Num_GT');

% GT Coverage - K=50 metrics (SEGMENT DATA, SAME COLUMN NAMES!)
seg_table = addvars(seg_table, seg_gt_expansion_50, 'Before', 'Num_GT', ...
    'NewVariableNames', 'ER@50_GT');
seg_table = addvars(seg_table, seg_gt_recall_50, 'Before', 'ER@50_GT', ...
    'NewVariableNames', 'R@50_GT');

% GT Coverage - K=10 metrics (SEGMENT DATA, SAME COLUMN NAMES!)
seg_table = addvars(seg_table, seg_gt_expansion_10, 'Before', 'R@50_GT', ...
    'NewVariableNames', 'ER@10_GT');
seg_table = addvars(seg_table, seg_gt_recall_10, 'Before', 'ER@10_GT', ...
    'NewVariableNames', 'R@10_GT');

% DTW vs Embedding Coverage - K=50 metrics (SEGMENT DATA, SAME COLUMN NAMES!)
seg_table = addvars(seg_table, seg_expansion_50, 'Before', 'R@10_GT', ...
    'NewVariableNames', 'ER@50_all');
seg_table = addvars(seg_table, seg_recall_cov_50, 'Before', 'ER@50_all', ...
    'NewVariableNames', 'R@50_all');

% DTW vs Embedding Coverage - K=10 metrics (SEGMENT DATA, SAME COLUMN NAMES!)
seg_table = addvars(seg_table, seg_expansion_10, 'Before', 'R@50_all', ...
    'NewVariableNames', 'ER@10_all');
seg_table = addvars(seg_table, seg_recall_cov_10, 'Before', 'ER@10_all', ...
    'NewVariableNames', 'R@10_all');

% Configuration columns (identical to trajectory table)
seg_table = addvars(seg_table, w_meta, 'Before', 'R@10_all', ...
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

% Level identifier
level_seg = repmat({'Segment'}, total_experiments, 1);
seg_table = addvars(seg_table, level_seg, 'Before', 'Embedding_Config', ...
    'NewVariableNames', 'Level');

% ========================================================================
%%  COMBINE TABLES
%  ========================================================================

combined_table = [traj_table; seg_table];

% ========================================================================
%%  ADD CONFIGURATION COLUMNS
%  ========================================================================
timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd''T''HHmmss'));

fprintf('=== Adding Configuration Columns ===\n');

% Create config arrays (same for all rows)
num_total_rows = height(combined_table);

timestamp_arr = repmat({timestamp}, num_total_rows, 1);
database_size_arr = repmat(base_config.database_sample_size, num_total_rows, 1);
top_k_arr = repmat(base_config.top_k_trajectories, num_total_rows, 1);
lb_kim_ratio_arr = repmat(dtw_config.lb_kim_keep_ratio, num_total_rows, 1);
lb_keogh_candidates_arr = repmat(dtw_config.lb_keogh_candidates, num_total_rows, 1);
dtw_window_arr = repmat(dtw_config.cdtw_window, num_total_rows, 1);
dtw_normalization_arr = repmat(dtw_config.normalize_dtw, num_total_rows, 1);
dtw_rotation_align_arr = repmat(dtw_config.use_rotation_alignment, num_total_rows, 1);
ground_truth_arr = repmat(base_config.has_ground_truth, num_total_rows, 1);

% Add columns at the BEGINNING (before 'Level')
combined_table = addvars(combined_table, ground_truth_arr, 'Before', 'Level', ...
    'NewVariableNames', 'Ground_Truth');
combined_table = addvars(combined_table, dtw_rotation_align_arr, 'Before', 'Ground_Truth', ...
    'NewVariableNames', 'DTW_RotationAlign');
combined_table = addvars(combined_table, dtw_normalization_arr, 'Before', 'DTW_RotationAlign', ...
    'NewVariableNames', 'DTW_Normalization');
combined_table = addvars(combined_table, dtw_window_arr, 'Before', 'DTW_Normalization', ...
    'NewVariableNames', 'DTW_Window');
combined_table = addvars(combined_table, lb_keogh_candidates_arr, 'Before', 'DTW_Window', ...
    'NewVariableNames', 'LB_Keogh_Candidates');
combined_table = addvars(combined_table, lb_kim_ratio_arr, 'Before', 'LB_Keogh_Candidates', ...
    'NewVariableNames', 'LB_Kim_Ratio');
combined_table = addvars(combined_table, top_k_arr, 'Before', 'LB_Kim_Ratio', ...
    'NewVariableNames', 'Top_K');
combined_table = addvars(combined_table, database_size_arr, 'Before', 'Top_K', ...
    'NewVariableNames', 'Database_Size');
combined_table = addvars(combined_table, timestamp_arr, 'Before', 'Database_Size', ...
    'NewVariableNames', 'Timestamp');

fprintf('✓ Added 9 configuration columns\n\n');

% ========================================================================
%%  SAVE RESULTS
%  ========================================================================

output_dir = fullfile(pwd, 'results');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

output_file = fullfile(output_dir, sprintf('embedding_validation_%s.csv', timestamp));
writetable(combined_table, output_file, 'WriteRowNames', false);

fprintf('✓ Results saved to: %s\n\n', output_file);

% ========================================================================
%%  ANALYSIS & SUMMARY
%  ========================================================================

fprintf('========================================\n');
fprintf('EXPERIMENT COMPLETED\n');
fprintf('========================================\n');
fprintf('Pre-loading time: %.1f minutes\n', preload_time/60);
fprintf('DTW pre-computation time: %.1f minutes\n', dtw_time/60);
fprintf('Embeddings pre-computation time: %.1f minutes\n', emb_time/60);
fprintf('Experiments runtime: %.1f minutes (%.2f hours)\n', total_time/60, total_time/3600);
fprintf('Total runtime: %.2f hours\n', (preload_time + dtw_time + emb_time + total_time)/3600);
fprintf('Total Experiments: %d\n', total_experiments);
fprintf('Average time per experiment: %.1f seconds (weights-only!)\n', total_time/total_experiments);
fprintf('Results file: %s\n\n', output_file);

% ========================================================================
%%  SAVE METADATA (CONFIGURATION PARAMETERS)
%  ========================================================================

fprintf('=== Saving Experiment Metadata ===\n');

% Build metadata struct
metadata = struct();

% General Configuration
metadata.timestamp = timestamp;
metadata.total_experiments = total_experiments;
metadata.database_sample_size = base_config.database_sample_size;
metadata.random_seed = base_config.random_seed;
metadata.top_k_trajectories = base_config.top_k_trajectories;

% DTW Configuration
metadata.lb_kim_keep_ratio = dtw_config.lb_kim_keep_ratio;
metadata.lb_keogh_candidates = dtw_config.lb_keogh_candidates;
metadata.cdtw_window = dtw_config.cdtw_window;
metadata.normalize_dtw = dtw_config.normalize_dtw;
metadata.use_rotation_alignment = dtw_config.use_rotation_alignment;

% Data Loading Configuration
metadata.chunk_size = chunk_size;
metadata.schema = schema;

% Timing
metadata.preload_time_minutes = preload_time / 60;
metadata.dtw_precompute_time_minutes = dtw_time / 60;
metadata.emb_precompute_time_minutes = emb_time / 60;
metadata.experiments_time_minutes = total_time / 60;
metadata.total_runtime_hours = (preload_time + dtw_time + emb_time + total_time) / 3600;
metadata.avg_time_per_experiment_seconds = total_time / total_experiments;

% Dimensions
metadata.num_embedding_configs = num_embeddings;
metadata.num_queries = num_queries;
metadata.num_weight_modes = num_weight_modes;

% Embedding Configurations (as comma-separated strings)
metadata.embedding_configs = strjoin(embedding_configs(:,1), '; ');

% Query IDs (as comma-separated string)
metadata.query_ids = strjoin(query_ids, '; ');

% Weight Modes (as comma-separated string)
metadata.weight_mode_names = strjoin(weight_mode_configs(:,1), '; ');

% Memory Usage
metadata.data_cache_size_MB = whos('data_cache').bytes / 1e6;
metadata.dtw_cache_size_MB = whos('dtw_cache').bytes / 1e6;
metadata.embeddings_cache_size_MB = whos('embeddings_cache').bytes / 1e6;

% Convert struct to table (transpose for vertical layout)
metadata_table = struct2table(metadata);

% Save as CSV
metadata_file = strrep(output_file, '.csv', '_metadata.csv');
writetable(metadata_table, metadata_file, 'WriteRowNames', false);

fprintf('✓ Metadata saved to: %s\n\n', metadata_file);

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  🎉 EXPERIMENT COMPLETED SUCCESSFULLY!                         ║\n');
fprintf('║                                                                ║\n');
fprintf('║  Total speedup with embeddings pre-computation:                ║\n');
fprintf('║  Before: ~4-5 hours (with embedding computation each time)     ║\n');
fprintf('║  After:  ~1 hour (embeddings computed once!)                   ║\n');
fprintf('║  Speedup: 4-5× FASTER! 🚀                                      ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');
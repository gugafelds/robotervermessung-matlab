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
base_config.database_sample_size = 1000;  % Fixed for fair comparison
base_config.random_seed = 42;
base_config.top_k_trajectories = 100;     % Fixed

% === DIMENSION 1: Embedding Architectures ===
embedding_configs = {
    % Name,                  n_coarse, n_fine, use_multi_scale
    'Single-Fine-75',        0,        75,     false;
    'Multi-Balanced-100',    25,       75,     true;
    'Single-Fine-150',       0,        150,    false;
    'Multi-Dense-200',       50,       150,    true;
};

% === DIMENSION 2: Query Trajectories ===
query_ids = {
    '1764766034'; % np = 2 / p = 843
    %'1764766001'; % np = 2 / p = 861
    %'1764765958'; % np = 2 / p = 1077
    %'1764765635'; % np = 2 / p = 1278
    %'1764765776'; % np = 2 / p = 1851
    %'1763567277'; % np = 3 / p = 1188
    %'1763567148'; % np = 3 / p = 1482  
    %'1763567026'; % np = 3 / p = 1434
    %'1764763889'; % np = 3 / p = 1392
    %'1764763510'; % np = 3 / p = 1554
    %'1764762584'; % np = 4 / p = 2034
    %'1764763238'; % np = 4 / p = 1341
    %'1764762971'; % np = 4 / p = 1353
    %'1764762831'; % np = 4 / p = 1326
    %'1764762655'; % np = 4 / p = 1539
    %'1764765476'; % np = 5 / p = 1398
    %'1764764632'; % np = 5 / p = 2430
    %'1764764821'; % np = 5 / p = 1944
    %'1764765396'; % np = 5 / p = 1692
    %'1764765286'; % np = 5 / p = 1377
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

if use_ground_truth && ~isempty(ground_truth_ids)
    query_id = query_ids{1};  % Erste Query
    query_field = sprintf('q_%s', strrep(query_id, '-', '_'));
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
dtw_config.lb_kim_keep_ratio = 0.5;
dtw_config.lb_keogh_candidates = 100;
dtw_config.cdtw_window = 0.05;
dtw_config.normalize_dtw = true;
dtw_config.use_rotation_alignment = false;
dtw_config.ground_truth_map = ground_truth_map;

dtw_cache = precomputeDTW(data_cache, query_ids, dtw_config);

dtw_time = toc(dtw_tic);

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  DTW PRE-COMPUTATION COMPLETED                                 ║\n');
fprintf('║  Time: %.1f minutes                                            ║\n', dtw_time/60);
fprintf('║  Memory: %.1f MB                                               ║\n', whos('dtw_cache').bytes / 1e6);
fprintf('║  All experiments will skip DTW computation!                    ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

% ========================================================================
%%  COMPUTE DTW GROUND TRUTH METRICS (BASELINE)
%  ========================================================================

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  COMPUTING DTW GT METRICS (BASELINE)                          ║\n');
fprintf('║  This establishes DTW performance for GT retrieval            ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

if use_ground_truth && ~isempty(ground_truth_ids)
    % Storage for DTW GT metrics (one per query × mode)
    dtw_gt_metrics = struct();
    
    for q_idx = 1:num_queries
        query_id = query_ids{q_idx};
        query_field = sprintf('q_%s', strrep(query_id, '-', '_'));
        
        % Get GT for this query
        if ~isfield(ground_truth_map, query_field)
            continue;
        end
        
        gt_ids_for_query = ground_truth_map.(query_field).trajectories;
        num_gt = length(gt_ids_for_query);
        
        fprintf('Query %s: %d GT trajectories\n', query_id, num_gt);
        
        % For each DTW mode
        dtw_modes_to_compute = {'position', 'joint_states'};
        
        for mode_idx = 1:length(dtw_modes_to_compute)
            dtw_mode = dtw_modes_to_compute{mode_idx};
            
            % Get DTW ranking for this query/mode
            if ~isfield(dtw_cache.(query_field), dtw_mode)
                continue;
            end
            
            % ============================================================
            % TRAJECTORY-LEVEL GT METRICS
            % ============================================================
            dtw_ranking = dtw_cache.(query_field).(dtw_mode).trajectory_ranking;
            
            % Find GT ranks in DTW ranking
            gt_ranks = zeros(num_gt, 1);
            gt_found = 0;
            
            for gt_idx = 1:num_gt
                gt_id = gt_ids_for_query{gt_idx};
                rank_idx = find(strcmp(dtw_ranking.bahn_id, gt_id), 1);
                
                if ~isempty(rank_idx)
                    gt_ranks(gt_idx) = rank_idx;
                    gt_found = gt_found + 1;
                else
                    gt_ranks(gt_idx) = inf;  % Not found
                end
            end
            
            % Calculate trajectory metrics
            valid_ranks = gt_ranks(gt_ranks < inf);
            
            if isempty(valid_ranks)
                % No GT found
                p_gt_dtw = inf;
                mean_gt_rank_dtw = inf;
                r50_gt_dtw = 0;
                r10_gt_dtw = 0;
                r5_gt_dtw = 0;
                r3_gt_dtw = 0;
                r1_gt_dtw = 0;
            else
                % P_GT: Maximum rank needed to cover all GT
                p_gt_dtw = max(valid_ranks);
                
                % Mean GT Rank
                mean_gt_rank_dtw = mean(valid_ranks);
                
                % R@K_GT: Fraction of GT in Top-K (CORRECTED!)
                % Normalize by min(K, num_gt) not by num_gt!
                r50_gt_dtw = sum(valid_ranks <= 50) / min(50, num_gt);
                r10_gt_dtw = sum(valid_ranks <= 10) / min(10, num_gt);
                r5_gt_dtw = sum(valid_ranks <= 5) / min(5, num_gt);
                r3_gt_dtw = sum(valid_ranks <= 3) / min(3, num_gt);
                r1_gt_dtw = sum(valid_ranks <= 1) / min(1, num_gt);
            end
            
            % ============================================================
            % SEGMENT-LEVEL GT METRICS
            % ============================================================
            segment_rankings = dtw_cache.(query_field).(dtw_mode).segment_rankings;
            num_segments = length(segment_rankings);
            
            % Aggregate segment ranks for GT trajectories
            seg_gt_ranks = [];
            
            for seg_idx = 1:num_segments
                if isempty(segment_rankings{seg_idx})
                    continue;
                end
                
                seg_ranking = segment_rankings{seg_idx};
                
                for gt_idx = 1:num_gt
                    gt_id = gt_ids_for_query{gt_idx};
                    
                    % Find where GT appears in this segment ranking
                    seg_gt_mask = strcmp(seg_ranking.bahn_id, gt_id);
                    
                    if any(seg_gt_mask)
                        seg_ranks = find(seg_gt_mask);
                        seg_gt_ranks = [seg_gt_ranks; seg_ranks];
                    end
                end
            end
            
            % Calculate segment metrics
            if isempty(seg_gt_ranks)
                % No GT found in segments
                seg_p_gt_dtw = inf;
                seg_mean_gt_rank_dtw = inf;
                seg_r50_gt_dtw = 0;
                seg_r10_gt_dtw = 0;
                seg_r5_gt_dtw = 0;
                seg_r3_gt_dtw = 0;
                seg_r1_gt_dtw = 0;
            else
                valid_seg_ranks = seg_gt_ranks(seg_gt_ranks < inf);
                
                if isempty(valid_seg_ranks)
                    seg_p_gt_dtw = inf;
                    seg_mean_gt_rank_dtw = inf;
                    seg_r50_gt_dtw = 0;
                    seg_r10_gt_dtw = 0;
                    seg_r5_gt_dtw = 0;
                    seg_r3_gt_dtw = 0;
                    seg_r1_gt_dtw = 0;
                else
                    % For segments: use average metrics across all segment instances
                    seg_p_gt_dtw = max(valid_seg_ranks);
                    seg_mean_gt_rank_dtw = mean(valid_seg_ranks);
                    
                    % Count unique GT-segment pairs in Top-K
                    total_gt_seg_pairs = num_gt * num_segments;
                    
                    % R@K_GT: (CORRECTED!)
                    % Normalize by min(K, total_gt_seg_pairs)
                    seg_r50_gt_dtw = sum(valid_seg_ranks <= 50) / (num_segments * min(50, num_gt));
                    seg_r10_gt_dtw = sum(valid_seg_ranks <= 10) / (num_segments * min(10, num_gt));
                    seg_r5_gt_dtw = sum(valid_seg_ranks <= 5) / (num_segments * min(5, num_gt));
                    seg_r3_gt_dtw = sum(valid_seg_ranks <= 3) / (num_segments * min(3, num_gt));
                    seg_r1_gt_dtw = sum(valid_seg_ranks <= 1) / (num_segments * min(1, num_gt));
                end
            end
            
            % ============================================================
            % Store all metrics
            % ============================================================
            dtw_gt_metrics.(query_field).(dtw_mode) = struct();
            
            % Trajectory metrics
            dtw_gt_metrics.(query_field).(dtw_mode).p_gt = p_gt_dtw;
            dtw_gt_metrics.(query_field).(dtw_mode).mean_gt_rank = mean_gt_rank_dtw;
            dtw_gt_metrics.(query_field).(dtw_mode).r50_gt = r50_gt_dtw;
            dtw_gt_metrics.(query_field).(dtw_mode).r10_gt = r10_gt_dtw;
            dtw_gt_metrics.(query_field).(dtw_mode).r5_gt = r5_gt_dtw;
            dtw_gt_metrics.(query_field).(dtw_mode).r3_gt = r3_gt_dtw;
            dtw_gt_metrics.(query_field).(dtw_mode).r1_gt = r1_gt_dtw;
            
            % Segment metrics
            dtw_gt_metrics.(query_field).(dtw_mode).seg_p_gt = seg_p_gt_dtw;
            dtw_gt_metrics.(query_field).(dtw_mode).seg_mean_gt_rank = seg_mean_gt_rank_dtw;
            dtw_gt_metrics.(query_field).(dtw_mode).seg_r50_gt = seg_r50_gt_dtw;
            dtw_gt_metrics.(query_field).(dtw_mode).seg_r10_gt = seg_r10_gt_dtw;
            dtw_gt_metrics.(query_field).(dtw_mode).seg_r5_gt = seg_r5_gt_dtw;
            dtw_gt_metrics.(query_field).(dtw_mode).seg_r3_gt = seg_r3_gt_dtw;
            dtw_gt_metrics.(query_field).(dtw_mode).seg_r1_gt = seg_r1_gt_dtw;
            
            % Common
            dtw_gt_metrics.(query_field).(dtw_mode).num_gt = num_gt;
            
            fprintf('  %s (Trajectory) - P_GT: %.0f, Mean_Rank: %.1f, R@50: %.2f, R@10: %.2f, R@5: %.2f, R@3: %.2f, R@1: %.2f\n', ...
                dtw_mode, p_gt_dtw, mean_gt_rank_dtw, r50_gt_dtw, r10_gt_dtw, r5_gt_dtw, r3_gt_dtw, r1_gt_dtw);
            fprintf('  %s (Segment)    - P_GT: %.0f, Mean_Rank: %.1f, R@50: %.2f, R@10: %.2f, R@5: %.2f, R@3: %.2f, R@1: %.2f\n', ...
                dtw_mode, seg_p_gt_dtw, seg_mean_gt_rank_dtw, seg_r50_gt_dtw, seg_r10_gt_dtw, seg_r5_gt_dtw, seg_r3_gt_dtw, seg_r1_gt_dtw);
        end
        
        fprintf('\n');
    end
    
    % Store in base_config for later use
    base_config.dtw_gt_metrics = dtw_gt_metrics;
    
    fprintf('✓ DTW GT metrics computed for all queries\n\n');
else
    fprintf('⚠ No ground truth - skipping DTW GT metrics\n\n');
end

% Close database connection (no longer needed!)
close(conn);
fprintf('✓ Database connection closed (no longer needed)\n\n');


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

% Preallocate matrices (MORE COLUMNS NOW!)
traj_results_matrix = zeros(total_experiments, 7);  % Spearman + 6 R@K metrics
seg_results_matrix = zeros(total_experiments, 7);   % Same for segments

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

% DTW vs Embedding Coverage metrics (trajectory level)
p_dtw_vs_eb_traj = zeros(total_experiments, 1);
recall_dtw_vs_eb_50 = zeros(total_experiments, 1);
recall_dtw_vs_eb_10 = zeros(total_experiments, 1);
recall_dtw_vs_eb_5 = zeros(total_experiments, 1);
recall_dtw_vs_eb_3 = zeros(total_experiments, 1);
recall_dtw_vs_eb_1 = zeros(total_experiments, 1);

% DTW vs Embedding Coverage metrics (segment level)
p_dtw_vs_eb_seg = zeros(total_experiments, 1);
seg_recall_dtw_vs_eb_50 = zeros(total_experiments, 1);
seg_recall_dtw_vs_eb_10 = zeros(total_experiments, 1);
seg_recall_dtw_vs_eb_5 = zeros(total_experiments, 1);
seg_recall_dtw_vs_eb_3 = zeros(total_experiments, 1);
seg_recall_dtw_vs_eb_1 = zeros(total_experiments, 1);

% GT vs Embedding metrics (trajectory level)
p_gt_vs_eb_traj = zeros(total_experiments, 1);
mean_gt_vs_eb_rank_traj = zeros(total_experiments, 1);
recall_gt_vs_eb_50_traj = zeros(total_experiments, 1);
recall_gt_vs_eb_10_traj = zeros(total_experiments, 1);
recall_gt_vs_eb_5_traj = zeros(total_experiments, 1);
recall_gt_vs_eb_3_traj = zeros(total_experiments, 1);
recall_gt_vs_eb_1_traj = zeros(total_experiments, 1);

% GT vs Embedding metrics (segment level)
p_gt_vs_eb_seg = zeros(total_experiments, 1);
mean_gt_vs_eb_rank_seg = zeros(total_experiments, 1);
recall_gt_vs_eb_50_seg = zeros(total_experiments, 1);
recall_gt_vs_eb_10_seg = zeros(total_experiments, 1);
recall_gt_vs_eb_5_seg = zeros(total_experiments, 1);
recall_gt_vs_eb_3_seg = zeros(total_experiments, 1);
recall_gt_vs_eb_1_seg = zeros(total_experiments, 1);

% GT vs DTW metrics (trajectory level) - NEW!
p_gt_vs_dtw_traj = zeros(total_experiments, 1);
mean_gt_vs_dtw_rank_traj = zeros(total_experiments, 1);
recall_gt_vs_dtw_50_traj = zeros(total_experiments, 1);
recall_gt_vs_dtw_10_traj = zeros(total_experiments, 1);
recall_gt_vs_dtw_5_traj = zeros(total_experiments, 1);
recall_gt_vs_dtw_3_traj = zeros(total_experiments, 1);
recall_gt_vs_dtw_1_traj = zeros(total_experiments, 1);

% GT vs DTW metrics (segment level) - NEW!
p_gt_vs_dtw_seg = zeros(total_experiments, 1);
mean_gt_vs_dtw_rank_seg = zeros(total_experiments, 1);
recall_gt_vs_dtw_50_seg = zeros(total_experiments, 1);
recall_gt_vs_dtw_10_seg = zeros(total_experiments, 1);
recall_gt_vs_dtw_5_seg = zeros(total_experiments, 1);
recall_gt_vs_dtw_3_seg = zeros(total_experiments, 1);
recall_gt_vs_dtw_1_seg = zeros(total_experiments, 1);

% Num_GT (same for all levels)
num_gt_arr = zeros(total_experiments, 1);

% Extract data from results
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
    multi_scale_arr(i) = (r.n_coarse > 0);
    runtimes(i) = r.exp_runtime / 60;
    
    % Weights
    w_pos(i) = r.weight_pos;
    w_joint(i) = r.weight_joint;
    w_orient(i) = r.weight_orient;
    w_vel(i) = r.weight_vel;
    w_meta(i) = r.weight_meta;
    
    % Trajectory-level: Spearman + DTW vs EB Coverage
    traj_results_matrix(i, :) = [
        r.spearman, r.p_at_k, r.p_at_10, r.p_at_5, r.p_at_3, r.p_at_1, r.p_at_1
    ];
    
    % Segment-level: Spearman + DTW vs EB Coverage
    seg_results_matrix(i, :) = [
        r.seg_spearman, r.seg_p_at_k, r.seg_p_at_10, ...
        r.seg_p_at_5, r.seg_p_at_3, r.seg_p_at_1, r.seg_p_at_1
    ];
    
    % ====================================================================
    % DTW vs Embedding Coverage (trajectory)
    % ====================================================================
    if isfield(r, 'p_all_traj')
        p_dtw_vs_eb_traj(i) = r.p_all_traj;
        recall_dtw_vs_eb_50(i) = r.recall_at_50;
        recall_dtw_vs_eb_10(i) = r.recall_at_10;
        recall_dtw_vs_eb_5(i) = r.recall_at_5;
        recall_dtw_vs_eb_3(i) = r.recall_at_3;
        recall_dtw_vs_eb_1(i) = r.recall_at_1;
    else
        p_dtw_vs_eb_traj(i) = NaN;
        recall_dtw_vs_eb_50(i) = NaN;
        recall_dtw_vs_eb_10(i) = NaN;
        recall_dtw_vs_eb_5(i) = NaN;
        recall_dtw_vs_eb_3(i) = NaN;
        recall_dtw_vs_eb_1(i) = NaN;
    end
    
    % DTW vs Embedding Coverage (segment)
    if isfield(r, 'p_all_seg')
        p_dtw_vs_eb_seg(i) = r.p_all_seg;
        seg_recall_dtw_vs_eb_50(i) = r.seg_recall_at_50;
        seg_recall_dtw_vs_eb_10(i) = r.seg_recall_at_10;
        seg_recall_dtw_vs_eb_5(i) = r.seg_recall_at_5;
        seg_recall_dtw_vs_eb_3(i) = r.seg_recall_at_3;
        seg_recall_dtw_vs_eb_1(i) = r.seg_recall_at_1;
    else
        p_dtw_vs_eb_seg(i) = NaN;
        seg_recall_dtw_vs_eb_50(i) = NaN;
        seg_recall_dtw_vs_eb_10(i) = NaN;
        seg_recall_dtw_vs_eb_5(i) = NaN;
        seg_recall_dtw_vs_eb_3(i) = NaN;
        seg_recall_dtw_vs_eb_1(i) = NaN;
    end
    
    % ====================================================================
    % GT vs Embedding (trajectory)
    % ====================================================================
    if isfield(r, 'p_gt_traj')
        p_gt_vs_eb_traj(i) = r.p_gt_traj;
        mean_gt_vs_eb_rank_traj(i) = r.mean_gt_rank;
        recall_gt_vs_eb_50_traj(i) = r.gt_recall_at_50;
        recall_gt_vs_eb_10_traj(i) = r.gt_recall_at_10;
        recall_gt_vs_eb_5_traj(i) = r.gt_recall_at_5;
        recall_gt_vs_eb_3_traj(i) = r.gt_recall_at_3;
        recall_gt_vs_eb_1_traj(i) = r.gt_recall_at_1;
    else
        p_gt_vs_eb_traj(i) = NaN;
        mean_gt_vs_eb_rank_traj(i) = NaN;
        recall_gt_vs_eb_50_traj(i) = NaN;
        recall_gt_vs_eb_10_traj(i) = NaN;
        recall_gt_vs_eb_5_traj(i) = NaN;
        recall_gt_vs_eb_3_traj(i) = NaN;
        recall_gt_vs_eb_1_traj(i) = NaN;
    end
    
    % GT vs Embedding (segment)
    if isfield(r, 'p_gt_seg')
        p_gt_vs_eb_seg(i) = r.p_gt_seg;
        mean_gt_vs_eb_rank_seg(i) = r.seg_mean_gt_rank;
        recall_gt_vs_eb_50_seg(i) = r.seg_gt_recall_at_50;
        recall_gt_vs_eb_10_seg(i) = r.seg_gt_recall_at_10;
        recall_gt_vs_eb_5_seg(i) = r.seg_gt_recall_at_5;
        recall_gt_vs_eb_3_seg(i) = r.seg_gt_recall_at_3;
        recall_gt_vs_eb_1_seg(i) = r.seg_gt_recall_at_1;
    else
        p_gt_vs_eb_seg(i) = NaN;
        mean_gt_vs_eb_rank_seg(i) = NaN;
        recall_gt_vs_eb_50_seg(i) = NaN;
        recall_gt_vs_eb_10_seg(i) = NaN;
        recall_gt_vs_eb_5_seg(i) = NaN;
        recall_gt_vs_eb_3_seg(i) = NaN;
        recall_gt_vs_eb_1_seg(i) = NaN;
    end
    
    % ====================================================================
    % GT vs DTW (from base_config.dtw_gt_metrics) - NEW!
    % ====================================================================
    if isfield(base_config, 'dtw_gt_metrics')
        query_field = sprintf('q_%s', strrep(r.query_id, '-', '_'));
        dtw_mode = r.dtw_mode;
        
        if isfield(base_config.dtw_gt_metrics, query_field) && ...
           isfield(base_config.dtw_gt_metrics.(query_field), dtw_mode)
            
            dtw_gt = base_config.dtw_gt_metrics.(query_field).(dtw_mode);
            
            % Trajectory level
            p_gt_vs_dtw_traj(i) = dtw_gt.p_gt;
            mean_gt_vs_dtw_rank_traj(i) = dtw_gt.mean_gt_rank;
            recall_gt_vs_dtw_50_traj(i) = dtw_gt.r50_gt;
            recall_gt_vs_dtw_10_traj(i) = dtw_gt.r10_gt;
            recall_gt_vs_dtw_5_traj(i) = dtw_gt.r5_gt;
            recall_gt_vs_dtw_3_traj(i) = dtw_gt.r3_gt;
            recall_gt_vs_dtw_1_traj(i) = dtw_gt.r1_gt;
            
            % Segment level (same as trajectory for DTW)
            p_gt_vs_dtw_seg(i) = dtw_gt.seg_p_gt;
            mean_gt_vs_dtw_rank_seg(i) = dtw_gt.seg_mean_gt_rank;
            recall_gt_vs_dtw_50_seg(i) = dtw_gt.seg_r50_gt;
            recall_gt_vs_dtw_10_seg(i) = dtw_gt.seg_r10_gt;
            recall_gt_vs_dtw_5_seg(i) = dtw_gt.seg_r5_gt;
            recall_gt_vs_dtw_3_seg(i) = dtw_gt.seg_r3_gt;
            recall_gt_vs_dtw_1_seg(i) = dtw_gt.seg_r1_gt;
        else
            % No DTW GT metrics available
            p_gt_vs_dtw_traj(i) = NaN;
            mean_gt_vs_dtw_rank_traj(i) = NaN;
            recall_gt_vs_dtw_50_traj(i) = NaN;
            recall_gt_vs_dtw_10_traj(i) = NaN;
            recall_gt_vs_dtw_5_traj(i) = NaN;
            recall_gt_vs_dtw_3_traj(i) = NaN;
            recall_gt_vs_dtw_1_traj(i) = NaN;
            
            p_gt_vs_dtw_seg(i) = NaN;
            mean_gt_vs_dtw_rank_seg(i) = NaN;
            recall_gt_vs_dtw_50_seg(i) = NaN;
            recall_gt_vs_dtw_10_seg(i) = NaN;
            recall_gt_vs_dtw_5_seg(i) = NaN;
            recall_gt_vs_dtw_3_seg(i) = NaN;
            recall_gt_vs_dtw_1_seg(i) = NaN;
        end
    else
        % No DTW GT metrics computed
        p_gt_vs_dtw_traj(i) = NaN;
        mean_gt_vs_dtw_rank_traj(i) = NaN;
        recall_gt_vs_dtw_50_traj(i) = NaN;
        recall_gt_vs_dtw_10_traj(i) = NaN;
        recall_gt_vs_dtw_5_traj(i) = NaN;
        recall_gt_vs_dtw_3_traj(i) = NaN;
        recall_gt_vs_dtw_1_traj(i) = NaN;
        
        p_gt_vs_dtw_seg(i) = NaN;
        mean_gt_vs_dtw_rank_seg(i) = NaN;
        recall_gt_vs_dtw_50_seg(i) = NaN;
        recall_gt_vs_dtw_10_seg(i) = NaN;
        recall_gt_vs_dtw_5_seg(i) = NaN;
        recall_gt_vs_dtw_3_seg(i) = NaN;
        recall_gt_vs_dtw_1_seg(i) = NaN;
    end
    
    % Num_GT (same for all)
    if isfield(r, 'num_gt')
        num_gt_arr(i) = r.num_gt;
    else
        num_gt_arr(i) = NaN;
    end
end

% ========================================================================
%%  BUILD TRAJECTORY-LEVEL TABLE
%  ========================================================================

traj_table = array2table(traj_results_matrix, ...
    'VariableNames', {'Spearman_DTWvsEB', ...
                      sprintf('R@K_DTWvsEB'), ...
                      'R@10_DTWvsEB', 'R@5_DTWvsEB', 'R@3_DTWvsEB', 'R@1_DTWvsEB', 'R@1_DTWvsEB_dup'});

% Remove duplicate R@1 column (was added twice in matrix)
traj_table = removevars(traj_table, 'R@1_DTWvsEB_dup');

% Add R@50_DTWvsEB and P_DTWvsEB
traj_table = addvars(traj_table, recall_dtw_vs_eb_50, 'After', sprintf('R@K_DTWvsEB'), ...
    'NewVariableNames', 'R@50_DTWvsEB');
traj_table = addvars(traj_table, p_dtw_vs_eb_traj, 'After', 'R@1_DTWvsEB', ...
    'NewVariableNames', 'P_DTWvsEB');

% Add GT vs Embedding metrics
traj_table = addvars(traj_table, recall_gt_vs_eb_50_traj, 'After', 'P_DTWvsEB', ...
    'NewVariableNames', 'R@50_GTvsEB');
traj_table = addvars(traj_table, recall_gt_vs_eb_10_traj, 'After', 'R@50_GTvsEB', ...
    'NewVariableNames', 'R@10_GTvsEB');
traj_table = addvars(traj_table, recall_gt_vs_eb_5_traj, 'After', 'R@10_GTvsEB', ...
    'NewVariableNames', 'R@5_GTvsEB');
traj_table = addvars(traj_table, recall_gt_vs_eb_3_traj, 'After', 'R@5_GTvsEB', ...
    'NewVariableNames', 'R@3_GTvsEB');
traj_table = addvars(traj_table, recall_gt_vs_eb_1_traj, 'After', 'R@3_GTvsEB', ...
    'NewVariableNames', 'R@1_GTvsEB');
traj_table = addvars(traj_table, p_gt_vs_eb_traj, 'After', 'R@1_GTvsEB', ...
    'NewVariableNames', 'P_GTvsEB');
traj_table = addvars(traj_table, mean_gt_vs_eb_rank_traj, 'After', 'P_GTvsEB', ...
    'NewVariableNames', 'Mean_GTvsEB_Rank');

% Add GT vs DTW metrics
traj_table = addvars(traj_table, recall_gt_vs_dtw_50_traj, 'After', 'Mean_GTvsEB_Rank', ...
    'NewVariableNames', 'R@50_GTvsDTW');
traj_table = addvars(traj_table, recall_gt_vs_dtw_10_traj, 'After', 'R@50_GTvsDTW', ...
    'NewVariableNames', 'R@10_GTvsDTW');
traj_table = addvars(traj_table, recall_gt_vs_dtw_5_traj, 'After', 'R@10_GTvsDTW', ...
    'NewVariableNames', 'R@5_GTvsDTW');
traj_table = addvars(traj_table, recall_gt_vs_dtw_3_traj, 'After', 'R@5_GTvsDTW', ...
    'NewVariableNames', 'R@3_GTvsDTW');
traj_table = addvars(traj_table, recall_gt_vs_dtw_1_traj, 'After', 'R@3_GTvsDTW', ...
    'NewVariableNames', 'R@1_GTvsDTW');
traj_table = addvars(traj_table, p_gt_vs_dtw_traj, 'After', 'R@1_GTvsDTW', ...
    'NewVariableNames', 'P_GTvsDTW');
traj_table = addvars(traj_table, mean_gt_vs_dtw_rank_traj, 'After', 'P_GTvsDTW', ...
    'NewVariableNames', 'Mean_GTvsDTW_Rank');

% Add Num_GT
traj_table = addvars(traj_table, num_gt_arr, 'After', 'Mean_GTvsDTW_Rank', ...
    'NewVariableNames', 'Num_GT');

% Add configuration columns BEFORE Spearman_DTWvsEB (from right to left)
traj_table = addvars(traj_table, dtw_modes, 'Before', 'Spearman_DTWvsEB', ...
    'NewVariableNames', 'DTW_Mode');
traj_table = addvars(traj_table, total_dims_arr, 'Before', 'DTW_Mode', ...
    'NewVariableNames', 'Total_Dims');
traj_table = addvars(traj_table, n_fine_arr, 'Before', 'Total_Dims', ...
    'NewVariableNames', 'N_Fine');
traj_table = addvars(traj_table, n_coarse_arr, 'Before', 'N_Fine', ...
    'NewVariableNames', 'N_Coarse');
traj_table = addvars(traj_table, weight_modes, 'Before', 'N_Coarse', ...
    'NewVariableNames', 'Weight_Mode');
traj_table = addvars(traj_table, embedding_names, 'Before', 'Weight_Mode', ...
    'NewVariableNames', 'Embedding_Config');

% Add level identifier
level_traj = repmat({'Trajectory'}, total_experiments, 1);
traj_table = addvars(traj_table, level_traj, 'Before', 'Embedding_Config', ...
    'NewVariableNames', 'Level');

% Add Query_Bahn_ID
traj_table = addvars(traj_table, query_bahn_ids, 'Before', 'Level', ...
    'NewVariableNames', 'Query_Bahn_ID');

% ========================================================================
%%  BUILD SEGMENT-LEVEL TABLE
%  ========================================================================

seg_table = array2table(seg_results_matrix, ...
    'VariableNames', {'Spearman_DTWvsEB', ...
                      sprintf('R@K_DTWvsEB'), ...
                      'R@10_DTWvsEB', 'R@5_DTWvsEB', 'R@3_DTWvsEB', 'R@1_DTWvsEB', 'R@1_DTWvsEB_dup'});

% Remove duplicate
seg_table = removevars(seg_table, 'R@1_DTWvsEB_dup');

% Add R@50_DTWvsEB and P_DTWvsEB
seg_table = addvars(seg_table, seg_recall_dtw_vs_eb_50, 'After', sprintf('R@K_DTWvsEB'), ...
    'NewVariableNames', 'R@50_DTWvsEB');
seg_table = addvars(seg_table, p_dtw_vs_eb_seg, 'After', 'R@1_DTWvsEB', ...
    'NewVariableNames', 'P_DTWvsEB');

% Add GT vs Embedding metrics
seg_table = addvars(seg_table, recall_gt_vs_eb_50_seg, 'After', 'P_DTWvsEB', ...
    'NewVariableNames', 'R@50_GTvsEB');
seg_table = addvars(seg_table, recall_gt_vs_eb_10_seg, 'After', 'R@50_GTvsEB', ...
    'NewVariableNames', 'R@10_GTvsEB');
seg_table = addvars(seg_table, recall_gt_vs_eb_5_seg, 'After', 'R@10_GTvsEB', ...
    'NewVariableNames', 'R@5_GTvsEB');
seg_table = addvars(seg_table, recall_gt_vs_eb_3_seg, 'After', 'R@5_GTvsEB', ...
    'NewVariableNames', 'R@3_GTvsEB');
seg_table = addvars(seg_table, recall_gt_vs_eb_1_seg, 'After', 'R@3_GTvsEB', ...
    'NewVariableNames', 'R@1_GTvsEB');
seg_table = addvars(seg_table, p_gt_vs_eb_seg, 'After', 'R@1_GTvsEB', ...
    'NewVariableNames', 'P_GTvsEB');
seg_table = addvars(seg_table, mean_gt_vs_eb_rank_seg, 'After', 'P_GTvsEB', ...
    'NewVariableNames', 'Mean_GTvsEB_Rank');

% Add GT vs DTW metrics
seg_table = addvars(seg_table, recall_gt_vs_dtw_50_seg, 'After', 'Mean_GTvsEB_Rank', ...
    'NewVariableNames', 'R@50_GTvsDTW');
seg_table = addvars(seg_table, recall_gt_vs_dtw_10_seg, 'After', 'R@50_GTvsDTW', ...
    'NewVariableNames', 'R@10_GTvsDTW');
seg_table = addvars(seg_table, recall_gt_vs_dtw_5_seg, 'After', 'R@10_GTvsDTW', ...
    'NewVariableNames', 'R@5_GTvsDTW');
seg_table = addvars(seg_table, recall_gt_vs_dtw_3_seg, 'After', 'R@5_GTvsDTW', ...
    'NewVariableNames', 'R@3_GTvsDTW');
seg_table = addvars(seg_table, recall_gt_vs_dtw_1_seg, 'After', 'R@3_GTvsDTW', ...
    'NewVariableNames', 'R@1_GTvsDTW');
seg_table = addvars(seg_table, p_gt_vs_dtw_seg, 'After', 'R@1_GTvsDTW', ...
    'NewVariableNames', 'P_GTvsDTW');
seg_table = addvars(seg_table, mean_gt_vs_dtw_rank_seg, 'After', 'P_GTvsDTW', ...
    'NewVariableNames', 'Mean_GTvsDTW_Rank');

% Add Num_GT
seg_table = addvars(seg_table, num_gt_arr, 'After', 'Mean_GTvsDTW_Rank', ...
    'NewVariableNames', 'Num_GT');

% Add configuration columns
seg_table = addvars(seg_table, dtw_modes, 'Before', 'Spearman_DTWvsEB', ...
    'NewVariableNames', 'DTW_Mode');
seg_table = addvars(seg_table, total_dims_arr, 'Before', 'DTW_Mode', ...
    'NewVariableNames', 'Total_Dims');
seg_table = addvars(seg_table, n_fine_arr, 'Before', 'Total_Dims', ...
    'NewVariableNames', 'N_Fine');
seg_table = addvars(seg_table, n_coarse_arr, 'Before', 'N_Fine', ...
    'NewVariableNames', 'N_Coarse');
seg_table = addvars(seg_table, weight_modes, 'Before', 'N_Coarse', ...
    'NewVariableNames', 'Weight_Mode');
seg_table = addvars(seg_table, embedding_names, 'Before', 'Weight_Mode', ...
    'NewVariableNames', 'Embedding_Config');

% Add level identifier
level_seg = repmat({'Segment'}, total_experiments, 1);
seg_table = addvars(seg_table, level_seg, 'Before', 'Embedding_Config', ...
    'NewVariableNames', 'Level');

% Add Query_Bahn_ID
seg_table = addvars(seg_table, query_bahn_ids, 'Before', 'Level', ...
    'NewVariableNames', 'Query_Bahn_ID');

% ========================================================================
%%  COMBINE TABLES
%  ========================================================================

combined_table = [traj_table; seg_table];

fprintf('✓ Tables created with new column structure\n\n');

% ========================================================================
%%  ADD CONFIGURATION COLUMNS
%  ========================================================================
timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd''T''HHmmss'));

fprintf('=== Adding Configuration Columns ===\n');

% Create config arrays (same for all rows)
num_total_rows = height(combined_table);

timestamp_arr = repmat({timestamp}, num_total_rows, 1);
database_size_arr = repmat(base_config.database_sample_size, num_total_rows, 1);

actual_num_trajs = length(candidate_ids);
actual_num_segs = height(data_cache.segments.metadata);  

num_trajs_arr = repmat(actual_num_trajs, num_total_rows, 1);
num_segs_arr = repmat(actual_num_segs, num_total_rows, 1);

top_k_arr = repmat(base_config.top_k_trajectories, num_total_rows, 1);
lb_kim_ratio_arr = repmat(dtw_config.lb_kim_keep_ratio, num_total_rows, 1);
lb_keogh_candidates_arr = repmat(dtw_config.lb_keogh_candidates, num_total_rows, 1);
dtw_normalization_arr = repmat(dtw_config.normalize_dtw, num_total_rows, 1);
dtw_rotation_align_arr = repmat(dtw_config.use_rotation_alignment, num_total_rows, 1);

% Add columns BEFORE Query_Bahn_ID (in reverse order to get correct sequence)
combined_table = addvars(combined_table, top_k_arr, 'Before', 'Query_Bahn_ID', ...
    'NewVariableNames', 'Top_K');
combined_table = addvars(combined_table, num_segs_arr, 'Before', 'Top_K', ...
    'NewVariableNames', 'Num_Segments');
combined_table = addvars(combined_table, num_trajs_arr, 'Before', 'Num_Segments', ...
    'NewVariableNames', 'Num_Trajectories');
combined_table = addvars(combined_table, database_size_arr, 'Before', 'Num_Trajectories', ...
    'NewVariableNames', 'Database_Size');
combined_table = addvars(combined_table, lb_keogh_candidates_arr, 'Before', 'Database_Size', ...
    'NewVariableNames', 'LB_Keogh_Candidates');
combined_table = addvars(combined_table, lb_kim_ratio_arr, 'Before', 'LB_Keogh_Candidates', ...
    'NewVariableNames', 'LB_Kim_Ratio');
combined_table = addvars(combined_table, dtw_rotation_align_arr, 'Before', 'LB_Kim_Ratio', ...
    'NewVariableNames', 'DTW_RotationAlign');
combined_table = addvars(combined_table, dtw_normalization_arr, 'Before', 'DTW_RotationAlign', ...
    'NewVariableNames', 'DTW_Normalization');
combined_table = addvars(combined_table, timestamp_arr, 'Before', 'DTW_Normalization', ...
    'NewVariableNames', 'Timestamp');

% Remove columns that are not in the desired order
cols_to_remove = {'Runtime_min', 'W_Pos', 'W_Joint', 'W_Orient', 'W_Vel', 'W_Meta', ...
                  'Multi_Scale', 'DTW_Window', 'Ground_Truth'};
for i = 1:length(cols_to_remove)
    if ismember(cols_to_remove{i}, combined_table.Properties.VariableNames)
        combined_table = removevars(combined_table, cols_to_remove{i});
    end
end

fprintf('✓ Added configuration columns and removed unused columns\n\n');

% ========================================================================
%%  SAVE RESULTS
%  ========================================================================

% Find robotervermessung-matlab root directory
current_path = mfilename('fullpath');
path_parts = strsplit(current_path, filesep);

% Find the index of 'robotervermessung-matlab' in path
matlab_project_idx = find(contains(path_parts, 'robotervermessung-matlab'), 1, 'last');

if ~isempty(matlab_project_idx)
    % Build path to project root
    % FIX: Don't use filesep as first argument on Windows!
    if ispc  % Windows
        project_root = fullfile(path_parts{1:matlab_project_idx});
    else  % Linux/Mac
        project_root = fullfile(filesep, path_parts{1:matlab_project_idx});
    end
    output_dir = fullfile(project_root, 'similarity');
else
    % Fallback: use current directory
    warning('Could not find robotervermessung-matlab in path, using current directory');
    output_dir = fullfile(pwd, 'similarity');
end

% Create directory if it doesn't exist
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
    fprintf('Created directory: %s\n', output_dir);
end

fprintf('Saving results to: %s\n', output_dir);

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

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  🎉 EXPERIMENT COMPLETED SUCCESSFULLY!                         ║\n');
fprintf('║                                                                ║\n');
fprintf('║  Total speedup with embeddings pre-computation:                ║\n');
fprintf('║  Before: ~4-5 hours (with embedding computation each time)     ║\n');
fprintf('║  After:  ~1 hour (embeddings computed once!)                   ║\n');
fprintf('║  Speedup: 4-5× FASTER! 🚀                                      ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

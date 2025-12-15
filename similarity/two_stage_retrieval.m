%  TWO-STAGE RETRIEVAL EVALUATION
%  ========================================================================
%  Goal: Evaluate embedding + DTW reranking efficiency vs. accuracy
%  
%  Pipeline:
%  Stage 1: Fast embedding-based retrieval → Top-K candidates
%  Stage 2: Precise DTW reranking on K candidates
%
%  Tests:
%  - K candidates: [10, 50, 100, 200, 500]
%  - Database sizes: [100, 500, 1000, 2000]
%  - 6 Query trajectories
%  - Top 2 embedding configs (from previous analysis)
%
%  Metrics:
%  - Recall@K (1, 5, 10, 50)
%  - NDCG@10, MRR
%  - Latency (Stage 1, Stage 2, Total)
%  - DTW calls saved
%  - Speedup vs. brute-force (DB ≤ 1000)
%
%  Total Experiments: 2 configs × 5 K-values × 4 DB-sizes × 6 queries = 240
%  Runtime: ~2-3 hours (depending on K and DB size)
%  Output: results/two_stage_retrieval_YYYY-MM-DD.csv
%  ========================================================================

clear; clc;

% Add paths
addpath(genpath(pwd));
addpath(genpath('../main'));
addpath(genpath('../lasertracker'));
addpath(genpath('../methods'));
addpath(genpath('results'));

fprintf('✓ Paths added successfully\n\n');

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║                                                                ║\n');
fprintf('║  TWO-STAGE RETRIEVAL EVALUATION                                ║\n');
fprintf('║  Embedding + DTW Reranking                                     ║\n');
fprintf('║                                                                ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

% ========================================================================
%% CONFIGURATION
% ========================================================================

% === Best Configurations (from previous analysis) ===
% TODO: Fill in your top 2 embedding configs from Figure 1
embedding_configs = {
    'Multi-Balanced-25',      5,  20,   true;   % Best overall
    'Single-Coarse-15',     0,    15,   false;  % Best single-scale
};

% Best weight mode (from Figure 2)
weight_mode_configs = {
    % Motion (Joint States)
    'Joint only',           'joint_states',  [0, 1, 0, 0, 0];  % Baseline
    'Joint + Pos',         'joint_states',  [1, 1, 0, 0, 0];  % Best (from Figure 2)
    
    % Space (Cartesian/Position)
    'Position only',        'position',      [1, 0, 0, 0, 0];  % Baseline
    'Pos + Joint',            'position',      [1, 1, 0, 0, 0];  % Best (from Figure 2)
};

% === Two-Stage Parameters ===
%k_candidates = [10, 50, 100, 200, 500];  % Stage 1 output sizes
%database_sizes = [100, 500, 1000, 2000]; % Database sizes to test
k_candidates = [50, 100, 200];  % Stage 1 output sizes
database_sizes = [5000]; % Database sizes to test
k_final = 10;                             % Final Top-K (fixed)

% === Query Trajectories ===
query_ids = {
    '1764766034'; % Good query
    '1765473159'; % Noisy
    '1765473166'; % Noisy
    '1765473097'; % Noisy
    '1765473008'; % Noisy
    '1765472956'; % Noisy
};

% === Base Configuration ===
base_config = struct();
base_config.random_seed = 42;
base_config.top_k_trajectories = max(k_candidates);  % Buffer for DTW !!!!!!!!!!!!!!!!!!!!!!!!! WICHTIG!!!!
use_ground_truth = false;

% === Calculate Total Experiments ===
num_embeddings = size(embedding_configs, 1);        % 2
num_weight_modes = size(weight_mode_configs, 1);    % 4
num_k_values = length(k_candidates);                % 5
num_db_sizes = length(database_sizes);              % 4
num_queries = length(query_ids);                    % 6

total_experiments = num_embeddings * num_weight_modes * num_k_values * num_db_sizes * num_queries;

fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('EXPERIMENT CONFIGURATION\n');
fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('  Embedding configs:  %d\n', num_embeddings);
fprintf('  Weight modes:       %d (Baseline + Best per DTW mode)\n', num_weight_modes);
fprintf('  K candidates:       %s\n', mat2str(k_candidates));
fprintf('  Database sizes:     %s\n', mat2str(database_sizes));
fprintf('  Queries:            %d\n', num_queries);
fprintf('  Total experiments:  %d\n', total_experiments);
fprintf('═══════════════════════════════════════════════════════════════\n\n');

% Print weight modes
fprintf('Weight Modes:\n');
for i = 1:size(weight_mode_configs, 1)
    fprintf('  %s (%s): %s\n', ...
        weight_mode_configs{i, 1}, ...
        weight_mode_configs{i, 2}, ...
        mat2str(weight_mode_configs{i, 3}));
end
fprintf('\n');

% ========================================================================
%% SECTION 2: DATABASE CONNECTION & SAMPLING
% ========================================================================

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  SECTION 2: DATABASE CONNECTION & SAMPLING                     ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

% === Connect to Database ===
fprintf('=== Connecting to Database ===\n');
conn = connectingToPostgres;

if isopen(conn)
    fprintf('✓ Database connection successful\n\n');
else
    error('✗ Database connection failed');
end

schema = 'bewegungsdaten';

% === Load Full Database Metadata ===
fprintf('=== Loading Database Metadata ===\n');
full_db_query = sprintf(...
    ['SELECT bahn_id FROM robotervermessung.%s.bahn_metadata ' ...
     'WHERE bahn_id = segment_id'], schema);

full_db_metadata = fetch(conn, full_db_query);
num_total_trajectories = height(full_db_metadata);

fprintf('  Total trajectories in database: %d\n', num_total_trajectories);
fprintf('  Maximum sample size needed: %d\n', max(database_sizes));
fprintf('  Random seed: %d\n\n', base_config.random_seed);

% === Sample Maximum Database Size ===
% We sample the maximum size once, then subsample for smaller DB sizes
max_db_size = max(database_sizes);
rng(base_config.random_seed);

sample_indices = randperm(num_total_trajectories, max_db_size);
sampled_metadata = full_db_metadata(sample_indices, :);
max_candidate_ids = sampled_metadata.bahn_id;

fprintf('✓ Sampled %d trajectories (maximum needed)\n', length(max_candidate_ids));
fprintf('  Smaller DB sizes will be subsampled from this pool\n\n');

% ========================================================================
%% GROUND TRUTH HANDLING
% ========================================================================

if use_ground_truth
    fprintf('=== Ground Truth Handling ===\n');
    
    % Get ground truth for all queries
    [ground_truth_ids, ground_truth_map] = getGTCandidates(conn, schema, query_ids);
    
    if ~isempty(ground_truth_ids)
        % Remove GT from sampled pool to avoid duplicates
        max_candidate_ids = setdiff(max_candidate_ids, ground_truth_ids);
        
        % Add GT to pool
        max_candidate_ids = [max_candidate_ids; ground_truth_ids];
        
        fprintf('  Ground truth trajectories: %d\n', length(ground_truth_ids));
        fprintf('  Random candidates: %d\n', length(max_candidate_ids) - length(ground_truth_ids));
        fprintf('  Total candidate pool: %d\n\n', length(max_candidate_ids));
        
        % Store for later use
        base_config.ground_truth_map = ground_truth_map;
        base_config.has_ground_truth = true;
        base_config.ground_truth_ids = ground_truth_ids;
    else
        fprintf('⚠ No ground truth found - continuing without GT evaluation\n\n');
        base_config.has_ground_truth = false;
    end
else
    fprintf('=== Ground Truth Disabled ===\n');
    fprintf('  Using only random samples\n\n');
    base_config.has_ground_truth = false;
end

% ========================================================================
%% CREATE SUBSAMPLED CANDIDATE POOLS FOR EACH DB SIZE
% ========================================================================

fprintf('=== Creating Subsampled Candidate Pools ===\n');

% Storage for candidate pools
candidate_pools = struct();

for i = 1:length(database_sizes)
    db_size = database_sizes(i);
    
    if db_size == max_db_size
        % Use full pool
        candidate_pools.(sprintf('db_%d', db_size)) = max_candidate_ids;
    else
        % Subsample from max pool (excluding GT)
        if base_config.has_ground_truth
            non_gt_candidates = setdiff(max_candidate_ids, ground_truth_ids);
            subsample_size = db_size - length(ground_truth_ids);
            
            % Take first N from non-GT candidates (deterministic)
            subsampled = non_gt_candidates(1:subsample_size);
            
            % Add GT back
            candidate_pools.(sprintf('db_%d', db_size)) = [subsampled; ground_truth_ids];
        else
            % Just take first N
            candidate_pools.(sprintf('db_%d', db_size)) = max_candidate_ids(1:db_size);
        end
    end
    
    fprintf('  DB size %4d: %d candidates', db_size, ...
        length(candidate_pools.(sprintf('db_%d', db_size))));
    
    if base_config.has_ground_truth
        fprintf(' (including %d GT)\n', length(ground_truth_ids));
    else
        fprintf('\n');
    end
end

fprintf('\n✓ All candidate pools created\n\n');

% Store in base_config
base_config.candidate_pools = candidate_pools;
base_config.max_candidate_ids = max_candidate_ids;

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  SECTION 2 COMPLETE                                            ║\n');
fprintf('║  Ready for Section 3: Data Preloading                          ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

% ========================================================================
%% SECTION 3: DATA PRELOADING
% ========================================================================

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  SECTION 3: DATA PRELOADING                                    ║\n');
fprintf('║  Loading all data ONCE for maximum DB size                     ║\n');
fprintf('║  Estimated time: 15-20 minutes                                 ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

preload_tic = tic;

% === Load Data for Maximum DB Size ===
chunk_size = 50;

fprintf('Loading data for %d trajectories...\n', length(max_candidate_ids));
data_cache = loadDataExperiment(conn, schema, max_candidate_ids, query_ids, chunk_size);

preload_time = toc(preload_tic);

fprintf('\n╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  DATA PRELOADING COMPLETED                                     ║\n');
fprintf('║  Time: %.1f minutes                                            ║\n', preload_time/60);
fprintf('║  Memory: %.1f MB                                               ║\n', whos('data_cache').bytes / 1e6);
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

% ========================================================================
%% CREATE SUBSAMPLED DATA CACHES FOR EACH DB SIZE
% ========================================================================

fprintf('=== Creating Subsampled Data Caches ===\n');

data_caches = struct();

for i = 1:length(database_sizes)
    db_size = database_sizes(i);
    pool_name = sprintf('db_%d', db_size);
    candidate_ids = candidate_pools.(pool_name);
    
    fprintf('  Creating cache for DB size %d...\n', db_size);
    
    % Filter data_cache to only include these candidates
    data_caches.(pool_name) = filterDataCache(data_cache, candidate_ids);
    
    fprintf('    ✓ Candidates: %d\n', length(data_caches.(pool_name).candidates.bahn_ids));
    fprintf('    ✓ Segments: %d\n', length(data_caches.(pool_name).segments.segment_ids));
end

fprintf('\n✓ All data caches created\n\n');

% Close database connection (no longer needed)
close(conn);
fprintf('✓ Database connection closed\n\n');

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  SECTION 3 COMPLETE                                            ║\n');
fprintf('║  Ready for Section 4: DTW Precomputation                       ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

% ========================================================================
%% SECTION 4: DTW PRECOMPUTATION (BASELINE)
% ========================================================================

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  SECTION 4: DTW PRECOMPUTATION (BASELINE)                      ║\n');
fprintf('║  Computing full DTW rankings for all DB sizes                  ║\n');
fprintf('║  Using LB_Kim + LB_Keogh for efficiency                        ║\n');
fprintf('║  Estimated time: 20-30 minutes (depends on DB sizes)           ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

dtw_tic = tic;

% === DTW Configuration ===
dtw_config = struct();
dtw_config.top_k_trajectories = base_config.top_k_trajectories;
dtw_config.lb_kim_keep_ratio = 0.9;      % Keep 80% after LB_Kim
dtw_config.lb_keogh_candidates = 400;     % Further reduce to 400 with LB_Keogh
dtw_config.cdtw_window = 0.2;
dtw_config.normalize_dtw = false;
dtw_config.use_rotation_alignment = false;

fprintf('DTW Configuration:\n');
fprintf('  LB_Kim keep ratio:      %.1f%% (filters to ~%d for DB=5000)\n', ...
    dtw_config.lb_kim_keep_ratio * 100, ...
    round(max(database_sizes) * dtw_config.lb_kim_keep_ratio));
fprintf('  LB_Keogh candidates:    %d\n', dtw_config.lb_keogh_candidates);
fprintf('  CDTW window:            %.1f\n', dtw_config.cdtw_window);
fprintf('  Normalization:          %d\n', dtw_config.normalize_dtw);
fprintf('  Rotation alignment:     %d\n\n', dtw_config.use_rotation_alignment);

% === Storage for DTW Baselines ===
dtw_baselines = struct();

% === Compute DTW ONCE per DB Size (computes all modes internally) ===
total_dtw_computations = length(database_sizes);
dtw_counter = 0;

for db_idx = 1:length(database_sizes)
    dtw_counter = dtw_counter + 1;
    db_size = database_sizes(db_idx);
    pool_name = sprintf('db_%d', db_size);
    
    fprintf('\n═══════════════════════════════════════════════════════════════\n');
    fprintf('[%d/%d] Computing DTW for DB size %d (all modes)\n', ...
        dtw_counter, total_dtw_computations, db_size);
    fprintf('═══════════════════════════════════════════════════════════════\n\n');
    
    % Get data cache for this DB size
    curr_data_cache = data_caches.(pool_name);
    
    % Compute DTW (computes position AND joint_states internally!)
    dtw_start = tic;
    dtw_cache_curr = precomputeDTW(curr_data_cache, query_ids, dtw_config);
    dtw_time_total = toc(dtw_start);
    
    fprintf('  ✓ Completed in %.1f minutes\n', dtw_time_total / 60);
    fprintf('  Memory: %.1f MB\n\n', whos('dtw_cache_curr').bytes / 1e6);
    
    % Extract times per mode from cache
    unique_dtw_modes = {'position', 'joint_states'};
    
    for mode_idx = 1:length(unique_dtw_modes)
        dtw_mode = unique_dtw_modes{mode_idx};
        
        baseline_key = sprintf('%s_%s', dtw_mode, pool_name);
        
        % Store cache
        dtw_baselines.(baseline_key) = dtw_cache_curr;
        
        % Extract time for this mode (average over all queries)
        total_time_for_mode = 0;
        
        for q_idx = 1:num_queries
            query_id = query_ids{q_idx};
            query_field = sprintf('q_%s', strrep(query_id, '-', '_'));
            
            % Get time from cache
            if isfield(dtw_cache_curr, query_field) && ...
               isfield(dtw_cache_curr.(query_field), dtw_mode) && ...
               isfield(dtw_cache_curr.(query_field).(dtw_mode), 'dtw_time')
                
                total_time_for_mode = total_time_for_mode + ...
                    dtw_cache_curr.(query_field).(dtw_mode).dtw_time;
            end
        end
        
        % Average time per query for this mode
        dtw_time_per_query = total_time_for_mode / num_queries;
        
        % Store timing
        time_key = [baseline_key '_time_per_query'];
        dtw_baselines.(time_key) = dtw_time_per_query;
        
        fprintf('    %s: %.2f sec/query\n', dtw_mode, dtw_time_per_query);
    end
    fprintf('\n');
end

dtw_total_time = toc(dtw_tic);

% After DTW computation, store as proxy GT:
if ~base_config.has_ground_truth
    fprintf('=== Using DTW Rankings as Proxy Ground Truth ===\n\n');
    
    % Create proxy GT map from DTW rankings
    proxy_gt_map = struct();
    
    for q_idx = 1:num_queries
        query_id = query_ids{q_idx};
        query_field = sprintf('q_%s', strrep(query_id, '-', '_'));
        
        % For each DTW mode
        for mode_idx = 1:length(unique_dtw_modes)
            dtw_mode = unique_dtw_modes{mode_idx};
            
            if isfield(dtw_cache_curr, query_field) && ...
               isfield(dtw_cache_curr.(query_field), dtw_mode)
                
                % === TRAJECTORY RANKING ===
                dtw_traj_ranking = dtw_cache_curr.(query_field).(dtw_mode).trajectory_ranking;
                
                % Take Top-50 trajectories as "proxy GT"
                proxy_k_traj = min(50, height(dtw_traj_ranking));
                proxy_gt_traj_ids = dtw_traj_ranking.bahn_id(1:proxy_k_traj);
                
                % === SEGMENT RANKINGS (per query segment) ===
                dtw_seg_rankings = dtw_cache_curr.(query_field).(dtw_mode).segment_rankings;
                num_query_segments = length(dtw_seg_rankings);
                
                % Store Top-50 for EACH query segment
                proxy_gt_seg_ids = cell(num_query_segments, 1);
                
                for seg_idx = 1:num_query_segments
                    seg_ranking_table = dtw_seg_rankings{seg_idx};
                    
                    % Take Top-50 from this segment ranking
                    proxy_k_seg = min(50, height(seg_ranking_table));
                    proxy_gt_seg_ids{seg_idx} = seg_ranking_table.segment_id(1:proxy_k_seg);
                end
                
                % Store both
                proxy_field = sprintf('%s_%s', query_field, dtw_mode);
                proxy_gt_map.(proxy_field).trajectories = proxy_gt_traj_ids;
                proxy_gt_map.(proxy_field).segments = proxy_gt_seg_ids;  % Cell array!
                proxy_gt_map.(proxy_field).num_trajectories = length(proxy_gt_traj_ids);
                proxy_gt_map.(proxy_field).num_query_segments = num_query_segments;
                
                fprintf('  %s - %s: Top-%d trajectories, %d query segments (Top-50 each)\n', ...
                    query_id, dtw_mode, proxy_k_traj, num_query_segments);
            end
        end
    end
    
    % Store proxy GT
    base_config.proxy_gt_map = proxy_gt_map;
    base_config.has_proxy_gt = true;
    
    fprintf('\n✓ Proxy GT created from DTW Top-50 (trajectories + segment rankings)\n\n');
else
    base_config.proxy_gt_map = [];
    base_config.has_proxy_gt = false;
end

% ========================================================================
%% COMPUTE DTW GT METRICS (BASELINE PERFORMANCE)
% ========================================================================

if base_config.has_ground_truth || base_config.has_proxy_gt
    fprintf('╔════════════════════════════════════════════════════════════════╗\n');
    fprintf('║  COMPUTING DTW GT METRICS (BASELINE)                          ║\n');
    
    if base_config.has_ground_truth
        fprintf('║  Using REAL ground truth                                       ║\n');
    else
        fprintf('║  Using PROXY GT (DTW Top-50) - metrics will be near-perfect   ║\n');
    end
    
    fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');
    
    dtw_gt_metrics = struct();
    
    % Loop through unique DTW modes
    for dtw_mode_idx = 1:length(unique_dtw_modes)
        dtw_mode_curr = unique_dtw_modes{dtw_mode_idx};
        
        fprintf('%s:\n', dtw_mode_curr);
        fprintf('────────────────────────────────────────────────────────────────\n');
        
        for db_idx = 1:length(database_sizes)
            db_size = database_sizes(db_idx);
            pool_name = sprintf('db_%d', db_size);
            
            % Get DTW baseline by mode
            baseline_key = sprintf('%s_%s', dtw_mode_curr, pool_name);
            dtw_cache_curr = dtw_baselines.(baseline_key);
            
            fprintf('  DB=%d:\n', db_size);
            
            % Compute GT metrics for each query
            for q_idx = 1:num_queries
                query_id = query_ids{q_idx};
                query_field = sprintf('q_%s', strrep(query_id, '-', '_'));
                
                % ========================================================
                % GET GT IDS (Real or Proxy)
                % ========================================================
                if base_config.has_ground_truth && isfield(ground_truth_map, query_field)
                    % Use real GT
                    gt_ids_for_query = ground_truth_map.(query_field).trajectories;
                    num_gt = length(gt_ids_for_query);
                elseif base_config.has_proxy_gt
                    % Use proxy GT (DTW Top-50)
                    proxy_field = sprintf('%s_%s', query_field, dtw_mode_curr);
                    if isfield(base_config.proxy_gt_map, proxy_field)
                        gt_ids_for_query = base_config.proxy_gt_map.(proxy_field).trajectories;
                        num_gt = length(gt_ids_for_query);
                    else
                        continue;  % No proxy GT for this query/mode
                    end
                else
                    continue;  % No GT at all
                end
                
                % Get DTW ranking
                if ~isfield(dtw_cache_curr, query_field) || ...
                   ~isfield(dtw_cache_curr.(query_field), dtw_mode_curr)
                    continue;
                end
                
                % ============================================================
                % TRAJECTORY-LEVEL GT METRICS
                % ============================================================
                dtw_ranking = dtw_cache_curr.(query_field).(dtw_mode_curr).trajectory_ranking;
                
                % Calculate GT metrics (Trajectory level)
                gt_ranks = zeros(num_gt, 1);
                for gt_idx = 1:num_gt
                    gt_id = gt_ids_for_query{gt_idx};
                    rank_idx = find(strcmp(dtw_ranking.bahn_id, gt_id), 1);
                    
                    if ~isempty(rank_idx)
                        gt_ranks(gt_idx) = rank_idx;
                    else
                        gt_ranks(gt_idx) = inf;
                    end
                end
                
                valid_ranks = gt_ranks(gt_ranks < inf);
                
                if ~isempty(valid_ranks)
                    r10_gt_dtw = sum(valid_ranks <= 10) / min(10, num_gt);
                    r50_gt_dtw = sum(valid_ranks <= 50) / min(50, num_gt);
                    r5_gt_dtw = sum(valid_ranks <= 5) / min(5, num_gt);
                    r3_gt_dtw = sum(valid_ranks <= 3) / min(3, num_gt);
                    r1_gt_dtw = sum(valid_ranks <= 1) / min(1, num_gt);
                    mean_gt_rank_dtw = mean(valid_ranks);
                    p_gt_dtw = max(valid_ranks);
                else
                    r10_gt_dtw = 0;
                    r50_gt_dtw = 0;
                    r5_gt_dtw = 0;
                    r3_gt_dtw = 0;
                    r1_gt_dtw = 0;
                    mean_gt_rank_dtw = inf;
                    p_gt_dtw = inf;
                end
                
                % ============================================================
                % SEGMENT-LEVEL GT METRICS
                % ============================================================
                if isfield(dtw_cache_curr.(query_field).(dtw_mode_curr), 'segment_rankings')
                    segment_rankings = dtw_cache_curr.(query_field).(dtw_mode_curr).segment_rankings;
                    num_query_segments = length(segment_rankings);
                    
                    % Get SEGMENT GT
                    if base_config.has_ground_truth && isfield(ground_truth_map, query_field)
                        if isfield(ground_truth_map.(query_field), 'segments')
                            gt_segs_raw = ground_truth_map.(query_field).segments;
                            
                            % Check structure type
                            if isstruct(gt_segs_raw)
                                % Struct format: seg_1, seg_2, etc. → convert to cell array
                                seg_fields = fieldnames(gt_segs_raw);
                                gt_seg_ids_for_query = cell(length(seg_fields), 1);
                                
                                for sf_idx = 1:length(seg_fields)
                                    field_name = seg_fields{sf_idx};
                                    gt_seg_ids_for_query{sf_idx} = gt_segs_raw.(field_name);
                                    
                                    % Ensure it's a cell array
                                    if ~iscell(gt_seg_ids_for_query{sf_idx})
                                        gt_seg_ids_for_query{sf_idx} = {gt_seg_ids_for_query{sf_idx}};
                                    end
                                end
                                has_segment_gt = true;
                                
                            elseif iscell(gt_segs_raw) && ~isempty(gt_segs_raw) && iscell(gt_segs_raw{1})
                                % Already cell array per segment
                                gt_seg_ids_for_query = gt_segs_raw;
                                has_segment_gt = true;
                                
                            elseif iscell(gt_segs_raw)
                                % Flat cell list
                                gt_seg_ids_for_query = repmat({gt_segs_raw}, num_query_segments, 1);
                                has_segment_gt = true;
                                
                            else
                                % Unknown format
                                has_segment_gt = false;
                            end
                        else
                            has_segment_gt = false;
                        end
                        
                        % FALLBACK if no proper segment GT
                        if ~has_segment_gt
                            gt_traj_ids = ground_truth_map.(query_field).trajectories;
                            
                            gt_seg_ids_all = {};
                            for traj_idx = 1:length(gt_traj_ids)
                                traj_id = gt_traj_ids{traj_idx};
                                seg_mask = strcmp(data_cache_curr.segments.bahn_ids, traj_id);
                                gt_seg_ids_all = [gt_seg_ids_all; data_cache_curr.segments.segment_ids(seg_mask)];
                            end
                            
                            gt_seg_ids_for_query = repmat({gt_seg_ids_all}, num_query_segments, 1);
                            has_segment_gt = true;
                        end
                        
                    elseif base_config.has_proxy_gt
                        % Proxy GT segments
                        proxy_field = sprintf('%s_%s', query_field, dtw_mode_curr);
                        if isfield(base_config.proxy_gt_map, proxy_field)
                            gt_seg_ids_for_query = base_config.proxy_gt_map.(proxy_field).segments;
                            has_segment_gt = true;
                        else
                            gt_seg_ids_for_query = {};
                            has_segment_gt = false;
                        end
                    else
                        gt_seg_ids_for_query = {};
                        has_segment_gt = false;
                    end
                    
                    if ~has_segment_gt || isempty(gt_seg_ids_for_query)
                        % No segment GT at all
                        seg_r1_gt_dtw = NaN;
                        seg_r3_gt_dtw = NaN;
                        seg_r5_gt_dtw = NaN;
                        seg_r10_gt_dtw = NaN;
                        seg_r50_gt_dtw = NaN;
                        seg_mean_gt_rank_dtw = NaN;
                        seg_p_gt_dtw = NaN;
                    else
                        % Calculate GT coverage for EACH segment separately
                        seg_gt_recalls = struct();
                        seg_gt_recalls.r1 = [];
                        seg_gt_recalls.r3 = [];
                        seg_gt_recalls.r5 = [];
                        seg_gt_recalls.r10 = [];
                        seg_gt_recalls.r50 = [];
                        seg_gt_recalls.mean_rank = [];
                        seg_gt_recalls.p_gt = [];
                        
                        for seg_idx = 1:num_query_segments
                            if isempty(segment_rankings{seg_idx})
                                continue;
                            end
                            
                            seg_ranking = segment_rankings{seg_idx};
                            
                            % Get GT for THIS specific query segment
                            if seg_idx <= length(gt_seg_ids_for_query)
                                gt_segs_curr = gt_seg_ids_for_query{seg_idx};
                                num_gt_seg_curr = length(gt_segs_curr);
                            else
                                continue;
                            end
                            
                            if isempty(gt_segs_curr)
                                continue;
                            end
                            
                            % Find GT ranks by SEGMENT_ID
                            seg_gt_ranks = zeros(num_gt_seg_curr, 1);
                            
                            for gt_idx = 1:num_gt_seg_curr
                                gt_seg_id = gt_segs_curr{gt_idx};
                                rank_idx = find(strcmp(seg_ranking.segment_id, gt_seg_id), 1);
                                
                                if ~isempty(rank_idx)
                                    seg_gt_ranks(gt_idx) = rank_idx;
                                else
                                    seg_gt_ranks(gt_idx) = inf;
                                end
                            end
                            
                            valid_seg_ranks = seg_gt_ranks(seg_gt_ranks < inf);
                            
                            if ~isempty(valid_seg_ranks)
                                seg_gt_recalls.r1 = [seg_gt_recalls.r1; sum(valid_seg_ranks <= 1) / min(1, num_gt_seg_curr)];
                                seg_gt_recalls.r3 = [seg_gt_recalls.r3; sum(valid_seg_ranks <= 3) / min(3, num_gt_seg_curr)];
                                seg_gt_recalls.r5 = [seg_gt_recalls.r5; sum(valid_seg_ranks <= 5) / min(5, num_gt_seg_curr)];
                                seg_gt_recalls.r10 = [seg_gt_recalls.r10; sum(valid_seg_ranks <= 10) / min(10, num_gt_seg_curr)];
                                seg_gt_recalls.r50 = [seg_gt_recalls.r50; sum(valid_seg_ranks <= 50) / min(50, num_gt_seg_curr)];
                                seg_gt_recalls.mean_rank = [seg_gt_recalls.mean_rank; mean(valid_seg_ranks)];
                                seg_gt_recalls.p_gt = [seg_gt_recalls.p_gt; max(valid_seg_ranks)];
                            end
                        end
                        
                        % Average across segments
                        if ~isempty(seg_gt_recalls.r1)
                            seg_r1_gt_dtw = mean(seg_gt_recalls.r1);
                            seg_r3_gt_dtw = mean(seg_gt_recalls.r3);
                            seg_r5_gt_dtw = mean(seg_gt_recalls.r5);
                            seg_r10_gt_dtw = mean(seg_gt_recalls.r10);
                            seg_r50_gt_dtw = mean(seg_gt_recalls.r50);
                            seg_mean_gt_rank_dtw = mean(seg_gt_recalls.mean_rank);
                            seg_p_gt_dtw = mean(seg_gt_recalls.p_gt);
                        else
                            seg_r1_gt_dtw = 0;
                            seg_r3_gt_dtw = 0;
                            seg_r5_gt_dtw = 0;
                            seg_r10_gt_dtw = 0;
                            seg_r50_gt_dtw = 0;
                            seg_mean_gt_rank_dtw = inf;
                            seg_p_gt_dtw = inf;
                        end
                    end
                else
                    % No segment rankings
                    seg_r1_gt_dtw = NaN;
                    seg_r3_gt_dtw = NaN;
                    seg_r5_gt_dtw = NaN;
                    seg_r10_gt_dtw = NaN;
                    seg_r50_gt_dtw = NaN;
                    seg_mean_gt_rank_dtw = NaN;
                    seg_p_gt_dtw = NaN;
                end
                
                % Store metrics (key by dtw_mode, not weight_mode!)
                metric_key = sprintf('%s_%s_q%s', dtw_mode_curr, pool_name, query_id);
                
                dtw_gt_metrics.(metric_key) = struct();
                
                % Trajectory metrics
                dtw_gt_metrics.(metric_key).r1_gt = r1_gt_dtw;
                dtw_gt_metrics.(metric_key).r3_gt = r3_gt_dtw;
                dtw_gt_metrics.(metric_key).r5_gt = r5_gt_dtw;
                dtw_gt_metrics.(metric_key).r10_gt = r10_gt_dtw;
                dtw_gt_metrics.(metric_key).r50_gt = r50_gt_dtw;
                dtw_gt_metrics.(metric_key).mean_rank = mean_gt_rank_dtw;
                dtw_gt_metrics.(metric_key).p_gt = p_gt_dtw;
                
                % Segment metrics
                dtw_gt_metrics.(metric_key).seg_r1_gt = seg_r1_gt_dtw;
                dtw_gt_metrics.(metric_key).seg_r3_gt = seg_r3_gt_dtw;
                dtw_gt_metrics.(metric_key).seg_r5_gt = seg_r5_gt_dtw;
                dtw_gt_metrics.(metric_key).seg_r10_gt = seg_r10_gt_dtw;
                dtw_gt_metrics.(metric_key).seg_r50_gt = seg_r50_gt_dtw;
                dtw_gt_metrics.(metric_key).seg_mean_rank = seg_mean_gt_rank_dtw;
                dtw_gt_metrics.(metric_key).seg_p_gt = seg_p_gt_dtw;
                
                % Common
                dtw_gt_metrics.(metric_key).num_gt = num_gt;
                
                fprintf('    Query %s:\n', query_id);
                fprintf('      Trajectory: R@10=%.2f, R@50=%.2f, MeanRank=%.1f\n', ...
                    r10_gt_dtw, r50_gt_dtw, mean_gt_rank_dtw);
                fprintf('      Segment:    R@10=%.2f, R@50=%.2f, MeanRank=%.1f\n', ...
                    seg_r10_gt_dtw, seg_r50_gt_dtw, seg_mean_gt_rank_dtw);
            end
            fprintf('\n');
        end
    end
    
    % Store in base_config
    base_config.dtw_gt_metrics = dtw_gt_metrics;
    
    fprintf('✓ DTW GT metrics computed for all configurations\n');
    
    if ~base_config.has_ground_truth
        fprintf('  Note: Metrics use Proxy GT (DTW Top-50), so values are near-perfect\n');
    end
    fprintf('\n');
end


fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  DTW PRECOMPUTATION COMPLETED                                  ║\n');
fprintf('║  Total time: %.1f minutes (%.2f hours)                         ║\n', ...
    dtw_total_time/60, dtw_total_time/3600);
fprintf('║  Total memory: %.1f MB                                         ║\n', ...
    whos('dtw_baselines').bytes / 1e6);
fprintf('║  Configurations computed: %d                                   ║\n', ...
    total_dtw_computations);
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

% ========================================================================
%% SECTION 5: EMBEDDINGS PRECOMPUTATION
% ========================================================================

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  SECTION 5: EMBEDDINGS PRECOMPUTATION                          ║\n');
fprintf('║  Computing embeddings for all DB sizes & configs               ║\n');
fprintf('║  Estimated time: 15-20 minutes                                 ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

emb_tic = tic;

% === Embedding Configuration ===
emb_config = struct();
emb_config.norm_strategy = 'max_extent';

fprintf('Embedding Configuration:\n');
fprintf('  Normalization strategy: %s\n', emb_config.norm_strategy);
fprintf('  Embedding configs:      %d\n', num_embeddings);
fprintf('  DB sizes:               %d\n\n', length(database_sizes));

% === Storage for Embeddings ===
embeddings_caches = struct();

% === Compute Embeddings for Each DB Size ===
total_emb_computations = num_embeddings * length(database_sizes);
emb_counter = 0;

for emb_idx = 1:num_embeddings
    emb_name = embedding_configs{emb_idx, 1};
    n_coarse = embedding_configs{emb_idx, 2};
    n_fine = embedding_configs{emb_idx, 3};
    multi_scale = embedding_configs{emb_idx, 4};
    
    fprintf('\n═══════════════════════════════════════════════════════════════\n');
    fprintf('Embedding Config: %s\n', emb_name);
    fprintf('  n_coarse=%d, n_fine=%d, multi_scale=%d\n', n_coarse, n_fine, multi_scale);
    fprintf('═══════════════════════════════════════════════════════════════\n\n');
    
    for db_idx = 1:length(database_sizes)
        emb_counter = emb_counter + 1;
        db_size = database_sizes(db_idx);
        pool_name = sprintf('db_%d', db_size);
        
        fprintf('[%d/%d] Computing embeddings for DB size %d...\n', ...
            emb_counter, total_emb_computations, db_size);
        
        % Get data cache for this DB size
        curr_data_cache = data_caches.(pool_name);
        
        % Compute embeddings
        emb_start = tic;
        emb_cache_curr = precomputeEmbeddings(curr_data_cache, query_ids, ...
        embedding_configs(emb_idx, :), emb_config);
        emb_time_curr = toc(emb_start);
        
        % Store in embeddings structure
        cache_key = sprintf('%s_%s', matlab.lang.makeValidName(emb_name), pool_name);

        embeddings_caches.(cache_key) = emb_cache_curr;
        
        fprintf('  ✓ Completed in %.1f seconds\n', emb_time_curr);
        fprintf('  Memory: %.1f MB\n\n', whos('emb_cache_curr').bytes / 1e6);
    end
end

emb_total_time = toc(emb_tic);

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  EMBEDDINGS PRECOMPUTATION COMPLETED                           ║\n');
fprintf('║  Total time: %.1f minutes                                      ║\n', ...
    emb_total_time/60);
fprintf('║  Total memory: %.1f MB                                         ║\n', ...
    whos('embeddings_caches').bytes / 1e6);
fprintf('║  Configurations computed: %d                                   ║\n', ...
    total_emb_computations);
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  SECTION 5 COMPLETE                                            ║\n');
fprintf('║  Ready for Section 6: Two-Stage Retrieval Experiments          ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

% ========================================================================
%% SECTION 6: TWO-STAGE RETRIEVAL EXPERIMENTS (TRAJECTORY + SEGMENT)
% ========================================================================

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  SECTION 6: TWO-STAGE RETRIEVAL EXPERIMENTS                    ║\n');
fprintf('║  Computing for BOTH Trajectory and Segment levels              ║\n');
fprintf('║  Stage 1: Embedding retrieval → Top-K candidates               ║\n');
fprintf('║  Stage 2: DTW reranking on K candidates (with adaptive LB)     ║\n');
fprintf('║  Total experiments: %d × 2 levels = %d                         ║\n', ...
    total_experiments, total_experiments * 2);
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

% === Storage for Results (2x because trajectory + segment) ===
all_results = cell(total_experiments * 2, 1);  % Double for both levels!
experiment_counter = 0;

% === Start Experiments ===
experiment_start_time = tic;

fprintf('Starting experiments...\n\n');

% === Loop: Embedding Config ===
for emb_idx = 1:num_embeddings
    emb_name = embedding_configs{emb_idx, 1};
    n_coarse = embedding_configs{emb_idx, 2};
    n_fine = embedding_configs{emb_idx, 3};
    multi_scale = embedding_configs{emb_idx, 4};
    total_dims = n_coarse + n_fine;
    
    % === Loop: Weight Mode ===
    for wm_idx = 1:num_weight_modes
        weight_mode_name = weight_mode_configs{wm_idx, 1};
        dtw_mode_curr = weight_mode_configs{wm_idx, 2};
        weights = weight_mode_configs{wm_idx, 3};
        
        % === Loop: K Candidates ===
        for k_idx = 1:num_k_values
            K = k_candidates(k_idx);
            
            % === Loop: DB Size ===
            for db_idx = 1:length(database_sizes)
                db_size = database_sizes(db_idx);
                pool_name = sprintf('db_%d', db_size);
                
                % === Loop: Query ===
                for q_idx = 1:num_queries
                    query_id = query_ids{q_idx};
                    query_field = sprintf('q_%s', strrep(query_id, '-', '_'));
                    
                    % ========================================================
                    % SHARED SETUP FOR BOTH LEVELS
                    % ========================================================
                    
                    % Get precomputed caches
                    emb_cache_key = sprintf('%s_%s', matlab.lang.makeValidName(emb_name), pool_name);
                    emb_cache = embeddings_caches.(emb_cache_key);
                    
                    dtw_baseline_key = sprintf('%s_%s', matlab.lang.makeValidName(dtw_mode_curr), pool_name);
                    dtw_baseline = dtw_baselines.(dtw_baseline_key);
                    
                    data_cache_curr = data_caches.(pool_name);
                    
                    % Get baseline DTW time
                    dtw_baseline_time_key = sprintf('%s_db_%d_time_per_query', dtw_mode_curr, db_size);
                    if isfield(dtw_baselines, dtw_baseline_time_key)
                        time_dtw_baseline = dtw_baselines.(dtw_baseline_time_key);
                    else
                        time_dtw_baseline = db_size * 0.001;
                    end
                    
                    % Check query exists
                    if ~isfield(emb_cache, query_field)
                        warning('Query %s not found - skipping', query_id);
                        continue;
                    end
                    
                    query_emb_data = emb_cache.(query_field);
                    emb_field_name = matlab.lang.makeValidName(emb_name);
                    
                    % Get query data
                    query_data = data_cache_curr.queries.(query_field);
                    
                    if strcmp(dtw_mode_curr, 'position')
                        query_seq = query_data.position;
                    else
                        query_seq = query_data.joint;
                    end
                    
                    % Get GT info (shared)
                    if base_config.has_ground_truth && isfield(ground_truth_map, query_field)
                        gt_traj_ids = ground_truth_map.(query_field).trajectories;
                        num_gt_traj = length(gt_traj_ids);
                        gt_type = 'real';
                    elseif base_config.has_proxy_gt
                        proxy_field = sprintf('%s_%s', query_field, dtw_mode_curr);
                        if isfield(base_config.proxy_gt_map, proxy_field)
                            gt_traj_ids = base_config.proxy_gt_map.(proxy_field).trajectories;
                            num_gt_traj = length(gt_traj_ids);
                            gt_type = 'proxy_dtw';
                        else
                            gt_traj_ids = {};
                            num_gt_traj = 0;
                            gt_type = 'none';
                        end
                    else
                        gt_traj_ids = {};
                        num_gt_traj = 0;
                        gt_type = 'none';
                    end
                    
                    % ====================================================
                    % LEVEL 1: TRAJECTORY-LEVEL TWO-STAGE
                    % ====================================================
                    
                    experiment_counter = experiment_counter + 1;
                    
                    
                    fprintf('\n╔════════════════════════════════════════════════════════════════╗\n');
                    fprintf('║  EXPERIMENT %4d/%4d (TRAJECTORY LEVEL)                      ║\n', ...
                        experiment_counter, total_experiments * 2);
                    fprintf('╠════════════════════════════════════════════════════════════════╣\n');
                    fprintf('║  Embedding:  %-48s ║\n', emb_name);
                    fprintf('║  Weight:     %-48s ║\n', weight_mode_name);
                    fprintf('║  K:          %-48d ║\n', K);
                    fprintf('║  DB Size:    %-48d ║\n', db_size);
                    fprintf('║  Query:      %-48s ║\n', query_id);
                    fprintf('╚════════════════════════════════════════════════════════════════╝\n');
                    
                    
                    % STAGE 1: Trajectory Embedding Retrieval
                    tic_stage1_traj = tic;
                    
                    query_embeddings_traj = query_emb_data.(emb_field_name).query_embeddings;
                    candidate_embeddings_traj = query_emb_data.(emb_field_name).candidate_embeddings;
                    
                    [top_k_ids_traj, ~, ~] = performRRFFusion(...
                        query_embeddings_traj, candidate_embeddings_traj, weights, K, ...
                        data_cache_curr.candidates.bahn_ids);
                    
                    time_stage1_traj = toc(tic_stage1_traj);
                    
                    % Embedding-Only Metrics (Trajectory)
                    [emb_only_metrics_traj] = computeEmbeddingOnlyMetrics(...
                        top_k_ids_traj, gt_traj_ids, num_gt_traj, K);
                    
                    % STAGE 2: Trajectory DTW Reranking
                    [dtw_reranked_ids_traj, time_stage2_traj, num_dtw_calls_traj] = ...
                        performDTWReranking(top_k_ids_traj, query_seq, data_cache_curr, ...
                        dtw_mode_curr, dtw_config, K);
                    
                    % Two-Stage Metrics (Trajectory)
                    [twostage_metrics_traj] = computeTwoStageMetrics(...
                        dtw_reranked_ids_traj, gt_traj_ids, num_gt_traj, K);
                    
                    % Get Baseline DTW GT Metrics (Trajectory)
                    metric_key = sprintf('%s_%s_q%s', dtw_mode_curr, pool_name, query_id);
                    if isfield(base_config.dtw_gt_metrics, metric_key)
                        baseline_metrics_traj = base_config.dtw_gt_metrics.(metric_key);
                    else
                        baseline_metrics_traj = struct();
                        baseline_metrics_traj.r10_gt = NaN;
                        baseline_metrics_traj.r50_gt = NaN;
                        baseline_metrics_traj.mean_rank = NaN;
                    end
                    
                    % Store Trajectory Result
                    result_traj = createResultStruct('Trajectory', emb_name, n_coarse, n_fine, ...
                        total_dims, multi_scale, weight_mode_name, dtw_mode_curr, K, db_size, ...
                        query_id, time_stage1_traj, time_stage2_traj, time_dtw_baseline, ...
                        num_dtw_calls_traj, emb_only_metrics_traj, twostage_metrics_traj, ...
                        baseline_metrics_traj, gt_type, num_gt_traj, 0, ...
                        top_k_ids_traj, dtw_reranked_ids_traj);
                    
                    all_results{experiment_counter} = result_traj;
                    
                    % ====================================================
                    % LEVEL 2: SEGMENT-LEVEL TWO-STAGE
                    % ====================================================
                    
                    experiment_counter = experiment_counter + 1;
                    
                    % Get segment embeddings
                    segment_embeddings = query_emb_data.(emb_field_name).segment_embeddings;
                    num_query_segments = length(segment_embeddings);
                    
                    if num_query_segments == 0
                        warning('No segments for query %s - skipping segment level', query_id);
                        continue;
                    end
                    
                    % Get GT segments
                    if base_config.has_ground_truth && isfield(ground_truth_map, query_field)
                    % Real GT segments
                    if isfield(ground_truth_map.(query_field), 'segments')
                        gt_segs_raw = ground_truth_map.(query_field).segments;
                        
                        % Convert struct to cell array (same as Section 4!)
                        if isstruct(gt_segs_raw)
                            seg_fields = fieldnames(gt_segs_raw);
                            gt_seg_ids = cell(length(seg_fields), 1);
                            
                            for sf_idx = 1:length(seg_fields)
                                field_name = seg_fields{sf_idx};
                                gt_seg_ids{sf_idx} = gt_segs_raw.(field_name);
                                
                                % Ensure it's a cell array
                                if ~iscell(gt_seg_ids{sf_idx})
                                    gt_seg_ids{sf_idx} = {gt_seg_ids{sf_idx}};
                                end
                            end
                            
                            num_gt_seg = length(gt_seg_ids{1});
                            gt_type_seg = 'real';
                        else
                            % Unknown format
                            gt_seg_ids = {};
                            num_gt_seg = 0;
                            gt_type_seg = 'none';
                        end
                    else
                        % No segments in ground truth
                        gt_seg_ids = {};
                        num_gt_seg = 0;
                        gt_type_seg = 'none';
                    end
                        
                    elseif base_config.has_proxy_gt
                        % Proxy GT segments (from DTW Top-50)
                        proxy_field = sprintf('%s_%s', query_field, dtw_mode_curr);
                        
                        if isfield(base_config.proxy_gt_map, proxy_field)
                            gt_seg_ids = base_config.proxy_gt_map.(proxy_field).segments;  % Cell array!
                            num_gt_seg = length(base_config.proxy_gt_map.(proxy_field).segments{1});
                            gt_type_seg = 'proxy_dtw';
                        else
                            gt_seg_ids = {};
                            num_gt_seg = 0;
                            gt_type_seg = 'none';
                        end
                        
                    else
                        % No GT at all
                        gt_seg_ids = {};
                        num_gt_seg = 0;
                        gt_type_seg = 'none';
                    end
                    
                    % Initialize aggregated metrics
                    seg_emb_only_metrics_all = [];
                    seg_twostage_metrics_all = [];
                    time_stage1_seg_total = 0;
                    time_stage2_seg_total = 0;
                    num_dtw_calls_seg_total = 0;
                    
                    % Process EACH query segment
                    for seg_idx = 1:num_query_segments
                        seg_emb = segment_embeddings{seg_idx};
                        
                        % STAGE 1: Segment Embedding Retrieval
                        tic_stage1_seg = tic;
                        
                        query_embeddings_seg = seg_emb.query;
                        candidate_embeddings_seg = seg_emb.candidates;
                        candidate_seg_ids = seg_emb.candidates.segment_ids;
                        
                        [top_k_ids_seg, ~, ~] = performRRFFusion(...
                            query_embeddings_seg, candidate_embeddings_seg, weights, K, ...
                            candidate_seg_ids);
                        
                        time_stage1_seg_total = time_stage1_seg_total + toc(tic_stage1_seg);
                        
                        % Get GT for this segment
                        if ~isempty(gt_seg_ids) && iscell(gt_seg_ids) && seg_idx <= length(gt_seg_ids)
                            gt_seg_ids_curr = gt_seg_ids{seg_idx};
                            num_gt_seg_curr = length(gt_seg_ids_curr);
                        else
                            gt_seg_ids_curr = {};
                            num_gt_seg_curr = 0;
                        end
                        
                        % Embedding-Only Metrics (Segment)
                        seg_emb_only_metrics = computeEmbeddingOnlyMetrics(...
                            top_k_ids_seg, gt_seg_ids_curr, num_gt_seg_curr, K);
                        seg_emb_only_metrics_all = [seg_emb_only_metrics_all; seg_emb_only_metrics];
                        
                        % STAGE 2: Segment DTW Reranking
                        % Get query segment data
                        query_seg_data = query_data.segments.position{seg_idx};
                        if strcmp(dtw_mode_curr, 'joint_states')
                            query_seg_seq = query_data.segments.joint{seg_idx};
                        else
                            query_seg_seq = query_seg_data;
                        end
                        
                        [dtw_reranked_ids_seg, time_stage2_seg, num_dtw_calls_seg] = ...
                            performSegmentDTWReranking(top_k_ids_seg, query_seg_seq, ...
                            data_cache_curr, dtw_mode_curr, dtw_config, K);
                        
                        time_stage2_seg_total = time_stage2_seg_total + time_stage2_seg;
                        num_dtw_calls_seg_total = num_dtw_calls_seg_total + num_dtw_calls_seg;
                        
                        % Two-Stage Metrics (Segment)
                        seg_twostage_metrics = computeTwoStageMetrics(...
                            dtw_reranked_ids_seg, gt_seg_ids_curr, num_gt_seg_curr, K);
                        seg_twostage_metrics_all = [seg_twostage_metrics_all; seg_twostage_metrics];
                    end
                    
                    % Average across segments
                    emb_only_metrics_seg = averageMetrics(seg_emb_only_metrics_all);
                    twostage_metrics_seg = averageMetrics(seg_twostage_metrics_all);
                    time_stage1_seg = time_stage1_seg_total / num_query_segments;
                    time_stage2_seg = time_stage2_seg_total / num_query_segments;
                    num_dtw_calls_seg = round(num_dtw_calls_seg_total / num_query_segments);
                    
                    % Baseline DTW Segment Metrics
                    if isfield(base_config.dtw_gt_metrics, metric_key)
                        baseline_metrics_seg = struct();
                        baseline_metrics_seg.r10_gt = base_config.dtw_gt_metrics.(metric_key).seg_r10_gt;
                        baseline_metrics_seg.r50_gt = base_config.dtw_gt_metrics.(metric_key).seg_r50_gt;
                        baseline_metrics_seg.mean_rank = base_config.dtw_gt_metrics.(metric_key).seg_mean_rank;
                    else
                        baseline_metrics_seg = struct();
                        baseline_metrics_seg.r10_gt = NaN;
                        baseline_metrics_seg.r50_gt = NaN;
                        baseline_metrics_seg.mean_rank = NaN;
                    end
                    
                    % Store Segment Result
                    result_seg = createResultStruct('Segment', emb_name, n_coarse, n_fine, ...
                        total_dims, multi_scale, weight_mode_name, dtw_mode_curr, K, db_size, ...
                        query_id, time_stage1_seg, time_stage2_seg, time_dtw_baseline, ...
                        num_dtw_calls_seg, emb_only_metrics_seg, twostage_metrics_seg, ...
                        baseline_metrics_seg, gt_type, 0, num_gt_seg, ...
                        {}, {});  % No rankings stored for segments (optional)
                    
                    all_results{experiment_counter} = result_seg;
                    
                    % Progress
                    if mod(experiment_counter, 100) == 0
                        elapsed = toc(experiment_start_time);
                        fprintf('  Progress: %d/%d (%.1f%%) | Elapsed: %.1f min\n', ...
                            experiment_counter, total_experiments * 2, ...
                            100*experiment_counter/(total_experiments*2), elapsed/60);
                    end
                    
                end % query loop
            end % db_size loop
        end % k_candidates loop
    end % weight_mode loop
end % embedding loop

fprintf('\n[DEBUG] Loop completed:\n');
fprintf('  experiment_counter = %d\n', experiment_counter);
fprintf('  Expected total = %d\n', total_experiments);
fprintf('  Results stored = %d\n', sum(~cellfun(@isempty, all_results)));

total_time = toc(experiment_start_time);

fprintf('\n╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  ALL EXPERIMENTS COMPLETED!                                    ║\n');
fprintf('║  Total time: %.1f minutes (%.2f hours)                         ║\n', ...
    total_time/60, total_time/3600);
fprintf('║  Total experiments: %d                                         ║\n', ...
    total_experiments);
fprintf('║  Average time per experiment: %.2f seconds                     ║\n', ...
    total_time/total_experiments);
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  SECTION 6 COMPLETE                                            ║\n');
fprintf('║  Ready for Section 7: Results Export                           ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

% ========================================================================
%% SECTION 7: RESULTS EXPORT (WITH LEVEL COLUMN)
% ========================================================================

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  SECTION 7: RESULTS EXPORT                                     ║\n');
fprintf('║  Exporting results with Level column (Trajectory + Segment)    ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

% === Create results directory ===
results_dir = 'results';
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

timestamp = datestr(now, 'yyyy-mm-dd_HH-MM');

% ========================================================================
%% CSV 1: MAIN RESULTS (All Metrics with Level)
% ========================================================================

fprintf('=== Exporting Main Results CSV ===\n');

% Count valid results
num_results = sum(~cellfun(@isempty, all_results));
fprintf('  Valid results: %d (including both levels)\n', num_results);

% Pre-allocate arrays
levels = cell(num_results, 1);
experiment_ids = zeros(num_results, 1);
embedding_configs_arr = cell(num_results, 1);
n_coarse_vals = zeros(num_results, 1);
n_fine_vals = zeros(num_results, 1);
total_dims_vals = zeros(num_results, 1);
multi_scale_vals = zeros(num_results, 1);
weight_modes = cell(num_results, 1);
dtw_modes = cell(num_results, 1);
k_candidates_vals = zeros(num_results, 1);
db_sizes_vals = zeros(num_results, 1);
query_ids_vals = cell(num_results, 1);

% Timing
time_stage1 = zeros(num_results, 1);
time_stage2 = zeros(num_results, 1);
time_total = zeros(num_results, 1);
time_baseline = zeros(num_results, 1);
speedup = zeros(num_results, 1);
dtw_calls_made = zeros(num_results, 1);
dtw_calls_saved_pct = zeros(num_results, 1);

% Lower bounds
used_lb_kim = zeros(num_results, 1);
used_lb_keogh = zeros(num_results, 1);

% Embedding-Only Accuracy
embedding_only_r1 = zeros(num_results, 1);
embedding_only_r5 = zeros(num_results, 1);
embedding_only_r10 = zeros(num_results, 1);
embedding_only_r50 = zeros(num_results, 1);
embedding_only_mrr = zeros(num_results, 1);
embedding_only_mean_rank = zeros(num_results, 1);
embedding_only_ndcg10 = zeros(num_results, 1);
embedding_only_ndcg50 = zeros(num_results, 1);

% Two-Stage Accuracy
twostage_r1 = zeros(num_results, 1);
twostage_r5 = zeros(num_results, 1);
twostage_r10 = zeros(num_results, 1);
twostage_r50 = zeros(num_results, 1);
twostage_mrr = zeros(num_results, 1);
twostage_mean_rank = zeros(num_results, 1);
twostage_ndcg10 = zeros(num_results, 1);
twostage_ndcg50 = zeros(num_results, 1);

% Baseline comparison
baseline_dtw_r10 = zeros(num_results, 1);
baseline_dtw_r50 = zeros(num_results, 1);
baseline_dtw_mean_rank = zeros(num_results, 1);

% GT info
gt_types = cell(num_results, 1);
num_gt_traj = zeros(num_results, 1);
num_gt_seg = zeros(num_results, 1);

% Fill arrays from ALL results (both levels!)
result_idx = 0;
for i = 1:length(all_results)
    if isempty(all_results{i})
        continue;
    end
    
    result_idx = result_idx + 1;
    r = all_results{i};
    
    % Configuration
    levels{result_idx} = r.level;  % 'Trajectory' or 'Segment'
    experiment_ids(result_idx) = ceil(i / 2);  % Original experiment ID
    embedding_configs_arr{result_idx} = r.embedding_config;
    n_coarse_vals(result_idx) = r.n_coarse;
    n_fine_vals(result_idx) = r.n_fine;
    total_dims_vals(result_idx) = r.total_dims;
    multi_scale_vals(result_idx) = r.multi_scale;
    weight_modes{result_idx} = r.weight_mode;
    dtw_modes{result_idx} = r.dtw_mode;
    k_candidates_vals(result_idx) = r.k_candidates;
    db_sizes_vals(result_idx) = r.db_size;
    query_ids_vals{result_idx} = r.query_id;
    
    % Timing
    time_stage1(result_idx) = r.time_stage1;
    time_stage2(result_idx) = r.time_stage2;
    time_total(result_idx) = r.time_total;
    time_baseline(result_idx) = r.time_dtw_baseline;
    speedup(result_idx) = r.speedup;
    dtw_calls_made(result_idx) = r.dtw_calls_made;
    dtw_calls_saved_pct(result_idx) = r.dtw_calls_saved_pct;
    
    % Lower bounds
    used_lb_kim(result_idx) = r.used_lb_kim;
    used_lb_keogh(result_idx) = r.used_lb_keogh;
    
    % Embedding-Only
    embedding_only_r1(result_idx) = r.embedding_only_r1;
    embedding_only_r5(result_idx) = r.embedding_only_r5;
    embedding_only_r10(result_idx) = r.embedding_only_r10;
    embedding_only_r50(result_idx) = r.embedding_only_r50;
    embedding_only_mrr(result_idx) = r.embedding_only_mrr;
    embedding_only_mean_rank(result_idx) = r.embedding_only_mean_rank;
    embedding_only_ndcg10(result_idx) = r.embedding_only_ndcg10;
    embedding_only_ndcg50(result_idx) = r.embedding_only_ndcg50;
    
    % Two-Stage
    twostage_r1(result_idx) = r.recall_at_1;
    twostage_r5(result_idx) = r.recall_at_5;
    twostage_r10(result_idx) = r.recall_at_10;
    twostage_r50(result_idx) = r.recall_at_50;
    twostage_mrr(result_idx) = r.mrr;
    twostage_mean_rank(result_idx) = r.mean_rank;
    twostage_ndcg10(result_idx) = r.ndcg_10;
    twostage_ndcg50(result_idx) = r.ndcg_50;
    
    % Baseline
    baseline_dtw_r10(result_idx) = r.baseline_r10;
    baseline_dtw_r50(result_idx) = r.baseline_r50;
    baseline_dtw_mean_rank(result_idx) = r.baseline_mean_rank;
    
    % GT
    gt_types{result_idx} = r.gt_type;
    num_gt_traj(result_idx) = r.num_gt_trajectories;
    num_gt_seg(result_idx) = r.num_gt_segments;
end

% Create table
results_table = table(...
    levels, ...
    experiment_ids, ...
    query_ids_vals, ...
    embedding_configs_arr, ...
    n_coarse_vals, ...
    n_fine_vals, ...
    total_dims_vals, ...
    multi_scale_vals, ...
    weight_modes, ...
    dtw_modes, ...
    k_candidates_vals, ...
    db_sizes_vals, ...
    time_stage1, ...
    time_stage2, ...
    time_total, ...
    time_baseline, ...
    speedup, ...
    dtw_calls_made, ...
    dtw_calls_saved_pct, ...
    used_lb_kim, ...
    used_lb_keogh, ...
    embedding_only_r1, ...
    embedding_only_r5, ...
    embedding_only_r10, ...
    embedding_only_r50, ...
    embedding_only_mrr, ...
    embedding_only_mean_rank, ...
    embedding_only_ndcg10, ...
    embedding_only_ndcg50, ...
    twostage_r1, ...
    twostage_r5, ...
    twostage_r10, ...
    twostage_r50, ...
    twostage_mrr, ...
    twostage_mean_rank, ...
    twostage_ndcg10, ...
    twostage_ndcg50, ...
    baseline_dtw_r10, ...
    baseline_dtw_r50, ...
    baseline_dtw_mean_rank, ...
    gt_types, ...
    num_gt_traj, ...
    num_gt_seg, ...
    'VariableNames', {...
    'Level', ...
    'Experiment_ID', ...
    'Query_ID', ...
    'Embedding_Config', ...
    'n_coarse', ...
    'n_fine', ...
    'total_dims', ...
    'multi_scale', ...
    'Weight_Mode', ...
    'DTW_Mode', ...
    'K_Candidates', ...
    'DB_Size', ...
    'T1_Embeddings', ...
    'T2_DTW', ...
    'T_2S_Total', ...
    'T_DTW_Baseline', ...
    'Speedup', ...
    'DTW_Calls_Made', ...
    'DTW_Calls_Saved_Pct', ...
    'Used_LB_Kim', ...
    'Used_LB_Keogh', ...
    'Emb_R1', ...
    'Emb_R5', ...
    'Emb_R10', ...
    'Emb_R50', ...
    'Emb_MRR', ...
    'Emb_Mean_Rank', ...
    'Emb_NDCG_10', ...
    'Emb_NDCG_50', ...
    'TwoStage_R1', ...
    'TwoStage_R5', ...
    'TwoStage_R10', ...
    'TwoStage_R50', ...
    'TwoStage_MRR', ...
    'TwoStage_Mean_Rank', ...
    'TwoStage_NDCG_10', ...
    'TwoStage_NDCG_50', ...
    'DTW_R10', ...
    'DTW_R50', ...
    'DTW_Mean_Rank', ...
    'GT_Type', ...
    'Num_GT_Trajectories', ...
    'Num_GT_Segments'});

% Save
main_csv_filename = fullfile(results_dir, sprintf('two_stage_results_%s.csv', timestamp));
writetable(results_table, main_csv_filename);

fprintf('  ✓ Main results saved: %s\n', main_csv_filename);
fprintf('  Rows: %d, Columns: %d\n\n', height(results_table), width(results_table));

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  SECTION 7 COMPLETE                                            ║\n');
fprintf('║  All results exported successfully!                            ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');
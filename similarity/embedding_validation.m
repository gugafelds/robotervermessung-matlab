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

% Or add specific folders:
addpath(genpath(pwd));
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
base_config.top_k_trajectories = 500;     % Fixed

% === DIMENSION 1: Embedding Architectures ===
embedding_configs = {
    % === UNTER SWEET SPOT (zeigen dass es schlechter ist) ===
    'Single-2',       0,    2,   false;   % 6-12 dims (zu wenig)
    'Single-5',       0,    5,   false;   % 15-30 dims (Grenze)
    
    % === SWEET SPOT (30-150 dims) - Hauptvergleich ===
    'Single-10',      0,   10,   false;   % 30-60 dims
    'Single-20',      0,   20,   false;   % 60-120 dims
    'Multi-15',       5,   10,   true;    % 45-90 dims (1:2 ratio)
    'Multi-25',       5,   20,   true;    % 75-150 dims (1:4 ratio)
    
    % === PLATEAU (zeigen dass mehr nichts bringt) ===
    'Single-50',      0,   50,   false;   % 150-300 dims
    'Single-100',     0,  100,   false;   % 300-600 dims
    
    % === OPTIONAL: Ein großer zum Absichern ===
    'Single-200',     0,  200,   false;   % 600-1200 dims
};

% === DIMENSION 2: Query Trajectories ===
query_ids = {
    %%%% NEW QUERIES %%%%
    %'1764766034'; % np = 2 / p = 843 / 9 GT

    %'1765987078'; % clean, np = 4 / 10 GT
    %'1765987155'; % clean, np = 4 / 20 GT
    %'1765987310'; % clean, np = 4 / 30 GT

    %% clean
    '1765989370'; % clean, np = 3 / 10 GT
    '1765989294'; % clean, np = 3 / 20 GT
    '1765988821'; % clean, np = 3 / 30 GT
    '1765988920'; % clean; np = 3 / 40 GT
    '1765989411'; % clean; np = 3 / 50 GT
    
    %% noisy - 2 mm

    %'1765990630'; % noisy; np = 3 / 10 GT
    %'1765990747'; % noisy; np = 3 / 20 GT
    %'1765990822'; % noisy; np = 3 / 30 GT
    %'1765991047'; % noisy; np = 3 / 40 GT
    %'1765991234'; % noisy; np = 3 / 50 GT
    
    %% noisy - 5 mm

    %'1765991190'; % noisy; np = 3 / 10 GT
    %'1765991445'; % noisy; np = 3 / 20 GT
    %'1765991515'; % noisy; np = 3 / 30 GT
    %'1765991949'; % noisy; np = 3 / 40 GT 
    %'1765991743'; % noisy; np = 3 / 50 GT
};

% === DIMENSION 3: DTW Mode + Weight Combinations ===
weight_mode_configs = {
    % Joint space
    'Joint only',           'joint_states',  [0, 1, 0, 0, 0];
    'Joint + Position',     'joint_states',  [1, 1, 0, 0, 0];
    'Joint + Orient',       'joint_states',  [0, 1, 1, 0, 0];
    'Joint + Velocity',     'joint_states',  [0, 1, 0, 1, 0];
    'Joint + Meta',         'joint_states',  [0, 1, 0, 0, 1];
    'Joint + All',          'joint_states',  [1, 1, 1, 1, 1];
    
    % Position space
    'Position only',        'position',      [1, 0, 0, 0, 0];
    'Pos + Joint',          'position',      [1, 1, 0, 0, 0];
    'Pos + Orient',         'position',      [1, 0, 1, 0, 0];
    'Pos + Velocity',       'position',      [1, 0, 0, 1, 0];
    'Pos + Meta',           'position',      [1, 0, 0, 0, 1];
    'Pos + All',            'position',      [1, 1, 1, 1, 1];
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

chunk_size = 100;  % For batch loading

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
dtw_config.lb_kim_keep_ratio = 1;
dtw_config.lb_keogh_candidates = 500;
dtw_config.cdtw_window = 0.2;
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
            % SEGMENT-LEVEL GT METRICS (CORRECTED - CONSISTENT WITH GTvsEB!)
            % ============================================================
            segment_rankings = dtw_cache.(query_field).(dtw_mode).segment_rankings;
            num_segments = length(segment_rankings);
            
            % Calculate GT coverage for EACH segment separately (like GTvsEB!)
            gt_coverage_segments = cell(num_segments, 1);
            
            for seg_idx = 1:num_segments
                if isempty(segment_rankings{seg_idx})
                    continue;
                end
                
                seg_ranking = segment_rankings{seg_idx};
                
                % Find GT ranks for this segment
                seg_gt_ranks = zeros(num_gt, 1);
                
                for gt_idx = 1:num_gt
                    gt_id = gt_ids_for_query{gt_idx};
                    rank_idx = find(strcmp(seg_ranking.bahn_id, gt_id), 1);
                    
                    if ~isempty(rank_idx)
                        seg_gt_ranks(gt_idx) = rank_idx;
                    else
                        seg_gt_ranks(gt_idx) = inf;
                    end
                end
                
                % Calculate metrics for THIS segment
                valid_seg_ranks = seg_gt_ranks(seg_gt_ranks < inf);
                
                if isempty(valid_seg_ranks)
                    % No GT found in this segment
                    seg_cov = struct();
                    seg_cov.k_values = [1, 3, 5, 10, 50];
                    seg_cov.recall_at_k = zeros(5, 1);
                    seg_cov.mean_gt_rank = inf;
                    seg_cov.overall_coverage_point = inf;
                else
                    % Calculate for this segment using SAME logic as GTvsEB!
                    seg_cov = struct();
                    seg_cov.k_values = [1, 3, 5, 10, 50];
                    seg_cov.recall_at_k = zeros(5, 1);
                    
                    k_vals = [1, 3, 5, 10, 50];
                    for k_idx = 1:5
                        K = k_vals(k_idx);
                        num_gt_in_top_k = sum(valid_seg_ranks <= K);
                        % SAME FORMULA AS calculateGTCoverage!
                        seg_cov.recall_at_k(k_idx) = num_gt_in_top_k / min(K, num_gt);
                    end
                    
                    seg_cov.mean_gt_rank = mean(valid_seg_ranks);
                    seg_cov.overall_coverage_point = max(valid_seg_ranks);
                end
                
                gt_coverage_segments{seg_idx} = seg_cov;
            end
            
            % Average across segments (using same logic as averageGTCoverage)
            valid_segments = gt_coverage_segments(~cellfun(@isempty, gt_coverage_segments));
            
            if isempty(valid_segments)
                seg_p_gt_dtw = inf;
                seg_mean_gt_rank_dtw = inf;
                seg_r50_gt_dtw = 0;
                seg_r10_gt_dtw = 0;
                seg_r5_gt_dtw = 0;
                seg_r3_gt_dtw = 0;
                seg_r1_gt_dtw = 0;
            else
                % Average recall values
                recall_sum = zeros(5, 1);
                count = 0;
                mean_rank_sum = 0;
                coverage_sum = 0;
                
                for i = 1:length(valid_segments)
                    seg = valid_segments{i};
                    recall_sum = recall_sum + seg.recall_at_k;
                    mean_rank_sum = mean_rank_sum + seg.mean_gt_rank;
                    coverage_sum = coverage_sum + seg.overall_coverage_point;
                    count = count + 1;
                end
                
                avg_recall = recall_sum / count;
                
                % Extract averaged values
                seg_r1_gt_dtw = avg_recall(1);   % K=1
                seg_r3_gt_dtw = avg_recall(2);   % K=3
                seg_r5_gt_dtw = avg_recall(3);   % K=5
                seg_r10_gt_dtw = avg_recall(4);  % K=10
                seg_r50_gt_dtw = avg_recall(5);  % K=50
                
                seg_mean_gt_rank_dtw = mean_rank_sum / count;
                seg_p_gt_dtw = coverage_sum / count;
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

            % ============================================================
            % NDCG METRICS FOR DTW BASELINE (Trajectory) - NEW!
            % ============================================================
            
            fprintf('  Computing DTW baseline NDCG metrics (trajectory)...\n');
            
            if num_gt > 0 && ~isempty(valid_ranks)
                % Extract GT similarities from DTW ranking
                gt_dtw_similarities = zeros(num_gt, 1);
                
                for gt_idx = 1:num_gt
                    gt_id = gt_ids_for_query{gt_idx};
                    rank_idx = find(strcmp(dtw_ranking.bahn_id, gt_id), 1);
                    
                    if ~isempty(rank_idx)
                        % Get DTW distance
                        dtw_distance = dtw_ranking.dtw_distance(rank_idx);
                        gt_dtw_similarities(gt_idx) = dtw_distance;
                    else
                        gt_dtw_similarities(gt_idx) = inf;
                    end
                end
                
                % Normalize distances to [0, 1] and convert to similarities
                % (Same logic as extractGTSimilarities)
                valid_dists = gt_dtw_similarities(gt_dtw_similarities < inf);
                
                if ~isempty(valid_dists)
                    max_dist = max(valid_dists);
                    if max_dist > 0
                        normalized_dists = gt_dtw_similarities / max_dist;
                    else
                        normalized_dists = gt_dtw_similarities;
                    end
                    gt_dtw_similarities = 1 ./ (1 + normalized_dists);
                    gt_dtw_similarities(gt_dtw_similarities == inf) = 0;
                else
                    gt_dtw_similarities = zeros(num_gt, 1);
                end
                
                % Compute NDCG metrics using DTW ranking
                top_k_ids_dtw = dtw_ranking.bahn_id;
                K = dtw_config.top_k_trajectories;
                
                dtw_ndcg_metrics = computeEmbeddingOnlyMetrics(...
                    top_k_ids_dtw, gt_ids_for_query, gt_dtw_similarities, num_gt, K);
                
                % Store DTW NDCG metrics
                dtw_gt_metrics.(query_field).(dtw_mode).ndcg_10 = dtw_ndcg_metrics.ndcg10;
                dtw_gt_metrics.(query_field).(dtw_mode).ndcg_50 = dtw_ndcg_metrics.ndcg50;
                dtw_gt_metrics.(query_field).(dtw_mode).mrr = dtw_ndcg_metrics.mrr;
                
                fprintf('    DTW NDCG@10: %.4f, NDCG@50: %.4f, MRR: %.4f\n', ...
                    dtw_ndcg_metrics.ndcg10, dtw_ndcg_metrics.ndcg50, dtw_ndcg_metrics.mrr);
            else
                % No GTs or not found
                dtw_gt_metrics.(query_field).(dtw_mode).ndcg_10 = NaN;
                dtw_gt_metrics.(query_field).(dtw_mode).ndcg_50 = NaN;
                dtw_gt_metrics.(query_field).(dtw_mode).mrr = NaN;
            end
            
            % Segment metrics
            dtw_gt_metrics.(query_field).(dtw_mode).seg_p_gt = seg_p_gt_dtw;
            dtw_gt_metrics.(query_field).(dtw_mode).seg_mean_gt_rank = seg_mean_gt_rank_dtw;
            dtw_gt_metrics.(query_field).(dtw_mode).seg_r50_gt = seg_r50_gt_dtw;
            dtw_gt_metrics.(query_field).(dtw_mode).seg_r10_gt = seg_r10_gt_dtw;
            dtw_gt_metrics.(query_field).(dtw_mode).seg_r5_gt = seg_r5_gt_dtw;
            dtw_gt_metrics.(query_field).(dtw_mode).seg_r3_gt = seg_r3_gt_dtw;
            dtw_gt_metrics.(query_field).(dtw_mode).seg_r1_gt = seg_r1_gt_dtw;

            % ============================================================
            % NDCG METRICS FOR DTW BASELINE (Segment) - NEW!
            % ============================================================
            
            fprintf('  Computing DTW baseline NDCG metrics (segments)...\n');
            
            % Calculate GT coverage for EACH segment separately
            segment_ndcg_results = cell(num_segments, 1);
            
            for seg_idx = 1:num_segments
                if isempty(segment_rankings{seg_idx})
                    continue;
                end
                
                seg_ranking = segment_rankings{seg_idx};
                
                % Extract GT similarities for this segment
                seg_gt_similarities = zeros(num_gt, 1);
                seg_gt_found = 0;
                
                for gt_idx = 1:num_gt
                    gt_id = gt_ids_for_query{gt_idx};
                    rank_idx = find(strcmp(seg_ranking.bahn_id, gt_id), 1);
                    
                    if ~isempty(rank_idx)
                        dtw_distance = seg_ranking.dtw_distance(rank_idx);
                        seg_gt_similarities(gt_idx) = dtw_distance;
                        seg_gt_found = seg_gt_found + 1;
                    else
                        seg_gt_similarities(gt_idx) = inf;
                    end
                end
                
                if seg_gt_found == 0
                    % No GTs found in this segment
                    segment_ndcg_results{seg_idx} = [];
                    continue;
                end
                
                % Normalize similarities
                valid_dists = seg_gt_similarities(seg_gt_similarities < inf);
                
                if ~isempty(valid_dists)
                    max_dist = max(valid_dists);
                    if max_dist > 0
                        normalized_dists = seg_gt_similarities / max_dist;
                    else
                        normalized_dists = seg_gt_similarities;
                    end
                    seg_gt_similarities = 1 ./ (1 + normalized_dists);
                    seg_gt_similarities(seg_gt_similarities == inf) = 0;
                end
                
                % Compute NDCG for this segment
                top_k_seg_ids = seg_ranking.bahn_id;
                K_seg = dtw_config.top_k_trajectories;
                
                seg_ndcg = computeEmbeddingOnlyMetrics(...
                    top_k_seg_ids, gt_ids_for_query, seg_gt_similarities, num_gt, K_seg);
                
                segment_ndcg_results{seg_idx} = seg_ndcg;
            end
            
            % Average across valid segments
            valid_seg_ndcg = segment_ndcg_results(~cellfun(@isempty, segment_ndcg_results));
            
            if isempty(valid_seg_ndcg)
                seg_ndcg_10_dtw = NaN;
                seg_ndcg_50_dtw = NaN;
                seg_mrr_dtw = NaN;
            else
                % Average NDCG values
                ndcg10_sum = 0;
                ndcg50_sum = 0;
                mrr_sum = 0;
                count = 0;
                
                for i = 1:length(valid_seg_ndcg)
                    seg = valid_seg_ndcg{i};
                    ndcg10_sum = ndcg10_sum + seg.ndcg10;
                    ndcg50_sum = ndcg50_sum + seg.ndcg50;
                    mrr_sum = mrr_sum + seg.mrr;
                    count = count + 1;
                end
                
                seg_ndcg_10_dtw = ndcg10_sum / count;
                seg_ndcg_50_dtw = ndcg50_sum / count;
                seg_mrr_dtw = mrr_sum / count;
            end
            
            % Store segment DTW NDCG metrics
            dtw_gt_metrics.(query_field).(dtw_mode).seg_ndcg_10 = seg_ndcg_10_dtw;
            dtw_gt_metrics.(query_field).(dtw_mode).seg_ndcg_50 = seg_ndcg_50_dtw;
            dtw_gt_metrics.(query_field).(dtw_mode).seg_mrr = seg_mrr_dtw;
            
            fprintf('    Segment Avg - NDCG@10: %.4f, NDCG@50: %.4f, MRR: %.4f\n', ...
                seg_ndcg_10_dtw, seg_ndcg_50_dtw, seg_mrr_dtw);
            
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
%% SIMPLIFIED CSV EXPORT (Two-Stage Style)
% ========================================================================

fprintf('\n=== Creating Results Tables (Simplified) ===\n');

% Count valid results (should be total_experiments)
num_results = length(all_results);
fprintf('  Total results: %d\n', num_results);

% ========================================================================
% PRE-ALLOCATE ARRAYS FOR ALL COLUMNS
% ========================================================================

% Configuration
levels = cell(num_results, 1);
query_ids_arr = cell(num_results, 1);
embedding_configs_arr = cell(num_results, 1);
n_coarse_vals = zeros(num_results, 1);
n_fine_vals = zeros(num_results, 1);
total_dims_vals = zeros(num_results, 1);
multi_scale_vals = zeros(num_results, 1);
weight_modes = cell(num_results, 1);
dtw_modes = cell(num_results, 1);
w_pos_vals = zeros(num_results, 1);
w_joint_vals = zeros(num_results, 1);
w_orient_vals = zeros(num_results, 1);
w_vel_vals = zeros(num_results, 1);
w_meta_vals = zeros(num_results, 1);

% Spearman & Coverage (DTW vs EB)
spearman_dtw_eb = zeros(num_results, 1);
rk_dtw_eb = zeros(num_results, 1);
r50_dtw_eb = zeros(num_results, 1);
r10_dtw_eb = zeros(num_results, 1);
r5_dtw_eb = zeros(num_results, 1);
r3_dtw_eb = zeros(num_results, 1);
r1_dtw_eb = zeros(num_results, 1);
p_dtw_eb = zeros(num_results, 1);
ndcg10_dtw_eb = zeros(num_results, 1);
ndcg50_dtw_eb = zeros(num_results, 1);

% GT vs Embedding
r50_gt_eb = zeros(num_results, 1);
r10_gt_eb = zeros(num_results, 1);
r5_gt_eb = zeros(num_results, 1);
r3_gt_eb = zeros(num_results, 1);
r1_gt_eb = zeros(num_results, 1);
p_gt_eb = zeros(num_results, 1);
mean_rank_gt_eb = zeros(num_results, 1);
ndcg10_gt_eb = zeros(num_results, 1);
ndcg50_gt_eb = zeros(num_results, 1);
mrr_gt_eb = zeros(num_results, 1);

% GT vs DTW (Baseline)
r50_gt_dtw = zeros(num_results, 1);
r10_gt_dtw = zeros(num_results, 1);
r5_gt_dtw = zeros(num_results, 1);
r3_gt_dtw = zeros(num_results, 1);
r1_gt_dtw = zeros(num_results, 1);
p_gt_dtw = zeros(num_results, 1);
mean_rank_gt_dtw = zeros(num_results, 1);
ndcg10_gt_dtw = zeros(num_results, 1);
ndcg50_gt_dtw = zeros(num_results, 1);
mrr_gt_dtw = zeros(num_results, 1);

% GT Info
num_gt_vals = zeros(num_results, 1);

% Timing
runtimes = zeros(num_results, 1);

% ========================================================================
% FILL ARRAYS FROM RESULTS (TRAJECTORY LEVEL)
% ========================================================================

fprintf('Extracting trajectory-level results...\n');

for i = 1:total_experiments
    r = all_results{i};
    
    % Configuration
    levels{i} = 'Trajectory';
    query_ids_arr{i} = r.query_id;
    embedding_configs_arr{i} = r.embedding_config;
    n_coarse_vals(i) = r.n_coarse;
    n_fine_vals(i) = r.n_fine;
    total_dims_vals(i) = r.total_dims;
    multi_scale_vals(i) = (r.n_coarse > 0);
    weight_modes{i} = r.weight_mode;
    dtw_modes{i} = r.dtw_mode;
    w_pos_vals(i) = r.weight_pos;
    w_joint_vals(i) = r.weight_joint;
    w_orient_vals(i) = r.weight_orient;
    w_vel_vals(i) = r.weight_vel;
    w_meta_vals(i) = r.weight_meta;
    
    % Spearman & Coverage (DTW vs EB)
    spearman_dtw_eb(i) = r.spearman;
    rk_dtw_eb(i) = r.p_at_k;
    r50_dtw_eb(i) = r.recall_at_50;
    r10_dtw_eb(i) = r.recall_at_10;
    r5_dtw_eb(i) = r.recall_at_5;
    r3_dtw_eb(i) = r.recall_at_3;
    r1_dtw_eb(i) = r.recall_at_1;
    p_dtw_eb(i) = r.p_all_traj;
    ndcg10_dtw_eb(i) = r.ndcg_10_dtw_eb;
    ndcg50_dtw_eb(i) = r.ndcg_50_dtw_eb;
    
    % GT vs Embedding
    if isfield(r, 'gt_recall_at_50')
        r50_gt_eb(i) = r.gt_recall_at_50;
        r10_gt_eb(i) = r.gt_recall_at_10;
        r5_gt_eb(i) = r.gt_recall_at_5;
        r3_gt_eb(i) = r.gt_recall_at_3;
        r1_gt_eb(i) = r.gt_recall_at_1;
        p_gt_eb(i) = r.p_gt_traj;
        mean_rank_gt_eb(i) = r.mean_gt_rank;
        ndcg10_gt_eb(i) = r.ndcg_10_traj;
        ndcg50_gt_eb(i) = r.ndcg_50_traj;
        mrr_gt_eb(i) = r.mrr_traj;
    else
        r50_gt_eb(i) = NaN;
        r10_gt_eb(i) = NaN;
        r5_gt_eb(i) = NaN;
        r3_gt_eb(i) = NaN;
        r1_gt_eb(i) = NaN;
        p_gt_eb(i) = NaN;
        mean_rank_gt_eb(i) = NaN;
        ndcg10_gt_eb(i) = NaN;
        ndcg50_gt_eb(i) = NaN;
        mrr_gt_eb(i) = NaN;
    end
    
    % GT vs DTW (Baseline)
    if isfield(base_config, 'dtw_gt_metrics')
        query_field = sprintf('q_%s', strrep(r.query_id, '-', '_'));
        dtw_mode = r.dtw_mode;
        
        if isfield(base_config.dtw_gt_metrics, query_field) && ...
           isfield(base_config.dtw_gt_metrics.(query_field), dtw_mode)
            
            dtw_gt = base_config.dtw_gt_metrics.(query_field).(dtw_mode);
            
            r50_gt_dtw(i) = dtw_gt.r50_gt;
            r10_gt_dtw(i) = dtw_gt.r10_gt;
            r5_gt_dtw(i) = dtw_gt.r5_gt;
            r3_gt_dtw(i) = dtw_gt.r3_gt;
            r1_gt_dtw(i) = dtw_gt.r1_gt;
            p_gt_dtw(i) = dtw_gt.p_gt;
            mean_rank_gt_dtw(i) = dtw_gt.mean_gt_rank;
            
            % NDCG DTW Baseline (NEW!)
            if isfield(dtw_gt, 'ndcg_10')
                ndcg10_gt_dtw(i) = dtw_gt.ndcg_10;
                ndcg50_gt_dtw(i) = dtw_gt.ndcg_50;
                mrr_gt_dtw(i) = dtw_gt.mrr;
            else
                ndcg10_gt_dtw(i) = NaN;
                ndcg50_gt_dtw(i) = NaN;
                mrr_gt_dtw(i) = NaN;
            end
        else
            r50_gt_dtw(i) = NaN;
            r10_gt_dtw(i) = NaN;
            r5_gt_dtw(i) = NaN;
            r3_gt_dtw(i) = NaN;
            r1_gt_dtw(i) = NaN;
            p_gt_dtw(i) = NaN;
            mean_rank_gt_dtw(i) = NaN;
            ndcg10_gt_dtw(i) = NaN;
            ndcg50_gt_dtw(i) = NaN;
            mrr_gt_dtw(i) = NaN;
        end
    else
        r50_gt_dtw(i) = NaN;
        r10_gt_dtw(i) = NaN;
        r5_gt_dtw(i) = NaN;
        r3_gt_dtw(i) = NaN;
        r1_gt_dtw(i) = NaN;
        p_gt_dtw(i) = NaN;
        mean_rank_gt_dtw(i) = NaN;
        ndcg10_gt_dtw(i) = NaN;
        ndcg50_gt_dtw(i) = NaN;
        mrr_gt_dtw(i) = NaN;
    end
    
    % GT Info
    if isfield(r, 'num_gt')
        num_gt_vals(i) = r.num_gt;
    else
        num_gt_vals(i) = NaN;
    end
    
    % Timing
    runtimes(i) = r.exp_runtime / 60;  % Convert to minutes
end

fprintf('✓ Trajectory-level: %d experiments\n', total_experiments);

% ========================================================================
% FILL ARRAYS FROM RESULTS (SEGMENT LEVEL)
% ========================================================================

fprintf('Extracting segment-level results...\n');

seg_offset = total_experiments;

for i = 1:total_experiments
    idx = seg_offset + i;
    r = all_results{i};
    
    % Configuration (same as trajectory)
    levels{idx} = 'Segment';
    query_ids_arr{idx} = r.query_id;
    embedding_configs_arr{idx} = r.embedding_config;
    n_coarse_vals(idx) = r.n_coarse;
    n_fine_vals(idx) = r.n_fine;
    total_dims_vals(idx) = r.total_dims;
    multi_scale_vals(idx) = (r.n_coarse > 0);
    weight_modes{idx} = r.weight_mode;
    dtw_modes{idx} = r.dtw_mode;
    w_pos_vals(idx) = r.weight_pos;
    w_joint_vals(idx) = r.weight_joint;
    w_orient_vals(idx) = r.weight_orient;
    w_vel_vals(idx) = r.weight_vel;
    w_meta_vals(idx) = r.weight_meta;
    
    % Spearman & Coverage (DTW vs EB) - SEGMENT
    spearman_dtw_eb(idx) = r.seg_spearman;
    rk_dtw_eb(idx) = r.seg_p_at_k;
    r50_dtw_eb(idx) = r.seg_recall_at_50;
    r10_dtw_eb(idx) = r.seg_recall_at_10;
    r5_dtw_eb(idx) = r.seg_recall_at_5;
    r3_dtw_eb(idx) = r.seg_recall_at_3;
    r1_dtw_eb(idx) = r.seg_recall_at_1;
    p_dtw_eb(idx) = r.p_all_seg;
    ndcg10_dtw_eb(idx) = r.seg_ndcg_10_dtw_eb;
    ndcg50_dtw_eb(idx) = r.seg_ndcg_50_dtw_eb;

    
    % GT vs Embedding - SEGMENT
    if isfield(r, 'seg_gt_recall_at_50')
        r50_gt_eb(idx) = r.seg_gt_recall_at_50;
        r10_gt_eb(idx) = r.seg_gt_recall_at_10;
        r5_gt_eb(idx) = r.seg_gt_recall_at_5;
        r3_gt_eb(idx) = r.seg_gt_recall_at_3;
        r1_gt_eb(idx) = r.seg_gt_recall_at_1;
        p_gt_eb(idx) = r.p_gt_seg;
        mean_rank_gt_eb(idx) = r.seg_mean_gt_rank;
        ndcg10_gt_eb(idx) = r.ndcg_10_seg;
        ndcg50_gt_eb(idx) = r.ndcg_50_seg;
        mrr_gt_eb(idx) = r.mrr_seg;
    else
        r50_gt_eb(idx) = NaN;
        r10_gt_eb(idx) = NaN;
        r5_gt_eb(idx) = NaN;
        r3_gt_eb(idx) = NaN;
        r1_gt_eb(idx) = NaN;
        p_gt_eb(idx) = NaN;
        mean_rank_gt_eb(idx) = NaN;
        ndcg10_gt_eb(idx) = NaN;
        ndcg50_gt_eb(idx) = NaN;
        mrr_gt_eb(idx) = NaN;
    end
    
    % GT vs DTW - SEGMENT
    if isfield(base_config, 'dtw_gt_metrics')
        query_field = sprintf('q_%s', strrep(r.query_id, '-', '_'));
        dtw_mode = r.dtw_mode;
        
        if isfield(base_config.dtw_gt_metrics, query_field) && ...
           isfield(base_config.dtw_gt_metrics.(query_field), dtw_mode)
            
            dtw_gt = base_config.dtw_gt_metrics.(query_field).(dtw_mode);
            
            r50_gt_dtw(idx) = dtw_gt.seg_r50_gt;
            r10_gt_dtw(idx) = dtw_gt.seg_r10_gt;
            r5_gt_dtw(idx) = dtw_gt.seg_r5_gt;
            r3_gt_dtw(idx) = dtw_gt.seg_r3_gt;
            r1_gt_dtw(idx) = dtw_gt.seg_r1_gt;
            p_gt_dtw(idx) = dtw_gt.seg_p_gt;
            mean_rank_gt_dtw(idx) = dtw_gt.seg_mean_gt_rank;
            
            % NDCG DTW Baseline - SEGMENT (NEW!)
            if isfield(dtw_gt, 'seg_ndcg_10')
                ndcg10_gt_dtw(idx) = dtw_gt.seg_ndcg_10;
                ndcg50_gt_dtw(idx) = dtw_gt.seg_ndcg_50;
                mrr_gt_dtw(idx) = dtw_gt.seg_mrr;
            else
                ndcg10_gt_dtw(idx) = NaN;
                ndcg50_gt_dtw(idx) = NaN;
                mrr_gt_dtw(idx) = NaN;
            end
        else
            r50_gt_dtw(idx) = NaN;
            r10_gt_dtw(idx) = NaN;
            r5_gt_dtw(idx) = NaN;
            r3_gt_dtw(idx) = NaN;
            r1_gt_dtw(idx) = NaN;
            p_gt_dtw(idx) = NaN;
            mean_rank_gt_dtw(idx) = NaN;
            ndcg10_gt_dtw(idx) = NaN;
            ndcg50_gt_dtw(idx) = NaN;
            mrr_gt_dtw(idx) = NaN;
        end
    else
        r50_gt_dtw(idx) = NaN;
        r10_gt_dtw(idx) = NaN;
        r5_gt_dtw(idx) = NaN;
        r3_gt_dtw(idx) = NaN;
        r1_gt_dtw(idx) = NaN;
        p_gt_dtw(idx) = NaN;
        mean_rank_gt_dtw(idx) = NaN;
        ndcg10_gt_dtw(idx) = NaN;
        ndcg50_gt_dtw(idx) = NaN;
        mrr_gt_dtw(idx) = NaN;
    end
    
    % GT Info (same as trajectory)
    num_gt_vals(idx) = num_gt_vals(i);
    
    % Timing (same as trajectory)
    runtimes(idx) = runtimes(i);
end

fprintf('✓ Segment-level: %d experiments\n', total_experiments);

% ========================================================================
% CREATE COMBINED TABLE
% ========================================================================

fprintf('Creating combined table...\n');

combined_table = table(...
    levels, ...
    query_ids_arr, ...
    embedding_configs_arr, ...
    n_coarse_vals, ...
    n_fine_vals, ...
    total_dims_vals, ...
    multi_scale_vals, ...
    weight_modes, ...
    dtw_modes, ...
    w_pos_vals, ...
    w_joint_vals, ...
    w_orient_vals, ...
    w_vel_vals, ...
    w_meta_vals, ...
    spearman_dtw_eb, ...
    rk_dtw_eb, ...
    r50_dtw_eb, ...
    r10_dtw_eb, ...
    r5_dtw_eb, ...
    r3_dtw_eb, ...
    r1_dtw_eb, ...
    p_dtw_eb, ...
    ndcg10_dtw_eb, ...
    ndcg50_dtw_eb, ...
    r50_gt_eb, ...
    r10_gt_eb, ...
    r5_gt_eb, ...
    r3_gt_eb, ...
    r1_gt_eb, ...
    p_gt_eb, ...
    mean_rank_gt_eb, ...
    ndcg10_gt_eb, ...
    ndcg50_gt_eb, ...
    mrr_gt_eb, ...
    r50_gt_dtw, ...
    r10_gt_dtw, ...
    r5_gt_dtw, ...
    r3_gt_dtw, ...
    r1_gt_dtw, ...
    p_gt_dtw, ...
    mean_rank_gt_dtw, ...
    ndcg10_gt_dtw, ...
    ndcg50_gt_dtw, ...
    mrr_gt_dtw, ...
    num_gt_vals, ...
    runtimes, ...
    'VariableNames', {...
    'Level', ...
    'Query_Bahn_ID', ...
    'Embedding_Config', ...
    'N_Coarse', ...
    'N_Fine', ...
    'Total_Dims', ...
    'Multi_Scale', ...
    'Weight_Mode', ...
    'DTW_Mode', ...
    'W_Pos', ...
    'W_Joint', ...
    'W_Orient', ...
    'W_Vel', ...
    'W_Meta', ...
    'Spearman_DTWvsEB', ...
    'R@K_DTWvsEB', ...
    'R@50_DTWvsEB', ...
    'R@10_DTWvsEB', ...
    'R@5_DTWvsEB', ...
    'R@3_DTWvsEB', ...
    'R@1_DTWvsEB', ...
    'P_DTWvsEB', ...
    'NDCG@10_DTWvsEB', ...
    'NDCG@50_DTWvsEB', ...
    'R@50_GTvsEB', ...
    'R@10_GTvsEB', ...
    'R@5_GTvsEB', ...
    'R@3_GTvsEB', ...
    'R@1_GTvsEB', ...
    'P_GTvsEB', ...
    'Mean_GTvsEB_Rank', ...
    'NDCG@10_GTvsEB', ...
    'NDCG@50_GTvsEB', ...
    'MRR_GTvsEB', ...
    'R@50_GTvsDTW', ...
    'R@10_GTvsDTW', ...
    'R@5_GTvsDTW', ...
    'R@3_GTvsDTW', ...
    'R@1_GTvsDTW', ...
    'P_GTvsDTW', ...
    'Mean_GTvsDTW_Rank', ...
    'NDCG@10_GTvsDTW', ...
    'NDCG@50_GTvsDTW', ...
    'MRR_GTvsDTW', ...
    'Num_GT', ...
    'Runtime_min'});

fprintf('✓ Combined table created: %d rows × %d columns\n', ...
    height(combined_table), width(combined_table));

fprintf('✓ Tables created successfully\n\n');

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
    output_dir = fullfile(project_root, 'similarity/results/');
else
    % Fallback: use current directory
    warning('Could not find robotervermessung-matlab in path, using current directory');
    output_dir = fullfile(pwd, 'similarity/results/');
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

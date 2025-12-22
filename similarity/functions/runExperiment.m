function results = runExperiment(config)
% RUNEXPERIMENT - Execute a single experiment with pre-loaded caches
%
%   This function runs dtw_baseline.m with the given configuration,
%   using pre-loaded data_cache, dtw_cache, and embeddings_cache.
%
%   INPUT:
%       config - Configuration struct with:
%           .query_bahn_id          - Query trajectory ID
%           .dtw_mode               - 'position' or 'joint_states'
%           .weights                - [pos, joint, orient, vel, meta]
%           .n_coarse, .n_fine      - Embedding dimensions
%           .use_multi_scale        - Boolean
%           .embedding_config_name  - Name for cache lookup
%           .data_cache             - Pre-loaded data cache
%           .dtw_cache              - Pre-loaded DTW cache
%           .embeddings_cache       - Pre-loaded embeddings cache (optional)
%           ... (other DTW/embedding parameters)
%
%   OUTPUT:
%       results - Struct with metrics:
%           .spearman, .p_at_k, .p_at_10, .p_at_5, .p_at_3, .p_at_1
%           .seg_spearman, .seg_p_at_k, etc.
%           .dtw_mode, .weights, .n_coarse, .n_fine, etc.
%
%   Author: Gustavo Barros
%   Date: 02.12.2025

%% ========================================================================
%  SETUP WORKSPACE FOR dtw_baseline.m
%  ========================================================================

% Clear only the variables that dtw_baseline.m will create
% Keep function scope clean
clearvars -except config

% ========================================================================
% EXTRACT CONFIGURATION PARAMETERS
% ========================================================================

% Query
query_bahn_id = config.query_bahn_id;

% DTW Mode & Weights
dtw_mode = config.dtw_mode;
weights = config.weights;

% Embedding Architecture
n_coarse = config.n_coarse;
n_fine = config.n_fine;
use_multi_scale = config.use_multi_scale;
embedding_config_name = config.embedding_config_name;

% Other parameters with defaults
if isfield(config, 'top_k_trajectories')
    top_k_trajectories = config.top_k_trajectories;
else
    top_k_trajectories = 50;
end

if isfield(config, 'rrf_k')
    rrf_k = config.rrf_k;
else
    rrf_k = 60;
end

if isfield(config, 'norm_strategy')
    norm_strategy = config.norm_strategy;
else
    norm_strategy = 'max_extent';
end

if isfield(config, 'cdtw_window')
    cdtw_window = config.cdtw_window;
else
    cdtw_window = 0.10;
end

if isfield(config, 'normalize_dtw')
    normalize_dtw = config.normalize_dtw;
else
    normalize_dtw = false;
end

if isfield(config, 'use_rotation_alignment')
    use_rotation_alignment = config.use_rotation_alignment;
else
    use_rotation_alignment = false;
end

% ========================================================================
% ⭐ CRITICAL: LOAD CACHES INTO WORKSPACE
% ========================================================================

if ~isfield(config, 'data_cache') || isempty(config.data_cache)
    error('data_cache is required in config');
end
data_cache = config.data_cache;

if ~isfield(config, 'dtw_cache') || isempty(config.dtw_cache)
    error('dtw_cache is required in config');
end
dtw_cache = config.dtw_cache;

% Embeddings cache is optional (for backward compatibility)
if isfield(config, 'embeddings_cache') && ~isempty(config.embeddings_cache)
    embeddings_cache = config.embeddings_cache;
else
    embeddings_cache = [];
end

%% ========================================================================
%  RUN dtw_baseline.m
%  ========================================================================

% dtw_baseline.m expects these variables in workspace:
% - query_bahn_id, dtw_mode, weights
% - n_coarse, n_fine, use_multi_scale, embedding_config_name
% - data_cache, dtw_cache, embeddings_cache
% - top_k_trajectories, rrf_k, norm_strategy, cdtw_window, etc.

dtw_baseline;

% ========================================================================
% 🔍 DEBUG: GT COVERAGE RESULTS FROM WORKSPACE
% ========================================================================

fprintf('\n╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║ 🔍 DEBUG: GT Coverage in Workspace (Query: %s)\n', query_bahn_id);
fprintf('╚════════════════════════════════════════════════════════════════╝\n');

% Check trajectory-level GT coverage
if exist('gt_coverage_traj', 'var')
    fprintf('✓ gt_coverage_traj EXISTS\n');
    if ~isempty(gt_coverage_traj)
        fprintf('  Fields: %s\n', strjoin(fieldnames(gt_coverage_traj), ', '));
        fprintf('  num_gt: %d\n', gt_coverage_traj.num_gt);
        fprintf('  num_gt_found: %d\n', gt_coverage_traj.num_gt_found);
        fprintf('  gt_ranks: [%s]\n', num2str(gt_coverage_traj.gt_ranks'));
        fprintf('  recall_at_k: [%s]\n', num2str(gt_coverage_traj.recall_at_k'));
        fprintf('  expansion_ratios: [%s]\n', num2str(gt_coverage_traj.expansion_ratios'));
        fprintf('  coverage_points: [%s]\n', num2str(gt_coverage_traj.coverage_points'));
        fprintf('  mean_gt_rank: %.2f\n', gt_coverage_traj.mean_gt_rank);
    else
        fprintf('⚠️ gt_coverage_traj is EMPTY\n');
    end
else
    fprintf('✗ gt_coverage_traj DOES NOT EXIST\n');
end

% Check segment-level GT coverage
if exist('gt_coverage_seg_avg', 'var')
    fprintf('\n✓ gt_coverage_seg_avg EXISTS\n');
    if ~isempty(gt_coverage_seg_avg)
        fprintf('  Fields: %s\n', strjoin(fieldnames(gt_coverage_seg_avg), ', '));
        fprintf('  recall_at_k: [%s]\n', num2str(gt_coverage_seg_avg.recall_at_k'));
        fprintf('  expansion_ratios: [%s]\n', num2str(gt_coverage_seg_avg.expansion_ratios'));
    else
        fprintf('⚠️ gt_coverage_seg_avg is EMPTY\n');
    end
else
    fprintf('✗ gt_coverage_seg_avg DOES NOT EXIST\n');
end

% Check ground_truth_map
if exist('config', 'var') && isfield(config, 'ground_truth_map')
    fprintf('\n✓ config.ground_truth_map EXISTS\n');
    query_field = ['q_' strrep(query_bahn_id, '-', '_')];
    if isfield(config.ground_truth_map, query_field)
        gt_data = config.ground_truth_map.(query_field);
        fprintf('  Query field: %s\n', query_field);
        if isfield(gt_data, 'trajectories')
            fprintf('  GT Trajectories (%d):\n', length(gt_data.trajectories));
            for i = 1:min(3, length(gt_data.trajectories))
                fprintf('    %d: %s\n', i, gt_data.trajectories{i});
            end
        end
    else
        fprintf('⚠️ Query field "%s" NOT FOUND in ground_truth_map\n', query_field);
        fprintf('  Available fields: %s\n', strjoin(fieldnames(config.ground_truth_map), ', '));
    end
else
    fprintf('✗ config.ground_truth_map DOES NOT EXIST\n');
end

% Check embedding_table
if exist('embedding_table', 'var')
    fprintf('\n✓ embedding_table EXISTS (%d rows)\n', height(embedding_table));
    fprintf('  First 5 bahn_ids:\n');
    for i = 1:min(5, height(embedding_table))
        fprintf('    %d: %s\n', i, embedding_table.bahn_id{i});
    end
else
    fprintf('✗ embedding_table DOES NOT EXIST\n');
end

fprintf('════════════════════════════════════════════════════════════════\n\n');

%% ========================================================================
%  EXTRACT RESULTS
%  ========================================================================

results = struct();

% ========================================================================
% Trajectory-Level Metrics
% ========================================================================

results.spearman = rho_spearman;
results.p_at_k = prec_k;

if exist('prec_10', 'var')
    results.p_at_10 = prec_10;
else
    results.p_at_10 = NaN;
end

if exist('prec_5', 'var')
    results.p_at_5 = prec_5;
else
    results.p_at_5 = NaN;
end

if exist('prec_3', 'var')
    results.p_at_3 = prec_3;
else
    results.p_at_3 = NaN;
end

results.p_at_1 = prec_1;

% ========================================================================
% Segment-Level Metrics
% ========================================================================

if num_query_segments > 0 && exist('seg_rho_all', 'var') && ~isempty(seg_rho_all)
    results.seg_spearman = mean(seg_rho_all);
    results.seg_p_at_k = mean(seg_prec_k_all);
    
    if ~isempty(seg_prec_10_all)
        results.seg_p_at_10 = mean(seg_prec_10_all);
    else
        results.seg_p_at_10 = NaN;
    end
    
    if ~isempty(seg_prec_5_all)
        results.seg_p_at_5 = mean(seg_prec_5_all);
    else
        results.seg_p_at_5 = NaN;
    end
    
    if ~isempty(seg_prec_3_all)
        results.seg_p_at_3 = mean(seg_prec_3_all);
    else
        results.seg_p_at_3 = NaN;
    end
    
    if ~isempty(seg_prec_1_all)
        results.seg_p_at_1 = mean(seg_prec_1_all);
    else
        results.seg_p_at_1 = NaN;
    end
else
    % No segments or empty results
    results.seg_spearman = NaN;
    results.seg_p_at_k = NaN;
    results.seg_p_at_10 = NaN;
    results.seg_p_at_5 = NaN;
    results.seg_p_at_3 = NaN;
    results.seg_p_at_1 = NaN;
end

% ========================================================================
% Configuration Metadata
% ========================================================================

results.dtw_mode = dtw_mode;
results.weight_pos = weights(1);
results.weight_joint = weights(2);
results.weight_orient = weights(3);
results.weight_vel = weights(4);
results.weight_meta = weights(5);

results.n_coarse = n_coarse;
results.n_fine = n_fine;
results.total_dims = n_coarse + n_fine;
results.use_multi_scale = use_multi_scale;

results.top_k = top_k_trajectories;
results.rrf_k = rrf_k;

% ========================================================================
% Additional Info
% ========================================================================

results.num_candidates = num_candidates;
results.num_query_segments = num_query_segments;

% ========================================================================
% Coverage Metrics (Trajectory Level) - EXISTING, keeping as is
% ========================================================================

if exist('coverage_traj', 'var') && ~isempty(coverage_traj)
    % Find indices for K=10, K=50, K=5, K=3, K=1
    idx_10 = find(coverage_traj.k_values == 10, 1);
    idx_50 = find(coverage_traj.k_values == 50, 1);
    idx_5 = find(coverage_traj.k_values == 5, 1);   % NEW
    idx_3 = find(coverage_traj.k_values == 3, 1);   % NEW
    idx_1 = find(coverage_traj.k_values == 1, 1);   
    
    % P_all: Maximum coverage point across all K values
    % This is the minimum embedding top-X needed to cover ALL DTW candidates
    results.p_all_traj = max(coverage_traj.coverage_points);
    
    if ~isempty(idx_10)
        results.recall_at_10 = coverage_traj.recall_at_k(idx_10);
    else
        results.recall_at_10 = NaN;
    end
    
    if ~isempty(idx_50)
        results.recall_at_50 = coverage_traj.recall_at_k(idx_50);
    else
        results.recall_at_50 = NaN;
    end
    
    % NEW: Add R@5, R@3, R@1
    if ~isempty(idx_5)
        results.recall_at_5 = coverage_traj.recall_at_k(idx_5);
    else
        results.recall_at_5 = NaN;
    end
    
    if ~isempty(idx_3)
        results.recall_at_3 = coverage_traj.recall_at_k(idx_3);
    else
        results.recall_at_3 = NaN;
    end
    
    if ~isempty(idx_1)
        results.recall_at_1 = coverage_traj.recall_at_k(idx_1);
    else
        results.recall_at_1 = NaN;
    end
else
    % No coverage data available
    results.p_all_traj = NaN;
    results.recall_at_10 = NaN;
    results.recall_at_50 = NaN;
    results.recall_at_5 = NaN;   
    results.recall_at_3 = NaN;   
    results.recall_at_1 = NaN;   
end

% ========================================================================
% Coverage Metrics (Segment Level) - EXISTING, keeping as is
% ========================================================================

if exist('coverage_seg_avg', 'var') && ~isempty(coverage_seg_avg)
    % Find indices for K=10, K=50, K=5, K=3, K=1
    idx_10 = find(coverage_seg_avg.k_values == 10, 1);
    idx_50 = find(coverage_seg_avg.k_values == 50, 1);
    idx_5 = find(coverage_seg_avg.k_values == 5, 1);  
    idx_3 = find(coverage_seg_avg.k_values == 3, 1);   
    idx_1 = find(coverage_seg_avg.k_values == 1, 1);   
    
    % P_all: Maximum coverage point across all K values
    results.p_all_seg = max(coverage_seg_avg.coverage_points);
    
    if ~isempty(idx_10)
        results.seg_recall_at_10 = coverage_seg_avg.recall_at_k(idx_10);
    else
        results.seg_recall_at_10 = NaN;
    end
    
    if ~isempty(idx_50)
        results.seg_recall_at_50 = coverage_seg_avg.recall_at_k(idx_50);
    else
        results.seg_recall_at_50 = NaN;
    end
    
    % NEW: Add R@5, R@3, R@1
    if ~isempty(idx_5)
        results.seg_recall_at_5 = coverage_seg_avg.recall_at_k(idx_5);
    else
        results.seg_recall_at_5 = NaN;
    end
    
    if ~isempty(idx_3)
        results.seg_recall_at_3 = coverage_seg_avg.recall_at_k(idx_3);
    else
        results.seg_recall_at_3 = NaN;
    end
    
    if ~isempty(idx_1)
        results.seg_recall_at_1 = coverage_seg_avg.recall_at_k(idx_1);
    else
        results.seg_recall_at_1 = NaN;
    end
else
    % No segment coverage data available
    results.p_all_seg = NaN;
    results.seg_recall_at_10 = NaN;
    results.seg_recall_at_50 = NaN;
    results.seg_recall_at_5 = NaN;   
    results.seg_recall_at_3 = NaN;   
    results.seg_recall_at_1 = NaN;   
end

% ========================================================================
% GT Coverage Metrics (Trajectory Level) - NEW!
% ========================================================================

if exist('gt_coverage_traj', 'var') && ~isempty(gt_coverage_traj)
    % Find indices for K=10, K=50, K=5, K=3, K=1
    idx_10 = find(gt_coverage_traj.k_values == 10, 1);
    idx_50 = find(gt_coverage_traj.k_values == 50, 1);
    idx_5 = find(gt_coverage_traj.k_values == 5, 1);  
    idx_3 = find(gt_coverage_traj.k_values == 3, 1);   
    idx_1 = find(gt_coverage_traj.k_values == 1, 1);   
    
    % P_GT: Overall coverage point for ALL GT trajectories
    results.p_gt_traj = gt_coverage_traj.overall_coverage_point;
    
    if ~isempty(idx_10)
        results.gt_recall_at_10 = gt_coverage_traj.recall_at_k(idx_10);
    else
        results.gt_recall_at_10 = NaN;
    end
    
    if ~isempty(idx_50)
        results.gt_recall_at_50 = gt_coverage_traj.recall_at_k(idx_50);
    else
        results.gt_recall_at_50 = NaN;
    end
    
    % NEW: Add R@5, R@3, R@1
    if ~isempty(idx_5)
        results.gt_recall_at_5 = gt_coverage_traj.recall_at_k(idx_5);
    else
        results.gt_recall_at_5 = NaN;
    end
    
    if ~isempty(idx_3)
        results.gt_recall_at_3 = gt_coverage_traj.recall_at_k(idx_3);
    else
        results.gt_recall_at_3 = NaN;
    end
    
    if ~isempty(idx_1)
        results.gt_recall_at_1 = gt_coverage_traj.recall_at_k(idx_1);
    else
        results.gt_recall_at_1 = NaN;
    end
    
    % Additional GT statistics
    results.num_gt = gt_coverage_traj.num_gt;
    results.num_gt_found = gt_coverage_traj.num_gt_found;
    results.mean_gt_rank = gt_coverage_traj.mean_gt_rank;
else
    % No GT coverage data available
    results.p_gt_traj = NaN;
    results.gt_recall_at_10 = NaN;
    results.gt_recall_at_50 = NaN;
    results.gt_recall_at_5 = NaN;   
    results.gt_recall_at_3 = NaN;   
    results.gt_recall_at_1 = NaN;   
    
    % Additional GT statistics
    results.num_gt = NaN;
    results.num_gt_found = NaN;
    results.mean_gt_rank = NaN;
end

% ========================================================================
% NDCG METRICS (Trajectory Level) - NEW!
% ========================================================================

if exist('gt_coverage_traj', 'var') && ~isempty(gt_coverage_traj) && gt_coverage_traj.num_gt > 0
    fprintf('  Computing NDCG metrics (trajectory)...\n');
    
    % Get GT IDs from ground_truth_map
    query_field = sprintf('q_%s', strrep(query_bahn_id, '-', '_'));
    
    if isfield(config, 'ground_truth_map') && isfield(config.ground_truth_map, query_field)
        gt_ids = config.ground_truth_map.(query_field).trajectories;
        
        % Extract GT similarities from DTW cache
        gt_similarities = extractGTSimilarities(...
            config.dtw_cache, query_bahn_id, config.dtw_mode, gt_ids);
        
        if ~isempty(gt_similarities)
            % Get top-K IDs from embedding ranking
            top_k_ids = fused_ranking.bahn_id;
            K = config.top_k_trajectories;  % Usually 100
            
            % Compute NDCG metrics using existing function
            embedding_metrics = computeEmbeddingOnlyMetrics(...
                top_k_ids, gt_ids, gt_similarities, gt_coverage_traj.num_gt, K);
            
            % Store in results
            results.ndcg_10_traj = embedding_metrics.ndcg10;
            results.ndcg_50_traj = embedding_metrics.ndcg50;
            results.mrr_traj = embedding_metrics.mrr;
            
            fprintf('    NDCG@10 (Traj): %.4f\n', results.ndcg_10_traj);
            fprintf('    NDCG@50 (Traj): %.4f\n', results.ndcg_50_traj);
            fprintf('    MRR (Traj):     %.4f\n', results.mrr_traj);
        else
            results.ndcg_10_traj = NaN;
            results.ndcg_50_traj = NaN;
            results.mrr_traj = NaN;
        end
    else
        results.ndcg_10_traj = NaN;
        results.ndcg_50_traj = NaN;
        results.mrr_traj = NaN;
    end
else
    results.ndcg_10_traj = NaN;
    results.ndcg_50_traj = NaN;
    results.mrr_traj = NaN;
end


% ========================================================================
% GT Coverage Metrics (Segment Level) - NEW!
% ========================================================================

if exist('gt_coverage_seg_avg', 'var') && ~isempty(gt_coverage_seg_avg)
    % Find indices for K=10, K=50, K=5, K=3, K=1
    idx_10 = find(gt_coverage_seg_avg.k_values == 10, 1);
    idx_50 = find(gt_coverage_seg_avg.k_values == 50, 1);
    idx_5 = find(gt_coverage_seg_avg.k_values == 5, 1);  
    idx_3 = find(gt_coverage_seg_avg.k_values == 3, 1);   
    idx_1 = find(gt_coverage_seg_avg.k_values == 1, 1);   
    
    % P_GT: Overall coverage point for ALL GT segments
    results.p_gt_seg = gt_coverage_seg_avg.overall_coverage_point;
    
    if ~isempty(idx_10)
        results.seg_gt_recall_at_10 = gt_coverage_seg_avg.recall_at_k(idx_10);
    else
        results.seg_gt_recall_at_10 = NaN;
    end
    
    if ~isempty(idx_50)
        results.seg_gt_recall_at_50 = gt_coverage_seg_avg.recall_at_k(idx_50);
    else
        results.seg_gt_recall_at_50 = NaN;
    end
    
    % NEW: Add R@5, R@3, R@1
    if ~isempty(idx_5)
        results.seg_gt_recall_at_5 = gt_coverage_seg_avg.recall_at_k(idx_5);
    else
        results.seg_gt_recall_at_5 = NaN;
    end
    
    if ~isempty(idx_3)
        results.seg_gt_recall_at_3 = gt_coverage_seg_avg.recall_at_k(idx_3);
    else
        results.seg_gt_recall_at_3 = NaN;
    end
    
    if ~isempty(idx_1)
        results.seg_gt_recall_at_1 = gt_coverage_seg_avg.recall_at_k(idx_1);
    else
        results.seg_gt_recall_at_1 = NaN;
    end
    
    % Additional segment GT statistics
    results.seg_mean_gt_rank = gt_coverage_seg_avg.mean_gt_rank;
else
    % No segment GT coverage data available
    results.p_gt_seg = NaN;
    results.seg_gt_recall_at_10 = NaN;
    results.seg_gt_recall_at_50 = NaN;
    results.seg_gt_recall_at_5 = NaN;  
    results.seg_gt_recall_at_3 = NaN;  
    results.seg_gt_recall_at_1 = NaN;  
    
    % Additional segment GT statistics
    results.seg_mean_gt_rank = NaN;
end

% ========================================================================
% NDCG METRICS (Segment Level) - NEW!
% ========================================================================

if exist('gt_coverage_seg_avg', 'var') && ~isempty(gt_coverage_seg_avg)
    fprintf('  Computing NDCG metrics (segments)...\n');
    
    query_field = sprintf('q_%s', strrep(query_bahn_id, '-', '_'));
    
    % Check if we have segment GT data
    if isfield(config, 'ground_truth_map') && ...
       isfield(config.ground_truth_map, query_field) && ...
       isfield(config.ground_truth_map.(query_field), 'segments')
        
        num_query_segments = length(segment_embedding_results);
        gt_seg_struct = config.ground_truth_map.(query_field).segments;
        seg_fields = fieldnames(gt_seg_struct);
        
        ndcg_10_segments = [];
        ndcg_50_segments = [];
        mrr_segments = [];
        
        for seg_idx = 1:num_query_segments
            % Get segment ranking
            seg_fused_ranking = segment_embedding_results{seg_idx};
            
            if isempty(seg_fused_ranking)
                continue;
            end
            
            % Get GT segment IDs for this query segment
            if seg_idx <= length(seg_fields)
                gt_seg_ids_curr = gt_seg_struct.(seg_fields{seg_idx});
                
                if ~iscell(gt_seg_ids_curr)
                    gt_seg_ids_curr = {gt_seg_ids_curr};
                end
                
                num_gt_seg = length(gt_seg_ids_curr);
                
                if num_gt_seg > 0
                    % Extract segment similarities from DTW cache
                    gt_seg_similarities = extractSegmentSimilarities(...
                        config.dtw_cache, query_bahn_id, config.dtw_mode, ...
                        seg_idx, gt_seg_ids_curr);
                    
                    if ~isempty(gt_seg_similarities)
                        % Compute NDCG for this segment
                        top_k_seg_ids = seg_fused_ranking.segment_id;
                        K_seg = config.top_k_trajectories;
                        
                        seg_metrics = computeEmbeddingOnlyMetrics(...
                            top_k_seg_ids, gt_seg_ids_curr, gt_seg_similarities, num_gt_seg, K_seg);
                        
                        ndcg_10_segments = [ndcg_10_segments; seg_metrics.ndcg10];
                        ndcg_50_segments = [ndcg_50_segments; seg_metrics.ndcg50];
                        mrr_segments = [mrr_segments; seg_metrics.mrr];
                    end
                end
            end
        end
        
        % Average across valid segments
        if ~isempty(ndcg_10_segments)
            results.ndcg_10_seg = mean(ndcg_10_segments);
            results.ndcg_50_seg = mean(ndcg_50_segments);
            results.mrr_seg = mean(mrr_segments);
            
            fprintf('    NDCG@10 (Seg Avg): %.4f\n', results.ndcg_10_seg);
            fprintf('    NDCG@50 (Seg Avg): %.4f\n', results.ndcg_50_seg);
            fprintf('    MRR (Seg Avg):     %.4f\n', results.mrr_seg);
        else
            results.ndcg_10_seg = NaN;
            results.ndcg_50_seg = NaN;
            results.mrr_seg = NaN;
        end
    else
        results.ndcg_10_seg = NaN;
        results.ndcg_50_seg = NaN;
        results.mrr_seg = NaN;
    end
else
    results.ndcg_10_seg = NaN;
    results.ndcg_50_seg = NaN;
    results.mrr_seg = NaN;
end

end
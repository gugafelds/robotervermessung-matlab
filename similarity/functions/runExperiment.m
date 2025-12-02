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

% Suppress output from dtw_baseline.m (we only want results)
% Use evalc to capture output but not display it
evalc('dtw_baseline');

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
% Coverage Metrics (Trajectory Level)
% ========================================================================

if exist('coverage_traj', 'var') && ~isempty(coverage_traj)
    % Find indices for K=10 and K=50
    idx_10 = find(coverage_traj.k_values == 10, 1);
    idx_50 = find(coverage_traj.k_values == 50, 1);
    
    if ~isempty(idx_10)
        results.coverage_at_10 = coverage_traj.coverage_points(idx_10);
        results.expansion_ratio_10 = coverage_traj.expansion_ratios(idx_10);
        results.recall_at_10 = coverage_traj.recall_at_k(idx_10);
    else
        results.coverage_at_10 = NaN;
        results.expansion_ratio_10 = NaN;
        results.recall_at_10 = NaN;
    end
    
    if ~isempty(idx_50)
        results.coverage_at_50 = coverage_traj.coverage_points(idx_50);
        results.expansion_ratio_50 = coverage_traj.expansion_ratios(idx_50);
        results.recall_at_50 = coverage_traj.recall_at_k(idx_50);
    else
        results.coverage_at_50 = NaN;
        results.expansion_ratio_50 = NaN;
        results.recall_at_50 = NaN;
    end
else
    % No coverage data available
    results.coverage_at_10 = NaN;
    results.expansion_ratio_10 = NaN;
    results.recall_at_10 = NaN;
    results.coverage_at_50 = NaN;
    results.expansion_ratio_50 = NaN;
    results.recall_at_50 = NaN;
end

% ========================================================================
% Coverage Metrics (Segment Level)
% ========================================================================

if exist('coverage_seg_avg', 'var') && ~isempty(coverage_seg_avg)
    % Find indices for K=10 and K=50
    idx_10 = find(coverage_seg_avg.k_values == 10, 1);
    idx_50 = find(coverage_seg_avg.k_values == 50, 1);
    
    if ~isempty(idx_10) && coverage_seg_avg.expansion_ratios(idx_10) > 0
        results.seg_coverage_at_10 = coverage_seg_avg.coverage_points(idx_10);
        results.seg_expansion_ratio_10 = coverage_seg_avg.expansion_ratios(idx_10);
        results.seg_recall_at_10 = coverage_seg_avg.recall_at_k(idx_10);
    else
        results.seg_coverage_at_10 = NaN;
        results.seg_expansion_ratio_10 = NaN;
        results.seg_recall_at_10 = NaN;
    end
    
    if ~isempty(idx_50) && coverage_seg_avg.expansion_ratios(idx_50) > 0
        results.seg_coverage_at_50 = coverage_seg_avg.coverage_points(idx_50);
        results.seg_expansion_ratio_50 = coverage_seg_avg.expansion_ratios(idx_50);
        results.seg_recall_at_50 = coverage_seg_avg.recall_at_k(idx_50);
    else
        results.seg_coverage_at_50 = NaN;
        results.seg_expansion_ratio_50 = NaN;
        results.seg_recall_at_50 = NaN;
    end
else
    % No segment coverage data available
    results.seg_coverage_at_10 = NaN;
    results.seg_expansion_ratio_10 = NaN;
    results.seg_recall_at_10 = NaN;
    results.seg_coverage_at_50 = NaN;
    results.seg_expansion_ratio_50 = NaN;
    results.seg_recall_at_50 = NaN;
end

end
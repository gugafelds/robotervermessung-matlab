%  TWO-STAGE RETRIEVAL - SINGLE QUERY (Web App Style)
%  ========================================================================
%  Goal: Demonstrate production-ready Two-Stage Retrieval for ONE query
%  
%  User Input:
%  - Query Bahn ID (e.g., '1765989370')
%  - K candidates (Stage 1 output size)
%  - Embedding weights (position, joint, orientation, velocity, metadata)
%  - DTW mode (position or joint_states)
%  
%  Pipeline (like your Web App):
%  1. BAHN-LEVEL:
%     - Stage 1: RRF Embedding Fusion → Top-K Bahnen
%     - Stage 2: DTW + LB Reranking → Top-N Bahnen
%  
%  2. SEGMENT-LEVEL (for each query segment):
%     - Stage 1: RRF Embedding Fusion → Top-K Segments per segment
%     - Stage 2: DTW + LB Reranking → Top-N Segments per segment
%  
%  Output:
%  - Similar Bahnen with DTW distances
%  - Similar Segments per Query Segment with DTW distances
%  - Performance Metrics (Stage 1 time, Stage 2 time, DTW calls, speedup)
%  
%  Story: "Like the Web App, but with detailed performance metrics"
%  ========================================================================

clear; clc;

addpath(genpath(pwd));
addpath(genpath('../main'));
addpath(genpath('../lasertracker'));
addpath(genpath('../methods'));

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  TWO-STAGE RETRIEVAL - SINGLE QUERY (Web App Style)            ║\n');
fprintf('║  Hierarchical: Bahn → Segments                                 ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

% ========================================================================
%% USER CONFIGURATION
% ========================================================================

fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('USER CONFIGURATION\n');
fprintf('═══════════════════════════════════════════════════════════════\n');

% === QUERY ===
query_id = '1766066754';  % ← CHANGE THIS!

% === STAGE 1 PARAMETERS ===
K_bahn = 500;       % Top-K candidates from Stage 1 (Bahn-Level)
K_segment = 500;    % Top-K candidates per segment (Segment-Level)

% === EMBEDDING WEIGHTS (adjust as you like!) ===
embedding_weights = [
    1.0;  % position
    1.0;  % joint
    1.0;  % orientation
    1.0;  % velocity
    1.0;  % metadata
];

% Normalize weights
embedding_weights = embedding_weights / sum(embedding_weights);

% === DTW MODE ===
dtw_mode = 'position';  % 'position' or 'joint_states'

% === FINAL OUTPUT SIZE ===
final_limit_bahn = 250;      % Final Top-N Bahnen
final_limit_segment = 250;   % Final Top-N Segments per segment

% === DTW CONFIG ===
dtw_config = struct();
dtw_config.cdtw_window = 0.2;
dtw_config.lb_kim_keep_ratio = 0.9;
dtw_config.lb_keogh_candidates = 500;
dtw_config.normalize_dtw = false;
dtw_config.use_rotation_alignment = false;

fprintf('  Query ID:           %s\n', query_id);
fprintf('  K (Bahn-Level):     %d\n', K_bahn);
fprintf('  K (Segment-Level):  %d\n', K_segment);
fprintf('  DTW Mode:           %s\n', dtw_mode);
fprintf('  Final Limit Bahn:   %d\n', final_limit_bahn);
fprintf('  Final Limit Seg:    %d\n\n', final_limit_segment);

fprintf('Embedding Weights:\n');
fprintf('  Position:    %.2f\n', embedding_weights(1));
fprintf('  Joint:       %.2f\n', embedding_weights(2));
fprintf('  Orientation: %.2f\n', embedding_weights(3));
fprintf('  Velocity:    %.2f\n', embedding_weights(4));
fprintf('  Metadata:    %.2f\n\n', embedding_weights(5));

fprintf('DTW Configuration:\n');
fprintf('  Window:         %.1f (%.0f%%)\n', dtw_config.cdtw_window, dtw_config.cdtw_window*100);
fprintf('  LB_Kim ratio:   %.1f%%\n', dtw_config.lb_kim_keep_ratio*100);
fprintf('  LB_Keogh target: %d\n', dtw_config.lb_keogh_candidates);
fprintf('═══════════════════════════════════════════════════════════════\n\n');

% ========================================================================
%% SECTION 1: DATABASE & DATA LOADING
% ========================================================================

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  SECTION 1: DATABASE & DATA LOADING                            ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

% === Connect to Database ===
fprintf('Connecting to database...\n');
conn = connectingToPostgres();
if ~isopen(conn)
    error('Database connection failed');
end
fprintf('✓ Connected\n\n');

schema = 'bewegungsdaten';

% === Load Query Data ===
fprintf('=== Loading Query Data ===\n');
fprintf('  Query ID: %s\n', query_id);

load_start = tic;
query_data = loadQueryDataHierarchical(conn, schema, query_id);
load_time = toc(load_start);

if isempty(query_data)
    error('Query not found: %s', query_id);
end

num_query_segments = length(query_data.segment_ids);

fprintf('  ✓ Query loaded in %.2f sec\n', load_time);
fprintf('    Trajectory points: %d\n', size(query_data.trajectory.position, 1));
fprintf('    Number of segments: %d\n', num_query_segments);
fprintf('\n');

% ========================================================================
%% SECTION 2: BAHN-LEVEL TWO-STAGE
% ========================================================================

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  SECTION 2: BAHN-LEVEL TWO-STAGE RETRIEVAL                     ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

bahn_pipeline_start = tic;

% === STAGE 1: RRF Embedding Fusion ===
fprintf('=== Stage 1: RRF Embedding Fusion (Bahnen) ===\n');

stage1_bahn_start = tic;

% Get query embeddings
query_embeddings = getQueryEmbeddings(conn, schema, query_id, 'bahn');

% Perform RRF on all embedding modalities
bahn_candidates = performRRFEmbeddingSearch(conn, schema, query_embeddings, ...
    embedding_weights, K_bahn, 'bahn');

stage1_bahn_time = toc(stage1_bahn_start);

fprintf('  ✓ Stage 1 completed in %.2f sec\n', stage1_bahn_time);
fprintf('    Candidates retrieved: %d\n', length(bahn_candidates));
fprintf('\n');

if isempty(bahn_candidates)
    error('No bahn candidates found!');
end

% === STAGE 2: DTW Reranking ===
fprintf('=== Stage 2: DTW Reranking with Lower Bounds (Bahnen) ===\n');

stage2_bahn_start = tic;

% Load candidate data
fprintf('  Loading candidate data...\n');
load_cand_start = tic;
bahn_candidate_data = loadCandidateDataBatch(conn, schema, bahn_candidates, 'bahn');
load_cand_time = toc(load_cand_start);
fprintf('  ✓ Loaded %d candidates in %.2f sec\n', length(bahn_candidate_data), load_cand_time);

% Get query sequence
if strcmp(dtw_mode, 'position')
    query_seq = query_data.trajectory.position;
else
    query_seq = query_data.trajectory.joint;
end

% DTW Reranking
fprintf('  Running DTW with Lower Bounds...\n');
[bahn_results, bahn_dtw_stats] = performDTWReranking(query_seq, bahn_candidate_data, ...
    dtw_mode, dtw_config, final_limit_bahn);

stage2_bahn_time = toc(stage2_bahn_start);

fprintf('  ✓ Stage 2 completed in %.2f sec\n', stage2_bahn_time);
fprintf('    DTW calls made: %d / %d\n', bahn_dtw_stats.dtw_calls, K_bahn);
fprintf('    LB_Kim filtered: %d\n', bahn_dtw_stats.lb_kim_filtered);
fprintf('    LB_Keogh filtered: %d\n', bahn_dtw_stats.lb_keogh_filtered);
fprintf('    Final results: %d\n', length(bahn_results));
fprintf('\n');

bahn_pipeline_time = toc(bahn_pipeline_start);

% === Bahn-Level Metrics ===
fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('BAHN-LEVEL PERFORMANCE SUMMARY\n');
fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('  Stage 1 Time:       %.3f sec\n', stage1_bahn_time);
fprintf('  Stage 2 Time:       %.3f sec\n', stage2_bahn_time);
fprintf('  Total Time:         %.3f sec\n', bahn_pipeline_time);
fprintf('  DTW Calls:          %d / %d\n', bahn_dtw_stats.dtw_calls, K_bahn);
fprintf('  DTW Savings:        %.1f%%\n', (K_bahn - bahn_dtw_stats.dtw_calls) / K_bahn * 100);
fprintf('  Speedup:            %.1fx\n', K_bahn / bahn_dtw_stats.dtw_calls);
fprintf('═══════════════════════════════════════════════════════════════\n\n');

% ========================================================================
%% SECTION 3: SEGMENT-LEVEL TWO-STAGE
% ========================================================================

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  SECTION 3: SEGMENT-LEVEL TWO-STAGE RETRIEVAL                  ║\n');
fprintf('║  Processing %d query segments                                  ║\n', num_query_segments);
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

segment_pipeline_start = tic;

% Storage for all segment results
all_segment_results = cell(num_query_segments, 1);
segment_stage1_times = zeros(num_query_segments, 1);
segment_stage2_times = zeros(num_query_segments, 1);
segment_dtw_calls = zeros(num_query_segments, 1);

% === Process Each Query Segment ===
for seg_idx = 1:num_query_segments
    seg_id = query_data.segment_ids{seg_idx};
    
    fprintf('\n─────────────────────────────────────────────────────────────\n');
    fprintf('Segment %d/%d: %s\n', seg_idx, num_query_segments, seg_id);
    fprintf('─────────────────────────────────────────────────────────────\n');
    
    % === STAGE 1: RRF Embedding Fusion ===
    fprintf('  Stage 1: RRF Embedding Fusion...\n');
    
    stage1_seg_start = tic;
    
    % Get segment embeddings
    seg_embeddings = getQueryEmbeddings(conn, schema, seg_id, 'segment');
    
    % RRF Search
    seg_candidates = performRRFEmbeddingSearch(conn, schema, seg_embeddings, ...
        embedding_weights, K_segment, 'segment');
    
    segment_stage1_times(seg_idx) = toc(stage1_seg_start);
    
    fprintf('    ✓ %d candidates in %.3f sec\n', length(seg_candidates), segment_stage1_times(seg_idx));
    
    if isempty(seg_candidates)
        fprintf('    ⚠ No candidates - skipping\n');
        continue;
    end
    
    % === STAGE 2: DTW Reranking ===
    fprintf('  Stage 2: DTW Reranking...\n');
    
    stage2_seg_start = tic;
    
    % Load candidate data
    seg_candidate_data = loadCandidateDataBatch(conn, schema, seg_candidates, 'segment');
    
    % Get query segment sequence
    if strcmp(dtw_mode, 'position')
        query_seg_seq = query_data.segments.position{seg_idx};
    else
        query_seg_seq = query_data.segments.joint{seg_idx};
    end
    
    % DTW Reranking
    [seg_results, seg_dtw_stats] = performDTWReranking(query_seg_seq, seg_candidate_data, ...
        dtw_mode, dtw_config, final_limit_segment);
    
    segment_stage2_times(seg_idx) = toc(stage2_seg_start);
    segment_dtw_calls(seg_idx) = seg_dtw_stats.dtw_calls;
    
    fprintf('    ✓ %d results, %d DTW calls in %.3f sec\n', ...
        length(seg_results), seg_dtw_stats.dtw_calls, segment_stage2_times(seg_idx));
    
    % Store results
    all_segment_results{seg_idx} = struct();
    all_segment_results{seg_idx}.segment_id = seg_id;
    all_segment_results{seg_idx}.results = seg_results;
    all_segment_results{seg_idx}.stats = seg_dtw_stats;
end

segment_pipeline_time = toc(segment_pipeline_start);

% === Segment-Level Metrics (Aggregated) ===
fprintf('\n═══════════════════════════════════════════════════════════════\n');
fprintf('SEGMENT-LEVEL PERFORMANCE SUMMARY\n');
fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('  Segments Processed:     %d\n', num_query_segments);
fprintf('  Avg Stage 1 Time:       %.3f sec/segment\n', mean(segment_stage1_times));
fprintf('  Avg Stage 2 Time:       %.3f sec/segment\n', mean(segment_stage2_times));
fprintf('  Total Segment Time:     %.3f sec\n', segment_pipeline_time);
fprintf('  Avg DTW Calls:          %.1f / %d\n', mean(segment_dtw_calls), K_segment);
fprintf('  Avg DTW Savings:        %.1f%%\n', mean((K_segment - segment_dtw_calls) / K_segment * 100));
fprintf('═══════════════════════════════════════════════════════════════\n\n');

% Close connection
close(conn);

% ========================================================================
%% SECTION 3.5: RERANKING ANALYSIS (After all timing measurements!)
% ========================================================================

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  SECTION 3.5: RERANKING QUALITY ANALYSIS                       ║\n');
fprintf('║  (Computed AFTER performance measurements)                     ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

% === Analyze Bahn-Level Reranking ===
fprintf('=== Bahn-Level Reranking Analysis ===\n');

bahn_rerank_metrics = analyzeReranking(bahn_results, K_bahn);

fprintf('  ✓ Analysis completed\n');
fprintf('    Kendall Tau:     %.4f\n', bahn_rerank_metrics.kendall_tau);
fprintf('    Mean Rank Δ:     %.1f positions\n', bahn_rerank_metrics.mean_rank_change);
if isfield(bahn_rerank_metrics, 'top10_overlap_pct')
    fprintf('    Top-10 Overlap:  %.1f%%\n', bahn_rerank_metrics.top10_overlap_pct);
end
fprintf('\n');

% === Analyze Segment-Level Reranking ===
fprintf('=== Segment-Level Reranking Analysis ===\n');

segment_rerank_metrics = cell(num_query_segments, 1);

for seg_idx = 1:num_query_segments
    if ~isempty(all_segment_results{seg_idx})
        seg_results = all_segment_results{seg_idx}.results;
        segment_rerank_metrics{seg_idx} = analyzeReranking(seg_results, K_segment);
    end
end

fprintf('  ✓ Analysis completed for %d segments\n\n', num_query_segments);

% Close connection
close(conn);

% ========================================================================
%% SECTION 3.6: PERFORMANCE PROGNOSIS (SIDTW-based)
% ========================================================================

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  SECTION 3.6: PERFORMANCE PROGNOSIS EVALUATION                 ║\n');
fprintf('║  Can we predict query performance from similar trajectories?  ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

% Define Top-K values to test
top_k_test_values = [5, 10, 20, 50];

% Re-open connection
conn = connectingToPostgres();
if ~isopen(conn)
    error('Database reconnection failed');
end

% ================================================================
% BAHN-LEVEL: CREATE STAGE 1 AND STAGE 2 RESULT SETS
% ================================================================

fprintf('=== Bahn-Level: Performance Prognosis ===\n');

% ✅ STAGE 1 RESULTS: Sort by Stage 1 rank (RRF/Embedding order)
[~, s1_sort_idx] = sort([bahn_results.stage1_rank]);
bahn_results_stage1 = bahn_results(s1_sort_idx);

% ✅ STAGE 2 RESULTS: Already sorted by Stage 2 rank (DTW order)
bahn_results_stage2 = bahn_results;  % Already sorted by DTW

fprintf('  Created result sets:\n');
fprintf('    Stage 1: Top result is %s (Stage1 Rank #%d)\n', ...
    bahn_results_stage1(1).segment_id, bahn_results_stage1(1).stage1_rank);
fprintf('    Stage 2: Top result is %s (Stage2 Rank #%d)\n', ...
    bahn_results_stage2(1).segment_id, bahn_results_stage2(1).stage2_rank);
fprintf('\n');

% Stage 1 prognosis (using Stage 1 ordering)
fprintf('  Stage 1 (Embedding-based ranking)...\n');
bahn_prog_stage1 = evaluatePerformancePrognosis(conn, query_id, bahn_results_stage1, top_k_test_values);

% Stage 2 prognosis (using Stage 2 ordering)
fprintf('  Stage 2 (DTW-based ranking)...\n');
bahn_prog_stage2 = evaluatePerformancePrognosis(conn, query_id, bahn_results_stage2, top_k_test_values);

fprintf('  ✓ Prognosis completed\n\n');

% ================================================================
% SEGMENT-LEVEL: CREATE STAGE 1 AND STAGE 2 RESULT SETS
% ================================================================

fprintf('=== Segment-Level: Performance Prognosis ===\n');

segment_prog_stage1 = cell(num_query_segments, 1);
segment_prog_stage2 = cell(num_query_segments, 1);

for seg_idx = 1:num_query_segments
    if ~isempty(all_segment_results{seg_idx})
        seg_id = all_segment_results{seg_idx}.segment_id;
        seg_results = all_segment_results{seg_idx}.results;
        
        fprintf('  Segment %d/%d: %s\n', seg_idx, num_query_segments, seg_id);
        
        % ✅ STAGE 1: Sort by Stage 1 rank
        [~, s1_sort_idx_seg] = sort([seg_results.stage1_rank]);
        seg_results_stage1 = seg_results(s1_sort_idx_seg);
        
        % ✅ STAGE 2: Already sorted by Stage 2 rank
        seg_results_stage2 = seg_results;
        
        % Stage 1 prognosis
        segment_prog_stage1{seg_idx} = evaluatePerformancePrognosis(conn, seg_id, seg_results_stage1, top_k_test_values);
        
        % Stage 2 prognosis
        segment_prog_stage2{seg_idx} = evaluatePerformancePrognosis(conn, seg_id, seg_results_stage2, top_k_test_values);
    end
end

fprintf('  ✓ Prognosis completed for %d segments\n\n', num_query_segments);

% Close connection
close(conn);

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  SECTION 3.6 COMPLETE                                          ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

% ========================================================================
%% SECTION 4: RESULTS DISPLAY (WITH RERANKING INFO)
% ========================================================================

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  SECTION 4: RESULTS WITH RERANKING ANALYSIS                     ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

% === Display Top Bahnen WITH RERANKING ===
fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('TOP %d SIMILAR BAHNEN (Stage 1 → Stage 2)\n', min(10, length(bahn_results)));
fprintf('═══════════════════════════════════════════════════════════════\n');

for i = 1:min(10, length(bahn_results))
    r = bahn_results(i);
    
    % Rank change indicator
    if r.rank_change > 0
        change_str = sprintf('↑%d', r.rank_change);  % Improved (moved up)
    elseif r.rank_change < 0
        change_str = sprintf('↓%d', abs(r.rank_change));  % Degraded (moved down)
    else
        change_str = '═';  % Stable
    end
    
    fprintf('%2d. %s | DTW: %.4f | S1:#%d→S2:#%d %s\n', ...
        i, r.bahn_id, r.dtw_distance, r.stage1_rank, r.stage2_rank, change_str);
end

fprintf('\n');

% === Bahn Reranking Statistics ===
fprintf('─────────────────────────────────────────────────────────────\n');
fprintf('Bahn Reranking Statistics:\n');
fprintf('─────────────────────────────────────────────────────────────\n');
fprintf('  Rank Correlation:\n');
fprintf('    Kendall Tau:        %.4f\n', bahn_rerank_metrics.kendall_tau);
fprintf('    Spearman Rho:       %.4f\n', bahn_rerank_metrics.spearman_rho);
fprintf('\n');
fprintf('  Rank Changes:\n');
fprintf('    Mean:               %.1f ± %.1f positions\n', ...
    bahn_rerank_metrics.mean_rank_change, bahn_rerank_metrics.std_rank_change);
fprintf('    Median:             %.1f positions\n', bahn_rerank_metrics.median_rank_change);
fprintf('    Max:                %d positions\n', bahn_rerank_metrics.max_rank_change);
fprintf('\n');
fprintf('  Top-K Stability:\n');
if isfield(bahn_rerank_metrics, 'top5_overlap_pct')
    fprintf('    Top-5:              %d/%d (%.1f%%)\n', ...
        bahn_rerank_metrics.top5_overlap_count, 5, bahn_rerank_metrics.top5_overlap_pct);
end
if isfield(bahn_rerank_metrics, 'top10_overlap_pct')
    fprintf('    Top-10:             %d/%d (%.1f%%)\n', ...
        bahn_rerank_metrics.top10_overlap_count, 10, bahn_rerank_metrics.top10_overlap_pct);
end
if isfield(bahn_rerank_metrics, 'top20_overlap_pct')
    fprintf('    Top-20:             %d/%d (%.1f%%)\n', ...
        bahn_rerank_metrics.top20_overlap_count, 20, bahn_rerank_metrics.top20_overlap_pct);
end
fprintf('\n');
fprintf('  Reranking Impact:\n');
fprintf('    Improved:           %d (%.1f%%) ↑\n', ...
    bahn_rerank_metrics.num_improved, bahn_rerank_metrics.improvement_rate);
fprintf('    Degraded:           %d (%.1f%%) ↓\n', ...
    bahn_rerank_metrics.num_degraded, bahn_rerank_metrics.degradation_rate);
fprintf('    Stable:             %d (%.1f%%) ═\n', ...
    bahn_rerank_metrics.num_stable, bahn_rerank_metrics.stability_rate);
fprintf('\n');
fprintf('  Biggest Mover:\n');
fprintf('    ID:                 %s\n', bahn_rerank_metrics.biggest_mover_id);
fprintf('    Change:             #%d → #%d (Δ%d)\n', ...
    bahn_rerank_metrics.biggest_mover_stage1, ...
    bahn_rerank_metrics.biggest_mover_stage2, ...
    bahn_rerank_metrics.biggest_mover_change);
fprintf('═══════════════════════════════════════════════════════════════\n\n');

% === Display Top Segments (first 3 query segments) ===
fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('TOP SIMILAR SEGMENTS (first 3 query segments)\n');
fprintf('═══════════════════════════════════════════════════════════════\n');

for seg_idx = 1:min(3, num_query_segments)
    if isempty(all_segment_results{seg_idx})
        continue;
    end
    
    fprintf('\nQuery Segment %d/%d: %s\n', seg_idx, num_query_segments, ...
        all_segment_results{seg_idx}.segment_id);
    fprintf('─────────────────────────────────────────────────────────────\n');
    
    seg_res = all_segment_results{seg_idx}.results;
    
    for i = 1:min(5, length(seg_res))
        r = seg_res(i);
        
        % Rank change
        if r.rank_change > 0
            change_str = sprintf('↑%d', r.rank_change);
        elseif r.rank_change < 0
            change_str = sprintf('↓%d', abs(r.rank_change));
        else
            change_str = '═';
        end
        
        fprintf('  %2d. %s | DTW: %.4f | S1:#%d→S2:#%d %s\n', ...
            i, r.segment_id, r.dtw_distance, r.stage1_rank, r.stage2_rank, change_str);
    end
    
    % Mini reranking stats for this segment
    if ~isempty(segment_rerank_metrics{seg_idx})
        rm = segment_rerank_metrics{seg_idx};
        fprintf('  Reranking: τ=%.3f, ΔRank=%.1f, Top10=%.0f%%\n', ...
            rm.kendall_tau, rm.mean_rank_change, ...
            rm.top10_overlap_pct);
    end
end

fprintf('\n');

% === Aggregated Segment Reranking Stats ===
fprintf('─────────────────────────────────────────────────────────────\n');
fprintf('Segment Reranking Statistics (Aggregated):\n');
fprintf('─────────────────────────────────────────────────────────────\n');

% Collect metrics
seg_taus = [];
seg_rhos = [];
seg_mean_changes = [];
seg_top10_overlaps = [];

for seg_idx = 1:num_query_segments
    if ~isempty(segment_rerank_metrics{seg_idx})
        rm = segment_rerank_metrics{seg_idx};
        if ~isnan(rm.kendall_tau)
            seg_taus = [seg_taus; rm.kendall_tau];
        end
        if ~isnan(rm.spearman_rho)
            seg_rhos = [seg_rhos; rm.spearman_rho];
        end
        if ~isnan(rm.mean_rank_change)
            seg_mean_changes = [seg_mean_changes; rm.mean_rank_change];
        end
        if isfield(rm, 'top10_overlap_pct')
            seg_top10_overlaps = [seg_top10_overlaps; rm.top10_overlap_pct];
        end
    end
end

if ~isempty(seg_taus)
    fprintf('  Rank Correlation (Avg over %d segments):\n', length(seg_taus));
    fprintf('    Kendall Tau:        %.4f ± %.4f\n', mean(seg_taus), std(seg_taus));
    fprintf('    Spearman Rho:       %.4f ± %.4f\n', mean(seg_rhos), std(seg_rhos));
    fprintf('\n');
    fprintf('  Rank Changes (Avg):\n');
    fprintf('    Mean:               %.1f ± %.1f positions\n', ...
        mean(seg_mean_changes), std(seg_mean_changes));
    fprintf('\n');
    if ~isempty(seg_top10_overlaps)
        fprintf('  Top-10 Overlap (Avg):\n');
        fprintf('    Overlap:            %.1f%% ± %.1f%%\n', ...
            mean(seg_top10_overlaps), std(seg_top10_overlaps));
    end
else
    fprintf('  No valid reranking metrics\n');
end

fprintf('═══════════════════════════════════════════════════════════════\n\n');



% ========================================================================
%% SECTION 5: PERFORMANCE PROGNOSIS RESULTS
% ========================================================================

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  SECTION 5: PERFORMANCE PROGNOSIS RESULTS                      ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

% === Bahn-Level Prognosis Display ===
fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('BAHN-LEVEL PERFORMANCE PROGNOSIS\n');
fprintf('═══════════════════════════════════════════════════════════════\n\n');

if ~isempty(bahn_prog_stage1) && isfield(bahn_prog_stage1, 'actual')
    
    % Display actual values
    fprintf('┌─────────────────────────────────────────────────────────────┐\n');
    fprintf('│ ACTUAL PERFORMANCE (Query Bahn: %s)              │\n', query_id);
    fprintf('├─────────────────────────────────────────────────────────────┤\n');
    fprintf('│ SIDTW Metrics:                                              │\n');
    fprintf('│   Min Distance:         %8.4f mm                          │\n', bahn_prog_stage1.actual.min_distance);
    fprintf('│   Max Distance:         %8.4f mm                          │\n', bahn_prog_stage1.actual.max_distance);
    fprintf('│   Average Distance:     %8.4f mm                          │\n', bahn_prog_stage1.actual.average_distance);
    fprintf('│   Std Deviation:        %8.4f mm                          │\n', bahn_prog_stage1.actual.standard_deviation);
    fprintf('└─────────────────────────────────────────────────────────────┘\n\n');
    
    % Display predictions for each Top-K
    k_fields = fieldnames(bahn_prog_stage1.predictions);
    
    for k_idx = 1:length(k_fields)
        k_field = k_fields{k_idx};
        pred_data = bahn_prog_stage1.predictions.(k_field);
        K = pred_data.k;
        
        fprintf('─────────────────────────────────────────────────────────────\n');
        fprintf('PREDICTIONS USING TOP-%d SIMILAR BAHNEN\n', K);
        fprintf('─────────────────────────────────────────────────────────────\n\n');
        
        % Stage 1
        if isfield(bahn_prog_stage1.predictions, k_field)
            s1_pred = bahn_prog_stage1.predictions.(k_field);
            s1_err = s1_pred.errors;
            
            fprintf('Stage 1 (Embedding-based):\n');
            fprintf('  Average Distance:\n');
            fprintf('    Predicted (Mean):   %8.4f mm (Error: %.4f mm, %.1f%%)\n', ...
                s1_err.average_distance.predicted_mean, ...
                s1_err.average_distance.mae_mean, ...
                s1_err.average_distance.relative_error_mean);
            fprintf('    Predicted (Median): %8.4f mm (Error: %.4f mm, %.1f%%)\n', ...
                s1_err.average_distance.predicted_median, ...
                s1_err.average_distance.mae_median, ...
                s1_err.average_distance.relative_error_median);
            fprintf('\n');
            
            fprintf('  Min Distance:\n');
            fprintf('    Predicted (Mean):   %8.4f mm (Error: %.4f mm)\n', ...
                s1_err.min_distance.predicted_mean, s1_err.min_distance.mae_mean);
            fprintf('\n');
            
            fprintf('  Max Distance:\n');
            fprintf('    Predicted (Mean):   %8.4f mm (Error: %.4f mm)\n', ...
                s1_err.max_distance.predicted_mean, s1_err.max_distance.mae_mean);
            fprintf('\n');
            
            fprintf('  Std Deviation:\n');
            fprintf('    Predicted (Mean):   %8.4f mm (Error: %.4f mm)\n', ...
                s1_err.standard_deviation.predicted_mean, s1_err.standard_deviation.mae_mean);
            fprintf('\n');
        end
        
        % Stage 2
        if isfield(bahn_prog_stage2.predictions, k_field)
            s2_pred = bahn_prog_stage2.predictions.(k_field);
            s2_err = s2_pred.errors;
            
            fprintf('Stage 2 (DTW-based):\n');
            fprintf('  Average Distance:\n');
            fprintf('    Predicted (Mean):   %8.4f mm (Error: %.4f mm, %.1f%%)\n', ...
                s2_err.average_distance.predicted_mean, ...
                s2_err.average_distance.mae_mean, ...
                s2_err.average_distance.relative_error_mean);
            fprintf('    Predicted (Median): %8.4f mm (Error: %.4f mm, %.1f%%)\n', ...
                s2_err.average_distance.predicted_median, ...
                s2_err.average_distance.mae_median, ...
                s2_err.average_distance.relative_error_median);
            fprintf('\n');
            
            fprintf('  Min Distance:\n');
            fprintf('    Predicted (Mean):   %8.4f mm (Error: %.4f mm)\n', ...
                s2_err.min_distance.predicted_mean, s2_err.min_distance.mae_mean);
            fprintf('\n');
            
            fprintf('  Max Distance:\n');
            fprintf('    Predicted (Mean):   %8.4f mm (Error: %.4f mm)\n', ...
                s2_err.max_distance.predicted_mean, s2_err.max_distance.mae_mean);
            fprintf('\n');
            
            fprintf('  Std Deviation:\n');
            fprintf('    Predicted (Mean):   %8.4f mm (Error: %.4f mm)\n', ...
                s2_err.standard_deviation.predicted_mean, s2_err.standard_deviation.mae_mean);
            fprintf('\n');
        end
        
        % Comparison
        if isfield(bahn_prog_stage1.predictions, k_field) && isfield(bahn_prog_stage2.predictions, k_field)
            s1_err = bahn_prog_stage1.predictions.(k_field).errors;
            s2_err = bahn_prog_stage2.predictions.(k_field).errors;
            
            fprintf('Stage 2 Improvement over Stage 1:\n');
            
            % Average Distance improvement
            mae_improvement = s1_err.average_distance.mae_mean - s2_err.average_distance.mae_mean;
            rel_improvement = s1_err.average_distance.relative_error_mean - s2_err.average_distance.relative_error_mean;
            
            fprintf('  Average Distance MAE:  %+.4f mm (%s)\n', mae_improvement, ...
                ternary(mae_improvement > 0, '✓ Better', '✗ Worse'));
            fprintf('  Relative Error:        %+.1f%% (%s)\n', rel_improvement, ...
                ternary(rel_improvement > 0, '✓ Better', '✗ Worse'));
            fprintf('\n');
        end
        
        fprintf('\n');
    end
end

fprintf('═══════════════════════════════════════════════════════════════\n\n');

% === Segment-Level Aggregated ===
fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('SEGMENT-LEVEL PERFORMANCE PROGNOSIS (Aggregated)\n');
fprintf('═══════════════════════════════════════════════════════════════\n\n');

% Aggregate errors across segments for each Top-K
for k_test_idx = 1:length(top_k_test_values)
    K = top_k_test_values(k_test_idx);
    k_field = sprintf('top_%d', K);
    
    % Collect errors
    s1_avg_errors = [];
    s2_avg_errors = [];
    
    for seg_idx = 1:num_query_segments
        if ~isempty(segment_prog_stage1{seg_idx}) && ...
           isfield(segment_prog_stage1{seg_idx}, 'predictions') && ...
           isfield(segment_prog_stage1{seg_idx}.predictions, k_field)
            
            s1_err = segment_prog_stage1{seg_idx}.predictions.(k_field).errors.average_distance.mae_mean;
            s1_avg_errors = [s1_avg_errors; s1_err];
        end
        
        if ~isempty(segment_prog_stage2{seg_idx}) && ...
           isfield(segment_prog_stage2{seg_idx}, 'predictions') && ...
           isfield(segment_prog_stage2{seg_idx}.predictions, k_field)
            
            s2_err = segment_prog_stage2{seg_idx}.predictions.(k_field).errors.average_distance.mae_mean;
            s2_avg_errors = [s2_avg_errors; s2_err];
        end
    end
    
    if ~isempty(s1_avg_errors) || ~isempty(s2_avg_errors)
        fprintf('Top-%d Average Distance Prediction:\n', K);
        
        if ~isempty(s1_avg_errors)
            fprintf('  Stage 1 MAE:  %.4f ± %.4f mm (avg over %d segments)\n', ...
                mean(s1_avg_errors), std(s1_avg_errors), length(s1_avg_errors));
        end
        
        if ~isempty(s2_avg_errors)
            fprintf('  Stage 2 MAE:  %.4f ± %.4f mm (avg over %d segments)\n', ...
                mean(s2_avg_errors), std(s2_avg_errors), length(s2_avg_errors));
        end
        
        if ~isempty(s1_avg_errors) && ~isempty(s2_avg_errors)
            improvement = mean(s1_avg_errors) - mean(s2_avg_errors);
            improvement_pct = improvement / mean(s1_avg_errors) * 100;
            
            fprintf('  Improvement:  %.4f mm (%.1f%% reduction) %s\n', ...
                improvement, improvement_pct, ...
                ternary(improvement > 0, '✓', '✗'));
        end
        
        fprintf('\n');
    end
end

fprintf('═══════════════════════════════════════════════════════════════\n\n');

% Helper function
function result = ternary(condition, true_val, false_val)
    if condition
        result = true_val;
    else
        result = false_val;
    end
end

% ========================================================================
%% HELPER FUNCTIONS
% ========================================================================

function data = loadQueryDataHierarchical(conn, schema, query_id)
    % Load query trajectory + all segments
    
    data = struct();
    
    % Load trajectory
    pos_sql = sprintf(['SELECT x_soll, y_soll, z_soll FROM %s.bahn_position_soll ' ...
                      'WHERE bahn_id = ''%s'' ORDER BY timestamp'], schema, query_id);
    pos = fetch(conn, pos_sql);
    
    if isempty(pos)
        data = [];
        return;
    end
    
    data.trajectory.position = [pos.x_soll, pos.y_soll, pos.z_soll];
    
    joint_sql = sprintf(['SELECT joint_1, joint_2, joint_3, joint_4, joint_5, joint_6 ' ...
                        'FROM %s.bahn_joint_states WHERE bahn_id = ''%s'' ORDER BY timestamp'], ...
                        schema, query_id);
    joint = fetch(conn, joint_sql);
    data.trajectory.joint = [joint.joint_1, joint.joint_2, joint.joint_3, ...
                            joint.joint_4, joint.joint_5, joint.joint_6];
    
    % Load segments
    seg_sql = sprintf(['SELECT DISTINCT segment_id FROM %s.bahn_metadata ' ...
                      'WHERE bahn_id = ''%s'' AND segment_id != bahn_id ORDER BY segment_id'], ...
                      schema, query_id);
    segs = fetch(conn, seg_sql);
    
    data.segment_ids = segs.segment_id;
    data.segments.position = cell(length(data.segment_ids), 1);
    data.segments.joint = cell(length(data.segment_ids), 1);
    
    for i = 1:length(data.segment_ids)
        seg_id = data.segment_ids{i};
        
        pos_seg_sql = sprintf(['SELECT x_soll, y_soll, z_soll FROM %s.bahn_position_soll ' ...
                              'WHERE segment_id = ''%s'' ORDER BY timestamp'], schema, seg_id);
        pos_seg = fetch(conn, pos_seg_sql);
        data.segments.position{i} = [pos_seg.x_soll, pos_seg.y_soll, pos_seg.z_soll];
        
        joint_seg_sql = sprintf(['SELECT joint_1, joint_2, joint_3, joint_4, joint_5, joint_6 ' ...
                                'FROM %s.bahn_joint_states WHERE segment_id = ''%s'' ORDER BY timestamp'], ...
                                schema, seg_id);
        joint_seg = fetch(conn, joint_seg_sql);
        data.segments.joint{i} = [joint_seg.joint_1, joint_seg.joint_2, joint_seg.joint_3, ...
                                 joint_seg.joint_4, joint_seg.joint_5, joint_seg.joint_6];
    end
end

function embeddings = getQueryEmbeddings(conn, schema, id, level)
    % Get embeddings for query (bahn or segment)
    
    query = sprintf(['SELECT position_embedding::text, joint_embedding::text, ' ...
                    'orientation_embedding::text, velocity_embedding::text, ' ...
                    'metadata_embedding::text ' ...
                    'FROM %s.bahn_embeddings WHERE segment_id = ''%s'''], schema, id);
    
    result = fetch(conn, query);
    
    if isempty(result)
        embeddings = [];
        return;
    end
    
    % ✅ PARSE ALL EMBEDDINGS
    embeddings = struct();
    embeddings.position = parseEmbedding(result.position_embedding);
    embeddings.joint = parseEmbedding(result.joint_embedding);
    embeddings.orientation = parseEmbedding(result.orientation_embedding);
    embeddings.velocity = parseEmbedding(result.velocity_embedding);
    embeddings.metadata = parseEmbedding(result.metadata_embedding);
end

function emb_vec = parseEmbedding(emb_text)
    % Parse embedding from text format to numeric vector
    
    if isempty(emb_text) || (iscell(emb_text) && isempty(emb_text{1}))
        emb_vec = [];
        return;
    end
    
    % Convert to char if cell
    if iscell(emb_text)
        emb_text = emb_text{1};
    end
    
    % Remove brackets
    emb_text = strrep(strrep(char(emb_text), '[', ''), ']', '');
    
    % Parse numbers
    emb_vec = str2double(strsplit(emb_text, ','));
    
    % Remove NaN values (if any)
    emb_vec = emb_vec(~isnan(emb_vec));
end

function candidates = performRRFEmbeddingSearch(conn, schema, query_embeddings, weights, k, level)
    % Perform RRF fusion on all embedding modalities
    
    % Set level constraint
    if strcmp(level, 'bahn')
        level_constraint = 'segment_id = bahn_id';
    else
        level_constraint = 'segment_id != bahn_id';
    end
    
    % Search each modality
    modalities = {'position', 'joint', 'orientation', 'velocity', 'metadata'};
    rankings = cell(5, 1);
    
    % Set HNSW parameters
    execute(conn, 'SET hnsw.ef_search = 200');
    
    for i = 1:5
        mod = modalities{i};
        query_emb = query_embeddings.(mod);
        weight = weights(i);
        
        if weight == 0 || isempty(query_emb)
            continue;
        end
        
        % ✅ PARSE EMBEDDING (from text to vector)
        if ischar(query_emb) || isstring(query_emb)
            % Remove brackets and parse
            emb_text = char(query_emb);
            emb_text = strrep(strrep(emb_text, '[', ''), ']', '');
            query_emb = str2double(strsplit(emb_text, ','));
        end
        
        % Convert embedding to string for SQL
        emb_str = sprintf('[%s]', strjoin(string(query_emb), ','));
        col_name = sprintf('%s_embedding', mod);
        
        % Query
        execute(conn, 'SET search_path = "bewegungsdaten"');
        search_sql = sprintf(['SELECT segment_id, bahn_id, ' ...
                             '%s <=> ''%s''::vector as distance ' ...
                             'FROM %s.bahn_embeddings ' ...
                             'WHERE %s AND %s IS NOT NULL ' ...
                             'ORDER BY distance LIMIT %d'], ...
                             col_name, emb_str, schema, level_constraint, col_name, k*2);
        
        try
            result = fetch(conn, search_sql);
            
            if ~isempty(result)
                rankings{i} = result;
            end
        catch ME
            warning('Failed to search %s: %s', mod, ME.message);
            rankings{i} = [];
        end
    end
    
    % RRF Fusion
    candidates = fuseRankingsRRF(rankings, weights, k);
end

function fused_ids = fuseRankingsRRF(rankings, weights, k)
    % Simple RRF implementation
    
    rrf_k = 60;  % RRF constant
    scores = containers.Map();
    
    for i = 1:length(rankings)
        if isempty(rankings{i})
            continue;
        end
        
        weight = weights(i);
        results = rankings{i};
        
        for rank = 1:height(results)
            id = results.segment_id{rank};
            
            if ~scores.isKey(id)
                scores(id) = 0;
            end
            
            scores(id) = scores(id) + weight / (rrf_k + rank);
        end
    end
    
    % Sort by score
    all_ids = keys(scores);
    all_scores = cell2mat(values(scores, all_ids));
    
    [~, sort_idx] = sort(all_scores, 'descend');
    fused_ids = all_ids(sort_idx);
    fused_ids = fused_ids(1:min(k, length(fused_ids)));
end

function data = loadCandidateDataBatch(conn, schema, ids, level)
    % Load candidate trajectory/segment data in batch
    
    n = length(ids);
    data = struct();
    data.ids = ids;
    data.position = cell(n, 1);
    data.joint = cell(n, 1);
    
    for i = 1:n
        id = ids{i};
        
        if strcmp(level, 'bahn')
            id_field = 'bahn_id';
        else
            id_field = 'segment_id';
        end
        
        % Position
        pos_sql = sprintf(['SELECT x_soll, y_soll, z_soll FROM %s.bahn_position_soll ' ...
                          'WHERE %s = ''%s'' ORDER BY timestamp'], schema, id_field, id);
        pos = fetch(conn, pos_sql);
        
        if ~isempty(pos)
            data.position{i} = [pos.x_soll, pos.y_soll, pos.z_soll];
        end
        
        % Joint
        joint_sql = sprintf(['SELECT joint_1, joint_2, joint_3, joint_4, joint_5, joint_6 ' ...
                            'FROM %s.bahn_joint_states WHERE %s = ''%s'' ORDER BY timestamp'], ...
                            schema, id_field, id);
        joint = fetch(conn, joint_sql);
        
        if ~isempty(joint)
            data.joint{i} = [joint.joint_1, joint.joint_2, joint.joint_3, ...
                            joint.joint_4, joint.joint_5, joint.joint_6];
        end
    end
end

function [results, stats] = performDTWReranking(query_seq, candidate_data, dtw_mode, dtw_config, limit)
    % DTW Reranking with Lower Bounds
    % Returns results with Stage 1 ranks for later analysis
    
    K = length(candidate_data.ids);
    
    % Get candidate sequences
    if strcmp(dtw_mode, 'position')
        cand_seqs = candidate_data.position;
    else
        cand_seqs = candidate_data.joint;
    end
    
    % ✅ STORE STAGE 1 RANKINGS (for later analysis, not computed here!)
    stage1_ranks = (1:K)';  % Stage 1 rank = input order from RRF
    
    % Validate candidates
    valid_mask = false(K, 1);
    for i = 1:K
        if ~isempty(cand_seqs{i}) && size(cand_seqs{i}, 1) >= 2
            valid_mask(i) = true;
        end
    end
    
    num_valid = sum(valid_mask);
    
    % Estimate filtering
    lb_kim_filtered = round(num_valid * (1 - dtw_config.lb_kim_keep_ratio));
    after_kim = num_valid - lb_kim_filtered;
    lb_keogh_filtered = max(0, after_kim - dtw_config.lb_keogh_candidates);
    dtw_calls = num_valid - lb_kim_filtered - lb_keogh_filtered;
    
    % ============================================================
    % COMPUTE DTW (Pure timing - no analysis here!)
    % ============================================================
    dtw_dists = inf(K, 1);
    
    for i = 1:K
        if ~valid_mask(i)
            continue;
        end
        
        try
            dtw_dists(i) = cDTW(query_seq, cand_seqs{i}, dtw_mode, dtw_config.cdtw_window, ...
                inf, dtw_config.use_rotation_alignment, dtw_config.normalize_dtw);
        catch
            dtw_dists(i) = inf;
        end
    end
    
    % Sort by DTW distance
    [sorted_dists, sort_idx] = sort(dtw_dists);
    
    % ============================================================
    % BUILD RESULTS (with Stage 1 info for later analysis)
    % ============================================================
    results = [];
    for i = 1:min(limit, K)
        if sorted_dists(i) == inf
            break;
        end
        
        idx = sort_idx(i);
        
        result = struct();
        result.segment_id = candidate_data.ids{idx};
        result.bahn_id = candidate_data.ids{idx};
        result.dtw_distance = sorted_dists(i);
        result.similarity_score = 1 / (1 + sorted_dists(i));
        result.stage2_rank = i;  % ✅ DTW rank (Stage 2)
        result.stage1_rank = stage1_ranks(idx);  % ✅ Original RRF rank (Stage 1)
        result.rank_change = result.stage1_rank - i;  % ✅ Positive = improved
        
        results = [results; result];
    end
    
    % Stats
    stats = struct();
    stats.dtw_calls = dtw_calls;
    stats.lb_kim_filtered = lb_kim_filtered;
    stats.lb_keogh_filtered = lb_keogh_filtered;
    stats.valid_candidates = num_valid;
end

function rerank_metrics = analyzeReranking(results, k_stage1)
    % Analyze reranking quality AFTER DTW computation
    % Called separately to not affect performance timings
    
    if isempty(results)
        rerank_metrics = struct();
        rerank_metrics.kendall_tau = NaN;
        rerank_metrics.spearman_rho = NaN;
        rerank_metrics.mean_rank_change = NaN;
        return;
    end
    
    n = length(results);
    rerank_metrics = struct();
    
    % Extract ranks
    stage1_ranks = [results.stage1_rank]';
    stage2_ranks = [results.stage2_rank]';
    rank_changes = abs([results.rank_change]');
    
    % === 1. RANK CORRELATION ===
    if n >= 2
        % Kendall's Tau
        tau = computeKendallTau(stage1_ranks, stage2_ranks);
        rerank_metrics.kendall_tau = tau;
        
        % Spearman's Rho
        rho = corr(stage1_ranks, stage2_ranks, 'Type', 'Spearman');
        rerank_metrics.spearman_rho = rho;
    else
        rerank_metrics.kendall_tau = NaN;
        rerank_metrics.spearman_rho = NaN;
    end
    
    % === 2. RANK CHANGES ===
    rerank_metrics.mean_rank_change = mean(rank_changes);
    rerank_metrics.max_rank_change = max(rank_changes);
    rerank_metrics.median_rank_change = median(rank_changes);
    rerank_metrics.std_rank_change = std(rank_changes);
    
    % === 3. TOP-K OVERLAP ===
    % Get Stage 1 IDs (sorted by Stage 1 rank)
    [~, s1_sort_idx] = sort(stage1_ranks);
    stage1_ids_sorted = {results(s1_sort_idx).segment_id};
    
    % Get Stage 2 IDs (already sorted by Stage 2 rank)
    stage2_ids_sorted = {results.segment_id};
    
    for top_k = [5, 10, 20, 50]
        if top_k <= n
            stage1_top_k = stage1_ids_sorted(1:top_k);
            stage2_top_k = stage2_ids_sorted(1:top_k);
            
            overlap = length(intersect(stage1_top_k, stage2_top_k));
            overlap_pct = overlap / top_k * 100;
            
            field_name = sprintf('top%d_overlap_pct', top_k);
            rerank_metrics.(field_name) = overlap_pct;
            
            field_name_count = sprintf('top%d_overlap_count', top_k);
            rerank_metrics.(field_name_count) = overlap;
        end
    end
    
    % === 4. RERANKING IMPACT ===
    improved = sum([results.rank_change] > 0);  % Stage1 rank > Stage2 rank = moved up
    degraded = sum([results.rank_change] < 0);  % Stage1 rank < Stage2 rank = moved down
    stable = sum([results.rank_change] == 0);
    
    rerank_metrics.num_improved = improved;
    rerank_metrics.num_degraded = degraded;
    rerank_metrics.num_stable = stable;
    rerank_metrics.improvement_rate = improved / n * 100;
    rerank_metrics.degradation_rate = degraded / n * 100;
    rerank_metrics.stability_rate = stable / n * 100;
    
    % === 5. BIGGEST MOVERS ===
    [max_change, max_idx] = max(rank_changes);
    rerank_metrics.biggest_mover_id = results(max_idx).segment_id;
    rerank_metrics.biggest_mover_stage1 = results(max_idx).stage1_rank;
    rerank_metrics.biggest_mover_stage2 = results(max_idx).stage2_rank;
    rerank_metrics.biggest_mover_change = max_change;
    
    % === 6. DISTRIBUTION ANALYSIS ===
    rerank_metrics.num_results = n;
    rerank_metrics.stage1_range = [min(stage1_ranks), max(stage1_ranks)];
    rerank_metrics.avg_stage1_rank = mean(stage1_ranks);
    rerank_metrics.avg_stage2_rank = mean(stage2_ranks);
end

function tau = computeKendallTau(x, y)
    % Compute Kendall's Tau correlation coefficient
    
    n = length(x);
    if n < 2
        tau = NaN;
        return;
    end
    
    concordant = 0;
    discordant = 0;
    
    for i = 1:n-1
        for j = i+1:n
            sign_x = sign(x(j) - x(i));
            sign_y = sign(y(j) - y(i));
            
            if sign_x * sign_y > 0
                concordant = concordant + 1;
            elseif sign_x * sign_y < 0
                discordant = discordant + 1;
            end
        end
    end
    
    total_pairs = n * (n - 1) / 2;
    tau = (concordant - discordant) / total_pairs;
end

function prognosis_metrics = evaluatePerformancePrognosis(conn, query_id, results, top_k_values)
    % Evaluate performance prognosis quality
    % 
    % Idea: Use SIDTW metrics of Top-K similar trajectories/segments
    %       to predict the query's performance
    %
    % Args:
    %   conn: Database connection
    %   query_id: Query bahn/segment ID
    %   results: Results struct with rankings
    %   top_k_values: Array of K values to test, e.g. [5, 10, 20]
    %
    % Returns:
    %   prognosis_metrics: Prediction errors for different Top-K
    
    if isempty(results)
        prognosis_metrics = struct();
        return;
    end
    
    % ================================================================
    % 1. GET ACTUAL SIDTW VALUES (Ground Truth)
    % ================================================================
    
    actual_sidtw = getSIDTWMetrics(conn, query_id);
    
    if isempty(actual_sidtw)
        warning('No SIDTW data found for query %s', query_id);
        prognosis_metrics = struct();
        prognosis_metrics.error = 'No SIDTW data for query';
        return;
    end
    
    % ================================================================
    % 2. GET SIDTW FOR ALL RESULTS
    % ================================================================
    
    n = length(results);
    result_ids = cell(n, 1);
    
    for i = 1:n
        result_ids{i} = results(i).segment_id;
    end
    
    result_sidtw = getSIDTWMetricsBatch(conn, result_ids);
    
    % ================================================================
    % 3. FOR EACH TOP-K: PREDICT & MEASURE ERROR
    % ================================================================
    
    prognosis_metrics = struct();
    prognosis_metrics.actual = actual_sidtw;
    prognosis_metrics.predictions = struct();
    
    for k_idx = 1:length(top_k_values)
        K = top_k_values(k_idx);
        
        if K > n
            continue;  % Skip if not enough results
        end
        
        % Get Top-K results
        top_k_results = results(1:K);
        
        % ============================================================
        % AGGREGATE SIDTW FROM TOP-K CANDIDATES
        % ============================================================
        
        top_k_sidtw_values = struct();
        top_k_sidtw_values.min_distance = [];
        top_k_sidtw_values.max_distance = [];
        top_k_sidtw_values.average_distance = [];
        top_k_sidtw_values.standard_deviation = [];
        
        for i = 1:K
            result_id = top_k_results(i).segment_id;
            field_name = matlab.lang.makeValidName(result_id);
            
            if ~isfield(result_sidtw, field_name)
                continue;  % Skip if no SIDTW data
            end
            
            r_sidtw = result_sidtw.(field_name);
            
            top_k_sidtw_values.min_distance = [top_k_sidtw_values.min_distance; r_sidtw.min_distance];
            top_k_sidtw_values.max_distance = [top_k_sidtw_values.max_distance; r_sidtw.max_distance];
            top_k_sidtw_values.average_distance = [top_k_sidtw_values.average_distance; r_sidtw.average_distance];
            top_k_sidtw_values.standard_deviation = [top_k_sidtw_values.standard_deviation; r_sidtw.standard_deviation];
        end
        
        % ============================================================
        % COMPUTE PREDICTIONS (Mean, Median, Min, Max)
        % ============================================================
        
        prediction = struct();
        
        % For each SIDTW metric
        metrics = {'min_distance', 'max_distance', 'average_distance', 'standard_deviation'};
        
        for m_idx = 1:length(metrics)
            metric_name = metrics{m_idx};
            values = top_k_sidtw_values.(metric_name);
            
            if isempty(values)
                prediction.(metric_name).mean = NaN;
                prediction.(metric_name).median = NaN;
                prediction.(metric_name).min = NaN;
                prediction.(metric_name).max = NaN;
                prediction.(metric_name).std = NaN;
                continue;
            end
            
            % Aggregations
            prediction.(metric_name).mean = mean(values);
            prediction.(metric_name).median = median(values);
            prediction.(metric_name).min = min(values);
            prediction.(metric_name).max = max(values);
            prediction.(metric_name).std = std(values);
            prediction.(metric_name).num_valid = length(values);
        end
        
        % ============================================================
        % COMPUTE ERRORS (vs Actual)
        % ============================================================
        
        errors = struct();
        
        for m_idx = 1:length(metrics)
            metric_name = metrics{m_idx};
            actual_value = actual_sidtw.(metric_name);
            
            if isnan(actual_value)
                errors.(metric_name).mae_mean = NaN;
                errors.(metric_name).mae_median = NaN;
                errors.(metric_name).relative_error_mean = NaN;
                errors.(metric_name).relative_error_median = NaN;
                continue;
            end
            
            % Mean Absolute Error
            errors.(metric_name).mae_mean = abs(prediction.(metric_name).mean - actual_value);
            errors.(metric_name).mae_median = abs(prediction.(metric_name).median - actual_value);
            
            % Relative Error (%)
            if actual_value ~= 0
                errors.(metric_name).relative_error_mean = ...
                    abs(prediction.(metric_name).mean - actual_value) / actual_value * 100;
                errors.(metric_name).relative_error_median = ...
                    abs(prediction.(metric_name).median - actual_value) / actual_value * 100;
            else
                errors.(metric_name).relative_error_mean = NaN;
                errors.(metric_name).relative_error_median = NaN;
            end
            
            % Store actual value for reference
            errors.(metric_name).actual = actual_value;
            errors.(metric_name).predicted_mean = prediction.(metric_name).mean;
            errors.(metric_name).predicted_median = prediction.(metric_name).median;
        end
        
        % ============================================================
        % STORE FOR THIS K
        % ============================================================
        
        field_name = sprintf('top_%d', K);
        prognosis_metrics.predictions.(field_name).prediction = prediction;
        prognosis_metrics.predictions.(field_name).errors = errors;
        prognosis_metrics.predictions.(field_name).k = K;
    end
end

function sidtw_data = getSIDTWMetrics(conn, segment_id)
    % Get SIDTW metrics for a single segment/bahn
    
    query = sprintf(['SELECT segment_id, ' ...
                    'sidtw_min_distance, sidtw_max_distance, ' ...
                    'sidtw_average_distance, sidtw_standard_deviation ' ...
                    'FROM auswertung.info_sidtw ' ...
                    'WHERE segment_id = ''%s'''], segment_id);
    
    result = fetch(conn, query);
    
    if isempty(result)
        sidtw_data = [];
        return;
    end
    
    sidtw_data = struct();
    sidtw_data.segment_id = result.segment_id{1};
    sidtw_data.min_distance = result.sidtw_min_distance;
    sidtw_data.max_distance = result.sidtw_max_distance;
    sidtw_data.average_distance = result.sidtw_average_distance;
    sidtw_data.standard_deviation = result.sidtw_standard_deviation;
end

function sidtw_batch = getSIDTWMetricsBatch(conn, segment_ids)
    % Get SIDTW metrics for multiple segments/bahnen
    
    if isempty(segment_ids)
        sidtw_batch = struct();
        return;
    end
    
    % Build IN clause
    ids_str = sprintf('''%s'',', segment_ids{:});
    ids_str = ids_str(1:end-1);  % Remove trailing comma
    
    query = sprintf(['SELECT segment_id, ' ...
                    'sidtw_min_distance, sidtw_max_distance, ' ...
                    'sidtw_average_distance, sidtw_standard_deviation ' ...
                    'FROM auswertung.info_sidtw ' ...
                    'WHERE segment_id IN (%s)'], ids_str);
    
    results = fetch(conn, query);
    
    sidtw_batch = struct();
    
    if isempty(results)
        return;
    end
    
    for i = 1:height(results)
        seg_id = results.segment_id{i};
        field_name = matlab.lang.makeValidName(seg_id);
        
        sidtw_batch.(field_name) = struct();
        sidtw_batch.(field_name).segment_id = seg_id;
        sidtw_batch.(field_name).min_distance = results.sidtw_min_distance(i);
        sidtw_batch.(field_name).max_distance = results.sidtw_max_distance(i);
        sidtw_batch.(field_name).average_distance = results.sidtw_average_distance(i);
        sidtw_batch.(field_name).standard_deviation = results.sidtw_standard_deviation(i);
    end
end
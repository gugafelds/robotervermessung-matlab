%  kNN METADATA BASELINE - WITH TIMING ANALYSIS
%  ========================================================================
%  Baseline comparison: kNN on metadata features (no embeddings, no DTW)
%  Includes timing measurements for deployment comparison
%  ========================================================================

clear; clc;

addpath(genpath(pwd));
addpath(genpath('../main'));

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  kNN METADATA BASELINE - TIMING ANALYSIS                       ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

% ========================================================================
%% SECTION 1: CONFIGURATION
% ========================================================================

queries_quantity = 500;
random_seed = 21;
all_k_values = [5, 10, 25, 50];

exclude_ids = {
    '1765989370'; '1765989294'; '1765988821'; '1765988920'; '1765989411';
    '1765990630'; '1765990747'; '1765990822'; '1765991047'; '1765991234';
    '1765991190'; '1765991445'; '1765991515'; '1765991949'; '1765991743';
    '1769770498'; '1769770684'; '1769770935'; '1769771107'; '1769771447';
    '1769773928'; '1769772060'; '1769772213'; '1769773985'; '1769774278';
    '1769772609'; '1769773593'; '1769772776'; '1769772900'; '1769773333';
    '1769774581';
};

fprintf('Queries:        %d\n', queries_quantity);
fprintf('K values:       [%s]\n', strjoin(string(all_k_values), ', '));
fprintf('Random seed:    %d\n\n', random_seed);

% ========================================================================
%% SECTION 2: LOAD ALL DATA (TIMED)
% ========================================================================

fprintf('═══ LOADING DATA ═══\n\n');

total_loading_start = tic;

conn = connectingToPostgres();
schema = 'bewegungsdaten';

% --- Bahn metadata ---
t1 = tic;
meta_sql = sprintf([...
    'SELECT m.bahn_id, m.movement_type, m.duration, m.weight, m.length, ' ...
    'm.min_twist_ist, m.max_twist_ist, m.mean_twist_ist, ' ...
    'm.median_twist_ist, m.std_twist_ist, ' ...
    'm.min_acceleration_ist, m.max_acceleration_ist, m.mean_acceleration_ist, ' ...
    'm.median_acceleration_ist, m.std_acceleration_ist, ' ...
    'm.position_x, m.position_y, m.position_z ' ...
    'FROM %s.bahn_metadata m ' ...
    'WHERE m.segment_id = m.bahn_id'], schema);
meta_all = fetch(conn, meta_sql);
time_load_bahn_meta = toc(t1);
fprintf('  ✓ Bahn metadata:    %d rows (%.2f sec)\n', height(meta_all), time_load_bahn_meta);

% --- Bahn ground truth ---
t1 = tic;
gt_sql = sprintf([...
    'SELECT segment_id, sidtw_average_distance ' ...
    'FROM auswertung.info_sidtw ' ...
    'WHERE segment_id IN (SELECT bahn_id FROM %s.bahn_metadata WHERE segment_id = bahn_id)'], schema);
gt_all = fetch(conn, gt_sql);
time_load_bahn_gt = toc(t1);
fprintf('  ✓ Bahn ground truth: %d rows (%.2f sec)\n', height(gt_all), time_load_bahn_gt);

% --- Segment metadata ---
t1 = tic;
seg_meta_sql = sprintf([...
    'SELECT m.bahn_id, m.segment_id, m.movement_type, m.duration, m.weight, m.length, ' ...
    'm.min_twist_ist, m.max_twist_ist, m.mean_twist_ist, ' ...
    'm.median_twist_ist, m.std_twist_ist, ' ...
    'm.min_acceleration_ist, m.max_acceleration_ist, m.mean_acceleration_ist, ' ...
    'm.median_acceleration_ist, m.std_acceleration_ist, ' ...
    'm.position_x, m.position_y, m.position_z ' ...
    'FROM %s.bahn_metadata m ' ...
    'WHERE m.segment_id != m.bahn_id'], schema);
seg_meta = fetch(conn, seg_meta_sql);
time_load_seg_meta = toc(t1);
fprintf('  ✓ Segment metadata:  %d rows (%.2f sec)\n', height(seg_meta), time_load_seg_meta);

% --- Segment ground truth ---
t1 = tic;
seg_gt_sql = sprintf([...
    'SELECT segment_id, sidtw_average_distance ' ...
    'FROM auswertung.info_sidtw ' ...
    'WHERE segment_id IN (SELECT segment_id FROM %s.bahn_metadata WHERE segment_id != bahn_id)'], schema);
seg_gt = fetch(conn, seg_gt_sql);
time_load_seg_gt = toc(t1);
fprintf('  ✓ Segment ground truth: %d rows (%.2f sec)\n', height(seg_gt), time_load_seg_gt);

% --- Segment lengths ---
t1 = tic;
seg_len_sql = sprintf([...
    'SELECT segment_id, COUNT(*) as num_points ' ...
    'FROM %s.bahn_position_soll ' ...
    'WHERE segment_id != bahn_id ' ...
    'GROUP BY segment_id'], schema);
seg_lengths_data = fetch(conn, seg_len_sql);
time_load_seg_len = toc(t1);
fprintf('  ✓ Segment lengths:   %d rows (%.2f sec)\n', height(seg_lengths_data), time_load_seg_len);

% --- Query IDs ---
ids_str_list = sprintf('''%s'',', exclude_ids{:});
ids_str_list = ids_str_list(1:end-1);

t1 = tic;
query_sql = sprintf([...
    'SELECT bahn_id FROM (' ...
    'SELECT DISTINCT b.bahn_id FROM %s.bahn_metadata b ' ...
    'INNER JOIN auswertung.info_sidtw s ON b.bahn_id = s.segment_id ' ...
    'WHERE b.segment_id = b.bahn_id ' ...
    'AND b.bahn_id NOT IN (' ...
        'SELECT t_all.bahn_id ' ...
        'FROM %s.bahn_info t_exclude ' ...
        'JOIN %s.bahn_info t_all ON t_exclude.record_filename = t_all.record_filename ' ...
        'WHERE t_exclude.bahn_id IN (%s) ' ...
    ') ' ...
    ') AS distinct_bahnen ' ...
    'ORDER BY md5(bahn_id || ''%d'') LIMIT %d'], ...
    schema, schema, schema, ids_str_list, random_seed, queries_quantity);
query_results = fetch(conn, query_sql);
query_ids = query_results.bahn_id;
time_load_queries = toc(t1);
fprintf('  ✓ Query IDs:         %d rows (%.2f sec)\n', length(query_ids), time_load_queries);

close(conn);

total_loading_time = toc(total_loading_start);
fprintf('\n  ══ TOTAL LOADING TIME: %.2f sec ══\n\n', total_loading_time);

% ========================================================================
%% SECTION 3: BUILD FEATURE MATRICES (TIMED)
% ========================================================================

fprintf('═══ BUILDING FEATURE MATRICES ═══\n\n');

t_build = tic;

% --- Bahn features ---
all_movement_types = meta_all.movement_type;
unique_types = unique(all_movement_types);
type_map = containers.Map();
for i = 1:length(unique_types)
    type_map(unique_types{i}) = i;
end

movement_numeric = zeros(height(meta_all), 1);
for i = 1:height(meta_all)
    movement_numeric(i) = type_map(all_movement_types{i});
end

feature_cols = [movement_numeric, ...
    meta_all.duration, meta_all.weight, meta_all.length, ...
    meta_all.min_twist_ist, meta_all.max_twist_ist, meta_all.mean_twist_ist, ...
    meta_all.median_twist_ist, meta_all.std_twist_ist, ...
    meta_all.min_acceleration_ist, meta_all.max_acceleration_ist, meta_all.mean_acceleration_ist, ...
    meta_all.median_acceleration_ist, meta_all.std_acceleration_ist, ...
    meta_all.position_x, meta_all.position_y, meta_all.position_z];

feature_cols(isnan(feature_cols)) = 0;
feat_min = min(feature_cols, [], 1);
feat_max = max(feature_cols, [], 1);
feat_range = feat_max - feat_min;
feat_range(feat_range == 0) = 1;
X = (feature_cols - feat_min) ./ feat_range;

all_bahn_ids = meta_all.bahn_id;
y = nan(height(meta_all), 1);
for i = 1:height(meta_all)
    idx = find(strcmp(gt_all.segment_id, all_bahn_ids{i}), 1);
    if ~isempty(idx)
        y(i) = gt_all.sidtw_average_distance(idx);
    end
end

% --- Segment features ---
seg_movement_numeric = zeros(height(seg_meta), 1);
for i = 1:height(seg_meta)
    mt = seg_meta.movement_type{i};
    if type_map.isKey(mt)
        seg_movement_numeric(i) = type_map(mt);
    else
        seg_movement_numeric(i) = 0;
    end
end

seg_feature_cols = [seg_movement_numeric, ...
    seg_meta.duration, seg_meta.weight, seg_meta.length, ...
    seg_meta.min_twist_ist, seg_meta.max_twist_ist, seg_meta.mean_twist_ist, ...
    seg_meta.median_twist_ist, seg_meta.std_twist_ist, ...
    seg_meta.min_acceleration_ist, seg_meta.max_acceleration_ist, seg_meta.mean_acceleration_ist, ...
    seg_meta.median_acceleration_ist, seg_meta.std_acceleration_ist, ...
    seg_meta.position_x, seg_meta.position_y, seg_meta.position_z];

seg_feature_cols(isnan(seg_feature_cols)) = 0;
seg_feat_min = min(seg_feature_cols, [], 1);
seg_feat_max = max(seg_feature_cols, [], 1);
seg_feat_range = seg_feat_max - seg_feat_min;
seg_feat_range(seg_feat_range == 0) = 1;
X_seg = (seg_feature_cols - seg_feat_min) ./ seg_feat_range;

all_seg_ids = seg_meta.segment_id;
all_seg_bahn_ids = seg_meta.bahn_id;
y_seg = nan(height(seg_meta), 1);
for i = 1:height(seg_meta)
    idx = find(strcmp(seg_gt.segment_id, all_seg_ids{i}), 1);
    if ~isempty(idx)
        y_seg(i) = seg_gt.sidtw_average_distance(idx);
    end
end

seg_len_map = containers.Map();
for i = 1:height(seg_lengths_data)
    seg_len_map(seg_lengths_data.segment_id{i}) = seg_lengths_data.num_points(i);
end

time_build = toc(t_build);
fprintf('  Feature matrices built in %.2f sec\n', time_build);
fprintf('  Bahn: %d x %d (valid GT: %d)\n', size(X), sum(~isnan(y)));
fprintf('  Segments: %d x %d (valid GT: %d)\n\n', size(X_seg), sum(~isnan(y_seg)));

% ========================================================================
%% SECTION 4: kNN PREDICTION WITH PER-QUERY TIMING
% ========================================================================

fprintf('═══ RUNNING kNN BASELINE (DIRECT + DECOMPOSED) ═══\n\n');

timestamp = datestr(now, 'yyyymmdd_HHMMSS');
if ~exist('results', 'dir'), mkdir('results'); end

all_results = [];
decomp_results = [];

% Timing arrays
query_times_direct = zeros(length(query_ids), 1);
query_times_decomp = zeros(length(query_ids), 1);
query_num_segments = zeros(length(query_ids), 1);

for q_idx = 1:length(query_ids)
    query_id = query_ids{q_idx};
    
    q_row = find(strcmp(all_bahn_ids, query_id), 1);
    if isempty(q_row) || isnan(y(q_row))
        continue;
    end
    
    gt_value = y(q_row);
    query_feat = X(q_row, :);
    
    % ════════════════════════════════════════
    % DIRECT kNN (timed)
    % ════════════════════════════════════════
    t_direct = tic;
    
    dists = vecnorm(X - query_feat, 2, 2);
    dists(q_row) = inf;
    dists(isnan(y)) = inf;
    [sorted_dists, sorted_idx] = sort(dists, 'ascend');
    
    for k_idx = 1:length(all_k_values)
        K = all_k_values(k_idx);
        neighbors = sorted_idx(1:K);
        neighbor_dists = sorted_dists(1:K);
        neighbor_values = y(neighbors);
        
        pred_simple = mean(neighbor_values);
        w = 1 ./ (neighbor_dists + 1e-6);
        w = w / sum(w);
        pred_weighted = sum(w .* neighbor_values);
        
        row = struct();
        row.query_id = query_id;
        row.K = K;
        row.ground_truth = gt_value;
        row.pred_simple = pred_simple;
        row.pred_weighted = pred_weighted;
        row.err_simple = abs(gt_value - pred_simple);
        row.err_weighted = abs(gt_value - pred_weighted);
        all_results = [all_results; row];
    end
    
    query_times_direct(q_idx) = toc(t_direct);
    
    % ════════════════════════════════════════
    % DECOMPOSED kNN (timed)
    % ════════════════════════════════════════
    t_decomp = tic;
    
    seg_mask = strcmp(all_seg_bahn_ids, query_id);
    query_seg_indices = find(seg_mask);
    query_num_segments(q_idx) = length(query_seg_indices);
    
    if isempty(query_seg_indices)
        query_times_decomp(q_idx) = toc(t_decomp);
        continue;
    end
    
    for k_idx = 1:length(all_k_values)
        K = all_k_values(k_idx);
        
        seg_preds_simple = zeros(length(query_seg_indices), 1);
        seg_preds_weighted = zeros(length(query_seg_indices), 1);
        seg_weights = zeros(length(query_seg_indices), 1);
        valid_segs = true(length(query_seg_indices), 1);
        
        for s = 1:length(query_seg_indices)
            si = query_seg_indices(s);
            seg_id = all_seg_ids{si};
            
            if isnan(y_seg(si))
                valid_segs(s) = false;
                continue;
            end
            
            query_seg_feat = X_seg(si, :);
            
            dists_seg = vecnorm(X_seg - query_seg_feat, 2, 2);
            dists_seg(si) = inf;
            
            same_bahn = strcmp(all_seg_bahn_ids, query_id);
            dists_seg(same_bahn) = inf;
            dists_seg(isnan(y_seg)) = inf;
            
            [sorted_d, sorted_i] = sort(dists_seg, 'ascend');
            
            if sorted_d(1) == inf
                valid_segs(s) = false;
                continue;
            end
            
            nb = sorted_i(1:K);
            nb_d = sorted_d(1:K);
            nb_y = y_seg(nb);
            
            seg_preds_simple(s) = mean(nb_y);
            w = 1 ./ (nb_d + 1e-6);
            w = w / sum(w);
            seg_preds_weighted(s) = sum(w .* nb_y);
            
            if seg_len_map.isKey(seg_id)
                seg_weights(s) = seg_len_map(seg_id);
            else
                seg_weights(s) = 1;
            end
        end
        
        if any(valid_segs)
            sw = seg_weights(valid_segs) / sum(seg_weights(valid_segs));
            agg_simple = sum(sw .* seg_preds_simple(valid_segs));
            agg_weighted = sum(sw .* seg_preds_weighted(valid_segs));
        else
            agg_simple = NaN;
            agg_weighted = NaN;
        end
        
        row = struct();
        row.query_id = query_id;
        row.K = K;
        row.ground_truth = gt_value;
        row.pred_simple = agg_simple;
        row.pred_weighted = agg_weighted;
        row.err_simple = abs(gt_value - agg_simple);
        row.err_weighted = abs(gt_value - agg_weighted);
        decomp_results = [decomp_results; row];
    end
    
    query_times_decomp(q_idx) = toc(t_decomp);
    
    if mod(q_idx, 50) == 0
        fprintf('  Progress: %d/%d (direct: %.4fs, decomp: %.4fs, %d segs)\n', ...
            q_idx, length(query_ids), ...
            query_times_direct(q_idx), query_times_decomp(q_idx), ...
            query_num_segments(q_idx));
    end
end

fprintf('\n  ✓ Completed all predictions\n\n');

% ========================================================================
%% SECTION 5: RESULTS SUMMARY
% ========================================================================

fprintf('═══ PREDICTION RESULTS (DIRECT) ═══\n\n');
fprintf('%-5s  %-12s %-12s %-12s %-12s\n', 'K', 'MAE_simple', 'MAE_weighted', 'RMSE_simple', 'RMSE_weighted');
fprintf('%s\n', repmat('-', 1, 57));

results_table = struct2table(all_results);
for k_idx = 1:length(all_k_values)
    K = all_k_values(k_idx);
    mask = results_table.K == K;
    fprintf('%-5d  %-12.4f %-12.4f %-12.4f %-12.4f\n', K, ...
        mean(results_table.err_simple(mask)), mean(results_table.err_weighted(mask)), ...
        sqrt(mean(results_table.err_simple(mask).^2)), sqrt(mean(results_table.err_weighted(mask).^2)));
end

valid_y = y(~isnan(y));
global_mean = mean(valid_y);
query_gt = results_table.ground_truth(results_table.K == all_k_values(1));
fprintf('\nGlobal mean baseline: %.4f mm\n', global_mean);
fprintf('Baseline MAE: %.4f mm\n', mean(abs(query_gt - global_mean)));

fprintf('\n═══ PREDICTION RESULTS (DECOMPOSED) ═══\n\n');
fprintf('%-5s  %-12s %-12s %-12s %-12s\n', 'K', 'MAE_simple', 'MAE_weighted', 'RMSE_simple', 'RMSE_weighted');
fprintf('%s\n', repmat('-', 1, 57));

decomp_table = struct2table(decomp_results);
for k_idx = 1:length(all_k_values)
    K = all_k_values(k_idx);
    mask = decomp_table.K == K;
    fprintf('%-5d  %-12.4f %-12.4f %-12.4f %-12.4f\n', K, ...
        mean(decomp_table.err_simple(mask)), mean(decomp_table.err_weighted(mask)), ...
        sqrt(mean(decomp_table.err_simple(mask).^2)), sqrt(mean(decomp_table.err_weighted(mask).^2)));
end

% ========================================================================
%% SECTION 6: TIMING ANALYSIS
% ========================================================================

fprintf('\n');
fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  TIMING ANALYSIS                                               ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

% Filter out zero times (skipped queries)
valid_direct = query_times_direct > 0;
valid_decomp = query_times_decomp > 0;

fprintf('═══ DATA LOADING (one-time cost) ═══\n\n');
fprintf('  Bahn metadata:       %.2f sec\n', time_load_bahn_meta);
fprintf('  Bahn ground truth:   %.2f sec\n', time_load_bahn_gt);
fprintf('  Segment metadata:    %.2f sec\n', time_load_seg_meta);
fprintf('  Segment ground truth:%.2f sec\n', time_load_seg_gt);
fprintf('  Segment lengths:     %.2f sec\n', time_load_seg_len);
fprintf('  Query selection:     %.2f sec\n', time_load_queries);
fprintf('  Feature matrix build:%.2f sec\n', time_build);
fprintf('  ─────────────────────────────\n');
fprintf('  TOTAL SETUP:         %.2f sec\n\n', total_loading_time + time_build);

fprintf('═══ PER-QUERY INFERENCE ═══\n\n');
fprintf('  DIRECT (Bahn-Level):\n');
fprintf('    Mean:    %.4f sec\n', mean(query_times_direct(valid_direct)));
fprintf('    Median:  %.4f sec\n', median(query_times_direct(valid_direct)));
fprintf('    Min:     %.4f sec\n', min(query_times_direct(valid_direct)));
fprintf('    Max:     %.4f sec\n', max(query_times_direct(valid_direct)));
fprintf('    Std:     %.4f sec\n\n', std(query_times_direct(valid_direct)));

fprintf('  DECOMPOSED (Segment-Level):\n');
fprintf('    Mean:    %.4f sec\n', mean(query_times_decomp(valid_decomp)));
fprintf('    Median:  %.4f sec\n', median(query_times_decomp(valid_decomp)));
fprintf('    Min:     %.4f sec\n', min(query_times_decomp(valid_decomp)));
fprintf('    Max:     %.4f sec\n', max(query_times_decomp(valid_decomp)));
fprintf('    Std:     %.4f sec\n', std(query_times_decomp(valid_decomp)));
fprintf('    Mean segments/query: %.1f\n\n', mean(query_num_segments(valid_decomp)));

fprintf('═══ DEPLOYMENT SCENARIOS ═══\n\n');

mean_direct = mean(query_times_direct(valid_direct));
mean_decomp = mean(query_times_decomp(valid_decomp));
setup_time = total_loading_time + time_build;

fprintf('  ┌─────────────────────┬──────────────┬──────────────┐\n');
fprintf('  │ Scenario            │ kNN Metadata │ Retrieval FW │\n');
fprintf('  ├─────────────────────┼──────────────┼──────────────┤\n');
fprintf('  │ Setup (one-time)    │ %8.2f sec │    0.00 sec  │\n', setup_time);
fprintf('  │ 1 query (direct)    │ %8.4f sec │   ~0.62 sec  │\n', mean_direct);
fprintf('  │ 1 query (decomp.)   │ %8.4f sec │   ~0.62 sec  │\n', mean_decomp);
fprintf('  │ 10 queries (decomp.)│ %8.2f sec │   ~6.2 sec   │\n', mean_decomp*10);
fprintf('  │ 100 queries (decomp)│ %8.2f sec │  ~62.0 sec   │\n', mean_decomp*100);
fprintf('  │ Total: 1 query+setup│ %8.2f sec │   ~0.62 sec  │\n', setup_time + mean_decomp);
fprintf('  │ Total:100 q.+setup  │ %8.2f sec │  ~62.0 sec   │\n', setup_time + mean_decomp*100);
fprintf('  └─────────────────────┴──────────────┴──────────────┘\n\n');

fprintf('  Note: Retrieval Framework times from paper (Table X).\n');
fprintf('        kNN requires loading entire feature matrix into memory.\n');
fprintf('        Retrieval Framework uses persistent HNSW index in PostgreSQL.\n');

% Save results
results_filename = sprintf('results/knn_baseline_timed_%s.csv', timestamp);
writetable(results_table, results_filename);

decomp_filename = sprintf('results/knn_decomposed_timed_%s.csv', timestamp);
writetable(decomp_table, decomp_filename);

fprintf('\n✓ Saved: %s\n', results_filename);
fprintf('✓ Saved: %s\n', decomp_filename);
fprintf('✓ Done!\n');
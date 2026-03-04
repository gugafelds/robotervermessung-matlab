%  kNN METADATA BASELINE - WITH TIMING ANALYSIS
%  ========================================================================
%  Baseline comparison: kNN on metadata features (no embeddings, no DTW)
%  Includes timing measurements for deployment comparison
%  ========================================================================

%clear; clc;

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

% === DTW SETTINGS ===
dtw_active = false;
dtw_mode = 'position';
dtw_window = 0.2;
normalize_dtw = false;
use_rotation_alignment = false;
lb_kim_keep_ratio = 1.0;
lb_keogh_candidates = 50;

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

movement_ratio = zeros(height(meta_all), 1);
for i = 1:height(meta_all)
    mt = lower(meta_all.movement_type{i});
    n_l = sum(mt == 'l');
    n_c = sum(mt == 'c');
    total = n_l + n_c;
    if total > 0
        movement_ratio(i) = n_c / total;
    else
        movement_ratio(i) = 0;
    end
end

feature_cols = [movement_ratio, ...
    %meta_all.duration, meta_all.weight, meta_all.length, ...
    meta_all.duration, meta_all.length, ...
    meta_all.min_twist_ist, meta_all.max_twist_ist, meta_all.mean_twist_ist, ...
    meta_all.median_twist_ist, meta_all.std_twist_ist, ...
    meta_all.min_acceleration_ist, meta_all.max_acceleration_ist, meta_all.mean_acceleration_ist, ...
    meta_all.median_acceleration_ist, meta_all.std_acceleration_ist, ...
    %meta_all.position_x, meta_all.position_y, meta_all.position_z];
    ];

feature_cols(isnan(feature_cols)) = 0;

vel_max = 3200.0;
accel_max = 10200.0;
max_length = 9000.0;
max_duration = 25.0;
max_payload = 60.0;
ws_min_x = 400; ws_max_x = 1900;
ws_min_y = -1100; ws_max_y = 1100;
ws_min_z = 400; ws_max_z = 2000;

X = zeros(size(feature_cols));
X(:,1) = feature_cols(:,1);                                            % movement_ratio [0,1]
X(:,2) = min(feature_cols(:,2) / max_duration, 1);                     % duration
%X(:,3) = min(feature_cols(:,3) / max_payload, 1);                      % weight
X(:,3) = min(feature_cols(:,3) / max_length, 1);                       % length
X(:,4) = min(max((feature_cols(:,4) + vel_max) / (2*vel_max), 0), 1);  % min_twist
X(:,5) = min(max((feature_cols(:,5) + vel_max) / (2*vel_max), 0), 1);  % max_twist
X(:,6) = min(max((feature_cols(:,6) + vel_max) / (2*vel_max), 0), 1);  % mean_twist
X(:,7) = min(max((feature_cols(:,7) + vel_max) / (2*vel_max), 0), 1);  % median_twist
X(:,8) = min(feature_cols(:,8) / vel_max, 1);                          % std_twist
X(:,9) = min(max((feature_cols(:,9) + accel_max) / (2*accel_max), 0), 1); % min_accel
X(:,10) = min(max((feature_cols(:,10) + accel_max) / (2*accel_max), 0), 1); % max_accel
X(:,11) = min(max((feature_cols(:,11) + accel_max) / (2*accel_max), 0), 1); % mean_accel
X(:,12) = min(max((feature_cols(:,12) + accel_max) / (2*accel_max), 0), 1); % median_accel
X(:,13) = min(feature_cols(:,13) / accel_max, 1);                      % std_accel
%X(:,15) = min(max((feature_cols(:,15) - ws_min_x) / (ws_max_x - ws_min_x), 0), 1); % pos_x
%X(:,16) = min(max((feature_cols(:,16) - ws_min_y) / (ws_max_y - ws_min_y), 0), 1); % pos_y
%X(:,17) = min(max((feature_cols(:,17) - ws_min_z) / (ws_max_z - ws_min_z), 0), 1); % pos_z

% X = zeros(size(feature_cols));
% X(:,1) = feature_cols(:,1);                                            % movement_ratio [0,1]
% X(:,2) = min(feature_cols(:,2) / max_duration, 1);                     % duration
% X(:,3) = min(feature_cols(:,3) / max_payload, 1);                      % weight
% X(:,4) = min(feature_cols(:,4) / max_length, 1);                       % length
% X(:,5) = min(max((feature_cols(:,5) + vel_max) / (2*vel_max), 0), 1);  % min_twist
% X(:,6) = min(max((feature_cols(:,6) + vel_max) / (2*vel_max), 0), 1);  % max_twist
% X(:,7) = min(max((feature_cols(:,7) + vel_max) / (2*vel_max), 0), 1);  % mean_twist
% X(:,8) = min(max((feature_cols(:,8) + vel_max) / (2*vel_max), 0), 1);  % median_twist
% X(:,9) = min(feature_cols(:,9) / vel_max, 1);                          % std_twist
% X(:,10) = min(max((feature_cols(:,10) + accel_max) / (2*accel_max), 0), 1); % min_accel
% X(:,11) = min(max((feature_cols(:,11) + accel_max) / (2*accel_max), 0), 1); % max_accel
% X(:,12) = min(max((feature_cols(:,12) + accel_max) / (2*accel_max), 0), 1); % mean_accel
% X(:,13) = min(max((feature_cols(:,13) + accel_max) / (2*accel_max), 0), 1); % median_accel
% X(:,14) = min(feature_cols(:,14) / accel_max, 1);                      % std_accel
% X(:,15) = min(max((feature_cols(:,15) - ws_min_x) / (ws_max_x - ws_min_x), 0), 1); % pos_x
% X(:,16) = min(max((feature_cols(:,16) - ws_min_y) / (ws_max_y - ws_min_y), 0), 1); % pos_y
% X(:,17) = min(max((feature_cols(:,17) - ws_min_z) / (ws_max_z - ws_min_z), 0), 1); % pos_z

all_bahn_ids = meta_all.bahn_id;
y = nan(height(meta_all), 1);
for i = 1:height(meta_all)
    idx = find(strcmp(gt_all.segment_id, all_bahn_ids{i}), 1);
    if ~isempty(idx)
        y(i) = gt_all.sidtw_average_distance(idx);
    end
end

% --- Segment features ---
seg_movement_ratio = zeros(height(seg_meta), 1);
for i = 1:height(seg_meta)
    mt = lower(seg_meta.movement_type{i});
    if contains(mt, 'circular')
        seg_movement_ratio(i) = 1.0;
    else
        seg_movement_ratio(i) = 0.0;
    end
end

seg_feature_cols = [seg_movement_ratio, ...
    %seg_meta.duration, seg_meta.weight, seg_meta.length, ...
    seg_meta.duration, seg_meta.length, ...
    seg_meta.min_twist_ist, seg_meta.max_twist_ist, seg_meta.mean_twist_ist, ...
    seg_meta.median_twist_ist, seg_meta.std_twist_ist, ...
    seg_meta.min_acceleration_ist, seg_meta.max_acceleration_ist, seg_meta.mean_acceleration_ist, ...
    seg_meta.median_acceleration_ist, seg_meta.std_acceleration_ist, ...
    %seg_meta.position_x, seg_meta.position_y, seg_meta.position_z];
    ];

seg_feature_cols(isnan(seg_feature_cols)) = 0;

X_seg = zeros(size(seg_feature_cols));
X_seg(:,1) = seg_feature_cols(:,1);
X_seg(:,2) = min(seg_feature_cols(:,2) / max_duration, 1);
%X_seg(:,3) = min(seg_feature_cols(:,3) / max_payload, 1);
X_seg(:,3) = min(seg_feature_cols(:,3) / max_length, 1);
X_seg(:,4) = min(max((seg_feature_cols(:,4) + vel_max) / (2*vel_max), 0), 1);
X_seg(:,5) = min(max((seg_feature_cols(:,5) + vel_max) / (2*vel_max), 0), 1);
X_seg(:,6) = min(max((seg_feature_cols(:,6) + vel_max) / (2*vel_max), 0), 1);
X_seg(:,7) = min(max((seg_feature_cols(:,7) + vel_max) / (2*vel_max), 0), 1);
X_seg(:,8) = min(seg_feature_cols(:,8) / vel_max, 1);
X_seg(:,9) = min(max((seg_feature_cols(:,9) + accel_max) / (2*accel_max), 0), 1);
X_seg(:,10) = min(max((seg_feature_cols(:,10) + accel_max) / (2*accel_max), 0), 1);
X_seg(:,11) = min(max((seg_feature_cols(:,11) + accel_max) / (2*accel_max), 0), 1);
X_seg(:,12) = min(max((seg_feature_cols(:,12) + accel_max) / (2*accel_max), 0), 1);
X_seg(:,13) = min(seg_feature_cols(:,13) / accel_max, 1);
%X_seg(:,15) = min(max((seg_feature_cols(:,15) - ws_min_x) / (ws_max_x - ws_min_x), 0), 1);
%X_seg(:,16) = min(max((seg_feature_cols(:,16) - ws_min_y) / (ws_max_y - ws_min_y), 0), 1);
%X_seg(:,17) = min(max((seg_feature_cols(:,17) - ws_min_z) / (ws_max_z - ws_min_z), 0), 1);

% X_seg = zeros(size(seg_feature_cols));
% X_seg(:,1) = seg_feature_cols(:,1);
% X_seg(:,2) = min(seg_feature_cols(:,2) / max_duration, 1);
% X_seg(:,3) = min(seg_feature_cols(:,3) / max_payload, 1);
% X_seg(:,4) = min(seg_feature_cols(:,4) / max_length, 1);
% X_seg(:,5) = min(max((seg_feature_cols(:,5) + vel_max) / (2*vel_max), 0), 1);
% X_seg(:,6) = min(max((seg_feature_cols(:,6) + vel_max) / (2*vel_max), 0), 1);
% X_seg(:,7) = min(max((seg_feature_cols(:,7) + vel_max) / (2*vel_max), 0), 1);
% X_seg(:,8) = min(max((seg_feature_cols(:,8) + vel_max) / (2*vel_max), 0), 1);
% X_seg(:,9) = min(seg_feature_cols(:,9) / vel_max, 1);
% X_seg(:,10) = min(max((seg_feature_cols(:,10) + accel_max) / (2*accel_max), 0), 1);
% X_seg(:,11) = min(max((seg_feature_cols(:,11) + accel_max) / (2*accel_max), 0), 1);
% X_seg(:,12) = min(max((seg_feature_cols(:,12) + accel_max) / (2*accel_max), 0), 1);
% X_seg(:,13) = min(max((seg_feature_cols(:,13) + accel_max) / (2*accel_max), 0), 1);
% X_seg(:,14) = min(seg_feature_cols(:,14) / accel_max, 1);
% X_seg(:,15) = min(max((seg_feature_cols(:,15) - ws_min_x) / (ws_max_x - ws_min_x), 0), 1);
% X_seg(:,16) = min(max((seg_feature_cols(:,16) - ws_min_y) / (ws_max_y - ws_min_y), 0), 1);
% X_seg(:,17) = min(max((seg_feature_cols(:,17) - ws_min_z) / (ws_max_z - ws_min_z), 0), 1);

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

% === Coverage analysis arrays ===
nearest_distances = nan(length(query_ids),1);
prediction_errors = nan(length(query_ids),1);

coverage_dists = [];
coverage_errors = [];

knn_dtw_results = [];
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

    % === Coverage metric ===
    nearest_distances(q_idx) = sorted_dists(1);
    
    for k_idx = 1:length(all_k_values)
        K = all_k_values(k_idx);
        neighbors = sorted_idx(1:K);
        neighbor_dists = sorted_dists(1:K);
        neighbor_values = y(neighbors);
        
        pred_simple = mean(neighbor_values);
        w = 1 ./ (neighbor_dists + 1e-6);
        w = w / sum(w);
        pred_weighted = sum(w .* neighbor_values);

        if K == 5
            prediction_errors(q_idx) = abs(gt_value - pred_weighted);
        end
        
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

    % === COVERAGE ANALYSIS (Decomposed) ===
    for s = 1:length(query_seg_indices)
        si = query_seg_indices(s);
        seg_id = all_seg_ids{si};
        
        if isnan(y_seg(si)), continue; end
        
        query_seg_feat = X_seg(si, :);
        dists_seg = vecnorm(X_seg - query_seg_feat, 2, 2);
        dists_seg(si) = inf;
        same_bahn = strcmp(all_seg_bahn_ids, query_id);
        dists_seg(same_bahn) = inf;
        dists_seg(isnan(y_seg)) = inf;
        
        [min_d, min_i] = min(dists_seg);
        if isinf(min_d), continue; end
        
        pred_err = abs(y_seg(si) - y_seg(min_i));
        coverage_dists = [coverage_dists; min_d];
        coverage_errors = [coverage_errors; pred_err];
    end
 
    % ════════════════════════════════════════
    % DECOMPOSED kNN + DTW RERANKING
    % ════════════════════════════════════════
    
    if dtw_active
        % For each segment, get kNN candidates, load sequences, run DTW
        K_pool = max(all_k_values); % Use largest K as candidate pool
        
        seg_dtw_preds_simple = zeros(length(query_seg_indices), length(all_k_values));
        seg_dtw_preds_weighted = zeros(length(query_seg_indices), length(all_k_values));
        seg_dtw_weights_len = zeros(length(query_seg_indices), 1);
        seg_dtw_valid = true(length(query_seg_indices), 1);
        
        for s = 1:length(query_seg_indices)
            si = query_seg_indices(s);
            seg_id = all_seg_ids{si};
            
            if isnan(y_seg(si))
                seg_dtw_valid(s) = false;
                continue;
            end
            
            % Load query sequence from DB
            conn_dtw = connectingToPostgres();
            q_pos_sql = sprintf([...
                'SELECT x_soll, y_soll, z_soll ' ...
                'FROM %s.bahn_position_soll ' ...
                'WHERE segment_id = ''%s'' ORDER BY timestamp'], schema, seg_id);
            q_pos = fetch(conn_dtw, q_pos_sql);
            if isempty(q_pos)
                seg_dtw_valid(s) = false;
                close(conn_dtw);
                continue;
            end
            query_seg_seq = [q_pos.x_soll, q_pos.y_soll, q_pos.z_soll];
            
            % kNN: find top K_pool candidates
            query_seg_feat = X_seg(si, :);
            dists_seg = vecnorm(X_seg - query_seg_feat, 2, 2);
            dists_seg(si) = inf;
            same_bahn = strcmp(all_seg_bahn_ids, query_id);
            dists_seg(same_bahn) = inf;
            dists_seg(isnan(y_seg)) = inf;
            [sorted_d, sorted_i] = sort(dists_seg, 'ascend');
            
            if sorted_d(1) == inf
                seg_dtw_valid(s) = false;
                continue;
            end
            
            nb_idx = sorted_i(1:K_pool);
            nb_ids = all_seg_ids(nb_idx);
            
            % Load candidate sequences from DB
            cand_ids_str = sprintf('''%s'',', nb_ids{:});
            cand_ids_str = cand_ids_str(1:end-1);
            cand_pos_sql = sprintf([...
                'SELECT segment_id, x_soll, y_soll, z_soll ' ...
                'FROM %s.bahn_position_soll ' ...
                'WHERE segment_id IN (%s) ' ...
                'ORDER BY segment_id, timestamp'], schema, cand_ids_str);
            cand_pos_data = fetch(conn_dtw, cand_pos_sql);
            close(conn_dtw);
            
            [Gc, gc_ids] = findgroups(cand_pos_data.segment_id);
            cand_vals = [cand_pos_data.x_soll, cand_pos_data.y_soll, cand_pos_data.z_soll];
            cand_map = containers.Map();
            for c = 1:length(gc_ids)
                cand_map(gc_ids{c}) = cand_vals(Gc == c, :);
            end
            
            cand_seqs = cell(K_pool, 1);
            for c = 1:K_pool
                if cand_map.isKey(nb_ids{c})
                    cand_seqs{c} = cand_map(nb_ids{c});
                else
                    cand_seqs{c} = [];
                end
            end
            
            % Build table for performReranking
            cand_table = table(nb_ids(:), sorted_d(1:K_pool), (1:K_pool)', cand_seqs, ...
                'VariableNames', {'segment_id', 'knn_dist', 'rank', 'sequence'});
            
            % Run DTW reranking
            cfg = struct();
            cfg.mode = dtw_mode;
            cfg.window = dtw_window;
            cfg.normalize = normalize_dtw;
            cfg.rot_align = use_rotation_alignment;
            cfg.lb_kim_ratio = lb_kim_keep_ratio;
            cfg.lb_keogh_n = lb_keogh_candidates;
            
            [reranked, ~] = performReranking(cand_table, query_seg_seq, cfg);
            
            % Get ground truth values for reranked candidates
            reranked_gt = nan(K_pool, 1);
            for c = 1:K_pool
                idx = find(strcmp(seg_gt.segment_id, reranked.segment_id{c}), 1);
                if ~isempty(idx)
                    reranked_gt(c) = seg_gt.sidtw_average_distance(idx);
                end
            end
            
            % Compute predictions for each K
            for k_idx = 1:length(all_k_values)
                K = all_k_values(k_idx);
                top_gt = reranked_gt(1:K);
                top_dtw = reranked.dtw_dist(1:K);
                
                valid_k = ~isnan(top_gt) & ~isinf(top_dtw);
                if ~any(valid_k)
                    seg_dtw_preds_simple(s, k_idx) = NaN;
                    seg_dtw_preds_weighted(s, k_idx) = NaN;
                    continue;
                end
                
                seg_dtw_preds_simple(s, k_idx) = mean(top_gt(valid_k));
                w = 1 ./ (top_dtw(valid_k) + 1e-6);
                w = w / sum(w);
                seg_dtw_preds_weighted(s, k_idx) = sum(w .* top_gt(valid_k));
            end
            
            if seg_len_map.isKey(seg_id)
                seg_dtw_weights_len(s) = seg_len_map(seg_id);
            else
                seg_dtw_weights_len(s) = 1;
            end
        end
        
        % Aggregate kNN+DTW predictions
        for k_idx = 1:length(all_k_values)
            K = all_k_values(k_idx);
            
            valid_s = seg_dtw_valid & ~isnan(seg_dtw_preds_simple(:, k_idx));
            if any(valid_s)
                sw = seg_dtw_weights_len(valid_s) / sum(seg_dtw_weights_len(valid_s));
                agg_simple = sum(sw .* seg_dtw_preds_simple(valid_s, k_idx));
                agg_weighted = sum(sw .* seg_dtw_preds_weighted(valid_s, k_idx));
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
            knn_dtw_results = [knn_dtw_results; row];
        end
        
        
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

if dtw_active
    fprintf('\n═══ PREDICTION RESULTS (kNN + DTW) ═══\n\n');
    fprintf('%-5s  %-12s %-12s %-12s %-12s\n', 'K', 'MAE_simple', 'MAE_weighted', 'RMSE_simple', 'RMSE_weighted');
    fprintf('%s\n', repmat('-', 1, 57));
    
    knn_dtw_table = struct2table(knn_dtw_results);
    for k_idx = 1:length(all_k_values)
        K = all_k_values(k_idx);
        mask = knn_dtw_table.K == K;
        fprintf('%-5d  %-12.4f %-12.4f %-12.4f %-12.4f\n', K, ...
            mean(knn_dtw_table.err_simple(mask)), mean(knn_dtw_table.err_weighted(mask)), ...
            sqrt(mean(knn_dtw_table.err_simple(mask).^2)), sqrt(mean(knn_dtw_table.err_weighted(mask).^2)));
    end
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

% ========================================================================
%% SECTION 7: DATABASE COVERAGE ANALYSIS
% ========================================================================

fprintf('\n═══ DATABASE COVERAGE ANALYSIS ═══\n\n');

valid_mask = ~isnan(nearest_distances) & ~isnan(prediction_errors);

distances = nearest_distances(valid_mask);
errors = prediction_errors(valid_mask);

fprintf('Nearest neighbor distance statistics:\n');
fprintf(' Mean: %.4f\n', mean(distances));
fprintf(' Median: %.4f\n', median(distances));
fprintf(' Min: %.4f\n', min(distances));
fprintf(' Max: %.4f\n\n', max(distances));

corr_val = corr(distances, errors);

fprintf('Correlation between NN distance and prediction error: %.3f\n\n', corr_val);

% === Plot 1: Histogram ===
figure;
histogram(distances, 40);
xlabel('Nearest Neighbor Distance (Euclidean)');
ylabel('Frequency');
title('Database Coverage: Distance to Nearest Neighbor');

% === Plot 2: Distance vs Prediction Error ===
figure;
scatter(distances, errors, 25, 'filled');
xlabel('Nearest Neighbor Distance');
ylabel('Prediction Error (mm)');
title('Prediction Error vs Database Coverage');

lsline
grid on

fprintf('\n═══ DATABASE COVERAGE ANALYSIS (Decomposed) ═══\n');
fprintf('Nearest neighbor distance statistics:\n');
fprintf('  Mean:   %.4f\n', mean(coverage_dists));
fprintf('  Median: %.4f\n', median(coverage_dists));
fprintf('  Min:    %.4f\n', min(coverage_dists));
fprintf('  Max:    %.4f\n', max(coverage_dists));
fprintf('Correlation between NN distance and prediction error: %.3f\n', ...
    corr(coverage_dists, coverage_errors));

figure;
subplot(1,2,1);
histogram(coverage_dists, 50);
xlabel('Nearest Neighbor Distance (Euclidean)');
ylabel('Frequency');
title('Database Coverage: Segment-Level');

subplot(1,2,2);
scatter(coverage_dists, coverage_errors, 20, 'filled', 'MarkerFaceAlpha', 0.5);
xlabel('Nearest Neighbor Distance');
ylabel('Prediction Error (mm)');
title(sprintf('Correlation: %.3f', corr(coverage_dists, coverage_errors)));
lsline;


function [sorted_table, stats] = performReranking(input_table, query_seq, cfg)
    n = height(input_table);
    dists = inf(n, 1);
    
    stats.lb_kim_calls = 0;
    stats.lb_keogh_calls = 0;
    stats.dtw_calls = 0;
    
    % === STAGE 2a: LB_Kim ===
    lb_kim_vals = inf(n, 1);
    for i = 1:n
        cand_seq = input_table.sequence{i};
        if isempty(cand_seq) || size(cand_seq,1) < 2, continue; end
        lb_kim_vals(i) = LB_Kim(query_seq, cand_seq, cfg.mode, cfg.rot_align, cfg.normalize);
        stats.lb_kim_calls = stats.lb_kim_calls + 1;
    end
    
    num_keep_kim = ceil(n * cfg.lb_kim_ratio);
    [~, kim_order] = sort(lb_kim_vals, 'ascend');
    kim_survivors = kim_order(1:min(num_keep_kim, n));
    
    % === STAGE 2b: LB_Keogh ===
    lb_keogh_vals = inf(n, 1);
    for i = kim_survivors'
        cand_seq = input_table.sequence{i};
        lb_keogh_vals(i) = LB_Keogh(query_seq, cand_seq, cfg.window, ...
                                    cfg.mode, cfg.rot_align, cfg.normalize);
        stats.lb_keogh_calls = stats.lb_keogh_calls + 1;
    end
    
    [~, keogh_order] = sort(lb_keogh_vals, 'ascend');
    num_keep_keogh = min(cfg.lb_keogh_n, length(kim_survivors));
    keogh_survivors = keogh_order(1:num_keep_keogh);
    
    % === STAGE 2c: DTW ===
    for i = keogh_survivors'
        cand_seq = input_table.sequence{i};
        dists(i) = cDTW(query_seq, cand_seq, cfg.mode, cfg.window, ...
                        Inf, cfg.rot_align, cfg.normalize);
        stats.dtw_calls = stats.dtw_calls + 1;
    end
    
    % === Sortieren & Ranks ===
    input_table.dtw_dist = dists;
    [~, sort_idx] = sort(dists, 'ascend');
    sorted_table = input_table(sort_idx, :);
    sorted_table.Properties.VariableNames{'rank'} = 'rank_s1';
    sorted_table.rank_s2 = (1:n)';
    sorted_table.rank_change = sorted_table.rank_s1 - sorted_table.rank_s2;
end



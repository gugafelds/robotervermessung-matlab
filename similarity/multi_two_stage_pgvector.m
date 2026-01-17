%  TWO-STAGE RETRIEVAL - MULTI-QUERY EVALUATION (CLEAN VERSION)
%  ========================================================================
%  Simple: Run pipeline for N queries, save clean CSVs
%  ✅ NOW USING REAL DTW FUNCTIONS FROM two_stage_retrieval.m
%  ========================================================================

clear; clc;

addpath(genpath(pwd));
addpath(genpath('../main'));
addpath(genpath('../lasertracker'));
addpath(genpath('../methods'));

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  TWO-STAGE RETRIEVAL - MULTI-QUERY EVALUATION                  ║\n');
fprintf('║  ✅ Using REAL DTW with Lower Bounds (same as two_stage)       ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

% ========================================================================
%% USER CONFIGURATION
% ========================================================================

queries_quantity = 250;
K_bahn = 60;
K_segment = 60;
embedding_weights = [1.0; 1.0; 1.0; 1.0; 1.0];
embedding_weights = embedding_weights / sum(embedding_weights);
dtw_mode = 'position';
final_limit_bahn = 50;
final_limit_segment = 50;

dtw_config = struct();
dtw_config.cdtw_window = 0.2;
dtw_config.lb_kim_keep_ratio = 1;  % ✅ CHANGED: Now filters 10% (was 1.0 = no filter)
dtw_config.lb_keogh_candidates = 51; % ✅ CHANGED: Reduce to 100 (was 500)
dtw_config.normalize_dtw = false;
dtw_config.use_rotation_alignment = false;

all_k_values = [5, 10, 25, 50];

fprintf('  Queries: %d\n', queries_quantity);
fprintf('  Top-K test: [%s]\n', strjoin(string(all_k_values), ', '));
fprintf('  ✅ LB_Kim keep ratio: %.1f%% (filters %.0f%%)\n', ...
    dtw_config.lb_kim_keep_ratio*100, (1-dtw_config.lb_kim_keep_ratio)*100);
fprintf('  ✅ LB_Keogh candidates: %d\n\n', dtw_config.lb_keogh_candidates);

% ========================================================================
%% SETUP
% ========================================================================

timestamp = datestr(now, 'yyyymmdd_HHMMSS');

% Create results folder if needed
if ~exist('results', 'dir')
    mkdir('results');
end

results_filename = sprintf('results/two_stage_pgvector_%s.csv', timestamp);

% Sample queries
conn = connectingToPostgres();
query_sql = sprintf(['SELECT bahn_id FROM (' ...
    'SELECT DISTINCT b.bahn_id FROM bewegungsdaten.bahn_metadata b ' ...
    'INNER JOIN auswertung.info_sidtw s ON b.bahn_id = s.segment_id ' ...
    'WHERE b.segment_id = b.bahn_id) AS distinct_bahnen ' ...
    'ORDER BY RANDOM() LIMIT %d'], queries_quantity);
query_results = fetch(conn, query_sql);
query_ids = query_results.bahn_id;
close(conn);

fprintf('✓ Sampled %d queries\n', length(query_ids));
fprintf('✓ Output: %s\n\n', results_filename);

% ========================================================================
%% PROCESS QUERIES
% ========================================================================

bahn_results_all = [];
segment_results_all = [];

overall_start = tic;

for q_idx = 1:length(query_ids)
    query_id = query_ids{q_idx};
    
    fprintf('Query %d/%d: %s\n', q_idx, length(query_ids), query_id);
    
    try
        conn = connectingToPostgres();
        schema = 'bewegungsdaten';
        
        % === LOAD DATA ===
        query_data = loadQueryDataHierarchical(conn, schema, query_id);
        if isempty(query_data)
            close(conn);
            continue;
        end
        num_query_segments = length(query_data.segment_ids);
        
        % === BAHN-LEVEL ===
        % Stage 1
        stage1_bahn_start = tic;
        query_embeddings = getQueryEmbeddings(conn, schema, query_id, 'bahn');
        bahn_candidates = performRRFEmbeddingSearch(conn, schema, query_embeddings, ...
            embedding_weights, K_bahn, 'bahn');
        
        % ✅ REMOVE QUERY FROM CANDIDATES (if present)
        bahn_candidates = bahn_candidates(~strcmp(bahn_candidates, query_id));
        
        stage1_bahn_time = toc(stage1_bahn_start);
        
        % Stage 2 - SEPARATE loading from DTW computation!
        % Loading (not counted in DTW time for fair comparison)
        load_start = tic;
        bahn_candidate_data = loadCandidateDataBatch(conn, schema, bahn_candidates, 'bahn');
        load_time_bahn = toc(load_start);
        
        if strcmp(dtw_mode, 'position')
            query_seq = query_data.trajectory.position;
        else
            query_seq = query_data.trajectory.joint;
        end
        
        % DTW computation only (✅ NOW WITH REAL LOWER BOUNDS!)
        dtw_start = tic;
        [bahn_results, bahn_dtw_stats] = performDTWReranking_REAL(query_seq, bahn_candidate_data, ...
            dtw_mode, dtw_config, final_limit_bahn);
        stage2_bahn_time = toc(dtw_start);  % ✅ ONLY DTW time!
        
        fprintf('  Bahn: S1=%.2fs, Load=%.2fs, DTW=%.2fs, DTW calls=%d/%d (saved %.1f%%)\n', ...
            stage1_bahn_time, load_time_bahn, stage2_bahn_time, ...
            bahn_dtw_stats.dtw_calls, K_bahn, ...
            (K_bahn - bahn_dtw_stats.dtw_calls)/K_bahn*100);
        
        % Prognosis
        [~, s1_idx] = sort([bahn_results.stage1_rank]);
        bahn_s1 = bahn_results(s1_idx);
        bahn_s2 = bahn_results;
        
        top_k_test = all_k_values(all_k_values <= length(bahn_results));
        prog_s1 = evaluatePerformancePrognosis(conn, query_id, bahn_s1, top_k_test);
        prog_s2 = evaluatePerformancePrognosis(conn, query_id, bahn_s2, top_k_test);
        
        % Extract BAHN metrics
        if ~isempty(prog_s1) && isfield(prog_s1, 'predictions')
            k_fields = fieldnames(prog_s1.predictions);
            for k_idx = 1:length(k_fields)
                k_field = k_fields{k_idx};
                K = prog_s1.predictions.(k_field).k;
                
                s1_err = prog_s1.predictions.(k_field).errors.average_distance;
                s2_err = prog_s2.predictions.(k_field).errors.average_distance;
                
                weighting_impr = (s1_err.mae_mean - s1_err.mae_weighted) / s1_err.mae_mean * 100;
                dtw_impr = (s1_err.mae_weighted - s2_err.mae_weighted) / s1_err.mae_weighted * 100;
                time_overhead = (stage2_bahn_time - stage1_bahn_time) / stage1_bahn_time * 100;
                efficiency = dtw_impr / time_overhead;
                
                row = struct();
                row.query_id = query_id;
                row.segment_id = query_id;
                row.level = 'bahn';
                row.length = size(query_data.trajectory.position, 1);
                row.K = K;
                row.actual_sidtw = prog_s1.actual.average_distance;
                row.s1_simple_mae = s1_err.mae_mean;
                row.s1_weighted_mae = s1_err.mae_weighted;
                row.s2_simple_mae = s2_err.mae_mean;
                row.s2_weighted_mae = s2_err.mae_weighted;
                row.stage1_time_sec = stage1_bahn_time;
                row.stage2_time_sec = stage2_bahn_time;
                row.total_time_sec = stage1_bahn_time + stage2_bahn_time;
                row.weighting_improvement_pct = weighting_impr;
                row.dtw_improvement_pct = dtw_impr;
                row.time_overhead_pct = time_overhead;
                row.efficiency = efficiency;
                row.dtw_calls_actual = bahn_dtw_stats.dtw_calls;
                row.lb_kim_filtered = bahn_dtw_stats.lb_kim_filtered;
                row.lb_keogh_filtered = bahn_dtw_stats.lb_keogh_filtered;
                
                bahn_results_all = [bahn_results_all; row];
            end
        end
        
        % === SEGMENT-LEVEL ===
        for seg_idx = 1:num_query_segments
            seg_id = query_data.segment_ids{seg_idx};
            
            % Stage 1
            stage1_seg_start = tic;
            seg_embeddings = getQueryEmbeddings(conn, schema, seg_id, 'segment');
            seg_candidates = performRRFEmbeddingSearch(conn, schema, seg_embeddings, ...
                embedding_weights, K_segment, 'segment');
            
            % ✅ REMOVE QUERY SEGMENT FROM CANDIDATES (if present)
            seg_candidates = seg_candidates(~strcmp(seg_candidates, seg_id));
            
            stage1_seg_time = toc(stage1_seg_start);
            
            if isempty(seg_candidates)
                continue;
            end
            
            % Stage 2 - SEPARATE loading from DTW!
            load_seg_start = tic;
            seg_candidate_data = loadCandidateDataBatch(conn, schema, seg_candidates, 'segment');
            load_time_seg = toc(load_seg_start);
            
            if strcmp(dtw_mode, 'position')
                query_seg_seq = query_data.segments.position{seg_idx};
            else
                query_seg_seq = query_data.segments.joint{seg_idx};
            end
            
            % DTW computation only (✅ NOW WITH REAL LOWER BOUNDS!)
            dtw_seg_start = tic;
            [seg_results, seg_dtw_stats] = performDTWReranking_REAL(query_seg_seq, seg_candidate_data, ...
                dtw_mode, dtw_config, final_limit_segment);
            stage2_seg_time = toc(dtw_seg_start);  % ✅ ONLY DTW time!
            
            % Prognosis
            [~, s1_idx_seg] = sort([seg_results.stage1_rank]);
            seg_s1 = seg_results(s1_idx_seg);
            seg_s2 = seg_results;
            
            prog_seg_s1 = evaluatePerformancePrognosis(conn, seg_id, seg_s1, top_k_test);
            prog_seg_s2 = evaluatePerformancePrognosis(conn, seg_id, seg_s2, top_k_test);
            
            % Extract SEGMENT metrics
            if ~isempty(prog_seg_s1) && isfield(prog_seg_s1, 'predictions')
                k_fields_seg = fieldnames(prog_seg_s1.predictions);
                for k_idx = 1:length(k_fields_seg)
                    k_field = k_fields_seg{k_idx};
                    K = prog_seg_s1.predictions.(k_field).k;
                    
                    s1_err = prog_seg_s1.predictions.(k_field).errors.average_distance;
                    s2_err = prog_seg_s2.predictions.(k_field).errors.average_distance;
                    
                    weighting_impr = (s1_err.mae_mean - s1_err.mae_weighted) / s1_err.mae_mean * 100;
                    dtw_impr = (s1_err.mae_weighted - s2_err.mae_weighted) / s1_err.mae_weighted * 100;
                    time_overhead = (stage2_seg_time - stage1_seg_time) / stage1_seg_time * 100;
                    efficiency = dtw_impr / time_overhead;
                    
                    row = struct();
                    row.query_id = query_id;
                    row.segment_id = seg_id;
                    row.level = 'segment';
                    row.length = size(query_data.segments.position{seg_idx}, 1);
                    row.K = K;
                    row.actual_sidtw = prog_seg_s1.actual.average_distance;
                    row.s1_simple_mae = s1_err.mae_mean;
                    row.s1_weighted_mae = s1_err.mae_weighted;
                    row.s2_simple_mae = s2_err.mae_mean;
                    row.s2_weighted_mae = s2_err.mae_weighted;
                    row.stage1_time_sec = stage1_seg_time;
                    row.stage2_time_sec = stage2_seg_time;
                    row.total_time_sec = stage1_seg_time + stage2_seg_time;
                    row.weighting_improvement_pct = weighting_impr;
                    row.dtw_improvement_pct = dtw_impr;
                    row.time_overhead_pct = time_overhead;
                    row.efficiency = efficiency;
                    row.dtw_calls_actual = seg_dtw_stats.dtw_calls;
                    row.lb_kim_filtered = seg_dtw_stats.lb_kim_filtered;
                    row.lb_keogh_filtered = seg_dtw_stats.lb_keogh_filtered;
                    
                    segment_results_all = [segment_results_all; row];
                end
            end
        end
        
        close(conn);
        fprintf('  ✓ Complete\n\n');
        
    catch ME
        warning('Query %s failed: %s', query_id, ME.message);
        try close(conn); catch; end
    end
end

overall_time = toc(overall_start);

% ========================================================================
%% SAVE RESULTS
% ========================================================================

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  SAVING RESULTS                                                ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

all_results = [];

if ~isempty(bahn_results_all)
    all_results = [all_results; bahn_results_all];
end

if ~isempty(segment_results_all)
    all_results = [all_results; segment_results_all];
end

if ~isempty(all_results)
    results_table = struct2table(all_results);
    writetable(results_table, results_filename);
    fprintf('✓ %s\n', results_filename);
    fprintf('  Total: %d rows\n', height(results_table));
    fprintf('  - Bahn:    %d rows\n', length(bahn_results_all));
    fprintf('  - Segment: %d rows\n', length(segment_results_all));
else
    fprintf('⚠ No results to save\n');
end

fprintf('\n');
fprintf('Time: %s\n', formatDuration(overall_time));
fprintf('\n✓ Done!\n');

% ========================================================================
%% HELPER FUNCTIONS
% ========================================================================

function str = formatDuration(seconds)
    hours = floor(seconds / 3600);
    minutes = floor(mod(seconds, 3600) / 60);
    secs = floor(mod(seconds, 60));
    if hours > 0
        str = sprintf('%dh %dm %ds', hours, minutes, secs);
    elseif minutes > 0
        str = sprintf('%dm %ds', minutes, secs);
    else
        str = sprintf('%.1fs', seconds);
    end
end

% ============================================================================
% ✅ REAL DTW FUNCTION (from two_stage_retrieval.m - EXACT COPY!)
% ============================================================================

function [results, stats] = performDTWReranking_REAL(query_seq, candidate_data, dtw_mode, dtw_config, limit)
    % REAL DTW Reranking with ACTUAL Lower Bounds
    % Same implementation as two_stage_retrieval.m
    
    K = length(candidate_data.ids);
    
    if strcmp(dtw_mode, 'position')
        cand_seqs = candidate_data.position;
    else
        cand_seqs = candidate_data.joint;
    end
    
    stage1_ranks = (1:K)';
    
    % Validate candidates
    valid_mask = false(K, 1);
    for i = 1:K
        if ~isempty(cand_seqs{i}) && size(cand_seqs{i}, 1) >= 2
            valid_mask(i) = true;
        end
    end
    
    valid_indices = find(valid_mask);
    num_valid = length(valid_indices);
    
    if num_valid == 0
        results = [];
        stats = struct('dtw_calls', 0, 'lb_kim_filtered', 0, 'lb_keogh_filtered', 0, 'valid_candidates', 0);
        return;
    end
    
    % ================================================================
    % ✅ STEP 1: LB_KIM FILTERING (REAL!)
    % ================================================================
    lb_kim_keep = round(num_valid * dtw_config.lb_kim_keep_ratio);
    lb_kim_dists = zeros(num_valid, 1);
    
    for i = 1:num_valid
        idx = valid_indices(i);
        seq = cand_seqs{idx};
        
        % LB_Kim: min distance using first/last points
        d1 = norm(query_seq(1,:) - seq(1,:));
        d2 = norm(query_seq(end,:) - seq(end,:));
        lb_kim_dists(i) = max(d1, d2);
    end
    
    [~, lb_kim_sort] = sort(lb_kim_dists);
    lb_kim_survivors = valid_indices(lb_kim_sort(1:lb_kim_keep));
    lb_kim_filtered = num_valid - lb_kim_keep;
    
    % ================================================================
    % ✅ STEP 2: LB_KEOGH FILTERING (REAL!)
    % ================================================================
    lb_keogh_target = min(dtw_config.lb_keogh_candidates, length(lb_kim_survivors));
    
    if length(lb_kim_survivors) > lb_keogh_target
        lb_keogh_dists = zeros(length(lb_kim_survivors), 1);
        
        % Compute envelope
        window = dtw_config.cdtw_window;
        w = max(1, round(size(query_seq, 1) * window));
        
        U = zeros(size(query_seq));
        L = zeros(size(query_seq));
        
        for d = 1:size(query_seq, 2)
            for t = 1:size(query_seq, 1)
                t_start = max(1, t - w);
                t_end = min(size(query_seq, 1), t + w);
                U(t, d) = max(query_seq(t_start:t_end, d));
                L(t, d) = min(query_seq(t_start:t_end, d));
            end
        end
        
        for i = 1:length(lb_kim_survivors)
            idx = lb_kim_survivors(i);
            cand = cand_seqs{idx};
            
            % Interpolate if needed
            if size(cand, 1) ~= size(query_seq, 1)
                cand_interp = zeros(size(query_seq));
                for d = 1:size(query_seq, 2)
                    cand_interp(:, d) = interp1(1:size(cand, 1), cand(:, d), ...
                        linspace(1, size(cand, 1), size(query_seq, 1)), 'linear');
                end
                cand = cand_interp;
            end
            
            % LB_Keogh
            lb_dist = 0;
            for t = 1:size(cand, 1)
                for d = 1:size(cand, 2)
                    if cand(t, d) > U(t, d)
                        lb_dist = lb_dist + (cand(t, d) - U(t, d))^2;
                    elseif cand(t, d) < L(t, d)
                        lb_dist = lb_dist + (L(t, d) - cand(t, d))^2;
                    end
                end
            end
            
            lb_keogh_dists(i) = sqrt(lb_dist);
        end
        
        [~, lb_keogh_sort] = sort(lb_keogh_dists);
        dtw_candidates = lb_kim_survivors(lb_keogh_sort(1:lb_keogh_target));
        lb_keogh_filtered = length(lb_kim_survivors) - lb_keogh_target;
    else
        dtw_candidates = lb_kim_survivors;
        lb_keogh_filtered = 0;
    end
    
    % ================================================================
    % ✅ STEP 3: FULL DTW (ONLY ON SURVIVORS!)
    % ================================================================
    dtw_dists = inf(K, 1);
    dtw_calls_actual = 0;
    
    for i = 1:length(dtw_candidates)
        idx = dtw_candidates(i);
        
        try
            dtw_dists(idx) = cDTW(query_seq, cand_seqs{idx}, dtw_mode, dtw_config.cdtw_window, ...
                inf, dtw_config.use_rotation_alignment, dtw_config.normalize_dtw);
            dtw_calls_actual = dtw_calls_actual + 1;
        catch
            dtw_dists(idx) = inf;
        end
    end
    
    [sorted_dists, sort_idx] = sort(dtw_dists);
    
    % Build results
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
        result.stage2_rank = i;
        result.stage1_rank = stage1_ranks(idx);
        result.rank_change = result.stage1_rank - i;
        
        results = [results; result];
    end
    
    % ✅ REAL STATS!
    stats = struct();
    stats.dtw_calls = dtw_calls_actual;
    stats.lb_kim_filtered = lb_kim_filtered;
    stats.lb_keogh_filtered = lb_keogh_filtered;
    stats.valid_candidates = num_valid;
end


% ========================================================================
%% HELPER FUNCTIONS FOR MULTI-QUERY MODE
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
    % Evaluate performance prognosis quality using k-NN Regression
    % ✅ EXCLUDES query from results to prevent overfitting!
    
    if isempty(results)
        prognosis_metrics = struct();
        return;
    end
    
    % ================================================================
    % 1. GET ACTUAL SIDTW VALUES (Ground Truth) - FIRST!
    % ================================================================
    
    actual_sidtw = getSIDTWMetrics(conn, query_id);
    
    if isempty(actual_sidtw)
        warning('No SIDTW data found for query %s', query_id);
        prognosis_metrics = struct();
        prognosis_metrics.error = 'No SIDTW data for query';
        return;
    end
    
    % ================================================================
    % 2. FILTER OUT QUERY FROM RESULTS
    % ================================================================
    
    % ✅ CRITICAL: Remove query itself from candidates!
    filtered_results = [];
    n_original = length(results);
    
    for i = 1:n_original
        result_id = results(i).segment_id;
        
        % Skip if this is the query itself
        if strcmp(result_id, query_id)
            fprintf('    ℹ Removed query %s from candidates (was at rank %d)\n', ...
                query_id, results(i).stage2_rank);
            continue;
        end
        
        filtered_results = [filtered_results; results(i)];
    end
    
    % Use filtered results
    n = length(filtered_results);
    
    if n == 0
        warning('No valid candidates after filtering query');
        prognosis_metrics = struct();
        prognosis_metrics.actual = actual_sidtw;  % ✅ Store actual even if no candidates
        return;
    end
    
    % ================================================================
    % 3. GET SIDTW FOR FILTERED RESULTS
    % ================================================================
    
    result_ids = cell(n, 1);
    for i = 1:n
        result_ids{i} = filtered_results(i).segment_id;
    end
    
    result_sidtw = getSIDTWMetricsBatch(conn, result_ids);
    
    % ================================================================
    % 4. FOR EACH TOP-K: PREDICT & MEASURE ERROR
    % ================================================================
    
    prognosis_metrics = struct();
    prognosis_metrics.actual = actual_sidtw;  % ✅ Store actual BEFORE predictions
    prognosis_metrics.predictions = struct();
    
    for k_idx = 1:length(top_k_values)
        K = top_k_values(k_idx);
        
        if K > n
            continue;  % Skip if not enough results
        end
        
        % Get Top-K results (from filtered!)
        top_k_results = filtered_results(1:K);
        
        % ============================================================
        % AGGREGATE SIDTW + DTW DISTANCES FROM TOP-K
        % ============================================================
        
        top_k_sidtw_values = struct();
        top_k_sidtw_values.min_distance = [];
        top_k_sidtw_values.max_distance = [];
        top_k_sidtw_values.average_distance = [];
        top_k_sidtw_values.standard_deviation = [];
        
        % Collect DTW distances
        top_k_dtw_dists = zeros(K, 1);
        
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
            
            % Get DTW distance
            top_k_dtw_dists(i) = top_k_results(i).dtw_distance;
        end
        
        % ============================================================
        % COMPUTE PREDICTIONS: SIMPLE + WEIGHTED
        % ============================================================
        
        prediction = struct();
        metrics = {'min_distance', 'max_distance', 'average_distance', 'standard_deviation'};
        
        for m_idx = 1:length(metrics)
            metric_name = metrics{m_idx};
            values = top_k_sidtw_values.(metric_name);
            
            if isempty(values)
                prediction.(metric_name).mean = NaN;
                prediction.(metric_name).weighted_mean = NaN;
                prediction.(metric_name).median = NaN;
                prediction.(metric_name).min = NaN;
                prediction.(metric_name).max = NaN;
                prediction.(metric_name).std = NaN;
                continue;
            end
            
            % ========================================================
            % METHOD 1: SIMPLE MEAN (Baseline)
            % ========================================================
            prediction.(metric_name).mean = mean(values);
            
            % ========================================================
            % METHOD 2: WEIGHTED MEAN (Dudani 1976)
            % ========================================================
            epsilon = 1e-6;
            
            % Get corresponding DTW distances (only for valid SIDTW values)
            valid_dtw_dists = top_k_dtw_dists(1:length(values));
            
            % Calculate weights: weight_i = 1 / (distance_i + epsilon)
            weights_dtw = 1 ./ (valid_dtw_dists + epsilon);
            
            % Normalize weights to sum to 1
            weights_dtw = weights_dtw / sum(weights_dtw);
            
            % Weighted average
            prediction.(metric_name).weighted_mean = sum(weights_dtw .* values);
            
            % ========================================================
            % ADDITIONAL STATISTICS
            % ========================================================
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
                errors.(metric_name).mae_weighted = NaN;
                errors.(metric_name).mae_median = NaN;
                errors.(metric_name).relative_error_mean = NaN;
                errors.(metric_name).relative_error_weighted = NaN;
                errors.(metric_name).relative_error_median = NaN;
                continue;
            end
            
            % === SIMPLE MEAN ERRORS ===
            errors.(metric_name).mae_mean = abs(prediction.(metric_name).mean - actual_value);
            
            if actual_value ~= 0
                errors.(metric_name).relative_error_mean = ...
                    abs(prediction.(metric_name).mean - actual_value) / actual_value * 100;
            else
                errors.(metric_name).relative_error_mean = NaN;
            end
            
            % === WEIGHTED MEAN ERRORS ===
            errors.(metric_name).mae_weighted = abs(prediction.(metric_name).weighted_mean - actual_value);
            
            if actual_value ~= 0
                errors.(metric_name).relative_error_weighted = ...
                    abs(prediction.(metric_name).weighted_mean - actual_value) / actual_value * 100;
            else
                errors.(metric_name).relative_error_weighted = NaN;
            end
            
            % === MEDIAN ERRORS ===
            errors.(metric_name).mae_median = abs(prediction.(metric_name).median - actual_value);
            
            if actual_value ~= 0
                errors.(metric_name).relative_error_median = ...
                    abs(prediction.(metric_name).median - actual_value) / actual_value * 100;
            else
                errors.(metric_name).relative_error_median = NaN;
            end
            
            % === STORE PREDICTIONS FOR REFERENCE ===
            errors.(metric_name).actual = actual_value;
            errors.(metric_name).predicted_mean = prediction.(metric_name).mean;
            errors.(metric_name).predicted_weighted = prediction.(metric_name).weighted_mean;
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
%  TWO-STAGE RETRIEVAL - MULTI-QUERY EVALUATION (CLEAN VERSION)
%  ========================================================================
%  Based on similarity_search.m - runs pipeline for N queries, saves CSV
%  ========================================================================

clear; clc;

addpath(genpath(pwd));
addpath(genpath('../main'));
addpath(genpath('../lasertracker'));
addpath(genpath('../methods'));

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  TWO-STAGE RETRIEVAL - MULTI-QUERY EVALUATION                  ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

% ========================================================================
%% SECTION 1: CONFIGURATION
% ========================================================================

fprintf('═══ CONFIGURATION ═══\n\n');

% === QUERY SETTINGS ===
queries_quantity = 100;
random_seed = 21;  % Für Reproduzierbarkeit

% === EMBEDDING SETTINGS ===
K = 50;  % Stage 1 candidates
rrf_k = 60;
weights_conf = {    [1.0; 0.0; 0.0; 0.0; 0.0];
                    [1.0; 1.0; 0.0; 0.0; 0.0];
                    [1.0; 0.0; 1.0; 0.0; 0.0];
                    [1.0; 0.0; 0.0; 1.0; 0.0];
                    [1.0; 0.0; 0.0; 0.0; 1.0];
                    [1.0; 1.0; 1.0; 1.0; 1.0];
                }; 

% === DTW SETTINGS ===
dtw_mode = 'position';  % 'position' oder 'joint_states'
dtw_window = 0.2;
normalize_dtw = false;
use_rotation_alignment = false;

% === LOWER BOUNDS ===
lb_kim_keep_ratio = 1.0;
lb_keogh_candidates = 50;

% === PROGNOSE SETTINGS ===
all_k_values = [5, 10, 25, 50];

fprintf('Queries:           %d\n', queries_quantity);
fprintf('K Candidates:      %d\n', K);
fprintf('Top-K test:        [%s]\n', strjoin(string(all_k_values), ', '));
fprintf('DTW Mode:          %s\n', dtw_mode);
fprintf('LB_Kim ratio:      %.0f%%\n', lb_kim_keep_ratio * 100);
fprintf('LB_Keogh target:   %d\n\n', lb_keogh_candidates);

% ========================================================================
%% SECTION 2: SETUP
% ========================================================================
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
if ~exist('results', 'dir')
    mkdir('results');
end
results_filename = sprintf('results/similarity_search_%s.csv', timestamp);

% === EXCLUSION LIST (GT / Noisy Trajectories) ===
% Diese IDs sollen NICHT im Random Sampling landen, um Data Leakage zu vermeiden
exclude_ids = {
    %% clean - 0 mm
    '1765989370'; % clean, np = 3 / 10 GT
    '1765989294'; % clean, np = 3 / 20 GT
    '1765988821'; % clean, np = 3 / 30 GT
    '1765988920'; % clean; np = 3 / 40 GT
    '1765989411'; % clean; np = 3 / 50 GT
    
    %% noisy - 2 mm

    '1765990630'; % noisy; np = 3 / 10 GT
    '1765990747'; % noisy; np = 3 / 20 GT
    '1765990822'; % noisy; np = 3 / 30 GT
    '1765991047'; % noisy; np = 3 / 40 GT
    '1765991234'; % noisy; np = 3 / 50 GT
    
    %% noisy - 5 mm

    '1765991190'; % noisy; np = 3 / 10 GT
    '1765991445'; % noisy; np = 3 / 20 GT
    '1765991515'; % noisy; np = 3 / 30 GT
    '1765991949'; % noisy; np = 3 / 40 GT 
    '1765991743'; % noisy; np = 3 / 50 GT

    %% clean - 0 mm
    '1769770498'; % clean, np = 3 / 10 GT
    '1769770684'; % clean, np = 3 / 20 GT
    '1769770935'; % clean, np = 3 / 30 GT
    '1769771107'; % clean; np = 3 / 40 GT
    '1769771447'; % clean; np = 3 / 50 GT
    
    %% noisy - 2 mm

    '1769773928'; % noisy; np = 3 / 10 GT
    '1769772060'; % noisy; np = 3 / 20 GT
    '1769772213'; % noisy; np = 3 / 30 GT
    '1769773985'; % noisy; np = 3 / 40 GT
    '1769774278'; % noisy; np = 3 / 50 GT
    
    %% noisy - 5 mm

    '1769772609'; % noisy; np = 3 / 10 GT
    '1769773593'; % noisy; np = 3 / 20 GT
    '1769772776'; % noisy; np = 3 / 30 GT
    '1769772900'; % noisy; np = 3 / 39 (40) GT 
    '1769773333'; % noisy; np = 3 / 50 GT

    %% noise - 10 mm

    '1769774581'; %noise; np = 3 / 10 GT
};

% Formatieren für SQL: 'id1', 'id2', 'id3' ...
exclude_str = sprintf('''%s'',', exclude_ids{:});
exclude_str = exclude_str(1:end-1); % Letztes Komma entfernen

% Sample queries
conn = connectingToPostgres();

% SQL Query mit NOT IN Filter
query_sql = sprintf(['SELECT bahn_id FROM (' ...
    'SELECT DISTINCT b.bahn_id FROM bewegungsdaten.bahn_metadata b ' ...
    'INNER JOIN auswertung.info_sidtw s ON b.bahn_id = s.segment_id ' ...
    'WHERE b.segment_id = b.bahn_id ' ...
    'AND b.bahn_id NOT IN (' ...
        'SELECT t_all.bahn_id ' ...
        'FROM bewegungsdaten.bahn_info t_exclude ' ...
        'JOIN bewegungsdaten.bahn_info t_all ON t_exclude.record_filename = t_all.record_filename ' ... % <--- HIER GEÄNDERT
        'WHERE t_exclude.bahn_id IN (%s) ' ...
    ') ' ...
    ') AS distinct_bahnen ' ...
    'ORDER BY md5(bahn_id || ''%d'') LIMIT %d'], ...
    exclude_str, random_seed, queries_quantity);

query_results = fetch(conn, query_sql);
query_ids = query_results.bahn_id;
close(conn);

fprintf('✓ Sampled %d queries (excluded %d GT/Noisy IDs)\n', ...
    length(query_ids), length(exclude_ids));
fprintf('✓ Output: %s\n\n', results_filename);

% ========================================================================
%% SECTION 3: MAIN LOOP
% ========================================================================

all_results = [];
schema = 'bewegungsdaten';

overall_start = tic;

for q_idx = 1:length(query_ids)
    query_id = query_ids{q_idx};
    fprintf('Query %d/%d: %s\n', q_idx, length(query_ids), query_id);
    
    try
        conn = connectingToPostgres();
        
        % ================================================================
        % LOAD QUERY DATA
        % ================================================================
        pos_sql = sprintf(['SELECT x_soll, y_soll, z_soll FROM %s.bahn_position_soll ' ...
                          'WHERE bahn_id = ''%s'' ORDER BY timestamp'], schema, query_id);
        pos = fetch(conn, pos_sql);
        
        if isempty(pos)
            close(conn);
            continue;
        end
        
        query_bahn_position = [pos.x_soll, pos.y_soll, pos.z_soll];
        
        joint_sql = sprintf(['SELECT joint_1, joint_2, joint_3, joint_4, joint_5, joint_6 ' ...
                            'FROM %s.bahn_joint_states WHERE bahn_id = ''%s'' ORDER BY timestamp'], ...
                            schema, query_id);
        joint = fetch(conn, joint_sql);
        query_bahn_joint = [joint.joint_1, joint.joint_2, joint.joint_3, ...
                           joint.joint_4, joint.joint_5, joint.joint_6];
        
        % Load segments
        seg_sql = sprintf(['SELECT DISTINCT segment_id FROM %s.bahn_metadata ' ...
                          'WHERE bahn_id = ''%s'' AND segment_id != bahn_id ORDER BY segment_id'], ...
                          schema, query_id);
        segs = fetch(conn, seg_sql);
        query_segment_ids = segs.segment_id;
        num_segments = length(query_segment_ids);
        
        query_segments_position = cell(num_segments, 1);
        query_segments_joint = cell(num_segments, 1);
        
        for i = 1:num_segments
            seg_id = query_segment_ids{i};
            pos_seg_sql = sprintf(['SELECT x_soll, y_soll, z_soll FROM %s.bahn_position_soll ' ...
                                  'WHERE segment_id = ''%s'' ORDER BY timestamp'], schema, seg_id);
            pos_seg = fetch(conn, pos_seg_sql);
            query_segments_position{i} = [pos_seg.x_soll, pos_seg.y_soll, pos_seg.z_soll];
            
            joint_seg_sql = sprintf(['SELECT joint_1, joint_2, joint_3, joint_4, joint_5, joint_6 ' ...
                                    'FROM %s.bahn_joint_states WHERE segment_id = ''%s'' ORDER BY timestamp'], ...
                                    schema, seg_id);
            joint_seg = fetch(conn, joint_seg_sql);
            query_segments_joint{i} = [joint_seg.joint_1, joint_seg.joint_2, joint_seg.joint_3, ...
                                      joint_seg.joint_4, joint_seg.joint_5, joint_seg.joint_6];
        end
        
        if strcmp(dtw_mode, 'position')
            query_bahn_seq = query_bahn_position;
            query_segment_seqs = query_segments_position;
        else
            query_bahn_seq = query_bahn_joint;
            query_segment_seqs = query_segments_joint;
        end
            
        for i = 1:length(weights_conf)
            weights = weights_conf{i};
            % ================================================================
            % BAHN-LEVEL: STAGE 1
            % ================================================================
            stage1_bahn_start = tic;
            
            emb_sql = sprintf(['SELECT position_embedding::text, joint_embedding::text, ' ...
                              'orientation_embedding::text, velocity_embedding::text, ' ...
                              'metadata_embedding::text ' ...
                              'FROM %s.bahn_embeddings WHERE segment_id = ''%s'''], schema, query_id);
            emb_result = fetch(conn, emb_sql);
            
            query_embeddings = {
                parseEmbedding(emb_result.position_embedding), ...
                parseEmbedding(emb_result.joint_embedding), ...
                parseEmbedding(emb_result.orientation_embedding), ...
                parseEmbedding(emb_result.velocity_embedding), ...
                parseEmbedding(emb_result.metadata_embedding)
            };
            
            modalities = {'position', 'joint', 'orientation', 'velocity', 'metadata'};
            
            execute(conn, 'SET hnsw.ef_search = 200');
            execute(conn, 'SET search_path = "bewegungsdaten"');
            
            rrf_scores = containers.Map();
                       
            for m = 1:5
                if weights(m) == 0 || isempty(query_embeddings{m}), continue; end
                
                emb_str = sprintf('[%s]', strjoin(string(query_embeddings{m}), ','));
                col_name = sprintf('%s_embedding', modalities{m});
                
                search_sql = sprintf(['SELECT segment_id, %s <=> ''%s''::vector as distance ' ...
                                     'FROM %s.bahn_embeddings ' ...
                                     'WHERE segment_id = bahn_id AND segment_id != ''%s'' AND %s IS NOT NULL ' ...
                                     'ORDER BY distance LIMIT %d'], ...
                                     col_name, emb_str, schema, query_id, col_name, K*10);
                
                result = fetch(conn, search_sql);
                
                for rank = 1:height(result)
                    id = result.segment_id{rank};
                    if ~rrf_scores.isKey(id), rrf_scores(id) = 0; end
                    rrf_scores(id) = rrf_scores(id) + weights(m) / (rrf_k + rank);
                end
            end
            
            all_ids = keys(rrf_scores);
            all_scores_vals = cell2mat(values(rrf_scores, all_ids));
            [sorted_scores, sort_idx] = sort(all_scores_vals, 'descend');
            sorted_ids = all_ids(sort_idx);
            
            num_keep = min(K, length(sorted_ids));
            stage1_bahn_results = table(sorted_ids(1:num_keep)', sorted_scores(1:num_keep)', (1:num_keep)', ...
                'VariableNames', {'bahn_id', 'rrf_score', 'rank'});
            
            stage1_bahn_time = toc(stage1_bahn_start);
            
            % ================================================================
            % BAHN-LEVEL: LOAD DATA
            % ================================================================
            loading_bahn_start = tic;
            
            ids_to_load = stage1_bahn_results.bahn_id;
            sequences = fetchBatchSequences(conn, schema, ids_to_load, 'bahn_id', dtw_mode);
            stage1_bahn_results.sequence = sequences;
            
            loading_bahn_time = toc(loading_bahn_start);
            
            % ================================================================
            % BAHN-LEVEL: STAGE 2
            % ================================================================
            config = struct();
            config.mode = dtw_mode;
            config.window = dtw_window;
            config.normalize = normalize_dtw;
            config.rot_align = use_rotation_alignment;
            config.lb_kim_ratio = lb_kim_keep_ratio;
            config.lb_keogh_n = lb_keogh_candidates;
            config.weights = strjoin(string(weights),',');
            
            stage2_bahn_start = tic;
            [stage2_bahn_results, bahn_stats] = performReranking(stage1_bahn_results, query_bahn_seq, config);
            stage2_bahn_time = toc(stage2_bahn_start);
            
            % ================================================================
            % SEGMENT-LEVEL: STAGE 1 + LOAD + STAGE 2
            % ================================================================
            stage1_seg_results = cell(num_segments, 1);
            stage2_seg_results = cell(num_segments, 1);
            seg_stats = cell(num_segments, 1);
            seg_times = struct('stage1', zeros(num_segments,1), 'loading', zeros(num_segments,1), 'stage2', zeros(num_segments,1));
            seg_lengths = zeros(num_segments, 1);
            
            for seg_idx = 1:num_segments
                seg_id = query_segment_ids{seg_idx};
                seg_lengths(seg_idx) = size(query_segment_seqs{seg_idx}, 1);
                
                % Stage 1
                seg_s1_start = tic;
                
                emb_sql = sprintf(['SELECT position_embedding::text, joint_embedding::text, ' ...
                                  'orientation_embedding::text, velocity_embedding::text, ' ...
                                  'metadata_embedding::text ' ...
                                  'FROM %s.bahn_embeddings WHERE segment_id = ''%s'''], schema, seg_id);
                emb_result = fetch(conn, emb_sql);
                
                seg_embeddings = {
                    parseEmbedding(emb_result.position_embedding), ...
                    parseEmbedding(emb_result.joint_embedding), ...
                    parseEmbedding(emb_result.orientation_embedding), ...
                    parseEmbedding(emb_result.velocity_embedding), ...
                    parseEmbedding(emb_result.metadata_embedding)
                };
                
                rrf_scores_seg = containers.Map();
                
                for m = 1:5
                    if weights(m) == 0 || isempty(seg_embeddings{m}), continue; end
                    
                    emb_str = sprintf('[%s]', strjoin(string(seg_embeddings{m}), ','));
                    col_name = sprintf('%s_embedding', modalities{m});
                    
                    search_sql = sprintf(['SELECT segment_id, %s <=> ''%s''::vector as distance ' ...
                     'FROM %s.bahn_embeddings ' ...
                     'WHERE segment_id != bahn_id ' ...
                     'AND bahn_id != ''%s'' ' ... 
                     'AND %s IS NOT NULL ' ...
                     'ORDER BY distance LIMIT %d'], ...
                     col_name, emb_str, schema, query_id, col_name, K*2);
                    
                    result = fetch(conn, search_sql);
                    
                    for rank = 1:height(result)
                        id = result.segment_id{rank};
                        if ~rrf_scores_seg.isKey(id), rrf_scores_seg(id) = 0; end
                        rrf_scores_seg(id) = rrf_scores_seg(id) + weights(m) / (rrf_k + rank);
                    end
                end
                
                if ~isempty(rrf_scores_seg)
                    all_ids_seg = keys(rrf_scores_seg);
                    all_scores_seg = cell2mat(values(rrf_scores_seg, all_ids_seg));
                    [sorted_scores_seg, sort_idx_seg] = sort(all_scores_seg, 'descend');
                    sorted_ids_seg = all_ids_seg(sort_idx_seg);
                    num_keep_seg = min(K, length(sorted_ids_seg));
                    
                    stage1_seg_results{seg_idx} = table(sorted_ids_seg(1:num_keep_seg)', ...
                        sorted_scores_seg(1:num_keep_seg)', (1:num_keep_seg)', ...
                        'VariableNames', {'segment_id', 'rrf_score', 'rank'});
                else
                    stage1_seg_results{seg_idx} = table();
                end
                
                seg_times.stage1(seg_idx) = toc(seg_s1_start);
                
                % Load
                seg_load_start = tic;
                if ~isempty(stage1_seg_results{seg_idx})
                    ids_to_load = stage1_seg_results{seg_idx}.segment_id;
                    sequences = fetchBatchSequences(conn, schema, ids_to_load, 'segment_id', dtw_mode);
                    stage1_seg_results{seg_idx}.sequence = sequences;
                end
                seg_times.loading(seg_idx) = toc(seg_load_start);
                
                % Stage 2
                seg_s2_start = tic;
                if ~isempty(stage1_seg_results{seg_idx})
                    [stage2_seg_results{seg_idx}, seg_stats{seg_idx}] = performReranking(...
                        stage1_seg_results{seg_idx}, query_segment_seqs{seg_idx}, config);
                else
                    stage2_seg_results{seg_idx} = table();
                    seg_stats{seg_idx} = struct('lb_kim_calls', 0, 'lb_keogh_calls', 0, 'dtw_calls', 0);
                end
                seg_times.stage2(seg_idx) = toc(seg_s2_start);
            end
            
            % ================================================================
            % PROGNOSE + COLLECT RESULTS
            % ================================================================
            
            % Ground Truth Bahn
            gt_sql = sprintf(['SELECT sidtw_average_distance FROM auswertung.info_sidtw ' ...
                              'WHERE segment_id = ''%s'''], query_id);
            gt_result = fetch(conn, gt_sql);
            gt_bahn = [];
            if ~isempty(gt_result)
                gt_bahn = gt_result.sidtw_average_distance(1);
            end
            
            % For each K value
            for k_idx = 1:length(all_k_values)
                prognose_top_n = all_k_values(k_idx);
                
                if prognose_top_n > height(stage2_bahn_results)
                    continue;
                end
                
                % === BAHN DIRECT ===
                [s1_simple, s1_weighted] = computePrognose(stage1_bahn_results, prognose_top_n, conn, 'rrf_score');
                [s2_simple, s2_weighted] = computePrognose(stage2_bahn_results, prognose_top_n, conn, 'dtw_dist');
                
                % === BAHN FROM SEGMENTS ===
                seg_prognose_s1_simple = zeros(num_segments, 1);
                seg_prognose_s1_weighted = zeros(num_segments, 1);
                seg_prognose_s2_simple = zeros(num_segments, 1);
                seg_prognose_s2_weighted = zeros(num_segments, 1);
                
                for seg_idx = 1:num_segments
                    if ~isempty(stage1_seg_results{seg_idx}) && height(stage1_seg_results{seg_idx}) >= prognose_top_n
                        [seg_prognose_s1_simple(seg_idx), seg_prognose_s1_weighted(seg_idx)] = ...
                            computePrognose(stage1_seg_results{seg_idx}, prognose_top_n, conn, 'rrf_score');
                    else
                        seg_prognose_s1_simple(seg_idx) = NaN;
                        seg_prognose_s1_weighted(seg_idx) = NaN;
                    end
                    
                    if ~isempty(stage2_seg_results{seg_idx}) && height(stage2_seg_results{seg_idx}) >= prognose_top_n
                        [seg_prognose_s2_simple(seg_idx), seg_prognose_s2_weighted(seg_idx)] = ...
                            computePrognose(stage2_seg_results{seg_idx}, prognose_top_n, conn, 'dtw_dist');
                    else
                        seg_prognose_s2_simple(seg_idx) = NaN;
                        seg_prognose_s2_weighted(seg_idx) = NaN;
                    end
                end
                
                % Aggregate from segments
                valid_s1 = ~isnan(seg_prognose_s1_simple);
                valid_s2 = ~isnan(seg_prognose_s2_simple);
                
                if any(valid_s1)
                    w = seg_lengths(valid_s1) / sum(seg_lengths(valid_s1));
                    agg_s1_simple = sum(w .* seg_prognose_s1_simple(valid_s1));
                    agg_s1_weighted = sum(w .* seg_prognose_s1_weighted(valid_s1));
                else
                    agg_s1_simple = NaN; agg_s1_weighted = NaN;
                end
                
                if any(valid_s2)
                    w = seg_lengths(valid_s2) / sum(seg_lengths(valid_s2));
                    agg_s2_simple = sum(w .* seg_prognose_s2_simple(valid_s2));
                    agg_s2_weighted = sum(w .* seg_prognose_s2_weighted(valid_s2));
                else
                    agg_s2_simple = NaN; agg_s2_weighted = NaN;
                end
                
                % === SAVE BAHN ROW ===
                row = struct();
                row.query_id = query_id;
                row.segment_id = query_id;  % Für Bahn = query_id
                row.level = 'bahn';
                row.K = prognose_top_n;
                row.length = size(query_bahn_position, 1);
                row.num_segments = num_segments;
                row.ground_truth = gt_bahn;
                
                % Direct prognose
                row.direct_s1_simple = s1_simple;
                row.direct_s1_weighted = s1_weighted;
                row.direct_s2_simple = s2_simple;
                row.direct_s2_weighted = s2_weighted;
                
                % From segments prognose
                row.fromsegs_s1_simple = agg_s1_simple;
                row.fromsegs_s1_weighted = agg_s1_weighted;
                row.fromsegs_s2_simple = agg_s2_simple;
                row.fromsegs_s2_weighted = agg_s2_weighted;
                
                % Errors
                if ~isempty(gt_bahn)
                    row.err_direct_s1_simple = abs(gt_bahn - s1_simple);
                    row.err_direct_s1_weighted = abs(gt_bahn - s1_weighted);
                    row.err_direct_s2_simple = abs(gt_bahn - s2_simple);
                    row.err_direct_s2_weighted = abs(gt_bahn - s2_weighted);
                    row.err_fromsegs_s1_simple = abs(gt_bahn - agg_s1_simple);
                    row.err_fromsegs_s1_weighted = abs(gt_bahn - agg_s1_weighted);
                    row.err_fromsegs_s2_simple = abs(gt_bahn - agg_s2_simple);
                    row.err_fromsegs_s2_weighted = abs(gt_bahn - agg_s2_weighted);
                else
                    row.err_direct_s1_simple = NaN;
                    row.err_direct_s1_weighted = NaN;
                    row.err_direct_s2_simple = NaN;
                    row.err_direct_s2_weighted = NaN;
                    row.err_fromsegs_s1_simple = NaN;
                    row.err_fromsegs_s1_weighted = NaN;
                    row.err_fromsegs_s2_simple = NaN;
                    row.err_fromsegs_s2_weighted = NaN;
                end
                
                % Times
                row.stage1_time_sec = stage1_bahn_time;
                row.loading_time_sec = loading_bahn_time;
                row.stage2_time_sec = stage2_bahn_time;
                row.total_time_sec = stage1_bahn_time + loading_bahn_time + stage2_bahn_time;
                
                % Stats
                row.lb_kim_calls = bahn_stats.lb_kim_calls;
                row.lb_keogh_calls = bahn_stats.lb_keogh_calls;
                row.dtw_calls = bahn_stats.dtw_calls;
                row.dtw_mode = dtw_mode;
                row.normalize_dtw = normalize_dtw;
                row.weights = strjoin(string(weights),',');
                
                all_results = [all_results; row];
            end
            
            % === SEGMENT ROWS ===
            for seg_idx = 1:num_segments
                seg_id = query_segment_ids{seg_idx};
                
                gt_seg_sql = sprintf(['SELECT sidtw_average_distance FROM auswertung.info_sidtw ' ...
                                      'WHERE segment_id = ''%s'''], seg_id);
                gt_seg_result = fetch(conn, gt_seg_sql);
                gt_seg = [];
                if ~isempty(gt_seg_result)
                    gt_seg = gt_seg_result.sidtw_average_distance(1);
                end
                
                for k_idx = 1:length(all_k_values)
                    prognose_top_n = all_k_values(k_idx);
                    
                    s1_tbl = stage1_seg_results{seg_idx};
                    s2_tbl = stage2_seg_results{seg_idx};
                    
                    if isempty(s2_tbl) || height(s2_tbl) < prognose_top_n
                        continue;
                    end
                    
                    [s1_simple, s1_weighted] = computePrognose(s1_tbl, prognose_top_n, conn, 'rrf_score');
                    [s2_simple, s2_weighted] = computePrognose(s2_tbl, prognose_top_n, conn, 'dtw_dist');
                    
                    row = struct();
                    row.query_id = query_id;
                    row.segment_id = seg_id;
                    row.level = 'segment';
                    row.K = prognose_top_n;
                    row.length = seg_lengths(seg_idx);
                    row.num_segments = 1;
                    row.ground_truth = gt_seg;
                    
                    row.direct_s1_simple = s1_simple;
                    row.direct_s1_weighted = s1_weighted;
                    row.direct_s2_simple = s2_simple;
                    row.direct_s2_weighted = s2_weighted;
                    
                    row.fromsegs_s1_simple = NaN;
                    row.fromsegs_s1_weighted = NaN;
                    row.fromsegs_s2_simple = NaN;
                    row.fromsegs_s2_weighted = NaN;
                    
                    if ~isempty(gt_seg)
                        row.err_direct_s1_simple = abs(gt_seg - s1_simple);
                        row.err_direct_s1_weighted = abs(gt_seg - s1_weighted);
                        row.err_direct_s2_simple = abs(gt_seg - s2_simple);
                        row.err_direct_s2_weighted = abs(gt_seg - s2_weighted);
                    else
                        row.err_direct_s1_simple = NaN;
                        row.err_direct_s1_weighted = NaN;
                        row.err_direct_s2_simple = NaN;
                        row.err_direct_s2_weighted = NaN;
                    end
                    
                    row.err_fromsegs_s1_simple = NaN;
                    row.err_fromsegs_s1_weighted = NaN;
                    row.err_fromsegs_s2_simple = NaN;
                    row.err_fromsegs_s2_weighted = NaN;
                    
                    row.stage1_time_sec = seg_times.stage1(seg_idx);
                    row.loading_time_sec = seg_times.loading(seg_idx);
                    row.stage2_time_sec = seg_times.stage2(seg_idx);
                    row.total_time_sec = seg_times.stage1(seg_idx) + seg_times.loading(seg_idx) + seg_times.stage2(seg_idx);
                    
                    if ~isempty(seg_stats{seg_idx})
                        row.lb_kim_calls = seg_stats{seg_idx}.lb_kim_calls;
                        row.lb_keogh_calls = seg_stats{seg_idx}.lb_keogh_calls;
                        row.dtw_calls = seg_stats{seg_idx}.dtw_calls;
                    else
                        row.lb_kim_calls = 0;
                        row.lb_keogh_calls = 0;
                        row.dtw_calls = 0;
                    end
    
                    row.dtw_mode = dtw_mode;
                    row.normalize_dtw = normalize_dtw;
                    row.weights = strjoin(string(weights),',');
                    
                    all_results = [all_results; row];
                end
            end
            
            
            fprintf('  ✓ Complete (Bahn: %.2fs, %d DTW calls)\n\n', ...
                stage1_bahn_time + loading_bahn_time + stage2_bahn_time, bahn_stats.dtw_calls);

        end

        close(conn);
            
        catch ME
            warning('Query %s failed: %s', query_id, ME.message);
            try close(conn); catch; end
        end
end

overall_time = toc(overall_start);

% ========================================================================
%% SECTION 4: SAVE RESULTS
% ========================================================================

fprintf('═══ SAVING RESULTS ═══\n\n');

if ~isempty(all_results)
    
    results_table = struct2table(all_results);
    writetable(results_table, results_filename);
    
    fprintf('✓ Saved: %s\n', results_filename);
    fprintf('  Total rows: %d\n', height(results_table));
    fprintf('  Bahn rows:  %d\n', sum(strcmp({all_results.level}, 'bahn')));
    fprintf('  Seg rows:   %d\n', sum(strcmp({all_results.level}, 'segment')));
else
    fprintf('⚠ No results to save\n');
end

fprintf('\nTotal time: %.1f min\n', overall_time / 60);
fprintf('✓ Done!\n');

% ========================================================================
%% HELPER FUNCTIONS
% ========================================================================

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

function sequences = fetchBatchSequences(conn, schema, ids, id_col_name, mode)
    sequences = cell(length(ids), 1);
    if isempty(ids), return; end
    
    id_list_str = sprintf('''%s'',', ids{:});
    id_list_str = id_list_str(1:end-1);
    
    if strcmp(mode, 'position')
        table_name = 'bahn_position_soll';
        cols = 'x_soll, y_soll, z_soll';
    else
        table_name = 'bahn_joint_states';
        cols = 'joint_1, joint_2, joint_3, joint_4, joint_5, joint_6';
    end
    
    sql = sprintf(['SELECT %s, %s FROM %s.%s ' ...
                   'WHERE %s IN (%s) ORDER BY timestamp'], ...
                   id_col_name, cols, schema, table_name, id_col_name, id_list_str);
    
    data = fetch(conn, sql);
    if isempty(data), return; end
    
    data_map = containers.Map();
    id_col_data = data.(id_col_name);
    
    if strcmp(mode, 'position')
        vals = [data.x_soll, data.y_soll, data.z_soll];
    else
        vals = [data.joint_1, data.joint_2, data.joint_3, ...
                data.joint_4, data.joint_5, data.joint_6];
    end
    
    [G, group_ids] = findgroups(id_col_data);
    
    for i = 1:length(group_ids)
        current_id = group_ids{i};
        current_seq = vals(G == i, :);
        data_map(current_id) = current_seq;
    end
    
    for i = 1:length(ids)
        req_id = ids{i};
        if isKey(data_map, req_id)
            sequences{i} = data_map(req_id);
        else
            sequences{i} = [];
        end
    end
end

function emb_vec = parseEmbedding(emb_text)
    if isempty(emb_text) || (iscell(emb_text) && isempty(emb_text{1}))
        emb_vec = [];
        return;
    end
    
    if iscell(emb_text)
        emb_text = emb_text{1};
    end
    
    emb_text = strrep(strrep(char(emb_text), '[', ''), ']', '');
    emb_vec = str2double(strsplit(emb_text, ','));
    emb_vec = emb_vec(~isnan(emb_vec));
end

function [simple_mean, weighted_mean] = computePrognose(results_table, top_n, conn, weight_col)
    n = min(top_n, height(results_table));
    if n == 0
        simple_mean = NaN;
        weighted_mean = NaN;
        return;
    end
    
    if ismember('bahn_id', results_table.Properties.VariableNames)
        ids = results_table.bahn_id(1:n);
    else
        ids = results_table.segment_id(1:n);
    end
    
    id_list = sprintf('''%s'',', ids{:});
    id_list = id_list(1:end-1);
    
    sql = sprintf(['SELECT segment_id, sidtw_average_distance FROM auswertung.info_sidtw ' ...
                   'WHERE segment_id IN (%s)'], id_list);
    data = fetch(conn, sql);
    
    if isempty(data)
        simple_mean = NaN;
        weighted_mean = NaN;
        return;
    end
    
    values = zeros(n, 1);
    for i = 1:n
        idx = find(strcmp(data.segment_id, ids{i}), 1);
        if ~isempty(idx)
            values(i) = data.sidtw_average_distance(idx);
        else
            values(i) = NaN;
        end
    end
    
    valid = ~isnan(values);
    values = values(valid);
    weights_vec = results_table.(weight_col)(1:n);
    weights_vec = weights_vec(valid);
    
    if isempty(values)
        simple_mean = NaN;
        weighted_mean = NaN;
        return;
    end
    
    simple_mean = mean(values);
    
    if strcmp(weight_col, 'dtw_dist')
        weights_vec = 1 ./ (weights_vec + 1e-6);
    end
    
    weights_vec = weights_vec / sum(weights_vec);
    weighted_mean = sum(values .* weights_vec);
end
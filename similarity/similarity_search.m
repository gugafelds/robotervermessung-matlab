%  TWO-STAGE RETRIEVAL - SIMPLIFIED & CLEAN
%  ========================================================================
%  Step-by-step trajectory similarity search
%  ========================================================================

clear; clc;

addpath(genpath(pwd));
addpath(genpath('../main'));
addpath(genpath('../lasertracker'));
addpath(genpath('../methods'));

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  TWO-STAGE RETRIEVAL - CLEAN VERSION                           ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

% ========================================================================
%% SECTION 1: SETUP & CONFIGURATION
% ========================================================================

fprintf('═══ SECTION 1: SETUP & CONFIGURATION ═══\n\n');

% === WAS SUCHEN WIR? ===
query_id = '1766066299';

% === WIE VIELE KANDIDATEN MIT EMBEDDINGS? ===
K = 250;  % Stage 1: Hole 50 Kandidaten mit Embeddings

% RRF Parameter (Paper Standard = 60, Python Default oft 50)
rrf_k = 60; 

% === EMBEDDING GEWICHTE ===
weights = [1.0; 1.0; 1.0; 1.0; 1.0];  % position, joint, orientation, velocity, metadata
weights = weights / sum(weights);     % Normalisieren

% === DTW EINSTELLUNGEN ===
dtw_mode = 'joint_states';  % 'position' oder 'joint_states'
dtw_window = 0.2;       % 20% window
normalize_dtw = true; %% besser für prognose bis jetzt
use_rotation_alignment = false;

% === PROGNOSE ===
prognose = true;
prognose_top_n = 10;  % Top-N für Prognose nutzen

% Lower Bounds (für Speedup)
lb_kim_keep_ratio = 1.0;      % 100% = alle durchlassen (kein LB_Kim)
lb_keogh_candidates = 500;     % Nach LB_Keogh: max 500 für DTW

% === WAS WOLLEN WIR AM ENDE ZEIGEN? ===
final_top_n = 5;  % Zeige nur Top-10

fprintf('Query ID:              %s\n', query_id);
fprintf('Stage 1 (Embeddings):  K = %d candidates\n', K);
fprintf('Stage 2 (DTW):         All K candidates, show top %d\n', final_top_n);
fprintf('\n');
fprintf('DTW Config:\n');
fprintf('  Mode:                %s\n', dtw_mode);
fprintf('  Window:              %.0f%%\n', dtw_window * 100);
fprintf('  Normalize:           %s\n', mat2str(normalize_dtw));
fprintf('  Rotation Alignment:  %s\n', mat2str(use_rotation_alignment));
fprintf('  LB_Kim keep ratio:   %.0f%%\n', lb_kim_keep_ratio * 100);
fprintf('  LB_Keogh target:     %d\n', lb_keogh_candidates);
fprintf('\n');

% === CONNECT TO DATABASE ===
fprintf('Connecting to database...\n');
conn = connectingToPostgres();
if ~isopen(conn)
    error('Database connection failed');
end
fprintf('✓ Connected\n\n');

fprintf('═══ SECTION 1 COMPLETE ═══\n\n');

% ========================================================================
%% SECTION 2: LOAD QUERY TRAJECTORY + SEGMENTS
% ========================================================================

fprintf('═══ SECTION 2: LOAD QUERY TRAJECTORY + SEGMENTS ═══\n\n');

schema = 'bewegungsdaten';

fprintf('Loading query: %s\n', query_id);

% === LADE BAHN (Full Trajectory) ===
pos_sql = sprintf(['SELECT x_soll, y_soll, z_soll FROM %s.bahn_position_soll ' ...
                  'WHERE bahn_id = ''%s'' ORDER BY timestamp'], schema, query_id);
pos = fetch(conn, pos_sql);

if isempty(pos)
    error('Query not found: %s', query_id);
end

query_bahn_position = [pos.x_soll, pos.y_soll, pos.z_soll];

joint_sql = sprintf(['SELECT joint_1, joint_2, joint_3, joint_4, joint_5, joint_6 ' ...
                    'FROM %s.bahn_joint_states WHERE bahn_id = ''%s'' ORDER BY timestamp'], ...
                    schema, query_id);
joint = fetch(conn, joint_sql);

query_bahn_joint = [joint.joint_1, joint.joint_2, joint.joint_3, ...
                   joint.joint_4, joint.joint_5, joint.joint_6];

% === LADE SEGMENTS ===
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
    
    % Position
    pos_seg_sql = sprintf(['SELECT x_soll, y_soll, z_soll FROM %s.bahn_position_soll ' ...
                          'WHERE segment_id = ''%s'' ORDER BY timestamp'], schema, seg_id);
    pos_seg = fetch(conn, pos_seg_sql);
    query_segments_position{i} = [pos_seg.x_soll, pos_seg.y_soll, pos_seg.z_soll];
    
    % Joint
    joint_seg_sql = sprintf(['SELECT joint_1, joint_2, joint_3, joint_4, joint_5, joint_6 ' ...
                            'FROM %s.bahn_joint_states WHERE segment_id = ''%s'' ORDER BY timestamp'], ...
                            schema, seg_id);
    joint_seg = fetch(conn, joint_seg_sql);
    query_segments_joint{i} = [joint_seg.joint_1, joint_seg.joint_2, joint_seg.joint_3, ...
                              joint_seg.joint_4, joint_seg.joint_5, joint_seg.joint_6];
end

% === WÄHLE SEQUENZEN FÜR DTW ===
if strcmp(dtw_mode, 'position')
    query_bahn_seq = query_bahn_position;
    query_segment_seqs = query_segments_position;
else
    query_bahn_seq = query_bahn_joint;
    query_segment_seqs = query_segments_joint;
end

fprintf('✓ Query loaded\n');
fprintf('  Bahn points:    %d\n', size(query_bahn_seq, 1));
fprintf('  Num segments:   %d\n', num_segments);
fprintf('  DTW mode:       %s\n', dtw_mode);
fprintf('\n');

fprintf('═══ SECTION 2 COMPLETE ═══\n\n');

% ========================================================================
%% SECTION 3: STAGE 1 - EMBEDDING SEARCH (RRF)
% ========================================================================
fprintf('═══ SECTION 3: STAGE 1 - EMBEDDING SEARCH ═══\n\n');

% ========================================================================
% BAHN-LEVEL
% ========================================================================
fprintf('--- BAHN-LEVEL ---\n');

% 1. Lade Query Embeddings
emb_sql = sprintf(['SELECT position_embedding::text, joint_embedding::text, ' ...
                  'orientation_embedding::text, velocity_embedding::text, ' ...
                  'metadata_embedding::text ' ...
                  'FROM %s.bahn_embeddings WHERE segment_id = ''%s'''], schema, query_id);
emb_result = fetch(conn, emb_sql);

% Parse Embeddings in Cell-Array für Loop
query_embeddings = {
    parseEmbedding(emb_result.position_embedding), ...
    parseEmbedding(emb_result.joint_embedding), ...
    parseEmbedding(emb_result.orientation_embedding), ...
    parseEmbedding(emb_result.velocity_embedding), ...
    parseEmbedding(emb_result.metadata_embedding)
};

modalities = {'position', 'joint', 'orientation', 'velocity', 'metadata'};

% 2. RRF Scores berechnen
fprintf('Searching with RRF (k=%d)...\n', rrf_k);
execute(conn, 'SET hnsw.ef_search = 200');
execute(conn, 'SET search_path = "bewegungsdaten"');

rrf_scores = containers.Map();

for m = 1:5
    if weights(m) == 0 || isempty(query_embeddings{m})
        continue;
    end
    
    % Embedding Vector -> String
    emb_str = sprintf('[%s]', strjoin(string(query_embeddings{m}), ','));
    col_name = sprintf('%s_embedding', modalities{m});
    
    % SQL: Vektor-Suche
    search_sql = sprintf(['SELECT segment_id, %s <=> ''%s''::vector as distance ' ...
                         'FROM %s.bahn_embeddings ' ...
                         'WHERE segment_id = bahn_id ' ...     % Nur ganze Bahnen
                         'AND segment_id != ''%s'' ' ...       % Exclude Query
                         'AND %s IS NOT NULL ' ...
                         'ORDER BY distance LIMIT %d'], ...
                         col_name, emb_str, schema, query_id, col_name, K);
    
    result = fetch(conn, search_sql);
    
    % RRF Akkumulation
    for rank = 1:height(result)
        id = result.segment_id{rank};
        
        if ~rrf_scores.isKey(id)
            rrf_scores(id) = 0;
        end
        
        % Formula: weight / (k + rank)
        rrf_scores(id) = rrf_scores(id) + weights(m) / (rrf_k + rank);
    end
end

% 3. Sortieren und Tabelle erstellen
all_ids = keys(rrf_scores);
all_scores = cell2mat(values(rrf_scores, all_ids));

[sorted_scores, sort_idx] = sort(all_scores, 'descend');
sorted_ids = all_ids(sort_idx);

% Cut to Top-K
num_keep = min(K, length(sorted_ids));
final_ids = sorted_ids(1:num_keep)';
final_scores = sorted_scores(1:num_keep)';
final_ranks = (1:num_keep)';

% ERGEBNIS TABELLE (BAHN)
stage1_bahn_results = table(final_ids, final_scores, final_ranks, ...
    'VariableNames', {'bahn_id', 'rrf_score', 'rank'});

fprintf('✓ Found %d bahn candidates\n', height(stage1_bahn_results));

% ========================================================================
% SEGMENT-LEVEL
% ========================================================================
fprintf('--- SEGMENT-LEVEL ---\n');

% Cell-Array für Tabellen initialisieren
stage1_seg_results = cell(num_segments, 1);

for seg_idx = 1:num_segments
    seg_id = query_segment_ids{seg_idx};
    fprintf('  Segment %d/%d: %s\n', seg_idx, num_segments, seg_id);
    
    % 1. Lade Segment Embeddings
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
    
    % 2. RRF Suche
    rrf_scores_seg = containers.Map();
    
    for m = 1:5
        if weights(m) == 0 || isempty(seg_embeddings{m})
            continue;
        end
        
        emb_str = sprintf('[%s]', strjoin(string(seg_embeddings{m}), ','));
        col_name = sprintf('%s_embedding', modalities{m});
        
        search_sql = sprintf(['SELECT segment_id, %s <=> ''%s''::vector as distance ' ...
                             'FROM %s.bahn_embeddings ' ...
                             'WHERE segment_id != bahn_id ' ...    % Nur Segmente
                             'AND segment_id != ''%s'' ' ...       % Exclude Query
                             'AND %s IS NOT NULL ' ...
                             'ORDER BY distance LIMIT %d'], ...
                             col_name, emb_str, schema, seg_id, col_name, K);
        
        result = fetch(conn, search_sql);
        
        for rank = 1:height(result)
            id = result.segment_id{rank};
            if ~rrf_scores_seg.isKey(id), rrf_scores_seg(id) = 0; end
            rrf_scores_seg(id) = rrf_scores_seg(id) + weights(m) / (rrf_k + rank);
        end
    end
    
    % 3. Tabelle erstellen
    if isempty(rrf_scores_seg)
        % Leere Tabelle falls keine Ergebnisse
        stage1_seg_results{seg_idx} = table({}, [], [], ...
            'VariableNames', {'segment_id', 'rrf_score', 'rank'});
        fprintf('    ⚠ No candidates found\n');
    else
        all_ids_seg = keys(rrf_scores_seg);
        all_scores_seg = cell2mat(values(rrf_scores_seg, all_ids_seg));
        
        [sorted_scores_seg, sort_idx_seg] = sort(all_scores_seg, 'descend');
        sorted_ids_seg = all_ids_seg(sort_idx_seg);
        
        num_keep_seg = min(K, length(sorted_ids_seg));
        
        % Daten für Tabelle
        t_ids = sorted_ids_seg(1:num_keep_seg)';
        t_scores = sorted_scores_seg(1:num_keep_seg)';
        t_ranks = (1:num_keep_seg)';
        
        % ERGEBNIS TABELLE (SEGMENT)
        stage1_seg_results{seg_idx} = table(t_ids, t_scores, t_ranks, ...
            'VariableNames', {'segment_id', 'rrf_score', 'rank'});
            
        fprintf('    ✓ %d candidates\n', height(stage1_seg_results{seg_idx}));
    end
end

fprintf('\n═══ STAGE 1 COMPLETE ═══\n\n');

clear m rank id emb_str col_name search_sql result
clear rrf_scores rrf_scores_seg all_ids all_scores sort_idx sorted_ids sorted_scores
clear all_ids_seg all_scores_seg sorted_ids_seg sorted_scores_seg sort_idx_seg
clear num_keep num_keep_seg final_ids final_scores final_ranks
clear t_ids t_scores t_ranks
clear emb_result emb_sql query_embeddings seg_embeddings modalities

% ========================================================================
%% SECTION 4: LOAD CANDIDATE DATA (BATCH)
% ========================================================================
fprintf('═══ SECTION 4: LOAD CANDIDATE DATA FOR DTW ═══\n\n');
fprintf('Loading mode: %s\n', dtw_mode);

load_data_start = tic;

% ========================================================================
% BAHN-LEVEL DATA LOADING
% ========================================================================
fprintf('--- BAHN-LEVEL (%d candidates) ---\n', height(stage1_bahn_results));

if ~isempty(stage1_bahn_results)
    % IDs aus der Tabelle holen
    ids_to_load = stage1_bahn_results.bahn_id;
    
    % Batch Load Helper aufrufen
    sequences = fetchBatchSequences(conn, schema, ids_to_load, 'bahn_id', dtw_mode);
    
    % Daten direkt in die Tabelle speichern
    stage1_bahn_results.sequence = sequences;
    
    fprintf('✓ Loaded %d sequences\n', length(sequences));
else
    fprintf('⚠ No candidates to load\n');
end
fprintf('\n');

% ========================================================================
% SEGMENT-LEVEL DATA LOADING
% ========================================================================
fprintf('--- SEGMENT-LEVEL ---\n');

for seg_idx = 1:num_segments
    current_table = stage1_seg_results{seg_idx};
    
    if ~isempty(current_table)
        % IDs holen
        ids_to_load = current_table.segment_id;
        
        % Batch Load
        sequences = fetchBatchSequences(conn, schema, ids_to_load, 'segment_id', dtw_mode);
        
        % In Tabelle speichern
        current_table.sequence = sequences;
        
        % Tabelle im Cell-Array updaten
        stage1_seg_results{seg_idx} = current_table;
        
        fprintf('  Segment %d/%d: Loaded %d sequences\n', ...
            seg_idx, num_segments, length(sequences));
    else
        fprintf('  Segment %d/%d: Skipped (no candidates)\n', seg_idx, num_segments);
    end
end

total_load_time = toc(load_data_start);
fprintf('\n✓ Data loading complete in %.3f sec\n', total_load_time);
fprintf('═══ SECTION 4 COMPLETE ═══\n\n');


% ========================================================================
%% SECTION 5: STAGE 2 - DTW RERANKING
% ========================================================================
fprintf('═══ SECTION 5: STAGE 2 - DTW RERANKING ═══\n\n');

dtw_start = tic;
total_dtw_calls = 0;

% DTW Configuration Struct für Helper
config = struct();
config.mode = dtw_mode;
config.window = dtw_window;
config.normalize = normalize_dtw;
config.rot_align = use_rotation_alignment;
config.lb_kim_ratio = lb_kim_keep_ratio;
config.lb_keogh_n = lb_keogh_candidates;

% ========================================================================
% BAHN-LEVEL RERANKING
% ========================================================================
fprintf('--- BAHN-LEVEL ---\n');

if ~isempty(stage1_bahn_results)
    % Helper aufrufen: Berechnet DTW, sortiert neu, fügt Ranks hinzu
    [stage2_bahn_results, stats] = performReranking(stage1_bahn_results, query_bahn_seq, config);
    
    total_dtw_calls = total_dtw_calls + stats.dtw_calls;
    
    % Zeige Top-3 Änderungen
    fprintf('✓ Reranked %d candidates. Top 3 results:\n', height(stage2_bahn_results));
else
    stage2_bahn_results = table();
    fprintf('⚠ No bahn candidates to rerank.\n');
end
fprintf('\n');


% ========================================================================
% SEGMENT-LEVEL RERANKING
% ========================================================================
fprintf('--- SEGMENT-LEVEL ---\n');

stage2_seg_results = cell(num_segments, 1);

for seg_idx = 1:num_segments
    % Tabelle und Query Sequenz holen
    current_table = stage1_seg_results{seg_idx};
    
    if strcmp(dtw_mode, 'position')
        q_seq = query_segments_position{seg_idx};
    else
        q_seq = query_segments_joint{seg_idx};
    end
    
    if ~isempty(current_table) && ~isempty(q_seq)
        % Reranking durchführen
        [reranked_table, stats] = performReranking(current_table, q_seq, config);
        
        stage2_seg_results{seg_idx} = reranked_table;
        total_dtw_calls = total_dtw_calls + stats.dtw_calls;
        
        % Kleine Statistik ausgeben
        best_change = max(reranked_table.rank_change);
        fprintf('  Segment %d: Top match ID: %s (Dist: %.2f) | Max promotion: +%d places\n', ...
            seg_idx, reranked_table.segment_id{1}, reranked_table.dtw_dist(1), best_change);
    else
        stage2_seg_results{seg_idx} = table();
        fprintf('  Segment %d: Skipped.\n', seg_idx);
    end
end

total_time = toc(dtw_start);
fprintf('\n✓ Stage 2 complete in %.3f sec (Total DTW calls: %d)\n', total_time, total_dtw_calls);
fprintf('═══ SECTION 5 COMPLETE ═══\n\n');

%% SECTION 6: VISUALIZATION
fprintf('═══ SECTION 6: VISUALIZATION ═══\n\n');
figure('Position', [100 100 1400 600]);

colors = lines(final_top_n);
lw_max = 3.0;
lw_min = 0.5;
line_widths = linspace(lw_max, lw_min, final_top_n);

if strcmp(dtw_mode, 'position')
    % === 3D PLOTS ===
    subplot(1,2,1); hold on;
    plot3(query_bahn_position(:,1), query_bahn_position(:,2), query_bahn_position(:,3), ...
          'k-', 'LineWidth', 4, 'DisplayName', 'Query');
    for i = 1:min(final_top_n, height(stage1_bahn_results))
        seq = stage1_bahn_results.sequence{i};
        if ~isempty(seq)
            plot3(seq(:,1), seq(:,2), seq(:,3), '-', 'Color', colors(i,:), ...
                  'LineWidth', line_widths(i), 'DisplayName', sprintf('#%d (RRF=%.4f)', i, stage1_bahn_results.rrf_score(i)));
        end
    end
    hold off; grid on; axis equal; view(3);
    title('Stage 1: Embedding Ranking'); legend('Location', 'bestoutside', 'FontSize', 7);
    
    subplot(1,2,2); hold on;
    plot3(query_bahn_position(:,1), query_bahn_position(:,2), query_bahn_position(:,3), ...
          'k-', 'LineWidth', 4, 'DisplayName', 'Query');
    for i = 1:min(final_top_n, height(stage2_bahn_results))
        seq = stage2_bahn_results.sequence{i};
        if ~isempty(seq)
            plot3(seq(:,1), seq(:,2), seq(:,3), '-', 'Color', colors(i,:), ...
                  'LineWidth', line_widths(i), 'DisplayName', sprintf('#%d (DTW=%.2f, Δ%+d)', i, stage2_bahn_results.dtw_dist(i), stage2_bahn_results.rank_change(i)));
        end
    end
    hold off; grid on; axis equal; view(3);
    title('Stage 2: DTW Reranked'); legend('Location', 'bestoutside', 'FontSize', 7);
    
else
    % === 2D JOINT PLOTS ===
    subplot(1,2,1); hold on;
    for j = 1:6
        plot(query_bahn_joint(:,j), 'k-', 'LineWidth', 2, 'HandleVisibility', 'off');
    end
    plot(NaN, NaN, 'k-', 'LineWidth', 2, 'DisplayName', 'Query');  % Dummy für Legende
    for i = 1:min(final_top_n, height(stage1_bahn_results))
        seq = stage1_bahn_results.sequence{i};
        if ~isempty(seq)
            for j = 1:6
                plot(seq(:,j), '-', 'Color', colors(i,:), 'LineWidth', line_widths(i), 'HandleVisibility', 'off');
            end
            plot(NaN, NaN, '-', 'Color', colors(i,:), 'LineWidth', line_widths(i), ...
                 'DisplayName', sprintf('#%d (RRF=%.4f)', i, stage1_bahn_results.rrf_score(i)));
        end
    end
    hold off; grid on; xlabel('Sample'); ylabel('Joint Angle [rad]');
    title('Stage 1: Embedding Ranking'); legend('Location', 'bestoutside', 'FontSize', 7);
    
    subplot(1,2,2); hold on;
    for j = 1:6
        plot(query_bahn_joint(:,j), 'k-', 'LineWidth', 2, 'HandleVisibility', 'off');
    end
    plot(NaN, NaN, 'k-', 'LineWidth', 2, 'DisplayName', 'Query');
    for i = 1:min(final_top_n, height(stage2_bahn_results))
        seq = stage2_bahn_results.sequence{i};
        if ~isempty(seq)
            for j = 1:6
                plot(seq(:,j), '-', 'Color', colors(i,:), 'LineWidth', line_widths(i), 'HandleVisibility', 'off');
            end
            plot(NaN, NaN, '-', 'Color', colors(i,:), 'LineWidth', line_widths(i), ...
                 'DisplayName', sprintf('#%d (DTW=%.2f, Δ%+d)', i, stage2_bahn_results.dtw_dist(i), stage2_bahn_results.rank_change(i)));
        end
    end
    hold off; grid on; xlabel('Sample'); ylabel('Joint Angle [rad]');
    title('Stage 2: DTW Reranked'); legend('Location', 'bestoutside', 'FontSize', 7);
end

sgtitle(sprintf('Query: %s | Mode: %s', query_id, dtw_mode));
fprintf('✓ Plots generated\n');

%% SECTION 7: PROGNOSE
if prognose
    fprintf('═══ SECTION 7: PROGNOSE ═══\n\n');
    
    % === BAHN-LEVEL ===
    fprintf('--- BAHN-LEVEL ---\n');
    
    [s1_simple, s1_weighted] = computePrognose(stage1_bahn_results, prognose_top_n, conn, 'rrf_score');
    fprintf('Stage 1 (Embedding):  Simple=%.4f | Weighted=%.4f\n', s1_simple, s1_weighted);
    
    [s2_simple, s2_weighted] = computePrognose(stage2_bahn_results, prognose_top_n, conn, 'dtw_dist');
    fprintf('Stage 2 (DTW):        Simple=%.4f | Weighted=%.4f\n', s2_simple, s2_weighted);
    
    gt_sql = sprintf(['SELECT sidtw_average_distance FROM auswertung.info_sidtw ' ...
                      'WHERE segment_id = ''%s'''], query_id);
    gt_result = fetch(conn, gt_sql);
    if ~isempty(gt_result)
        gt_value = gt_result.sidtw_average_distance(1);
        fprintf('\nGround Truth:         %.4f\n', gt_value);
        fprintf('Error Stage 1:        Simple=%.4f | Weighted=%.4f\n', ...
                abs(gt_value - s1_simple), abs(gt_value - s1_weighted));
        fprintf('Error Stage 2:        Simple=%.4f | Weighted=%.4f\n', ...
                abs(gt_value - s2_simple), abs(gt_value - s2_weighted));
    end
    
    % === SEGMENT-LEVEL ===
    fprintf('\n--- SEGMENT-LEVEL ---\n');
    
    for seg_idx = 1:num_segments
        seg_id = query_segment_ids{seg_idx};
        fprintf('\nSegment %d: %s\n', seg_idx, seg_id);
        
        s1_table = stage1_seg_results{seg_idx};
        s2_table = stage2_seg_results{seg_idx};
        
        if ~isempty(s1_table)
            [s1_simple, s1_weighted] = computePrognose(s1_table, prognose_top_n, conn, 'rrf_score');
            fprintf('  Stage 1:  Simple=%.4f | Weighted=%.4f\n', s1_simple, s1_weighted);
        end
        
        if ~isempty(s2_table)
            [s2_simple, s2_weighted] = computePrognose(s2_table, prognose_top_n, conn, 'dtw_dist');
            fprintf('  Stage 2:  Simple=%.4f | Weighted=%.4f\n', s2_simple, s2_weighted);
        end
        
        % Ground Truth für Segment
        gt_sql = sprintf(['SELECT sidtw_average_distance FROM auswertung.info_sidtw ' ...
                          'WHERE segment_id = ''%s'''], seg_id);
        gt_result = fetch(conn, gt_sql);
        if ~isempty(gt_result)
            gt_val = gt_result.sidtw_average_distance(1);
            fprintf('  GT=%.4f | Err S1: %.4f / %.4f | Err S2: %.4f / %.4f\n', ...
                    gt_val, abs(gt_val - s1_simple), abs(gt_val - s1_weighted), ...
                    abs(gt_val - s2_simple), abs(gt_val - s2_weighted));
        end
    end
    
    fprintf('\n═══ SECTION 7 COMPLETE ═══\n\n');
end



% ========================================================================
% HELPER FUNCTION: RERANKING LOGIC
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
    % Lädt viele Sequenzen auf einmal (WHERE id IN (...))
    % Garantiert die gleiche Reihenfolge wie 'ids' Input!
    
    sequences = cell(length(ids), 1);
    if isempty(ids), return; end
    
    % 1. SQL Query bauen
    % IDs für SQL formatieren: 'id1','id2','id3'
    id_list_str = sprintf('''%s'',', ids{:});
    id_list_str = id_list_str(1:end-1); % Letztes Komma weg
    
    if strcmp(mode, 'position')
        table_name = 'bahn_position_soll';
        cols = 'x_soll, y_soll, z_soll';
    else
        table_name = 'bahn_joint_states';
        cols = 'joint_1, joint_2, joint_3, joint_4, joint_5, joint_6';
    end
    
    % Wir holen auch die ID, um die Daten zuzuordnen (SQL sortiert zufällig!)
    sql = sprintf(['SELECT %s, %s FROM %s.%s ' ...
                   'WHERE %s IN (%s) ORDER BY timestamp'], ...
                   id_col_name, cols, schema, table_name, id_col_name, id_list_str);
                   
    % 2. Abfrage ausführen
    data = fetch(conn, sql);
    
    
    if isempty(data), return; end
    
    % 3. Daten in Map speichern für schnelle Zuordnung
    % Key: ID, Value: Matrix
    data_map = containers.Map();
    
    % MATLAB fetch gibt Columns zurück, wir müssen gruppieren.
    % Da 'fetch' alle Zeilen aller IDs flach zurückgibt, müssen wir splitten.
    % Trick: 'findgroups' auf die ID Spalte
    
    % Hole die ID Spalte als Cell Array of Strings oder Categorical
    id_col_data = data.(id_col_name);
    
    % Konvertiere Daten zu Matrix
    if strcmp(mode, 'position')
        vals = [data.x_soll, data.y_soll, data.z_soll];
    else
        vals = [data.joint_1, data.joint_2, data.joint_3, ...
                data.joint_4, data.joint_5, data.joint_6];
    end
    
    % Gruppieren nach ID
    [G, group_ids] = findgroups(id_col_data);
    
    for i = 1:length(group_ids)
        current_id = group_ids{i};
        % Extrahiere Zeilen für diese ID
        current_seq = vals(G == i, :);
        data_map(current_id) = current_seq;
    end
    
    % 4. In der ursprünglichen Reihenfolge zurückgeben
    for i = 1:length(ids)
        req_id = ids{i};
        if isKey(data_map, req_id)
            sequences{i} = data_map(req_id);
        else
            sequences{i} = []; % Sollte nicht passieren
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
    % Berechnet Simple und Weighted Mean für sidtw_average_distance
    
    n = min(top_n, height(results_table));
    if n == 0
        simple_mean = NaN;
        weighted_mean = NaN;
        return;
    end
    
    % IDs holen (bahn_id oder segment_id)
    if ismember('bahn_id', results_table.Properties.VariableNames)
        ids = results_table.bahn_id(1:n);
    else
        ids = results_table.segment_id(1:n);
    end
    
    % Werte aus DB holen
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
    
    % Werte in richtige Reihenfolge bringen
    values = zeros(n, 1);
    for i = 1:n
        idx = find(strcmp(data.segment_id, ids{i}), 1);
        if ~isempty(idx)
            values(i) = data.sidtw_average_distance(idx);
        else
            values(i) = NaN;
        end
    end
    
    % NaN entfernen
    valid = ~isnan(values);
    values = values(valid);
    weights = results_table.(weight_col)(1:n);
    weights = weights(valid);
    
    if isempty(values)
        simple_mean = NaN;
        weighted_mean = NaN;
        return;
    end
    
    % Simple Mean
    simple_mean = mean(values);
    
    % Weighted Mean
    if strcmp(weight_col, 'dtw_dist')
        % DTW: kleiner = besser -> invertieren
        weights = 1 ./ (weights + 1e-6);
    end
    % RRF: größer = besser -> direkt nutzen
    
    weights = weights / sum(weights);  % Normalisieren
    weighted_mean = sum(values .* weights);
end
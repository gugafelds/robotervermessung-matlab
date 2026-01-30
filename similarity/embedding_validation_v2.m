%% SECTION 1: CONFIG
clear; clc;

addpath(genpath(pwd));
addpath(genpath('../main'));
addpath(genpath('../lasertracker'));
addpath(genpath('../methods'));

% === Experiment Settings ===
database_sample_size = 7500;
random_seed = 21;
top_k = 50;

% === Query IDs ===
query_ids = {
    % linear
    %% clean
    %'1765989370'; % clean, np = 3 / 10 GT
    %'1765989294'; % clean, np = 3 / 20 GT
    %'1765988821'; % clean, np = 3 / 30 GT
    %'1765988920'; % clean; np = 3 / 40 GT
    %'1765989411'; % clean; np = 3 / 50 GT
    
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

    % zirkular
    %% clean
    '1769770498'; % clean, np = 3 / 10 GT
    '1769770684'; % clean, np = 3 / 20 GT
    '1769770935'; % clean, np = 3 / 30 GT
    '1769771107'; % clean; np = 3 / 40 GT
    '1769771447'; % clean; np = 3 / 50 GT
    
    %% noisy - 2 mm

    %'1769773928'; % noisy; np = 3 / 10 GT
    %'1769772060'; % noisy; np = 3 / 20 GT
    %'1769772213'; % noisy; np = 3 / 30 GT
    %'1769773985'; % noisy; np = 3 / 40 GT
    %'1769774278'; % noisy; np = 3 / 50 GT
    
    %% noisy - 5 mm

    %'1769772609'; % noisy; np = 3 / 10 GT
    %'1769773593'; % noisy; np = 3 / 20 GT
    %'1769772776'; % noisy; np = 3 / 30 GT
    %'1769772900'; % noisy; np = 3 / 39 (40) GT 
    %'1769773333'; % noisy; np = 3 / 50 GT

    %% noise - 10 mm

    %'1769774581'; %noise; np = 3 / 10 GT
    };

% === Embedding Configs: {name, n_coarse, n_fine, multi_scale} ===
embedding_configs = {
    'Single-2',    0,   2, false;
    'Single-5',    0,   5, false;
    'Single-10',   0,  10, false;
    'Single-20',   0,  20, false;
    'Multi-15',    5,  10, true;
    'Multi-25',    5,  20, true;
    'Single-50',   0,  50, false;
    'Single-100',  0, 100, false;
};

% === Weight Modes: {name, dtw_mode, [pos,joint,orient,vel,meta]} ===
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

% === DTW Settings ===
dtw_window = 0.2;
dtw_normalize = false;
dtw_rot_align = false;
lb_kim_ratio = 1.0;
lb_keogh_n = top_k * 10;

% === Counts ===
num_queries = length(query_ids);
num_embeddings = size(embedding_configs, 1);
num_weight_modes = size(weight_mode_configs, 1);

fprintf('Queries: %d | Embeddings: %d | Weight Modes: %d\n', ...
    num_queries, num_embeddings, num_weight_modes);

%% SECTION 2: DATABASE + SAMPLING + GT

conn = connectingToPostgres();
schema = 'bewegungsdaten';

% === Sample Candidates ===
full_db_query = sprintf(...
    'SELECT bahn_id FROM robotervermessung.%s.bahn_metadata WHERE bahn_id = segment_id', schema);
full_db = fetch(conn, full_db_query);

rng(random_seed);
sample_idx = randperm(height(full_db), database_sample_size);
candidate_ids = full_db.bahn_id(sample_idx);

fprintf('Sampled: %d trajectories\n', length(candidate_ids));

% === Ground Truth ===
[ground_truth_ids, ground_truth_map] = getGTCandidates(conn, schema, query_ids);

if ~isempty(ground_truth_ids)
    candidate_ids = setdiff(candidate_ids, ground_truth_ids);
    candidate_ids = [candidate_ids; ground_truth_ids];
    fprintf('GT added: %d trajectories\n', length(ground_truth_ids));
else
    ground_truth_map = struct();
    fprintf('No GT found\n');
end

fprintf('Total candidates: %d\n', length(candidate_ids));

%% SECTION 3: PRELOAD DATA

fprintf('Loading data...\n');
tic;

chunk_size = 100;
data_cache = loadDataExperiment(conn, schema, candidate_ids, query_ids, chunk_size);

fprintf('Data loaded: %.1f min\n', toc/60);

%% SECTION 4: PRECOMPUTE DTW

fprintf('Precomputing DTW...\n');
tic;

dtw_cfg = struct('window', dtw_window, 'normalize', dtw_normalize, ...
                 'rot_align', dtw_rot_align, 'lb_kim_ratio', lb_kim_ratio, ...
                 'lb_keogh_n', lb_keogh_n);

dtw_cache = precomputeDTW_v2(data_cache, query_ids, dtw_cfg);

fprintf('DTW precomputed: %.1f min\n', toc/60);

%% SECTION 5: PRECOMPUTE EMBEDDINGS

fprintf('Precomputing embeddings...\n');
tic;

emb_config = struct();
emb_config.norm_strategy = 'max_extent';

emb_cache = precomputeEmbeddings_v2(data_cache, query_ids, embedding_configs);

fprintf('Embeddings precomputed: %.1f s\n', toc);

% === DB nicht mehr nötig ===
close(conn);
fprintf('DB connection closed\n');

%% SECTION 6: RUN EXPERIMENTS

fprintf('Running experiments...\n');
tic;

all_results = [];

timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd''T''HHmmss'));

for emb_idx = 1:num_embeddings
    emb_name = embedding_configs{emb_idx, 1};
    ef = matlab.lang.makeValidName(emb_name);
    
    for wm_idx = 1:num_weight_modes
        wm_name = weight_mode_configs{wm_idx, 1};
        dtw_mode = weight_mode_configs{wm_idx, 2};
        weights = weight_mode_configs{wm_idx, 3};
        weights = weights(:) / sum(weights);
        
        for q_idx = 1:num_queries
            query_id = query_ids{q_idx};
            qf = sprintf('q_%s', strrep(query_id, '-', '_'));
            
            % === Get caches ===
            query_emb = emb_cache.(qf).(ef);
            query_dtw = dtw_cache.(qf).(dtw_mode);
            
            % === GT ===
            if isfield(ground_truth_map, qf)
                gt_traj_ids = ground_truth_map.(qf).trajectories;
                gt_seg_struct = ground_truth_map.(qf).segments;
            else
                gt_traj_ids = {};
                gt_seg_struct = struct();
            end
            
            % ============================================================
            % BAHN-LEVEL
            % ============================================================
            
            % Embedding ranking (RRF)
            bahn_emb_ranking = fuseRRF(query_emb.query, query_emb.candidates, weights, top_k);
            
            % DTW ranking
            bahn_dtw_ranking = query_dtw.trajectory_ranking;
            
            % Metriken
            [dtw_eb, gt_eb, gt_dtw] = computeAllMetrics(...
                bahn_emb_ranking, bahn_dtw_ranking, gt_traj_ids, 'bahn_id');

            % Im Loop, vor createRow:
            exp_config = struct();
            exp_config.n_coarse = embedding_configs{emb_idx, 2};
            exp_config.n_fine = embedding_configs{emb_idx, 3};
            exp_config.top_k = top_k;
            exp_config.database_size = database_sample_size;
            exp_config.num_trajectories = length(candidate_ids);
            exp_config.num_segments = height(data_cache.segments.metadata);
            exp_config.lb_kim_ratio = lb_kim_ratio;
            exp_config.lb_keogh_n = lb_keogh_n;
            exp_config.dtw_normalize = dtw_normalize;
            exp_config.dtw_rot_align = dtw_rot_align;
            exp_config.timestamp = timestamp;
            
            % Zeile speichern
            row = createRow('bahn', query_id, query_id, emb_name, wm_name, dtw_mode, ...
                dtw_eb, gt_eb, gt_dtw, length(gt_traj_ids), exp_config);
            all_results = [all_results; row];
            
            % ============================================================
            % SEGMENT-LEVEL
            % ============================================================
            
            num_segs = length(query_emb.segments);
            seg_fields = fieldnames(gt_seg_struct);
            
            for s = 1:num_segs
                seg_id = sprintf('%s_%d', query_id, s);
                
                % Embedding ranking
                seg_emb_data = query_emb.segments{s};
                if isempty(seg_emb_data.candidates.segment_ids), continue; end
                
                seg_emb_ranking = fuseRRF(seg_emb_data.query, seg_emb_data.candidates, weights, top_k);
                
                % DTW ranking
                if s <= length(query_dtw.segment_rankings)
                    seg_dtw_ranking = query_dtw.segment_rankings{s};
                else
                    continue;
                end
                
                % GT für Segment
                if s <= length(seg_fields)
                    gt_seg_ids = gt_seg_struct.(seg_fields{s});
                    if ~iscell(gt_seg_ids), gt_seg_ids = {gt_seg_ids}; end
                else
                    gt_seg_ids = {};
                end
                
                % Metriken
                [seg_dtw_eb, seg_gt_eb, seg_gt_dtw] = computeAllMetrics(...
                    seg_emb_ranking, seg_dtw_ranking, gt_seg_ids, 'segment_id');
                
                % Zeile speichern
                row = createRow('segment', query_id, seg_id, emb_name, wm_name, dtw_mode, ...
                seg_dtw_eb, seg_gt_eb, seg_gt_dtw, length(gt_seg_ids), exp_config);
                all_results = [all_results; row];
            end
        end
    end
    
    fprintf('  %s done\n', emb_name);
end

fprintf('Experiments done: %.1f min\n', toc/60);

%% SECTION 7: SAVE RESULTS

fprintf('Saving results...\n');

% Results folder
if ~exist('results', 'dir')
    mkdir('results');
end

% Convert to table
results_table = struct2table(all_results);

% Filename
filename = sprintf('results/embedding_validation_%s.csv', timestamp);

% Save
writetable(results_table, filename);

fprintf('Saved: %s\n', filename);
fprintf('  Rows: %d\n', height(results_table));
fprintf('  Bahn: %d\n', sum(strcmp({all_results.level}, 'bahn')));
fprintf('  Segment: %d\n', sum(strcmp({all_results.level}, 'segment')));
fprintf('\nDone!\n');

%% === HELPER FUNCTIONS ===

function ranking = fuseRRF(query_emb, cand_emb, weights, top_k)
    % Erstellt Rankings pro Modalität, dann nutzt fuseRankingsRRF
    
    modalities = {'position', 'joint', 'orientation', 'velocity', 'metadata'};
    
    % ID-Feld bestimmen
    if isfield(cand_emb, 'bahn_ids')
        id_field = 'bahn_id';
        ids = cand_emb.bahn_ids;
    else
        id_field = 'segment_id';
        ids = cand_emb.segment_ids;
    end
    
    n = length(ids);
    rankings = struct();
    
    for m = 1:5
        mod_name = modalities{m};
        
        q_vec = query_emb.(mod_name);
        c_mat = cand_emb.(mod_name);
        
        if isempty(q_vec) || isempty(c_mat) || weights(m) == 0
            rankings.(mod_name) = [];
            continue;
        end
        
        % Cosine distance
        dists = 1 - (c_mat * q_vec') ./ (vecnorm(c_mat,2,2) * norm(q_vec) + 1e-10);
        
        % Sortieren
        [sorted_dists, sort_idx] = sort(dists, 'ascend');
        sorted_ids = ids(sort_idx);
        
        % Table erstellen
        rankings.(mod_name) = table(sorted_ids, sorted_dists, ...
            'VariableNames', {id_field, 'distance'});
    end
    
    % Deine Funktion nutzen
    rrf_k = 60;
    ranking = fuseRankingsRRF(rankings, weights, rrf_k, id_field);
    
    % Top-K limitieren
    if height(ranking) > top_k
        ranking = ranking(1:top_k, :);
    end
end

function [dtw_eb, gt_eb, gt_dtw] = computeAllMetrics(emb_ranking, dtw_ranking, gt_ids, id_field)
    % DTW vs Embedding
    dtw_eb = computeDTWvsEB(emb_ranking, dtw_ranking, id_field);
    
    % GT vs Embedding
    gt_eb = computeGTMetrics(emb_ranking, gt_ids, dtw_ranking, id_field);
    
    % GT vs DTW
    gt_dtw = computeGTMetrics(dtw_ranking, gt_ids, dtw_ranking, id_field);
end

function m = computeDTWvsEB(emb_ranking, dtw_ranking, id_field)
    % Join auf gemeinsame IDs
    [common_ids, emb_idx, dtw_idx] = intersect(emb_ranking.(id_field), dtw_ranking.(id_field), 'stable');
    
    m = struct();
    n_common = length(common_ids);
    
    if n_common < 2
        m.spearman = NaN; m.ndcg_10 = NaN; m.ndcg_50 = NaN;
        m.r1 = NaN; m.r5 = NaN; m.r10 = NaN; m.r50 = NaN;
        m.p_at_50 = NaN;
        return;
    end
    
    % Ranks
    emb_ranks = (1:height(emb_ranking))';
    dtw_ranks = (1:height(dtw_ranking))';
    
    emb_ranks_common = emb_ranks(emb_idx);
    dtw_ranks_common = dtw_ranks(dtw_idx);
    
    % Spearman
    m.spearman = corr(emb_ranks_common, dtw_ranks_common, 'Type', 'Spearman');
    
    % Recall@K: Wie viele der DTW Top-K sind in Embedding Top-K?
    for k = [1, 5, 10, 50]
        dtw_top_k = dtw_ranking.(id_field)(1:min(k, height(dtw_ranking)));
        emb_top_k = emb_ranking.(id_field)(1:min(k, height(emb_ranking)));
        overlap = length(intersect(dtw_top_k, emb_top_k));
        
        switch k
            case 1,  m.r1 = overlap / min(k, length(dtw_top_k));
            case 5,  m.r5 = overlap / min(k, length(dtw_top_k));
            case 10, m.r10 = overlap / min(k, length(dtw_top_k));
            case 50, m.r50 = overlap / min(k, length(dtw_top_k));
        end
    end
    
    % P@K (Coverage Point) - auf gemeinsamer Schnittmenge
    m.p_at_50 = NaN;
    dtw_top_50 = dtw_ranking.(id_field)(1:min(50, height(dtw_ranking)));
    dtw_top_50_common = intersect(dtw_top_50, common_ids, 'stable');
    
    if ~isempty(dtw_top_50_common)
        for k = 1:height(emb_ranking)
            emb_top_k = emb_ranking.(id_field)(1:k);
            if all(ismember(dtw_top_50_common, emb_top_k))
                m.p_at_50 = k;
                break;
            end
        end
    end
    
    % NDCG: DTW-Distanz als Relevanz
    dtw_dists = dtw_ranking.dtw_distance(dtw_idx);
    rel = 1 ./ (1 + dtw_dists);
    
    m.ndcg_10 = computeNDCG(rel, emb_idx, 10);
    m.ndcg_50 = computeNDCG(rel, emb_idx, 50);
end

function m = computeGTMetrics(ranking, gt_ids, dtw_ranking, id_field)
    m = struct();
    m.ndcg_10 = NaN; m.ndcg_50 = NaN; m.mrr = NaN;
    m.r1 = NaN; m.r5 = NaN; m.r10 = NaN; m.r50 = NaN;
    m.mean_rank = NaN;
    
    if isempty(gt_ids), return; end
    
    ranking_ids = ranking.(id_field);
    n_gt = length(gt_ids);
    
    % Finde GT Ränge
    gt_ranks = zeros(n_gt, 1);
    for i = 1:n_gt
        idx = find(strcmp(ranking_ids, gt_ids{i}), 1);
        if ~isempty(idx)
            gt_ranks(i) = idx;
        else
            gt_ranks(i) = inf;
        end
    end
    
    valid_ranks = gt_ranks(gt_ranks < inf);
    if isempty(valid_ranks), return; end
    
    % Recall@K
    m.r1 = sum(valid_ranks <= 1) / n_gt;
    m.r5 = sum(valid_ranks <= 5) / n_gt;
    m.r10 = sum(valid_ranks <= 10) / n_gt;
    m.r50 = sum(valid_ranks <= 50) / n_gt;
    
    % MRR
    m.mrr = mean(1 ./ valid_ranks);
    
    % Mean Rank
    m.mean_rank = mean(valid_ranks);
    
    % NDCG mit DTW-Similarity als Relevanz
    gt_sims = zeros(n_gt, 1);
    for i = 1:n_gt
        dtw_idx = find(strcmp(dtw_ranking.(id_field), gt_ids{i}), 1);
        if ~isempty(dtw_idx)
            gt_sims(i) = 1 / (1 + dtw_ranking.dtw_distance(dtw_idx));
        end
    end
    
    m.ndcg_10 = computeNDCGfromGT(gt_ranks, gt_sims, 10);
    m.ndcg_50 = computeNDCGfromGT(gt_ranks, gt_sims, 50);
end

function ndcg = computeNDCG(rel, idx, k)
    n = min(k, length(idx));
    if n == 0, ndcg = NaN; return; end
    
    % DCG
    dcg = sum(rel(1:n) ./ log2((1:n)' + 1));
    
    % IDCG
    ideal_rel = sort(rel, 'descend');
    idcg = sum(ideal_rel(1:n) ./ log2((1:n)' + 1));
    
    if idcg > 0
        ndcg = dcg / idcg;
    else
        ndcg = 0;
    end
end

function ndcg = computeNDCGfromGT(gt_ranks, gt_sims, k)
    valid = gt_ranks <= k & gt_ranks < inf;
    if ~any(valid), ndcg = 0; return; end
    
    % DCG
    dcg = sum(gt_sims(valid) ./ log2(gt_ranks(valid) + 1));
    
    % IDCG (perfect ranking)
    sorted_sims = sort(gt_sims, 'descend');
    n = min(k, length(sorted_sims));
    idcg = sum(sorted_sims(1:n) ./ log2((1:n)' + 1));
    
    if idcg > 0
        ndcg = dcg / idcg;
    else
        ndcg = 0;
    end
end

function row = createRow(level, query_id, segment_id, emb_name, wm_name, dtw_mode, ...
                         dtw_eb, gt_eb, gt_dtw, num_gt, config)
    row = struct();
    row.level = level;
    row.query_id = query_id;
    row.segment_id = segment_id;
    row.emb_config = emb_name;
    row.weight_mode = wm_name;
    row.dtw_mode = dtw_mode;
    
    % DTW vs EB
    row.spearman_dtw_eb = dtw_eb.spearman;
    row.ndcg_10_dtw_eb = dtw_eb.ndcg_10;
    row.ndcg_50_dtw_eb = dtw_eb.ndcg_50;
    row.r1_dtw_eb = dtw_eb.r1;
    row.r5_dtw_eb = dtw_eb.r5;
    row.r10_dtw_eb = dtw_eb.r10;
    row.r50_dtw_eb = dtw_eb.r50;
    row.p_at_50_dtw_eb = dtw_eb.p_at_50;
    
    % GT vs EB
    row.ndcg_10_gt_eb = gt_eb.ndcg_10;
    row.ndcg_50_gt_eb = gt_eb.ndcg_50;
    row.mrr_gt_eb = gt_eb.mrr;
    row.r1_gt_eb = gt_eb.r1;
    row.r5_gt_eb = gt_eb.r5;
    row.r10_gt_eb = gt_eb.r10;
    row.r50_gt_eb = gt_eb.r50;
    row.mean_rank_gt_eb = gt_eb.mean_rank;
    
    % GT vs DTW
    row.ndcg_10_gt_dtw = gt_dtw.ndcg_10;
    row.ndcg_50_gt_dtw = gt_dtw.ndcg_50;
    row.mrr_gt_dtw = gt_dtw.mrr;
    row.r1_gt_dtw = gt_dtw.r1;
    row.r5_gt_dtw = gt_dtw.r5;
    row.r10_gt_dtw = gt_dtw.r10;
    row.r50_gt_dtw = gt_dtw.r50;
    row.mean_rank_gt_dtw = gt_dtw.mean_rank;
    
    % Meta
    row.num_gt = num_gt;

    % === Config Info ===
    row.n_coarse = config.n_coarse;
    row.n_fine = config.n_fine;
    row.top_k = config.top_k;
    row.database_size = config.database_size;
    row.num_trajectories = config.num_trajectories;
    row.num_segments = config.num_segments;
    row.lb_kim_ratio = config.lb_kim_ratio;
    row.lb_keogh_n = config.lb_keogh_n;
    row.dtw_normalize = config.dtw_normalize;
    row.dtw_rot_align = config.dtw_rot_align;
    row.timestamp = config.timestamp;
end

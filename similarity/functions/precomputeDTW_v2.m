function dtw_cache = precomputeDTW_v2(data_cache, query_ids, config)
% Berechnet DTW für alle Queries × Modes (position, joint_states)

fprintf('Computing DTW...\n');
tic;

dtw_modes = {'position', 'joint_states'};
dtw_cache = struct();

for q_idx = 1:length(query_ids)
    query_id = query_ids{q_idx};
    qf = sprintf('q_%s', strrep(query_id, '-', '_'));
    
    if ~isfield(data_cache.queries, qf), continue; end
    query = data_cache.queries.(qf);
    
    dtw_cache.(qf) = struct();
    
    for m = 1:length(dtw_modes)
        mode = dtw_modes{m};
        
        % === Trajectory-Level ===
        if strcmp(mode, 'position')
            query_seq = query.position;
            cand_seqs = data_cache.candidates.position;
        else
            query_seq = query.joint;
            cand_seqs = data_cache.candidates.joint;
        end
        
        traj_ranking = computeDTWRanking(query_seq, cand_seqs, ...
            data_cache.candidates.bahn_ids, 'bahn_id', mode, config);
        
        % === Segment-Level ===
        num_segs = query.metadata.num_segments;
        seg_rankings = cell(num_segs, 1);
        
        for s = 1:num_segs
            own_seg_id = sprintf('%s_%d', query_id, s);
            mask = ~strcmp(data_cache.segments.segment_ids, own_seg_id);
            
            if strcmp(mode, 'position')
                q_seg = query.segments.position{s};
                c_segs = data_cache.segments.position(mask);
            else
                q_seg = query.segments.joint{s};
                c_segs = data_cache.segments.joint(mask);
            end
            
            seg_ids = data_cache.segments.segment_ids(mask);
            seg_rankings{s} = computeDTWRanking(q_seg, c_segs, seg_ids, 'segment_id', mode, config);
        end
        
        dtw_cache.(qf).(mode).trajectory_ranking = traj_ranking;
        dtw_cache.(qf).(mode).segment_rankings = seg_rankings;
    end
    
    fprintf('  Query %d/%d done\n', q_idx, length(query_ids));
end

fprintf('DTW done: %.1f min\n', toc/60);
end

%% === HELPER ===

function ranking = computeDTWRanking(query_seq, cand_seqs, cand_ids, id_field, mode, cfg)
    n = length(cand_seqs);
    
    % LB_Kim
    lb_kim = zeros(n, 1);
    for i = 1:n
        lb_kim(i) = LB_Kim(query_seq, cand_seqs{i}, mode, cfg.rot_align, cfg.normalize);
    end
    [~, kim_ord] = sort(lb_kim);
    kim_n = round(n * cfg.lb_kim_ratio);
    kim_idx = kim_ord(1:kim_n);
    
    % LB_Keogh
    lb_keogh = inf(n, 1);
    for i = kim_idx'
        lb_keogh(i) = LB_Keogh(query_seq, cand_seqs{i}, cfg.window, mode, cfg.rot_align, cfg.normalize);
    end
    [~, keogh_ord] = sort(lb_keogh);
    keogh_n = min(cfg.lb_keogh_n, kim_n);
    keogh_idx = keogh_ord(1:keogh_n);
    
    % DTW
    dtw_dist = inf(n, 1);
    for i = keogh_idx'
        dtw_dist(i) = cDTW(query_seq, cand_seqs{i}, mode, cfg.window, inf, cfg.rot_align, cfg.normalize);
    end
    
    % Build table
    ranking = table(cand_ids(keogh_idx), dtw_dist(keogh_idx), 'VariableNames', {id_field, 'dtw_distance'});
    ranking = sortrows(ranking, 'dtw_distance');
end
function [dtw_reranked_ids, time_stage2, num_dtw_calls] = performDTWReranking(...
    top_k_ids, query_seq, data_cache, dtw_mode, dtw_config, K)
    % Perform DTW reranking on trajectory-level candidates
    
    tic_stage2 = tic;
    
    % Adaptive LB Strategy
    if K <= 50
        use_lb_kim = false;
        use_lb_keogh = false;
        lb_kim_keep = K;
        lb_keogh_keep = K;
    elseif K <= 200
        use_lb_kim = true;
        use_lb_keogh = false;
        lb_kim_keep = round(K * 0.6);
        lb_keogh_keep = K;
    else
        use_lb_kim = true;
        use_lb_keogh = true;
        lb_kim_keep = round(K * 0.7);
        lb_keogh_keep = round(K * 0.3);
    end
    
    % PHASE 1: LB_Kim
    if use_lb_kim
        lb_kim_scores = zeros(K, 1);
        for i = 1:K
            candidate_id = top_k_ids{i};
            cand_idx = find(strcmp(data_cache.candidates.bahn_ids, candidate_id), 1);
            
            if ~isempty(cand_idx)
                if strcmp(dtw_mode, 'position')
                    cand_seq = data_cache.candidates.position{cand_idx};
                else
                    cand_seq = data_cache.candidates.joint{cand_idx};
                end
                lb_kim_scores(i) = LB_Kim(query_seq, cand_seq);
            else
                lb_kim_scores(i) = inf;
            end
        end
        [~, kim_sort_idx] = sort(lb_kim_scores);
        candidates_after_kim = kim_sort_idx(1:lb_kim_keep);
    else
        candidates_after_kim = 1:K;
    end
    
    % PHASE 2: LB_Keogh
    if use_lb_keogh
        lb_keogh_scores = zeros(length(candidates_after_kim), 1);
        for i = 1:length(candidates_after_kim)
            cand_original_idx = candidates_after_kim(i);
            candidate_id = top_k_ids{cand_original_idx};
            cand_idx = find(strcmp(data_cache.candidates.bahn_ids, candidate_id), 1);
            
            if ~isempty(cand_idx)
                if strcmp(dtw_mode, 'position')
                    cand_seq = data_cache.candidates.position{cand_idx};
                else
                    cand_seq = data_cache.candidates.joint{cand_idx};
                end
                lb_keogh_scores(i) = LB_Keogh(query_seq, cand_seq, dtw_config.cdtw_window);
            else
                lb_keogh_scores(i) = inf;
            end
        end
        [~, keogh_sort_idx] = sort(lb_keogh_scores);
        candidates_after_keogh_local = keogh_sort_idx(1:min(lb_keogh_keep, length(keogh_sort_idx)));
        final_candidates = candidates_after_kim(candidates_after_keogh_local);
    else
        final_candidates = candidates_after_kim;
    end
    
    % PHASE 3: Full DTW
    num_dtw_calls = length(final_candidates);
    dtw_scores = inf(K, 1);
    
    for i = 1:length(final_candidates)
        cand_idx_in_topk = final_candidates(i);
        candidate_id = top_k_ids{cand_idx_in_topk};
        cand_idx = find(strcmp(data_cache.candidates.bahn_ids, candidate_id), 1);
        
        if ~isempty(cand_idx)
            if strcmp(dtw_mode, 'position')
                cand_seq = data_cache.candidates.position{cand_idx};
            else
                cand_seq = data_cache.candidates.joint{cand_idx};
            end
            
            dtw_dist = cDTW(query_seq, cand_seq, dtw_mode, ...
                dtw_config.cdtw_window, inf, ...
                dtw_config.use_rotation_alignment, ...
                dtw_config.normalize_dtw);
            
            dtw_scores(cand_idx_in_topk) = dtw_dist;
        end
    end
    
    % Sort by DTW score
    [~, sort_idx] = sort(dtw_scores);
    dtw_reranked_ids = top_k_ids(sort_idx);
    
    time_stage2 = toc(tic_stage2);
end
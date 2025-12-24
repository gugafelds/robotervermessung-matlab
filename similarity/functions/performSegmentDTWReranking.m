function [dtw_reranked_ids, time_stage2, num_dtw_calls, used_lb_kim, used_lb_keogh] = performSegmentDTWReranking(...
    top_k_seg_ids, query_seg_seq, data_cache, dtw_mode, dtw_config, K)
    % Perform DTW reranking on segment-level candidates
    %
    % UPDATED: Now returns used_lb_kim and used_lb_keogh flags
    
    tic_stage2 = tic;
    
    % Adaptive LB Strategy (same as trajectory level)
    if K <= 500
        use_lb_kim = false;
        use_lb_keogh = false;
        lb_kim_keep = K;
        lb_keogh_keep = K;
    else
        use_lb_kim = true;
        use_lb_keogh = true;
        lb_kim_keep = round(K * 0.9);
        lb_keogh_keep = 500;
    end
    
    % PHASE 1: LB_Kim
    if use_lb_kim
        lb_kim_scores = zeros(K, 1);
        for i = 1:K
            segment_id = top_k_seg_ids{i};
            seg_idx = find(strcmp(data_cache.segments.segment_ids, segment_id), 1);
            
            if ~isempty(seg_idx)
                if strcmp(dtw_mode, 'position')
                    cand_seq = data_cache.segments.position{seg_idx};
                else
                    cand_seq = data_cache.segments.joint{seg_idx};
                end
                lb_kim_scores(i) = LB_Kim(query_seg_seq, cand_seq);
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
            segment_id = top_k_seg_ids{cand_original_idx};
            seg_idx = find(strcmp(data_cache.segments.segment_ids, segment_id), 1);
            
            if ~isempty(seg_idx)
                if strcmp(dtw_mode, 'position')
                    cand_seq = data_cache.segments.position{seg_idx};
                else
                    cand_seq = data_cache.segments.joint{seg_idx};
                end
                lb_keogh_scores(i) = LB_Keogh(query_seg_seq, cand_seq, dtw_config.cdtw_window);
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
    
    % PHASE 3: Full DTW on segments
    num_dtw_calls = length(final_candidates);
    dtw_scores = inf(K, 1);
    
    for i = 1:length(final_candidates)
        cand_idx_in_topk = final_candidates(i);
        segment_id = top_k_seg_ids{cand_idx_in_topk};
        seg_idx = find(strcmp(data_cache.segments.segment_ids, segment_id), 1);
        
        if ~isempty(seg_idx)
            if strcmp(dtw_mode, 'position')
                cand_seq = data_cache.segments.position{seg_idx};
            else
                cand_seq = data_cache.segments.joint{seg_idx};
            end
            
            dtw_dist = cDTW(query_seg_seq, cand_seq, dtw_mode, ...
                dtw_config.cdtw_window, inf, ...
                dtw_config.use_rotation_alignment, ...
                dtw_config.normalize_dtw);
            
            dtw_scores(cand_idx_in_topk) = dtw_dist;
        end
    end
    
    % Sort by DTW score
    [~, sort_idx] = sort(dtw_scores);
    dtw_reranked_ids = top_k_seg_ids(sort_idx);
    
    time_stage2 = toc(tic_stage2);
    
    % Return lower bound usage flags
    used_lb_kim = use_lb_kim;
    used_lb_keogh = use_lb_keogh;
end
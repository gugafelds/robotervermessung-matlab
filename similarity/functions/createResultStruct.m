function result = createResultStruct(level, emb_name, n_coarse, n_fine, ...
    total_dims, multi_scale, weight_mode, dtw_mode, K, db_size, query_id, ...
    time_stage1, time_stage2, time_dtw_baseline, num_dtw_calls, ...
    emb_only_metrics, twostage_metrics, baseline_metrics, gt_type, ...
    num_gt_traj, num_gt_seg, embedding_only_ranking, twostage_final_ranking, ...
    used_lb_kim, used_lb_keogh)
    % Create standardized result struct
    %
    % UPDATED: Now accepts used_lb_kim and used_lb_keogh as parameters
    %          instead of computing them internally
    
    result = struct();
    
    % Configuration
    result.level = level;
    result.embedding_config = emb_name;
    result.n_coarse = n_coarse;
    result.n_fine = n_fine;
    result.total_dims = total_dims;
    result.multi_scale = multi_scale;
    result.weight_mode = weight_mode;
    result.dtw_mode = dtw_mode;
    result.k_candidates = K;
    result.db_size = db_size;
    result.query_id = query_id;
    
    % Timing
    result.time_stage1 = time_stage1;
    result.time_stage2 = time_stage2;
    result.time_total = time_stage1 + time_stage2;
    result.time_dtw_baseline = time_dtw_baseline;
    
    if ~isnan(time_dtw_baseline) && time_dtw_baseline > 0
        result.speedup = time_dtw_baseline / result.time_total;
    else
        result.speedup = NaN;
    end
    
    result.dtw_calls_made = num_dtw_calls;
    result.dtw_calls_saved_pct = (db_size - num_dtw_calls) / db_size * 100;
    
    % Embedding-Only Metrics
    result.embedding_only_r1 = emb_only_metrics.r1;
    result.embedding_only_r5 = emb_only_metrics.r5;
    result.embedding_only_r10 = emb_only_metrics.r10;
    result.embedding_only_r50 = emb_only_metrics.r50;
    result.embedding_only_rk = emb_only_metrics.rk;
    result.embedding_only_mrr = emb_only_metrics.mrr;
    result.embedding_only_mean_rank = emb_only_metrics.mean_rank;
    result.embedding_only_ndcg10 = emb_only_metrics.ndcg10;
    result.embedding_only_ndcg50 = emb_only_metrics.ndcg50;
    
    % Two-Stage Metrics
    result.recall_at_1 = twostage_metrics.r1;
    result.recall_at_5 = twostage_metrics.r5;
    result.recall_at_10 = twostage_metrics.r10;
    result.recall_at_50 = twostage_metrics.r50;
    result.recall_at_k = twostage_metrics.rk;
    result.mrr = twostage_metrics.mrr;
    result.mean_rank = twostage_metrics.mean_rank;
    result.ndcg_10 = twostage_metrics.ndcg10;
    result.ndcg_50 = twostage_metrics.ndcg50;
    
    % Baseline Metrics
    result.baseline_r1 = baseline_metrics.r1_gt;
    result.baseline_r5 = baseline_metrics.r5_gt;
    result.baseline_r10 = baseline_metrics.r10_gt;
    result.baseline_r50 = baseline_metrics.r50_gt;
    result.baseline_rk = baseline_metrics.rk_gt;
    result.baseline_mean_rank = baseline_metrics.mean_rank;
    
    % GT Info
    result.gt_type = gt_type;
    result.num_gt_trajectories = num_gt_traj;
    result.num_gt_segments = num_gt_seg;
    
    % Rankings (optional, only for trajectory level)
    if ~isempty(embedding_only_ranking)
        result.embedding_only_ranking = embedding_only_ranking;
        result.twostage_final_ranking = twostage_final_ranking;
    end
    
    % Lower bounds usage (now from actual DTW reranking functions!)
    result.used_lb_kim = used_lb_kim;
    result.used_lb_keogh = used_lb_keogh;
end
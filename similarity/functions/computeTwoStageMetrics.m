function metrics = computeTwoStageMetrics(dtw_reranked_ids, gt_ids, gt_similarities, num_gt, K)
    % Same as computeEmbeddingOnlyMetrics but for two-stage results
    metrics = computeEmbeddingOnlyMetrics(dtw_reranked_ids, gt_ids, gt_similarities, num_gt, K);
end
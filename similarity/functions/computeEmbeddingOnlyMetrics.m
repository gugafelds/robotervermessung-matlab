function metrics = computeEmbeddingOnlyMetrics(top_k_ids, gt_ids, num_gt, K)
    % Compute metrics for embedding-only retrieval
    
    if num_gt == 0 || isempty(gt_ids)
        % No GT available
        metrics.r1 = NaN;
        metrics.r5 = NaN;
        metrics.r10 = NaN;
        metrics.r50 = NaN;
        metrics.mrr = NaN;
        metrics.mean_rank = NaN;
        metrics.ndcg10 = NaN;
        metrics.ndcg50 = NaN;
        return;
    end
    
    % Find GT ranks
    gt_ranks = inf(num_gt, 1);
    for gt_idx = 1:num_gt
        gt_id = gt_ids{gt_idx};
        rank = find(strcmp(top_k_ids, gt_id), 1);
        if ~isempty(rank)
            gt_ranks(gt_idx) = rank;
        end
    end
    
    valid_ranks = gt_ranks(gt_ranks < inf);
    
    if isempty(valid_ranks)
        % No GT found
        metrics.r1 = 0;
        metrics.r5 = 0;
        metrics.r10 = 0;
        metrics.r50 = 0;
        metrics.mrr = 0;
        metrics.mean_rank = inf;
        metrics.ndcg10 = 0;
        metrics.ndcg50 = NaN;
        return;
    end
    
    % Recall@K
    metrics.r1 = sum(valid_ranks <= 1) / min(1, num_gt);
    metrics.r5 = sum(valid_ranks <= 5) / min(5, num_gt);
    metrics.r10 = sum(valid_ranks <= 10) / min(10, num_gt);
    metrics.r50 = sum(valid_ranks <= 50) / min(50, num_gt);
    
    % MRR
    metrics.mrr = mean(1 ./ valid_ranks);
    
    % Mean Rank
    metrics.mean_rank = mean(valid_ranks);
    
    % NDCG@10
    dcg_10 = 0;
    for i = 1:length(valid_ranks)
        if valid_ranks(i) <= 10
            dcg_10 = dcg_10 + 1 / log2(valid_ranks(i) + 1);
        end
    end
    idcg_10 = sum(1 ./ log2((1:min(10, num_gt)) + 1));
    if idcg_10 > 0
        metrics.ndcg10 = dcg_10 / idcg_10;
    else
        metrics.ndcg10 = 0;
    end
    
    % NDCG@50 (only if K >= 50)
    if K >= 50
        dcg_50 = 0;
        for i = 1:length(valid_ranks)
            if valid_ranks(i) <= 50
                dcg_50 = dcg_50 + 1 / log2(valid_ranks(i) + 1);
            end
        end
        idcg_50 = sum(1 ./ log2((1:min(50, num_gt)) + 1));
        if idcg_50 > 0
            metrics.ndcg50 = dcg_50 / idcg_50;
        else
            metrics.ndcg50 = 0;
        end
    else
        metrics.ndcg50 = NaN;
    end
end
function metrics = computeEmbeddingOnlyMetrics(top_k_ids, gt_ids, gt_similarities, num_gt, K)
    % Compute metrics for embedding-only retrieval with GRADED relevance
    % gt_similarities: DTW-based relevance scores for each GT
    
    % Handle empty cases
    if num_gt == 0
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
    
    if isempty(gt_ids)
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
        metrics.ndcg50 = 0;
        return;
    end
    
    % Recall@K (unchanged - binary)
    metrics.r1 = sum(valid_ranks <= 1) / min(1, num_gt);
    metrics.r5 = sum(valid_ranks <= 5) / min(5, num_gt);
    metrics.r10 = sum(valid_ranks <= 10) / min(10, num_gt);
    metrics.r50 = sum(valid_ranks <= 50) / min(50, num_gt);
    
    % MRR (unchanged - binary)
    metrics.mrr = mean(1 ./ valid_ranks);
    
    % Mean Rank (unchanged - binary)
    metrics.mean_rank = mean(valid_ranks);
    
    % =========================================================================
    % NDCG@10 with GRADED relevance (DTW-based)
    % =========================================================================
    dcg_10 = 0;
    for i = 1:num_gt
        if gt_ranks(i) <= 10
            relevance = gt_similarities(i);  % ✅ Use DTW similarity as relevance!
            dcg_10 = dcg_10 + relevance / log2(gt_ranks(i) + 1);
        end
    end
    
    % IDCG@10 (perfect ranking = sorted by similarity DESC)
    [sorted_sims, ~] = sort(gt_similarities, 'descend');
    idcg_10 = 0;
    for i = 1:min(10, num_gt)
        idcg_10 = idcg_10 + sorted_sims(i) / log2(i + 1);
    end
    
    if idcg_10 > 0
        metrics.ndcg10 = dcg_10 / idcg_10;
    else
        metrics.ndcg10 = 0;
    end
    
    % =========================================================================
    % NDCG@50 with GRADED relevance (only if K >= 50)
    % =========================================================================
    if K >= 50
        dcg_50 = 0;
        for i = 1:num_gt
            if gt_ranks(i) <= 50
                relevance = gt_similarities(i);  % ✅ Use DTW similarity!
                dcg_50 = dcg_50 + relevance / log2(gt_ranks(i) + 1);
            end
        end
        
        % IDCG@50
        idcg_50 = 0;
        for i = 1:min(50, num_gt)
            idcg_50 = idcg_50 + sorted_sims(i) / log2(i + 1);
        end
        
        if idcg_50 > 0
            metrics.ndcg50 = dcg_50 / idcg_50;
        else
            metrics.ndcg50 = 0;
        end
    else
        metrics.ndcg50 = NaN;
    end
end
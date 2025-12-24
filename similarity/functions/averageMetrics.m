function avg_metrics = averageMetrics(metrics_array)
    % Average metrics across multiple segments
    
    if isempty(metrics_array)
        avg_metrics.r1 = NaN;
        avg_metrics.r5 = NaN;
        avg_metrics.r10 = NaN;
        avg_metrics.r50 = NaN;
        avg_metrics.rk = NaN;
        avg_metrics.mrr = NaN;
        avg_metrics.mean_rank = NaN;
        avg_metrics.ndcg10 = NaN;
        avg_metrics.ndcg50 = NaN;
        return;
    end
    
    % Extract and average each field
    avg_metrics.r1 = mean([metrics_array.r1], 'omitnan');
    avg_metrics.r5 = mean([metrics_array.r5], 'omitnan');
    avg_metrics.r10 = mean([metrics_array.r10], 'omitnan');
    avg_metrics.r50 = mean([metrics_array.r50], 'omitnan');
    avg_metrics.rk = mean([metrics_array.rk], 'omitnan');
    avg_metrics.mrr = mean([metrics_array.mrr], 'omitnan');
    avg_metrics.mean_rank = mean([metrics_array.mean_rank], 'omitnan');
    avg_metrics.ndcg10 = mean([metrics_array.ndcg10], 'omitnan');
    avg_metrics.ndcg50 = mean([metrics_array.ndcg50], 'omitnan');
end

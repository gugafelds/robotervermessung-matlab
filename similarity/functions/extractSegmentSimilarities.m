function gt_seg_similarities = extractSegmentSimilarities(dtw_cache, query_id, dtw_mode, seg_idx, gt_seg_ids)
    % Extract DTW similarities for ground truth segments
    %
    % INPUTS:
    %   dtw_cache    - DTW cache with segment rankings
    %   query_id     - Query ID (e.g., '1765473159')
    %   dtw_mode     - DTW mode ('position' or 'joint_states')
    %   seg_idx      - Query segment index (1, 2, 3, ...)
    %   gt_seg_ids   - Cell array of GT segment IDs
    %
    % OUTPUT:
    %   gt_seg_similarities - Array of DTW similarity scores [num_gt_seg × 1]
    
    query_field = sprintf('q_%s', strrep(query_id, '-', '_'));
    
    % Check if query exists in cache
    if ~isfield(dtw_cache, query_field) || ~isfield(dtw_cache.(query_field), dtw_mode)
        warning('Query %s with mode %s not found in DTW cache', query_id, dtw_mode);
        gt_seg_similarities = [];
        return;
    end
    
    % Check if segment rankings exist
    if ~isfield(dtw_cache.(query_field).(dtw_mode), 'segment_rankings')
        warning('No segment rankings found for query %s', query_id);
        gt_seg_similarities = [];
        return;
    end
    
    segment_rankings = dtw_cache.(query_field).(dtw_mode).segment_rankings;
    
    if seg_idx > length(segment_rankings) || isempty(segment_rankings{seg_idx})
        warning('Segment %d not found in rankings', seg_idx);
        gt_seg_similarities = [];
        return;
    end
    
    % Get DTW ranking for this segment
    seg_dtw_ranking = segment_rankings{seg_idx};
    
    % Extract similarities for each GT segment
    num_gt_seg = length(gt_seg_ids);
    gt_seg_similarities = zeros(num_gt_seg, 1);
    
    for gt_idx = 1:num_gt_seg
        gt_seg_id = gt_seg_ids{gt_idx};
        
        % Find GT segment in DTW ranking
        rank_idx = find(strcmp(seg_dtw_ranking.segment_id, gt_seg_id), 1);
        
        if ~isempty(rank_idx)
            % Get DTW distance
            dtw_distance = seg_dtw_ranking.dtw_distance(rank_idx);
            gt_seg_similarities(gt_idx) = dtw_distance;
        else
            % GT segment not found
            gt_seg_similarities(gt_idx) = inf;
        end
    end
    
    % Normalize distances to [0, 1] and convert to similarities
    valid_dists = gt_seg_similarities(gt_seg_similarities < inf);
    
    if ~isempty(valid_dists)
        max_dist = max(valid_dists);
        if max_dist > 0
            normalized_dists = gt_seg_similarities / max_dist;
        else
            normalized_dists = gt_seg_similarities;
        end
        gt_seg_similarities = 1 ./ (1 + normalized_dists);
        gt_seg_similarities(gt_seg_similarities == inf) = 0;  % ✅ KORRIGIERT!
    else
        gt_seg_similarities = zeros(num_gt_seg, 1);
    end
end
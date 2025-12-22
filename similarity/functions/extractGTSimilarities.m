function gt_similarities = extractGTSimilarities(dtw_cache, query_id, dtw_mode, gt_ids)
    % Extract DTW similarities for ground truth trajectories
    %
    % INPUTS:
    %   dtw_cache  - DTW cache with trajectory rankings
    %   query_id   - Query ID (e.g., '1765473159')
    %   dtw_mode   - DTW mode ('position' or 'joint_states')
    %   gt_ids     - Cell array of GT trajectory IDs
    %
    % OUTPUT:
    %   gt_similarities - Array of DTW similarity scores [num_gt × 1]
    
    query_field = sprintf('q_%s', strrep(query_id, '-', '_'));
    
    % Check if query exists in cache
    if ~isfield(dtw_cache, query_field) || ~isfield(dtw_cache.(query_field), dtw_mode)
        error('Query %s with mode %s not found in DTW cache', query_id, dtw_mode);
    end
    
    % Get DTW ranking
    dtw_ranking = dtw_cache.(query_field).(dtw_mode).trajectory_ranking;
    
    % Extract similarities for each GT
    num_gt = length(gt_ids);
    gt_similarities = zeros(num_gt, 1);
    
    for gt_idx = 1:num_gt
        gt_id = gt_ids{gt_idx};
        
        % Find GT in DTW ranking
        rank_idx = find(strcmp(dtw_ranking.bahn_id, gt_id), 1);
        
        if ~isempty(rank_idx)
            % Get DTW distance
            dtw_distance = dtw_ranking.dtw_distance(rank_idx);
            
            % Convert to similarity (same as in two_stage_retrieval.m)
            % Similarity = 1 / (1 + normalized_distance)
            % For single GT, we use raw distance and normalize later
            gt_similarities(gt_idx) = dtw_distance;
        else
            % GT not found in ranking (should not happen)
            gt_similarities(gt_idx) = inf;
        end
    end
    
    % Normalize distances to [0, 1] and convert to similarities
    valid_dists = gt_similarities(gt_similarities < inf);
    
    if ~isempty(valid_dists)
        max_dist = max(valid_dists);
        if max_dist > 0
            normalized_dists = gt_similarities / max_dist;
        else
            normalized_dists = gt_similarities;
        end
        gt_similarities = 1 ./ (1 + normalized_dists);
        gt_similarities(gt_similarities == inf) = 0;
    else
        gt_similarities = zeros(num_gt, 1);
    end
end
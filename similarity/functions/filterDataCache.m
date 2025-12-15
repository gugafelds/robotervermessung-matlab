function filtered_cache = filterDataCache(data_cache, candidate_ids)
    % FILTERDATACACHE - Filter data cache to only include specific candidates
    %
    % Inputs:
    %   data_cache - Original data cache (from loadDataExperiment)
    %   candidate_ids - Cell array of bahn_ids to keep
    %
    % Output:
    %   filtered_cache - Filtered data cache
    
    filtered_cache = struct();
    
    % === Filter Candidates (Trajectories) ===
    if isfield(data_cache, 'candidates')
        candidates = data_cache.candidates;
        
        % Find indices to keep
        [~, keep_idx] = ismember(candidate_ids, candidates.bahn_ids);
        keep_idx = keep_idx(keep_idx > 0);  % Remove zeros (not found)
        
        % Filter all candidate fields
        filtered_cache.candidates = struct();
        filtered_cache.candidates.bahn_ids = candidates.bahn_ids(keep_idx);
        filtered_cache.candidates.metadata = candidates.metadata(keep_idx, :);
        
        if isfield(candidates, 'position')
            filtered_cache.candidates.position = candidates.position(keep_idx);
        end
        if isfield(candidates, 'joint')
            filtered_cache.candidates.joint = candidates.joint(keep_idx);
        end
        if isfield(candidates, 'orientation')
            filtered_cache.candidates.orientation = candidates.orientation(keep_idx);
        end
    end
    
    % === Filter Segments ===
    if isfield(data_cache, 'segments')
        segments = data_cache.segments;
        
        % Keep segments whose parent trajectory (bahn_id) is in candidate_ids
        seg_keep_mask = ismember(segments.bahn_ids, candidate_ids);
        
        % Filter all segment fields
        filtered_cache.segments = struct();
        filtered_cache.segments.segment_ids = segments.segment_ids(seg_keep_mask);
        filtered_cache.segments.bahn_ids = segments.bahn_ids(seg_keep_mask);
        filtered_cache.segments.metadata = segments.metadata(seg_keep_mask, :);
        
        if isfield(segments, 'position')
            filtered_cache.segments.position = segments.position(seg_keep_mask);
        end
        if isfield(segments, 'joint')
            filtered_cache.segments.joint = segments.joint(seg_keep_mask);
        end
        if isfield(segments, 'orientation')
            filtered_cache.segments.orientation = segments.orientation(seg_keep_mask);
        end
    end
    
    % === Copy Queries (unchanged) ===
    if isfield(data_cache, 'queries')
        filtered_cache.queries = data_cache.queries;
    end
    
end
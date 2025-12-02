function [ground_truth_ids, ground_truth_map] = getGTCandidates(conn, schema, query_ids)
% GETGTCANDIDATES - Find GT trajectories AND segments from same recordings
%
%   Finds all trajectories AND their segments that share the same 
%   record_filename as the query trajectories. This gets everything in
%   one efficient query.
%
%   INPUTS:
%       conn       - Database connection object
%       schema     - Database schema name (e.g., 'bewegungsdaten')
%       query_ids  - Cell array of query trajectory IDs
%
%   OUTPUTS:
%       ground_truth_ids - Cell array of ALL GT IDs (trajectories + segments)
%       ground_truth_map - Struct mapping each query to its GT:
%           .q_QUERYID.trajectories - Cell array of GT trajectory IDs
%           .q_QUERYID.segments     - Struct with segment IDs per segment index
%               .seg_1 - Cell array of GT segment IDs for segment 1
%               .seg_2 - Cell array of GT segment IDs for segment 2
%               ...
%
%   EXAMPLE OUTPUT:
%       ground_truth_map.q_1763567277.trajectories = {'1763567278', ..., '1763567286'};
%       ground_truth_map.q_1763567277.segments.seg_1 = {'1763567278_1', ..., '1763567286_1'};
%       ground_truth_map.q_1763567277.segments.seg_2 = {'1763567278_2', ..., '1763567286_2'};
%
%   Author: Gustavo Barros
%   Date: 03.12.2025

fprintf('\n=== Finding Ground Truth Trajectories & Segments ===\n');

num_queries = length(query_ids);
ground_truth_map = struct();
all_ground_truth_ids = {};

for q_idx = 1:num_queries
    query_id = query_ids{q_idx};
    
    fprintf('\nProcessing Query %d/%d: %s\n', q_idx, num_queries, query_id);
    
    % ====================================================================
    % STEP 1: Get record_filename for this query
    % ====================================================================
    
    recording_query = sprintf(...
        ['SELECT record_filename ' ...
         'FROM robotervermessung.%s.bahn_info ' ...
         'WHERE bahn_id = ''%s'''], ...
        schema, query_id);
    
    recording_result = fetch(conn, recording_query);
    
    if isempty(recording_result) || isempty(recording_result.record_filename)
        fprintf('  ⚠ No record_filename found - skipping\n');
        
        % Initialize empty GT for this query
        query_field = sprintf('q_%s', strrep(query_id, '-', '_'));
        ground_truth_map.(query_field).trajectories = {};
        ground_truth_map.(query_field).segments = struct();
        continue;
    end
    
    record_filename = recording_result.record_filename{1};
    fprintf('  Record: %s\n', record_filename);
    
    % ====================================================================
    % STEP 2: Get ALL bahn_ids + segment_ids from same recording
    % ====================================================================
    % This gets EVERYTHING (trajectories + all segments) in ONE query!
    
    gt_query = sprintf(...
        ['SELECT bm.bahn_id, bm.segment_id ' ...
         'FROM robotervermessung.%s.bahn_info bi ' ...
         'JOIN robotervermessung.%s.bahn_metadata bm ' ...
         'ON bi.bahn_id = bm.bahn_id ' ...
         'WHERE bi.record_filename = ''%s'' ' ...
         'AND bi.bahn_id != ''%s'' ' ...  % Exclude query itself
         'ORDER BY bm.bahn_id, bm.segment_id'], ...
        schema, schema, record_filename, query_id);
    
    gt_result = fetch(conn, gt_query);
    
    if isempty(gt_result)
        fprintf('  ⚠ No ground truth found\n');
        
        % Initialize empty GT for this query
        query_field = sprintf('q_%s', strrep(query_id, '-', '_'));
        ground_truth_map.(query_field).trajectories = {};
        ground_truth_map.(query_field).segments = struct();
        continue;
    end
    
    % ====================================================================
    % STEP 3: Separate trajectories from segments
    % ====================================================================
    % Rule: segment_id == bahn_id → Trajectory
    %       segment_id != bahn_id → Segment
    
    bahn_ids = gt_result.bahn_id;
    segment_ids = gt_result.segment_id;

    if isstring(bahn_ids)
        bahn_ids = cellstr(bahn_ids);
    end
    if isstring(segment_ids)
        segment_ids = cellstr(segment_ids);
    end
    
    % Find trajectories (where segment_id == bahn_id)
    is_trajectory = strcmp(bahn_ids, segment_ids);
    gt_trajectories = segment_ids(is_trajectory);
    
    % Find segments (where segment_id != bahn_id)
    is_segment = ~is_trajectory;
    gt_segments_all = segment_ids(is_segment);
    
    num_gt_traj = length(gt_trajectories);
    num_gt_seg = length(gt_segments_all);
    
    fprintf('  Found %d GT trajectories\n', num_gt_traj);
    fprintf('  Found %d GT segments\n', num_gt_seg);
    
    % ====================================================================
    % STEP 4: Organize segments by segment index
    % ====================================================================
    % Parse segment IDs to group by segment index
    % E.g., '1763567278_1', '1763567279_1' → seg_1
    %       '1763567278_2', '1763567279_2' → seg_2
    
    segments_by_index = struct();
    
    for seg_id_idx = 1:length(gt_segments_all)
        seg_id = gt_segments_all{seg_id_idx};
        
        % Extract segment index from ID (e.g., '1763567278_2' → '2')
        parts = strsplit(seg_id, '_');
        if length(parts) >= 2
            seg_index_str = parts{end};  % Last part is segment index
            seg_field = sprintf('seg_%s', seg_index_str);
            
            % Add to struct
            if ~isfield(segments_by_index, seg_field)
                segments_by_index.(seg_field) = {};
            end
            segments_by_index.(seg_field){end+1} = seg_id;
        end
    end
    
    % Convert cell arrays to column vectors
    seg_fields = fieldnames(segments_by_index);
    for sf_idx = 1:length(seg_fields)
        sf = seg_fields{sf_idx};
        segments_by_index.(sf) = segments_by_index.(sf)';  % Convert to column
    end
    
    % Display segment organization
    if ~isempty(seg_fields)
        fprintf('  Segments organized into %d groups:\n', length(seg_fields));
        for sf_idx = 1:length(seg_fields)
            sf = seg_fields{sf_idx};
            fprintf('    %s: %d segments\n', sf, length(segments_by_index.(sf)));
        end
    end
    
    % ====================================================================
    % STEP 5: Store in ground_truth_map
    % ====================================================================
    
    query_field = sprintf('q_%s', strrep(query_id, '-', '_'));
    ground_truth_map.(query_field).trajectories = gt_trajectories;
    ground_truth_map.(query_field).segments = segments_by_index;
    
    % ====================================================================
    % STEP 6: Add to global list (all IDs)
    % ====================================================================
    
    % ⭐ CRITICAL: Only add TRAJECTORIES to global list!
    % Segments are already included when their parent trajectory is loaded
    all_ground_truth_ids = [all_ground_truth_ids; gt_trajectories];
    % NOT: all_ground_truth_ids = [all_ground_truth_ids; gt_trajectories; gt_segments_all];
end

% ========================================================================
% FINAL SUMMARY
% ========================================================================

% Remove duplicates
ground_truth_ids = unique(all_ground_truth_ids);

fprintf('\n=== Ground Truth Summary ===\n');
fprintf('Total unique GT trajectories (for data loading): %d\n', length(ground_truth_ids));
fprintf('Total GT segments (stored in map): %d\n', ...
    sum(cellfun(@(x) sum(structfun(@length, x.segments)), struct2cell(ground_truth_map))));
fprintf('Queries with ground truth: %d/%d\n', ...
    sum(structfun(@(x) ~isempty(x.trajectories), ground_truth_map)), num_queries);

if ~isempty(ground_truth_ids)
    fprintf('\n✓ Ground truth ready for experiments\n');
    fprintf('  → Trajectories will be loaded via loadDataExperiment()\n');
    fprintf('  → Segments will be used from ground_truth_map for coverage analysis\n');
else
    fprintf('\n⚠ No ground truth found - experiments will use only random sampling\n');
end

fprintf('========================================\n\n');

end
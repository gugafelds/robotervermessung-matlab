function [ground_truth_ids, ground_truth_map] = getGroundTruthCandidates(conn, schema, query_ids)
% GETGROUNDTRUTHCANDIDATES - Find trajectories from same recordings as queries
%
%   Finds all trajectories that share the same record_filename as the
%   query trajectories. These serve as ground truth for evaluation since
%   they were recorded under identical conditions.
%
%   INPUTS:
%       conn       - Database connection object
%       schema     - Database schema name (e.g., 'bewegungsdaten')
%       query_ids  - Cell array of query trajectory IDs
%
%   OUTPUTS:
%       ground_truth_ids - Cell array of all ground truth trajectory IDs
%       ground_truth_map - Struct mapping each query to its ground truth:
%           .(query_id) - Cell array of ground truth IDs for this query
%
%   USAGE:
%       [gt_ids, gt_map] = getGroundTruthCandidates(conn, schema, query_ids);
%       candidate_ids = [random_sample; gt_ids];  % Mix random + ground truth
%
%   Author: Gustavo Barros
%   Date: 29.11.2025

fprintf('\n=== Finding Ground Truth Trajectories ===\n');

num_queries = length(query_ids);
ground_truth_map = struct();
all_ground_truth_ids = {};

for q_idx = 1:num_queries
    query_id = query_ids{q_idx};
    
    % Get record_filename for this query
    recording_query = sprintf(...
        ['SELECT record_filename ' ...
         'FROM robotervermessung.%s.bahn_info ' ...
         'WHERE bahn_id = ''%s'''], ...
        schema, query_id);
    
    recording_result = fetch(conn, recording_query);
    
    if isempty(recording_result) || isempty(recording_result.record_filename)
        fprintf('  Query %s: No record_filename found - skipping\n', query_id);
        ground_truth_map.(sprintf('q_%s', strrep(query_id, '-', '_'))) = {};
        continue;
    end
    
    record_filename = recording_result.record_filename{1};
    
    % Find all trajectories with same record_filename (excluding the query itself)
    gt_query = sprintf(...
        ['SELECT bi.bahn_id ' ...
         'FROM robotervermessung.%s.bahn_info bi ' ...
         'JOIN robotervermessung.%s.bahn_metadata bm ' ...
         'ON bi.bahn_id = bm.bahn_id ' ...
         'WHERE bi.record_filename = ''%s'' ' ...
         'AND bi.bahn_id != ''%s'' ' ...
         'AND bm.bahn_id = bm.segment_id'], ...  
        schema, schema, record_filename, query_id);
    
    gt_result = fetch(conn, gt_query);
    
    if isempty(gt_result)
        fprintf('  Query %s (%s): No ground truth found\n', query_id, record_filename);
        ground_truth_map.(sprintf('q_%s', strrep(query_id, '-', '_'))) = {};
    else
        gt_ids = gt_result.bahn_id;
        num_gt = length(gt_ids);
        
        fprintf('  Query %s (%s): Found %d ground truth trajectories\n', ...
            query_id, record_filename, num_gt);
        
        % Store in map
        query_field = sprintf('q_%s', strrep(query_id, '-', '_'));
        ground_truth_map.(query_field) = gt_ids;
        
        % Add to global list
        all_ground_truth_ids = [all_ground_truth_ids; gt_ids];
    end
end

% Remove duplicates (in case multiple queries share recordings)
ground_truth_ids = unique(all_ground_truth_ids);

fprintf('\n--- Ground Truth Summary ---\n');
fprintf('  Total unique ground truth trajectories: %d\n', length(ground_truth_ids));
fprintf('  Queries with ground truth: %d/%d\n', ...
    sum(structfun(@(x) ~isempty(x), ground_truth_map)), num_queries);

if ~isempty(ground_truth_ids)
    fprintf('\n✓ Ground truth candidates ready for mixing with random sample\n');
else
    fprintf('\n⚠ No ground truth found - experiments will use only random sampling\n');
end

fprintf('========================================\n\n');

end
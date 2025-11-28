function data_cache = loadDataExperiment(conn, schema, candidate_ids, query_ids, chunk_size)
% LOADDATAEXPERIMENT - Pre-load all data for multiple experiments
%
%   This function loads ALL necessary data ONCE before running experiments,
%   avoiding redundant database queries and massively improving performance.
%
%   INPUTS:
%       conn          - Database connection object
%       schema        - Database schema name (e.g., 'bewegungsdaten')
%       candidate_ids - Cell array of candidate trajectory IDs (e.g., 1100 sampled)
%       query_ids     - Cell array of query trajectory IDs (e.g., 4 queries)
%       chunk_size    - Chunk size for batch loading (default: 100)
%
%   OUTPUT:
%       data_cache - Struct containing all pre-loaded data:
%           .candidates - Candidate trajectory data
%               .bahn_ids       - Cell array of IDs
%               .metadata       - Table with metadata (length, duration, twist, accel, etc.)
%               .position       - Cell array of position trajectories
%               .joint          - Cell array of joint state trajectories
%               .orientation    - Cell array of orientation trajectories
%           .segments - Segment data for all candidates
%               .segment_ids    - Cell array of segment IDs
%               .bahn_ids       - Parent trajectory IDs
%               .metadata       - Table with segment metadata
%               .position       - Cell array of position segments
%               .joint          - Cell array of joint state segments
%               .orientation    - Cell array of orientation segments
%           .queries - Query trajectory data (one entry per query)
%               .(query_id) - Struct for each query
%                   .metadata       - Query metadata struct
%                   .position       - Position trajectory
%                   .joint          - Joint state trajectory
%                   .orientation    - Orientation trajectory
%                   .segments       - Cell array of query segments
%                       .metadata   - Segment metadata
%                       .position   - Position segments
%                       .joint      - Joint state segments
%                       .orientation - Orientation segments
%
%   USAGE:
%       % Before experiments
%       data_cache = loadDataExperiment(conn, schema, candidate_ids, query_ids, 100);
%       
%       % Pass to experiments via config
%       config.data_cache = data_cache;
%       results = runExperiment(config);
%
%   PERFORMANCE:
%       Loading 1100 candidates + 4 queries + segments:
%       - Time: ~15-20 minutes (one-time cost)
%       - Memory: ~500-800 MB
%       - Saves: ~3-4 hours over 128 experiments (40-50% reduction!)
%
%   Author: Gustavo Barros
%   Date: 28.11.2025

%% ========================================================================
%  VALIDATION & SETUP
%  ========================================================================

if nargin < 5
    chunk_size = 100;
end

fprintf('\n========================================\n');
fprintf('PRE-LOADING EXPERIMENT DATA\n');
fprintf('========================================\n\n');

num_candidates = length(candidate_ids);
num_queries = length(query_ids);

fprintf('Configuration:\n');
fprintf('  Candidates: %d trajectories\n', num_candidates);
fprintf('  Queries: %d trajectories\n', num_queries);
fprintf('  Chunk size: %d\n', chunk_size);
fprintf('  Schema: %s\n\n', schema);

data_cache = struct();
overall_tic = tic;

%% ========================================================================
%  PART 1: LOAD CANDIDATE TRAJECTORIES
%  ========================================================================

fprintf('=== PART 1/3: Loading Candidate Trajectories ===\n');

candidates_tic = tic;

% ------------------------------------------------------------------------
% STEP 1.1: Load Candidate Metadata
% ------------------------------------------------------------------------
fprintf('  Step 1/4: Loading metadata...\n');

bahn_id_list = sprintf('''%s''', strjoin(candidate_ids, ''','''));

metadata_query = sprintf(...
    ['SELECT bahn_id, length, duration, movement_type, ' ...
     'min_twist_ist, max_twist_ist, mean_twist_ist, median_twist_ist, std_twist_ist, ' ...
     'min_acceleration_ist, max_acceleration_ist, mean_acceleration_ist, ' ...
     'median_acceleration_ist, std_acceleration_ist ' ...
     'FROM robotervermessung.%s.bahn_metadata ' ...
     'WHERE bahn_id = segment_id ' ...
     'AND bahn_id IN (%s)'], ...
    schema, bahn_id_list);

candidate_metadata = fetch(conn, metadata_query);

fprintf('    ✓ Loaded metadata for %d trajectories\n', height(candidate_metadata));

% ------------------------------------------------------------------------
% STEP 1.2: Load Position Data
% ------------------------------------------------------------------------
fprintf('  Step 2/4: Loading position data...\n');

position_data = batchLoadTrajectoriesChunked(conn, schema, candidate_ids, 'position', chunk_size);

fprintf('    ✓ Loaded position data for %d trajectories\n', length(position_data));

% ------------------------------------------------------------------------
% STEP 1.3: Load Joint State Data
% ------------------------------------------------------------------------
fprintf('  Step 3/4: Loading joint state data...\n');

joint_data = batchLoadTrajectoriesChunked(conn, schema, candidate_ids, 'joint_states', chunk_size);

fprintf('    ✓ Loaded joint state data for %d trajectories\n', length(joint_data));

% ------------------------------------------------------------------------
% STEP 1.4: Load Orientation Data
% ------------------------------------------------------------------------
fprintf('  Step 4/4: Loading orientation data...\n');

orientation_data = batchLoadTrajectoriesChunked(conn, schema, candidate_ids, 'orientation', chunk_size);

fprintf('    ✓ Loaded orientation data for %d trajectories\n', length(orientation_data));

% ------------------------------------------------------------------------
% Store Candidate Data
% ------------------------------------------------------------------------
data_cache.candidates = struct();
data_cache.candidates.bahn_ids = candidate_ids;
data_cache.candidates.metadata = candidate_metadata;
data_cache.candidates.position = position_data;
data_cache.candidates.joint = joint_data;
data_cache.candidates.orientation = orientation_data;

candidates_time = toc(candidates_tic);
fprintf('  ✓ Part 1 completed in %.1f minutes\n\n', candidates_time/60);

%% ========================================================================
%  PART 2: LOAD SEGMENT DATA (for all candidates)
%  ========================================================================

fprintf('=== PART 2/3: Loading Segment Data ===\n');

segments_tic = tic;

% ------------------------------------------------------------------------
% STEP 2.1: Load Segment Metadata
% ------------------------------------------------------------------------
fprintf('  Step 1/4: Loading segment metadata...\n');

segment_metadata_query = sprintf(...
    ['SELECT segment_id, bahn_id, length, duration, movement_type, ' ...
     'min_twist_ist, max_twist_ist, mean_twist_ist, median_twist_ist, std_twist_ist, ' ...
     'min_acceleration_ist, max_acceleration_ist, mean_acceleration_ist, ' ...
     'median_acceleration_ist, std_acceleration_ist ' ...
     'FROM robotervermessung.%s.bahn_metadata ' ...
     'WHERE segment_id != bahn_id ' ...
     'AND bahn_id IN (%s)'], ...
    schema, bahn_id_list);

segment_metadata = fetch(conn, segment_metadata_query);
num_segments = height(segment_metadata);

fprintf('    ✓ Found %d segments from %d trajectories\n', num_segments, num_candidates);

if num_segments > 0
    segment_ids = segment_metadata.segment_id;
    
    % ------------------------------------------------------------------------
    % STEP 2.2: Load Segment Position Data
    % ------------------------------------------------------------------------
    fprintf('  Step 2/4: Loading segment position data...\n');
    
    segment_position_data = batchLoadSegmentsChunked(conn, schema, segment_ids, 'position', chunk_size);
    
    fprintf('    ✓ Loaded position data for %d segments\n', length(segment_position_data));
    
    % ------------------------------------------------------------------------
    % STEP 2.3: Load Segment Joint State Data
    % ------------------------------------------------------------------------
    fprintf('  Step 3/4: Loading segment joint state data...\n');
    
    segment_joint_data = batchLoadSegmentsChunked(conn, schema, segment_ids, 'joint_states', chunk_size);
    
    fprintf('    ✓ Loaded joint state data for %d segments\n', length(segment_joint_data));
    
    % ------------------------------------------------------------------------
    % STEP 2.4: Load Segment Orientation Data
    % ------------------------------------------------------------------------
    fprintf('  Step 4/4: Loading segment orientation data...\n');
    
    segment_orientation_data = batchLoadSegmentsChunked(conn, schema, segment_ids, 'orientation', chunk_size);
    
    fprintf('    ✓ Loaded orientation data for %d segments\n', length(segment_orientation_data));
    
    % ------------------------------------------------------------------------
    % Store Segment Data
    % ------------------------------------------------------------------------
    data_cache.segments = struct();
    data_cache.segments.segment_ids = segment_ids;
    data_cache.segments.bahn_ids = segment_metadata.bahn_id;
    data_cache.segments.metadata = segment_metadata;
    data_cache.segments.position = segment_position_data;
    data_cache.segments.joint = segment_joint_data;
    data_cache.segments.orientation = segment_orientation_data;
else
    fprintf('    ⚠ No segments found\n');
    data_cache.segments = struct();
    data_cache.segments.segment_ids = {};
    data_cache.segments.bahn_ids = {};
    data_cache.segments.metadata = table();
    data_cache.segments.position = {};
    data_cache.segments.joint = {};
    data_cache.segments.orientation = {};
end

segments_time = toc(segments_tic);
fprintf('  ✓ Part 2 completed in %.1f minutes\n\n', segments_time/60);

%% ========================================================================
%  PART 3: LOAD QUERY TRAJECTORIES
%  ========================================================================

fprintf('=== PART 3/3: Loading Query Trajectories ===\n');

queries_tic = tic;

data_cache.queries = struct();

for q_idx = 1:num_queries
    query_id = query_ids{q_idx};
    
    fprintf('  Query %d/%d: %s\n', q_idx, num_queries, query_id);
    
    % ------------------------------------------------------------------------
    % STEP 3.1: Load Query Metadata
    % ------------------------------------------------------------------------
    fprintf('    Step 1/5: Loading metadata...\n');
    
    query_meta_query = sprintf(...
        ['SELECT movement_type, min_twist_ist, max_twist_ist, mean_twist_ist, ' ...
         'median_twist_ist, std_twist_ist, min_acceleration_ist, max_acceleration_ist, ' ...
         'mean_acceleration_ist, median_acceleration_ist, std_acceleration_ist ' ...
         'FROM robotervermessung.%s.bahn_metadata ' ...
         'WHERE bahn_id = ''%s'' AND segment_id = ''%s'''], ...
        schema, query_id, query_id);
    
    query_meta_result = fetch(conn, query_meta_query);
    
    % Also load basic metadata (length, duration, num_segments)
    query_metadata = loadTrajectoryMetadata(conn, schema, query_id);
    
    % Add extended metadata
    query_metadata.movement_type = query_meta_result.movement_type{1};
    query_metadata.min_twist = query_meta_result.min_twist_ist;
    query_metadata.max_twist = query_meta_result.max_twist_ist;
    query_metadata.mean_twist = query_meta_result.mean_twist_ist;
    query_metadata.median_twist = query_meta_result.median_twist_ist;
    query_metadata.std_twist = query_meta_result.std_twist_ist;
    query_metadata.min_acceleration = query_meta_result.min_acceleration_ist;
    query_metadata.max_acceleration = query_meta_result.max_acceleration_ist;
    query_metadata.mean_acceleration = query_meta_result.mean_acceleration_ist;
    query_metadata.median_acceleration = query_meta_result.median_acceleration_ist;
    query_metadata.std_acceleration = query_meta_result.std_acceleration_ist;
    
    fprintf('      ✓ Metadata: length=%.1fmm, duration=%.2fs, segments=%d\n', ...
        query_metadata.length, query_metadata.duration, query_metadata.num_segments);
    
    % ------------------------------------------------------------------------
    % STEP 3.2: Load Query Position Data
    % ------------------------------------------------------------------------
    fprintf('    Step 2/5: Loading position data...\n');
    
    query_position = loadTrajectoryPosition(conn, schema, query_id);
    
    fprintf('      ✓ Position: %d points\n', size(query_position, 1));
    
    % ------------------------------------------------------------------------
    % STEP 3.3: Load Query Joint State Data
    % ------------------------------------------------------------------------
    fprintf('    Step 3/5: Loading joint state data...\n');
    
    query_joint = loadTrajectoryJointStates(conn, schema, query_id);
    
    fprintf('      ✓ Joint states: %d points\n', size(query_joint, 1));
    
    % ------------------------------------------------------------------------
    % STEP 3.4: Load Query Orientation Data
    % ------------------------------------------------------------------------
    fprintf('    Step 4/5: Loading orientation data...\n');
    
    query_orientation = loadTrajectoryOrientation(conn, schema, query_id);
    
    fprintf('      ✓ Orientation: %d points\n', size(query_orientation, 1));
    
    % ------------------------------------------------------------------------
    % STEP 3.5: Load Query Segments
    % ------------------------------------------------------------------------
    fprintf('    Step 5/5: Loading segments...\n');
    
    num_query_segments = query_metadata.num_segments;
    
    query_segments_data = struct();
    query_segments_data.metadata = cell(num_query_segments, 1);
    query_segments_data.position = cell(num_query_segments, 1);
    query_segments_data.joint = cell(num_query_segments, 1);
    query_segments_data.orientation = cell(num_query_segments, 1);
    
    for seg_idx = 1:num_query_segments
        query_segment_id = sprintf('%s_%d', query_id, seg_idx);
        
        % Segment metadata
        seg_meta_query = sprintf(...
            ['SELECT movement_type, min_twist_ist, max_twist_ist, mean_twist_ist, ' ...
             'median_twist_ist, std_twist_ist, min_acceleration_ist, max_acceleration_ist, ' ...
             'mean_acceleration_ist, median_acceleration_ist, std_acceleration_ist ' ...
             'FROM robotervermessung.%s.bahn_metadata ' ...
             'WHERE segment_id = ''%s'''], ...
            schema, query_segment_id);
        seg_meta_result = fetch(conn, seg_meta_query);
        
        seg_metadata = loadSegmentMetadata(conn, schema, query_segment_id);
        seg_metadata.movement_type = seg_meta_result.movement_type{1};
        seg_metadata.min_twist = seg_meta_result.min_twist_ist;
        seg_metadata.max_twist = seg_meta_result.max_twist_ist;
        seg_metadata.mean_twist = seg_meta_result.mean_twist_ist;
        seg_metadata.median_twist = seg_meta_result.median_twist_ist;
        seg_metadata.std_twist = seg_meta_result.std_twist_ist;
        seg_metadata.min_acceleration = seg_meta_result.min_acceleration_ist;
        seg_metadata.max_acceleration = seg_meta_result.max_acceleration_ist;
        seg_metadata.mean_acceleration = seg_meta_result.mean_acceleration_ist;
        seg_metadata.median_acceleration = seg_meta_result.median_acceleration_ist;
        seg_metadata.std_acceleration = seg_meta_result.std_acceleration_ist;
        
        query_segments_data.metadata{seg_idx} = seg_metadata;
        
        % Segment trajectory data
        query_segments_data.position{seg_idx} = loadTrajectoryPosition(conn, schema, query_segment_id);
        query_segments_data.joint{seg_idx} = loadTrajectoryJointStates(conn, schema, query_segment_id);
        query_segments_data.orientation{seg_idx} = loadTrajectoryOrientation(conn, schema, query_segment_id);
    end
    
    fprintf('      ✓ Loaded %d segments\n', num_query_segments);
    
    % ------------------------------------------------------------------------
    % Store Query Data (use valid field name - replace dashes/underscores)
    % ------------------------------------------------------------------------
    query_field_name = ['q_' strrep(query_id, '-', '_')];
    
    data_cache.queries.(query_field_name) = struct();
    data_cache.queries.(query_field_name).bahn_id = query_id;
    data_cache.queries.(query_field_name).metadata = query_metadata;
    data_cache.queries.(query_field_name).position = query_position;
    data_cache.queries.(query_field_name).joint = query_joint;
    data_cache.queries.(query_field_name).orientation = query_orientation;
    data_cache.queries.(query_field_name).segments = query_segments_data;
    
    fprintf('    ✓ Query %d/%d completed\n\n', q_idx, num_queries);
end

queries_time = toc(queries_tic);
fprintf('  ✓ Part 3 completed in %.1f minutes\n\n', queries_time/60);

%% ========================================================================
%  SUMMARY
%  ========================================================================

total_time = toc(overall_tic);

fprintf('========================================\n');
fprintf('DATA LOADING COMPLETED\n');
fprintf('========================================\n\n');

fprintf('Total Time: %.1f minutes (%.2f seconds)\n\n', total_time/60, total_time);

fprintf('--- Breakdown ---\n');
fprintf('  Candidates:  %.1f min (%.1f%%)\n', candidates_time/60, 100*candidates_time/total_time);
fprintf('  Segments:    %.1f min (%.1f%%)\n', segments_time/60, 100*segments_time/total_time);
fprintf('  Queries:     %.1f min (%.1f%%)\n\n', queries_time/60, 100*queries_time/total_time);

fprintf('--- Data Summary ---\n');
fprintf('  Candidates:\n');
fprintf('    - Trajectories: %d\n', num_candidates);
fprintf('    - Modalities: Position, Joint States, Orientation\n');
fprintf('  Segments:\n');
fprintf('    - Total segments: %d\n', num_segments);
fprintf('    - Modalities: Position, Joint States, Orientation\n');
fprintf('  Queries:\n');
fprintf('    - Trajectories: %d\n', num_queries);
fprintf('  - Total query segments: %d\n', sum(cellfun(@(x) x.metadata.num_segments, struct2cell(data_cache.queries))));

% Memory usage
cache_info = whos('data_cache');
fprintf('\n--- Memory Usage ---\n');
fprintf('  Cache size: %.1f MB\n', cache_info.bytes / 1e6);

fprintf('\n✓ Data cache ready for experiments!\n');
fprintf('========================================\n\n');

end
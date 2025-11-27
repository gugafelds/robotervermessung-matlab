%% ========================================================================
%  EMBEDDING-BASED TRAJECTORY SIMILARITY SEARCH
%  ========================================================================
%  Three-stage retrieval pipeline:
%  1. pgvector HNSW embedding search
%  2. Cascading lower bounds: LB_Kim + LB_Keogh
%  3. Constrained DTW refinement
%
%  Author: Gustavo Barros
%  Date: 25.11.2025
%  ========================================================================

%% SECTION 1: CONFIGURATION & PARAMETERS
%  ========================================================================

fprintf('\n=== CONFIGURATION ===\n');

% === Query Parameters ===
query_bahn_id = '1763738777';        % ID of the query trajectory

% === Retrieval Pipeline Parameters ===
% Stage 1: pgvector HNSW retrieval
pgvector_top_k = 250;                % Retrieve top-250 per modality

% Stage 2: LB_Keogh pre-filtering
lb_keogh_keep = 100;                 % Keep top-120 after LB_Keogh

% Stage 3: DTW refinement
dtw_top_n = 50;                     % Final top-100 results

% === RRF Fusion Parameters ===
rrf_k = 60;                          % RRF constant (lower = more weight on top ranks)

% RRF weights for each modality
weight_position = 0.0;               % Geometric shape
weight_joint = 0.0;                  % Joint configuration
weight_orientation = 0.0;            % End-effector orientation
weight_velocity = 0.0;               % Velocity profile
weight_metadata = 1.0;               % Metadata features (disabled for now)

% === DTW Parameters ===
dtw_mode = 'joint_states';               % 'position' or 'joint_states'
cdtw_window = 0.05;                  % Constrained DTW window (5% of sequence length)
use_rotation_alignment = false;      % Apply Procrustes alignment before DTW
normalize_dtw = false;

chunk_size = 100;

% === Database Configuration ===
schema = 'bewegungsdaten';
db_name = 'robotervermessung';

% ========================================================================
% Display Configuration Summary
% ========================================================================

fprintf('\n--- Pipeline Configuration ---\n');
fprintf('Query Trajectory: %s\n', query_bahn_id);
fprintf('\n');
fprintf('Stage 1 - pgvector HNSW:  Top-%d per modality\n', pgvector_top_k);
fprintf('Stage 2 - LB_Keogh:       Keep top-%d\n', lb_keogh_keep);
fprintf('Stage 3 - DTW:            Final top-%d\n', dtw_top_n);
fprintf('\n');
fprintf('RRF Weights:\n');
fprintf('  Position:    %.2f\n', weight_position);
fprintf('  Joint:       %.2f\n', weight_joint);
fprintf('  Orientation: %.2f\n', weight_orientation);
fprintf('  Velocity:    %.2f\n', weight_velocity);
fprintf('  Metadata:    %.2f\n', weight_metadata);
fprintf('\n');
fprintf('DTW Mode: %s\n', dtw_mode);
fprintf('DTW Window: %.1f%%\n', cdtw_window * 100);
fprintf('\n');

%% SECTION 2: DATABASE CONNECTION
%  ========================================================================

fprintf('=== DATABASE CONNECTION ===\n');

conn = connectingToPostgres;

if isopen(conn)
    fprintf('✓ Connected to database: %s.%s\n', db_name, schema);

    execute(conn, 'SET search_path = "bewegungsdaten";');
else
    error('✗ Database connection failed');
end

fprintf('\n');

%% SECTION 3: LOAD QUERY EMBEDDINGS FROM DATABASE
%  ========================================================================

fprintf('=== LOADING QUERY EMBEDDINGS ===\n');

% ========================================================================
% 3.1: Load Trajectory-Level Query Embedding
% ========================================================================

fprintf('\n--- Trajectory Level ---\n');

query_trajectory_embedding_query = sprintf(...
    ['SELECT bahn_id, segment_id, ' ...
     'position_embedding::text as position_embedding, ' ...
     'joint_embedding::text as joint_embedding, ' ...
     'orientation_embedding::text as orientation_embedding, ' ...
     'velocity_embedding::text as velocity_embedding, ' ...
     'metadata_embedding::text as metadata_embedding ' ...
     'FROM %s.%s.bahn_embeddings ' ...
     'WHERE bahn_id = ''%s'' AND segment_id = ''%s'''], ...
    db_name, schema, query_bahn_id, query_bahn_id);

query_traj_result = fetch(conn, query_trajectory_embedding_query);

if isempty(query_traj_result) || height(query_traj_result) == 0
    error('Query trajectory embedding not found in database: %s', query_bahn_id);
end

% Parse embeddings (stored as text strings: '[1.2, 3.4, ...]')
query_position_emb = parseEmbeddingString(query_traj_result.position_embedding{1});
query_joint_emb = parseEmbeddingString(query_traj_result.joint_embedding{1});
query_orientation_emb = parseEmbeddingString(query_traj_result.orientation_embedding{1});
query_velocity_emb = parseEmbeddingString(query_traj_result.velocity_embedding{1});
query_metadata_emb = parseEmbeddingString(query_traj_result.metadata_embedding{1});

fprintf('✓ Query trajectory embedding loaded: %s\n', query_bahn_id);
fprintf('  Position:    %d dimensions\n', length(query_position_emb));
fprintf('  Joint:       %d dimensions\n', length(query_joint_emb));
fprintf('  Orientation: %d dimensions\n', length(query_orientation_emb));
fprintf('  Velocity:    %d dimensions\n', length(query_velocity_emb));
fprintf('  Metadata:    %d dimensions\n', length(query_metadata_emb));

% ========================================================================
% 3.2: Check which embeddings are available
% ========================================================================

fprintf('\nAvailable Embeddings:\n');

has_position = ~isempty(query_position_emb) && length(query_position_emb) > 0;
has_joint = ~isempty(query_joint_emb) && length(query_joint_emb) > 0;
has_orientation = ~isempty(query_orientation_emb) && length(query_orientation_emb) > 0;
has_velocity = ~isempty(query_velocity_emb) && length(query_velocity_emb) > 0;
has_metadata = ~isempty(query_metadata_emb) && length(query_metadata_emb) > 0;

fprintf('  Position:    %s\n', iif(has_position, '✓', '✗'));
fprintf('  Joint:       %s\n', iif(has_joint, '✓', '✗'));
fprintf('  Orientation: %s\n', iif(has_orientation, '✓', '✗'));
fprintf('  Velocity:    %s\n', iif(has_velocity, '✓', '✗'));
fprintf('  Metadata:    %s\n', iif(has_metadata, '✓', '✗'));

% Warn if weights are set but embeddings are missing
if weight_position > 0 && ~has_position
    warning('Position weight > 0 but embedding not available!');
end
if weight_joint > 0 && ~has_joint
    warning('Joint weight > 0 but embedding not available!');
end
if weight_orientation > 0 && ~has_orientation
    warning('Orientation weight > 0 but embedding not available!');
end
if weight_velocity > 0 && ~has_velocity
    warning('Velocity weight > 0 but embedding not available!');
end
if weight_metadata > 0 && ~has_metadata
    warning('Metadata weight > 0 but embedding not available!');
end

% ========================================================================
% 3.3: Load Segment-Level Query Embeddings
% ========================================================================

fprintf('\n--- Segment Level ---\n');

query_segments_embedding_query = sprintf(...
    ['SELECT segment_id, bahn_id, ' ...
     'position_embedding::text as position_embedding, ' ...
     'joint_embedding::text as joint_embedding, ' ...
     'orientation_embedding::text as orientation_embedding, ' ...
     'velocity_embedding::text as velocity_embedding, ' ...
     'metadata_embedding::text as metadata_embedding ' ...
     'FROM %s.%s.bahn_embeddings ' ...
     'WHERE bahn_id = ''%s'' AND segment_id != ''%s'' ' ...
     'ORDER BY segment_id'], ...
    db_name, schema, query_bahn_id, query_bahn_id);

query_seg_result = fetch(conn, query_segments_embedding_query);

num_query_segments = height(query_seg_result);

if num_query_segments == 0
    fprintf('⚠ No segments found for query trajectory\n');
    query_segment_embeddings = [];
else
    % Store segment embeddings in cell array
    query_segment_embeddings = cell(num_query_segments, 1);
    
    for seg_idx = 1:num_query_segments
        segment_id = query_seg_result.segment_id{seg_idx};
        
        query_segment_embeddings{seg_idx} = struct();
        query_segment_embeddings{seg_idx}.segment_id = segment_id;
        query_segment_embeddings{seg_idx}.bahn_id = query_seg_result.bahn_id{seg_idx};
        query_segment_embeddings{seg_idx}.position = parseEmbeddingString(query_seg_result.position_embedding{seg_idx});
        query_segment_embeddings{seg_idx}.joint = parseEmbeddingString(query_seg_result.joint_embedding{seg_idx});
        query_segment_embeddings{seg_idx}.orientation = parseEmbeddingString(query_seg_result.orientation_embedding{seg_idx});
        query_segment_embeddings{seg_idx}.velocity = parseEmbeddingString(query_seg_result.velocity_embedding{seg_idx});
        query_segment_embeddings{seg_idx}.metadata = parseEmbeddingString(query_seg_result.metadata_embedding{seg_idx});
    end
    
    fprintf('✓ Loaded %d query segment embeddings\n', num_query_segments);
    for seg_idx = 1:num_query_segments
        fprintf('  Segment %d: %s\n', seg_idx, query_segment_embeddings{seg_idx}.segment_id);
    end
end

fprintf('\n');

%% SECTION 4: TRAJECTORY-LEVEL SEARCH
%  ========================================================================

fprintf('=== TRAJECTORY-LEVEL SEARCH ===\n');

% ========================================================================
% 4.1: pgvector HNSW Retrieval (Top-K per modality)
% ========================================================================

fprintf('\n--- Stage 1: pgvector HNSW Retrieval ---\n');

tic;

% Initialize rankings for each modality
rankings = struct();

% ------------------------------------------------------------------------
% Position Embedding Search
% ------------------------------------------------------------------------

if weight_position > 0 && has_position
    fprintf('  Searching Position embeddings...\n');
    
    % Convert embedding to PostgreSQL vector string format
    position_vector_str = sprintf('[%s]', strjoin(string(query_position_emb), ','));
    
    position_query = sprintf(...
        ['SELECT ' ...
         'e.segment_id as bahn_id, ' ...
         'm.bahn_id as parent_bahn_id, ' ...
         'e.position_embedding <=> ''%s''::vector as distance ' ...
         'FROM %s.%s.bahn_embeddings e ' ...
         'JOIN %s.%s.bahn_metadata m ON e.segment_id = m.segment_id ' ...
         'WHERE e.segment_id != ''%s'' ' ...
         'AND e.segment_id = e.bahn_id ' ...
         'AND e.position_embedding IS NOT NULL ' ...
         'ORDER BY distance ' ...
         'LIMIT %d'], ...
        position_vector_str, ...
        db_name, schema, db_name, schema, ...
        query_bahn_id, pgvector_top_k);
    
    position_results = fetch(conn, position_query);
    
    % Rename column for consistency
    position_results.Properties.VariableNames{'bahn_id'} = 'bahn_id';
    
    % Create ranking table
    position_results.rank = (1:height(position_results))';
    rankings.position = position_results;
    
    fprintf('    ✓ Retrieved %d candidates (best: %s, dist: %.6f)\n', ...
        height(position_results), position_results.bahn_id{1}, position_results.distance(1));
end

% ------------------------------------------------------------------------
% Joint Embedding Search
% ------------------------------------------------------------------------

if weight_joint > 0 && has_joint
    fprintf('  Searching Joint embeddings...\n');
    
    joint_vector_str = sprintf('[%s]', strjoin(string(query_joint_emb), ','));
    
    joint_query = sprintf(...
        ['SELECT ' ...
         'e.segment_id as bahn_id, ' ...
         'm.bahn_id as parent_bahn_id, ' ...
         'e.joint_embedding <=> ''%s''::vector as distance ' ...
         'FROM %s.%s.bahn_embeddings e ' ...
         'JOIN %s.%s.bahn_metadata m ON e.segment_id = m.segment_id ' ...
         'WHERE e.segment_id != ''%s'' ' ...
         'AND e.segment_id = e.bahn_id ' ...
         'AND e.joint_embedding IS NOT NULL ' ...
         'ORDER BY distance ' ...
         'LIMIT %d'], ...
        joint_vector_str, ...
        db_name, schema, db_name, schema, ...
        query_bahn_id, pgvector_top_k);
    
    joint_results = fetch(conn, joint_query);
    joint_results.Properties.VariableNames{'bahn_id'} = 'bahn_id';
    joint_results.rank = (1:height(joint_results))';
    rankings.joint = joint_results;
    
    fprintf('    ✓ Retrieved %d candidates (best: %s, dist: %.6f)\n', ...
        height(joint_results), joint_results.bahn_id{1}, joint_results.distance(1));
end

% ------------------------------------------------------------------------
% Orientation Embedding Search
% ------------------------------------------------------------------------

if weight_orientation > 0 && has_orientation
    fprintf('  Searching Orientation embeddings...\n');
    
    orientation_vector_str = sprintf('[%s]', strjoin(string(query_orientation_emb), ','));
    
    orientation_query = sprintf(...
        ['SELECT ' ...
         'e.segment_id as bahn_id, ' ...
         'm.bahn_id as parent_bahn_id, ' ...
         'e.orientation_embedding <=> ''%s''::vector as distance ' ...
         'FROM %s.%s.bahn_embeddings e ' ...
         'JOIN %s.%s.bahn_metadata m ON e.segment_id = m.segment_id ' ...
         'WHERE e.segment_id != ''%s'' ' ...
         'AND e.segment_id = e.bahn_id ' ...
         'AND e.orientation_embedding IS NOT NULL ' ...
         'ORDER BY distance ' ...
         'LIMIT %d'], ...
        orientation_vector_str, ...
        db_name, schema, db_name, schema, ...
        query_bahn_id, pgvector_top_k);
    
    orientation_results = fetch(conn, orientation_query);
    orientation_results.Properties.VariableNames{'bahn_id'} = 'bahn_id';
    orientation_results.rank = (1:height(orientation_results))';
    rankings.orientation = orientation_results;
    
    fprintf('    ✓ Retrieved %d candidates (best: %s, dist: %.6f)\n', ...
        height(orientation_results), orientation_results.bahn_id{1}, orientation_results.distance(1));
end

% ------------------------------------------------------------------------
% Velocity Embedding Search
% ------------------------------------------------------------------------

if weight_velocity > 0 && has_velocity
    fprintf('  Searching Velocity embeddings...\n');
    
    velocity_vector_str = sprintf('[%s]', strjoin(string(query_velocity_emb), ','));
    
    velocity_query = sprintf(...
        ['SELECT ' ...
         'e.segment_id as bahn_id, ' ...
         'm.bahn_id as parent_bahn_id, ' ...
         'e.velocity_embedding <=> ''%s''::vector as distance ' ...
         'FROM %s.%s.bahn_embeddings e ' ...
         'JOIN %s.%s.bahn_metadata m ON e.segment_id = m.segment_id ' ...
         'WHERE e.segment_id != ''%s'' ' ...
         'AND e.segment_id = e.bahn_id ' ...
         'AND e.velocity_embedding IS NOT NULL ' ...
         'ORDER BY distance ' ...
         'LIMIT %d'], ...
        velocity_vector_str, ...
        db_name, schema, db_name, schema, ...
        query_bahn_id, pgvector_top_k);
    
    velocity_results = fetch(conn, velocity_query);
    velocity_results.Properties.VariableNames{'bahn_id'} = 'bahn_id';
    velocity_results.rank = (1:height(velocity_results))';
    rankings.velocity = velocity_results;
    
    fprintf('    ✓ Retrieved %d candidates (best: %s, dist: %.6f)\n', ...
        height(velocity_results), velocity_results.bahn_id{1}, velocity_results.distance(1));
end

% ------------------------------------------------------------------------
% Metadata Embedding Search
% ------------------------------------------------------------------------

if weight_metadata > 0 && has_metadata
    fprintf('  Searching Metadata embeddings...\n');
    
    metadata_vector_str = sprintf('[%s]', strjoin(string(query_metadata_emb), ','));
    
    metadata_query = sprintf(...
        ['SELECT ' ...
         'e.segment_id as bahn_id, ' ...
         'm.bahn_id as parent_bahn_id, ' ...
         'e.metadata_embedding <=> ''%s''::vector as distance ' ...
         'FROM %s.%s.bahn_embeddings e ' ...
         'JOIN %s.%s.bahn_metadata m ON e.segment_id = m.segment_id ' ...
         'WHERE e.segment_id != ''%s'' ' ...
         'AND e.segment_id = e.bahn_id ' ...
         'AND e.metadata_embedding IS NOT NULL ' ...
         'ORDER BY distance ' ...
         'LIMIT %d'], ...
        metadata_vector_str, ...
        db_name, schema, db_name, schema, ...
        query_bahn_id, pgvector_top_k);
    
    metadata_results = fetch(conn, metadata_query);
    metadata_results.Properties.VariableNames{'bahn_id'} = 'bahn_id';
    metadata_results.rank = (1:height(metadata_results))';
    rankings.metadata = metadata_results;
    
    fprintf('    ✓ Retrieved %d candidates (best: %s, dist: %.6f)\n', ...
        height(metadata_results), metadata_results.bahn_id{1}, metadata_results.distance(1));
end

pgvector_time = toc;
fprintf('\n  ⏱ pgvector retrieval time: %.4f seconds\n', pgvector_time);

% ========================================================================
% 4.2: RRF Fusion
% ========================================================================

fprintf('\n--- Stage 2: RRF Fusion ---\n');

% Check if we have any rankings
if isempty(fieldnames(rankings))
    error('No embeddings available for search! Check weights and embedding availability.');
end

% Prepare weights struct (only include active modalities)
weights = struct();
if isfield(rankings, 'position'), weights.position = weight_position; end
if isfield(rankings, 'joint'), weights.joint = weight_joint; end
if isfield(rankings, 'orientation'), weights.orientation = weight_orientation; end
if isfield(rankings, 'velocity'), weights.velocity = weight_velocity; end
if isfield(rankings, 'metadata'), weights.metadata = weight_metadata; end

% Perform RRF fusion
fused_ranking = fuseRankingsRRF(rankings, weights, rrf_k, 'bahn_id');

fprintf('✓ RRF fusion completed: %d unique candidates\n', height(fused_ranking));
fprintf('  Top-3:\n');
for i = 1:min(3, height(fused_ranking))
    fprintf('    %d. %s (RRF score: %.6f)\n', ...
        i, fused_ranking.bahn_id{i}, fused_ranking.rrf_score(i));
end

% Take top-K for next stage
fused_top_pgvector = fused_ranking(1:min(pgvector_top_k, height(fused_ranking)), :);

% ========================================================================
% 4.3: Cascading Lower Bounds Pre-filtering
% ========================================================================

fprintf('\n--- Stage 3: Cascading Lower Bounds Pre-filtering ---\n');

tic;

% Load query trajectory data
if strcmp(dtw_mode, 'position')
    query_data = loadTrajectoryPosition(conn, schema, query_bahn_id);
else
    query_data = loadTrajectoryJointStates(conn, schema, query_bahn_id);
end

% Normalize to start point
query_start = query_data(1, :);
query_data = query_data - query_start;

% ------------------------------------------------------------------------
% ⭐ Batch Load ALL candidates at once (1 query instead of 250!)
% ------------------------------------------------------------------------

fprintf('  Loading all %d candidates in batch...\n', height(fused_top_pgvector));

candidate_ids = fused_top_pgvector.bahn_id;
candidate_trajectories = batchLoadTrajectoriesChunked(conn, schema, candidate_ids, dtw_mode, chunk_size);

fprintf('    ✓ Batch loading completed\n');

% ------------------------------------------------------------------------
% Step 1: LB_Kim (Ultra-fast, no more DB queries!)
% ------------------------------------------------------------------------

fprintf('\n  Step 1: LB_Kim screening...\n');

kim_distances = zeros(height(fused_top_pgvector), 1);

for i = 1:height(fused_top_pgvector)
    % ✅ Data already in memory!
    candidate_data = candidate_trajectories{i};
    
    % Normalize to start point
    if ~isempty(candidate_data)
        candidate_data = candidate_data - candidate_data(1, :);
    end
    
    % Compute LB_Kim (ultra-fast)
    kim_distances(i) = LB_Kim(query_data, candidate_data, dtw_mode);
end

fprintf('    ✓ LB_Kim computed for %d candidates\n', height(fused_top_pgvector));

% Sort by Kim distance and keep top 2× lb_keogh_keep
[~, kim_order] = sort(kim_distances);
kim_keep = min(lb_keogh_keep * 2, height(fused_top_pgvector));
candidates_after_kim = fused_top_pgvector(kim_order(1:kim_keep), :);

fprintf('    ✓ LB_Kim: %d → %d candidates (%.1f%% pruned)\n', ...
    height(fused_top_pgvector), height(candidates_after_kim), ...
    100 * (1 - height(candidates_after_kim) / height(fused_top_pgvector)));

% Keep corresponding trajectory data
candidate_trajectories_after_kim = candidate_trajectories(kim_order(1:kim_keep));

% ------------------------------------------------------------------------
% Step 2: LB_Keogh (Precise filtering, data already in memory!)
% ------------------------------------------------------------------------

fprintf('\n  Step 2: LB_Keogh refinement...\n');

keogh_distances = zeros(height(candidates_after_kim), 1);

for i = 1:height(candidates_after_kim)
    % ✅ Data already in memory!
    candidate_data = candidate_trajectories_after_kim{i};
    
    % Normalize to start point
    if ~isempty(candidate_data)
        candidate_data = candidate_data - candidate_data(1, :);
    end
    
    % Compute LB_Keogh
    keogh_distances(i) = LB_Keogh(query_data, candidate_data, cdtw_window, dtw_mode, ...
                                   use_rotation_alignment, normalize_dtw);
end

% Add Keogh distances and sort
candidates_after_kim.lb_distance = keogh_distances;
candidates_after_kim = sortrows(candidates_after_kim, 'lb_distance', 'ascend');
lb_top_trajs = candidates_after_kim(1:min(lb_keogh_keep, height(candidates_after_kim)), :);

lb_keogh_time = toc;

fprintf('    ✓ LB_Keogh: %d → %d candidates\n', ...
    height(candidates_after_kim), height(lb_top_trajs));
fprintf('\n  ⏱ Total cascading LB time: %.4f seconds\n', lb_keogh_time);

% ========================================================================
% 4.4: DTW Refinement
% ========================================================================

fprintf('\n--- Stage 4: DTW Refinement ---\n');

tic;

% ------------------------------------------------------------------------
% ⭐ Batch Load DTW candidates (if not already loaded)
% ------------------------------------------------------------------------

% Check if we need to reload data (if Kim pruned some candidates)
% Get the indices of lb_top_trajs in the original candidates_after_kim
[~, dtw_indices] = ismember(lb_top_trajs.bahn_id, candidates_after_kim.bahn_id);

% Extract corresponding trajectory data
candidate_trajectories_for_dtw = candidate_trajectories_after_kim(dtw_indices);

fprintf('  Computing DTW distances for %d candidates (data already in memory)...\n', ...
    height(lb_top_trajs));

% ------------------------------------------------------------------------
% DTW Computation (no DB queries!)
% ------------------------------------------------------------------------

dtw_distances = zeros(height(lb_top_trajs), 1);

for i = 1:height(lb_top_trajs)
    % ✅ Data already in memory!
    candidate_data = candidate_trajectories_for_dtw{i};
    
    % Normalize to start point
    if ~isempty(candidate_data)
        candidate_data = candidate_data - candidate_data(1, :);
    end
    
    % Early abandoning: use k-th best distance so far
    if i <= dtw_top_n
        best_so_far = inf;
    else
        sorted_dists = sort(dtw_distances(1:i-1));
        best_so_far = sorted_dists(dtw_top_n);
    end
    
    % Compute DTW with early abandoning
    dtw_distances(i) = cDTW(query_data, candidate_data, dtw_mode, cdtw_window, ...
                            best_so_far, use_rotation_alignment, normalize_dtw);
    
    if mod(i, 20) == 0
        fprintf('    Progress: %d/%d\n', i, height(lb_top_trajs));
    end
end

% Add DTW distances to table
lb_top_trajs.dtw_distance = dtw_distances;

% Sort by DTW distance and keep top-N
lb_top_trajs = sortrows(lb_top_trajs, 'dtw_distance', 'ascend');
final_top_k = lb_top_trajs(1:min(dtw_top_n, height(lb_top_trajs)), :);

dtw_time = toc;

fprintf('  ✓ DTW completed: Final top-%d results\n', height(final_top_k));
fprintf('  ⏱ DTW time: %.4f seconds\n', dtw_time);

% ========================================================================
% Summary
% ========================================================================

fprintf('\n--- Trajectory-Level Summary ---\n');
fprintf('  Total time: %.4f seconds\n', pgvector_time + lb_keogh_time + dtw_time);
fprintf('  Pipeline: %d → %d → %d → %d → %d\n', ...
    pgvector_top_k, height(fused_top_pgvector), height(candidates_after_kim), ...
    height(lb_top_trajs), height(final_top_k));
fprintf('\n  Top-10 Results:\n');
for i = 1:min(10, height(final_top_k))
    fprintf('    %d. %s (DTW: %.2f, RRF: %.6f)\n', ...
        i, final_top_k.bahn_id{i}, final_top_k.dtw_distance(i), final_top_k.rrf_score(i));
end

fprintf('\n');

%% SECTION 5: SEGMENT-LEVEL SEARCH
%  ========================================================================

if num_query_segments == 0
    fprintf('=== SEGMENT-LEVEL SEARCH ===\n');
    fprintf('⚠ No segments to process - skipping\n\n');
    segment_results = [];
else
    fprintf('=== SEGMENT-LEVEL SEARCH ===\n');
    fprintf('Processing %d query segments...\n\n', num_query_segments);
    
    % Initialize storage for segment results
    segment_results = cell(num_query_segments, 1);
    
    % ========================================================================
    % Process each query segment independently
    % ========================================================================
    
    for seg_idx = 1:num_query_segments
        query_segment_id = query_segment_embeddings{seg_idx}.segment_id;
        
        fprintf('--- Query Segment %d/%d: %s ---\n', seg_idx, num_query_segments, query_segment_id);
        
        % Get segment embeddings
        seg_position_emb = query_segment_embeddings{seg_idx}.position;
        seg_joint_emb = query_segment_embeddings{seg_idx}.joint;
        seg_orientation_emb = query_segment_embeddings{seg_idx}.orientation;
        seg_velocity_emb = query_segment_embeddings{seg_idx}.velocity;
        seg_metadata_emb = query_segment_embeddings{seg_idx}.metadata;
        
        % Check availability
        seg_has_position = ~isempty(seg_position_emb) && length(seg_position_emb) > 0;
        seg_has_joint = ~isempty(seg_joint_emb) && length(seg_joint_emb) > 0;
        seg_has_orientation = ~isempty(seg_orientation_emb) && length(seg_orientation_emb) > 0;
        seg_has_velocity = ~isempty(seg_velocity_emb) && length(seg_velocity_emb) > 0;
        seg_has_metadata = ~isempty(seg_metadata_emb) && length(seg_metadata_emb) > 0;
        
        % ====================================================================
        % Stage 1: pgvector HNSW Retrieval
        % ====================================================================
        
        fprintf('\n  Stage 1: pgvector HNSW Retrieval\n');
        
        tic;
        
        seg_rankings = struct();
        
        % --------------------------------------------------------------------
        % Position Embedding Search
        % --------------------------------------------------------------------
        
        if weight_position > 0 && seg_has_position
            fprintf('    Searching Position embeddings...\n');
            
            position_vector_str = sprintf('[%s]', strjoin(string(seg_position_emb), ','));
            
            seg_position_query = sprintf(...
                ['SELECT ' ...
                 'e.segment_id, ' ...
                 'm.bahn_id, ' ...
                 'e.position_embedding <=> ''%s''::vector as distance ' ...
                 'FROM %s.%s.bahn_embeddings e ' ...
                 'JOIN %s.%s.bahn_metadata m ON e.segment_id = m.segment_id ' ...
                 'WHERE e.segment_id != ''%s'' ' ...
                 'AND e.segment_id != e.bahn_id ' ...
                 'AND e.position_embedding IS NOT NULL ' ...
                 'ORDER BY distance ' ...
                 'LIMIT %d'], ...
                position_vector_str, ...
                db_name, schema, db_name, schema, ...
                query_segment_id, pgvector_top_k);
            
            seg_position_results = fetch(conn, seg_position_query);
            
            if ~isempty(seg_position_results)
                seg_position_results.rank = (1:height(seg_position_results))';
                seg_rankings.position = seg_position_results;
                fprintf('      ✓ Retrieved %d candidates\n', height(seg_position_results));
            end
        end
        
        % --------------------------------------------------------------------
        % Joint Embedding Search
        % --------------------------------------------------------------------
        
        if weight_joint > 0 && seg_has_joint
            fprintf('    Searching Joint embeddings...\n');
            
            joint_vector_str = sprintf('[%s]', strjoin(string(seg_joint_emb), ','));
            
            seg_joint_query = sprintf(...
                ['SELECT ' ...
                 'e.segment_id, ' ...
                 'm.bahn_id, ' ...
                 'e.joint_embedding <=> ''%s''::vector as distance ' ...
                 'FROM %s.%s.bahn_embeddings e ' ...
                 'JOIN %s.%s.bahn_metadata m ON e.segment_id = m.segment_id ' ...
                 'WHERE e.segment_id != ''%s'' ' ...
                 'AND e.segment_id != e.bahn_id ' ...
                 'AND e.joint_embedding IS NOT NULL ' ...
                 'ORDER BY distance ' ...
                 'LIMIT %d'], ...
                joint_vector_str, ...
                db_name, schema, db_name, schema, ...
                query_segment_id, pgvector_top_k);
            
            seg_joint_results = fetch(conn, seg_joint_query);
            
            if ~isempty(seg_joint_results)
                seg_joint_results.rank = (1:height(seg_joint_results))';
                seg_rankings.joint = seg_joint_results;
                fprintf('      ✓ Retrieved %d candidates\n', height(seg_joint_results));
            end
        end
        
        % --------------------------------------------------------------------
        % Orientation Embedding Search
        % --------------------------------------------------------------------
        
        if weight_orientation > 0 && seg_has_orientation
            fprintf('    Searching Orientation embeddings...\n');
            
            orientation_vector_str = sprintf('[%s]', strjoin(string(seg_orientation_emb), ','));
            
            seg_orientation_query = sprintf(...
                ['SELECT ' ...
                 'e.segment_id, ' ...
                 'm.bahn_id, ' ...
                 'e.orientation_embedding <=> ''%s''::vector as distance ' ...
                 'FROM %s.%s.bahn_embeddings e ' ...
                 'JOIN %s.%s.bahn_metadata m ON e.segment_id = m.segment_id ' ...
                 'WHERE e.segment_id != ''%s'' ' ...
                 'AND e.segment_id != e.bahn_id ' ...
                 'AND e.orientation_embedding IS NOT NULL ' ...
                 'ORDER BY distance ' ...
                 'LIMIT %d'], ...
                orientation_vector_str, ...
                db_name, schema, db_name, schema, ...
                query_segment_id, pgvector_top_k);
            
            seg_orientation_results = fetch(conn, seg_orientation_query);
            
            if ~isempty(seg_orientation_results)
                seg_orientation_results.rank = (1:height(seg_orientation_results))';
                seg_rankings.orientation = seg_orientation_results;
                fprintf('      ✓ Retrieved %d candidates\n', height(seg_orientation_results));
            end
        end
        
        % --------------------------------------------------------------------
        % Velocity Embedding Search
        % --------------------------------------------------------------------
        
        if weight_velocity > 0 && seg_has_velocity
            fprintf('    Searching Velocity embeddings...\n');
            
            velocity_vector_str = sprintf('[%s]', strjoin(string(seg_velocity_emb), ','));
            
            seg_velocity_query = sprintf(...
                ['SELECT ' ...
                 'e.segment_id, ' ...
                 'm.bahn_id, ' ...
                 'e.velocity_embedding <=> ''%s''::vector as distance ' ...
                 'FROM %s.%s.bahn_embeddings e ' ...
                 'JOIN %s.%s.bahn_metadata m ON e.segment_id = m.segment_id ' ...
                 'WHERE e.segment_id != ''%s'' ' ...
                 'AND e.segment_id != e.bahn_id ' ...
                 'AND e.velocity_embedding IS NOT NULL ' ...
                 'ORDER BY distance ' ...
                 'LIMIT %d'], ...
                velocity_vector_str, ...
                db_name, schema, db_name, schema, ...
                query_segment_id, pgvector_top_k);
            
            seg_velocity_results = fetch(conn, seg_velocity_query);
            
            if ~isempty(seg_velocity_results)
                seg_velocity_results.rank = (1:height(seg_velocity_results))';
                seg_rankings.velocity = seg_velocity_results;
                fprintf('      ✓ Retrieved %d candidates\n', height(seg_velocity_results));
            end
        end
        
        % --------------------------------------------------------------------
        % Metadata Embedding Search
        % --------------------------------------------------------------------
        
        if weight_metadata > 0 && seg_has_metadata
            fprintf('    Searching Metadata embeddings...\n');
            
            metadata_vector_str = sprintf('[%s]', strjoin(string(seg_metadata_emb), ','));
            
            seg_metadata_query = sprintf(...
                ['SELECT ' ...
                 'e.segment_id, ' ...
                 'm.bahn_id, ' ...
                 'e.metadata_embedding <=> ''%s''::vector as distance ' ...
                 'FROM %s.%s.bahn_embeddings e ' ...
                 'JOIN %s.%s.bahn_metadata m ON e.segment_id = m.segment_id ' ...
                 'WHERE e.segment_id != ''%s'' ' ...
                 'AND e.segment_id != e.bahn_id ' ...
                 'AND e.metadata_embedding IS NOT NULL ' ...
                 'ORDER BY distance ' ...
                 'LIMIT %d'], ...
                metadata_vector_str, ...
                db_name, schema, db_name, schema, ...
                query_segment_id, pgvector_top_k);
            
            seg_metadata_results = fetch(conn, seg_metadata_query);
            
            if ~isempty(seg_metadata_results)
                seg_metadata_results.rank = (1:height(seg_metadata_results))';
                seg_rankings.metadata = seg_metadata_results;
                fprintf('      ✓ Retrieved %d candidates\n', height(seg_metadata_results));
            end
        end
        
        seg_pgvector_time = toc;
        fprintf('    ⏱ pgvector time: %.4f seconds\n', seg_pgvector_time);
        
        % ====================================================================
        % Stage 2: RRF Fusion
        % ====================================================================
        
        fprintf('\n  Stage 2: RRF Fusion\n');
        
        % Check if we have any rankings
        if isempty(fieldnames(seg_rankings))
            fprintf('  ⚠ No embeddings/results found for segment %s - skipping.\n', query_segment_id);
            segment_results{seg_idx} = [];
            continue;
        end
        
        % Prepare weights struct
        seg_weights = struct();
        if isfield(seg_rankings, 'position'), seg_weights.position = weight_position; end
        if isfield(seg_rankings, 'joint'), seg_weights.joint = weight_joint; end
        if isfield(seg_rankings, 'orientation'), seg_weights.orientation = weight_orientation; end
        if isfield(seg_rankings, 'velocity'), seg_weights.velocity = weight_velocity; end
        if isfield(seg_rankings, 'metadata'), seg_weights.metadata = weight_metadata; end
        
        % Perform RRF fusion
        seg_fused_ranking = fuseRankingsRRF(seg_rankings, seg_weights, rrf_k, 'segment_id');
        
        fprintf('    ✓ RRF fusion: %d unique candidates\n', height(seg_fused_ranking));
        
        % Take top-K for next stage
        seg_fused_top = seg_fused_ranking(1:min(pgvector_top_k, height(seg_fused_ranking)), :);
        
        % ====================================================================
        % Stage 3: Cascading Lower Bounds Pre-filtering
        % ====================================================================
        
        fprintf('\n  Stage 3: Cascading Lower Bounds Pre-filtering\n');
        
        tic;
        
        % Load query segment data
        if strcmp(dtw_mode, 'position')
            query_seg_data = loadSegmentData(conn, schema, query_segment_id, 'position');
        else
            query_seg_data = loadSegmentData(conn, schema, query_segment_id, 'joint_states');
        end
        
        % Normalize to start point
        if ~isempty(query_seg_data)
            query_seg_start = query_seg_data(1, :);
            query_seg_data = query_seg_data - query_seg_start;
        end
        
        % --------------------------------------------------------------------
        % ⭐ Batch Load ALL segment candidates at once
        % --------------------------------------------------------------------
        
        fprintf('    Loading all %d segment candidates in batch...\n', height(seg_fused_top));
        
        seg_candidate_ids = seg_fused_top.segment_id;
        seg_candidate_segments = batchLoadSegmentsChunked(conn, schema, seg_candidate_ids, dtw_mode, chunk_size);
        
        fprintf('      ✓ Batch loading completed\n');
        
        % --------------------------------------------------------------------
        % Step 1: LB_Kim (Ultra-fast, no DB queries!)
        % --------------------------------------------------------------------
        
        fprintf('    Step 1: LB_Kim screening...\n');
        
        seg_kim_distances = zeros(height(seg_fused_top), 1);
        
        for i = 1:height(seg_fused_top)
            % ✅ Data already in memory!
            candidate_seg_data = seg_candidate_segments{i};
            
            % Normalize to start point
            if ~isempty(candidate_seg_data)
                candidate_seg_data = candidate_seg_data - candidate_seg_data(1, :);
            end
            
            % Compute LB_Kim
            seg_kim_distances(i) = LB_Kim(query_seg_data, candidate_seg_data, dtw_mode);
        end
        
        fprintf('      ✓ LB_Kim computed for %d candidates\n', height(seg_fused_top));
        
        % Sort and keep top 2× lb_keogh_keep
        [~, seg_kim_order] = sort(seg_kim_distances);
        seg_kim_keep = min(lb_keogh_keep * 2, height(seg_fused_top));
        seg_candidates_after_kim = seg_fused_top(seg_kim_order(1:seg_kim_keep), :);
        
        fprintf('      ✓ LB_Kim: %d → %d candidates (%.1f%% pruned)\n', ...
            height(seg_fused_top), height(seg_candidates_after_kim), ...
            100 * (1 - height(seg_candidates_after_kim) / height(seg_fused_top)));
        
        % Keep corresponding segment data
        seg_candidate_segments_after_kim = seg_candidate_segments(seg_kim_order(1:seg_kim_keep));
        
        % --------------------------------------------------------------------
        % Step 2: LB_Keogh (Precise, data already in memory!)
        % --------------------------------------------------------------------
        
        fprintf('    Step 2: LB_Keogh refinement...\n');
        
        seg_keogh_distances = zeros(height(seg_candidates_after_kim), 1);
        
        for i = 1:height(seg_candidates_after_kim)
            % ✅ Data already in memory!
            candidate_seg_data = seg_candidate_segments_after_kim{i};
            
            % Normalize to start point
            if ~isempty(candidate_seg_data)
                candidate_seg_data = candidate_seg_data - candidate_seg_data(1, :);
            end
            
            % Compute LB_Keogh
            seg_keogh_distances(i) = LB_Keogh(query_seg_data, candidate_seg_data, ...
                                              cdtw_window, dtw_mode, use_rotation_alignment, normalize_dtw);
        end
        
        % Add Keogh distances and sort
        seg_candidates_after_kim.lb_distance = seg_keogh_distances;
        seg_candidates_after_kim = sortrows(seg_candidates_after_kim, 'lb_distance', 'ascend');
        seg_lb_top = seg_candidates_after_kim(1:min(lb_keogh_keep, height(seg_candidates_after_kim)), :);
        
        seg_lb_time = toc;
        
        fprintf('      ✓ LB_Keogh: %d → %d candidates\n', ...
            height(seg_candidates_after_kim), height(seg_lb_top));
        fprintf('    ⏱ Cascading LB time: %.4f seconds\n', seg_lb_time);
        
        % ====================================================================
        % Stage 4: DTW Refinement
        % ====================================================================
        
        fprintf('\n  Stage 4: DTW Refinement\n');
        
        tic;
        
        % --------------------------------------------------------------------
        % ⭐ Get DTW candidates from memory (no DB queries!)
        % --------------------------------------------------------------------
        
        % Map lb_top segments to their data in seg_candidates_after_kim
        [~, seg_dtw_indices] = ismember(seg_lb_top.segment_id, seg_candidates_after_kim.segment_id);
        seg_candidate_segments_for_dtw = seg_candidate_segments_after_kim(seg_dtw_indices);
        
        fprintf('    Computing DTW for %d candidates (data already in memory)...\n', ...
            height(seg_lb_top));
        
        seg_dtw_distances = zeros(height(seg_lb_top), 1);
        
        for i = 1:height(seg_lb_top)
            % ✅ Data already in memory!
            candidate_seg_data = seg_candidate_segments_for_dtw{i};
            
            % Normalize to start point
            if ~isempty(candidate_seg_data)
                candidate_seg_data = candidate_seg_data - candidate_seg_data(1, :);
            end
            
            % Early abandoning
            if i <= dtw_top_n
                best_so_far = inf;
            else
                sorted_seg_dists = sort(seg_dtw_distances(1:i-1));
                best_so_far = sorted_seg_dists(dtw_top_n);
            end
            
            % Compute DTW
            seg_dtw_distances(i) = cDTW(query_seg_data, candidate_seg_data, dtw_mode, ...
                                        cdtw_window, best_so_far, use_rotation_alignment, normalize_dtw);
            
            if mod(i, 20) == 0
                fprintf('      Progress: %d/%d\n', i, height(seg_lb_top));
            end
        end
        
        % Add DTW distances and sort
        seg_lb_top.dtw_distance = seg_dtw_distances;
        seg_lb_top = sortrows(seg_lb_top, 'dtw_distance', 'ascend');
        seg_final_top = seg_lb_top(1:min(dtw_top_n, height(seg_lb_top)), :);
        
        seg_dtw_time = toc;
        
        fprintf('    ✓ DTW: Final top-%d results\n', height(seg_final_top));
        fprintf('    ⏱ DTW time: %.4f seconds\n', seg_dtw_time);
        
        % ====================================================================
        % Store results for this segment
        % ====================================================================
        
        segment_results{seg_idx} = seg_final_top;
        
        % ====================================================================
        % Segment Summary
        % ====================================================================
        
        fprintf('\n  Segment %d Summary:\n', seg_idx);
        fprintf('    Total time: %.4f seconds\n', ...
            seg_pgvector_time + seg_lb_time + seg_dtw_time);
        fprintf('    Pipeline: %d → %d → %d → %d → %d\n', ...
            pgvector_top_k, height(seg_fused_top), height(seg_candidates_after_kim), ...
            height(seg_lb_top), height(seg_final_top));
        
        if ~isempty(seg_final_top)
            fprintf('\n    Top-3 Results:\n');
            for i = 1:min(3, height(seg_final_top))
                fprintf('      %d. %s (DTW: %.2f)\n', ...
                    i, seg_final_top.segment_id{i}, seg_final_top.dtw_distance(i));
            end
        end
        
        fprintf('\n');
    end
    
    fprintf('✓ Segment-level search completed for all %d segments\n\n', num_query_segments);
end

%% SECTION 6: RESULTS & VISUALIZATION
%  ========================================================================

fprintf('=== RESULTS SUMMARY ===\n\n');

% ========================================================================
% 6.1: Trajectory-Level Results
% ========================================================================

fprintf('--- TRAJECTORY-LEVEL RESULTS ---\n\n');

fprintf('Query: %s\n', query_bahn_id);
fprintf('Pipeline: pgvector (%d) → LB_Keogh (%d) → DTW (%d) → Final (%d)\n\n', ...
    pgvector_top_k, lb_keogh_keep, dtw_top_n, height(final_top_k));

fprintf('Top-10 Similar Trajectories:\n');
fprintf('Rank | Bahn ID      | DTW Distance | RRF Score  | LB Distance\n');
fprintf('-----|--------------|--------------|------------|-------------\n');

for i = 1:min(10, height(final_top_k))
    fprintf('%4d | %12s | %12.2f | %10.6f | %11.2f\n', ...
        i, ...
        final_top_k.bahn_id{i}, ...
        final_top_k.dtw_distance(i), ...
        final_top_k.rrf_score(i), ...
        final_top_k.lb_distance(i));
end

fprintf('\n');

% ========================================================================
% 6.2: Segment-Level Results
% ========================================================================
if ~isempty(segment_results)
    fprintf('--- SEGMENT-LEVEL RESULTS ---\n\n');
    
    for seg_idx = 1:num_query_segments
        % ✅ FIX: Zugriff auf die ID aus query_segment_embeddings (cell array)
        query_segment_id = query_segment_embeddings{seg_idx}.segment_id;
        
        seg_final = segment_results{seg_idx};
        
        fprintf('Query Segment %d: %s\n', seg_idx, query_segment_id);
        
        if isempty(seg_final)
            fprintf('  No similar segments found.\n\n');
            continue;
        end
        
        fprintf('Top-5 Similar Segments:\n');
        fprintf('Rank | Segment ID        | Parent Bahn  | DTW Distance | RRF Score\n');
        fprintf('-----|-------------------|--------------|--------------|----------\n');
        
        for i = 1:min(5, height(seg_final))
            fprintf('%4d | %17s | %12s | %12.2f | %10.6f\n', ...
                i, ...
                seg_final.segment_id{i}, ...
                seg_final.bahn_id{i}, ...
                seg_final.dtw_distance(i), ...
                seg_final.rrf_score(i));
        end
        
        fprintf('\n');
    end
end

% ========================================================================
% 6.3: Timing Summary
% ========================================================================
fprintf('--- PERFORMANCE SUMMARY ---\n\n');
fprintf('Trajectory-Level:\n');
fprintf('  pgvector:  %.4f seconds\n', pgvector_time);
fprintf('  LB_Keogh:  %.4f seconds\n', lb_keogh_time);
fprintf('  DTW:       %.4f seconds\n', dtw_time);
fprintf('  Total:     %.4f seconds\n', pgvector_time + lb_keogh_time + dtw_time);
fprintf('\n');

if ~isempty(segment_results)
    fprintf('Segment-Level:\n');
    fprintf('  Processed %d query segments.\n', num_query_segments);
    fprintf('  (See console output above for detailed per-segment timing)\n');
    fprintf('\n');
end

% ========================================================================
% 6.4: Visualization - Trajectory Level
% ========================================================================

fprintf('--- CREATING VISUALIZATIONS ---\n\n');

% Limit visualization to top-10 for clarity
viz_k = min(10, height(final_top_k));

fprintf('Creating trajectory-level visualization (Top-%d)...\n', viz_k);

% Load query trajectory position data for visualization
query_vis_data = loadTrajectoryPosition(conn, schema, query_bahn_id);
query_vis_start = query_vis_data(1, :);
query_vis_data = query_vis_data - query_vis_start;

% Create figure
fig_traj = figure('Position', [100, 100, 1600, 700], ...
                  'Name', 'Trajectory-Level Results', ...
                  'NumberTitle', 'off');

% ------------------------------------------------------------------------
% Left plot: Top-5 by RRF (Embedding-based)
% ------------------------------------------------------------------------

subplot(1, 2, 1);
hold on;

% Plot query (thick black line)
plot3(query_vis_data(:,1), query_vis_data(:,2), query_vis_data(:,3), ...
    'k-', 'LineWidth', 3, 'DisplayName', sprintf('Query: %s', query_bahn_id));

% Plot top-5 by RRF score
colors_rrf = autumn(5);
for i = 1:min(5, height(final_top_k))
    candidate_id = final_top_k.bahn_id{i};
    candidate_vis_data = loadTrajectoryPosition(conn, schema, candidate_id);
    candidate_vis_data = candidate_vis_data - candidate_vis_data(1, :);
    
    plot3(candidate_vis_data(:,1), candidate_vis_data(:,2), candidate_vis_data(:,3), ...
        '-', 'Color', colors_rrf(i,:), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('#%d: %s (RRF: %.4f)', i, candidate_id, final_top_k.rrf_score(i)));
end

grid on;
axis equal;
view(3);
xlabel('X (mm)', 'FontSize', 11);
ylabel('Y (mm)', 'FontSize', 11);
zlabel('Z (mm)', 'FontSize', 11);
title('Top-5 by RRF Score (Embedding-based)', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'bestoutside', 'FontSize', 9);

% ------------------------------------------------------------------------
% Right plot: Top-5 by DTW (Final ranking)
% ------------------------------------------------------------------------

subplot(1, 2, 2);
hold on;

% Plot query (thick black line)
plot3(query_vis_data(:,1), query_vis_data(:,2), query_vis_data(:,3), ...
    'k-', 'LineWidth', 3, 'DisplayName', sprintf('Query: %s', query_bahn_id));

% Plot top-5 by DTW distance
colors_dtw = lines(5);
for i = 1:min(5, height(final_top_k))
    candidate_id = final_top_k.bahn_id{i};
    candidate_vis_data = loadTrajectoryPosition(conn, schema, candidate_id);
    candidate_vis_data = candidate_vis_data - candidate_vis_data(1, :);
    
    plot3(candidate_vis_data(:,1), candidate_vis_data(:,2), candidate_vis_data(:,3), ...
        '-', 'Color', colors_dtw(i,:), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('#%d: %s (DTW: %.1f)', i, candidate_id, final_top_k.dtw_distance(i)));
end

grid on;
axis equal;
view(3);
xlabel('X (mm)', 'FontSize', 11);
ylabel('Y (mm)', 'FontSize', 11);
zlabel('Z (mm)', 'FontSize', 11);
title('Top-5 by DTW Distance (Final Ranking)', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'bestoutside', 'FontSize', 9);

% Super title
sgtitle(sprintf('Query: %s | Pipeline: pgvector → LB_Keogh → DTW', query_bahn_id), ...
    'FontSize', 14, 'FontWeight', 'bold');

drawnow;

fprintf('  ✓ Trajectory-level visualization created\n');

% ========================================================================
% 6.5: Visualization - Segment Level
% ========================================================================
if ~isempty(segment_results)
    fprintf('\nCreating segment-level visualizations...\n');
    
    for seg_idx = 1:num_query_segments
        % UPDATED: Zugriff über die Tabelle statt das alte struct array
        query_segment_id = query_segment_embeddings{seg_idx}.segment_id;
        
        seg_final = segment_results{seg_idx};
        
        if isempty(seg_final)
            fprintf('  Segment %d: No results - skipping\n', seg_idx);
            continue;
        end
        
        fprintf('  Creating visualization for Segment %d/%d (%s)...\n', ...
            seg_idx, num_query_segments, query_segment_id);
        
        % Load query segment position data
        query_seg_vis_data = loadSegmentData(conn, schema, query_segment_id, 'position');
        if ~isempty(query_seg_vis_data)
            query_seg_vis_start = query_seg_vis_data(1, :);
            query_seg_vis_data = query_seg_vis_data - query_seg_vis_start;
        end
        
        % Create figure for this segment
        fig_seg = figure('Position', [100, 100, 1600, 700], ...
                         'Name', sprintf('Segment %d Results', seg_idx), ...
                         'NumberTitle', 'off');
        
        % --------------------------------------------------------------------
        % Left plot: Top-5 by RRF (Embedding-based)
        % --------------------------------------------------------------------
        subplot(1, 2, 1);
        hold on;
        
        % Plot query segment
        plot3(query_seg_vis_data(:,1), query_seg_vis_data(:,2), query_seg_vis_data(:,3), ...
            'k-', 'LineWidth', 3, 'DisplayName', sprintf('Query: %s', query_segment_id));
        
        % IMPROVEMENT: Sort by RRF score specifically for this plot
        seg_viz_rrf = sortrows(seg_final, 'rrf_score', 'descend');
        colors_seg_rrf = autumn(5);
        num_viz = min(5, height(seg_viz_rrf));
        
        for i = 1:num_viz
            candidate_segment_id = seg_viz_rrf.segment_id{i};
            candidate_seg_vis_data = loadSegmentData(conn, schema, candidate_segment_id, 'position');
            
            if ~isempty(candidate_seg_vis_data)
                candidate_seg_vis_data = candidate_seg_vis_data - candidate_seg_vis_data(1, :);
                
                plot3(candidate_seg_vis_data(:,1), candidate_seg_vis_data(:,2), candidate_seg_vis_data(:,3), ...
                    '-', 'Color', colors_seg_rrf(i,:), 'LineWidth', 1.5, ...
                    'DisplayName', sprintf('#%d: %s (RRF: %.4f)', i, candidate_segment_id, seg_viz_rrf.rrf_score(i)));
            end
        end
        
        grid on;
        axis equal;
        view(3);
        xlabel('X (mm)', 'FontSize', 11);
        ylabel('Y (mm)', 'FontSize', 11);
        zlabel('Z (mm)', 'FontSize', 11);
        title(sprintf('Segment %d - Top-5 by RRF Score', seg_idx), 'FontSize', 13, 'FontWeight', 'bold');
        legend('Location', 'bestoutside', 'FontSize', 9);
        
        % --------------------------------------------------------------------
        % Right plot: Top-5 by DTW (Refined)
        % --------------------------------------------------------------------
        subplot(1, 2, 2);
        hold on;
        
        % Plot query segment
        plot3(query_seg_vis_data(:,1), query_seg_vis_data(:,2), query_seg_vis_data(:,3), ...
            'k-', 'LineWidth', 3, 'DisplayName', sprintf('Query: %s', query_segment_id));
        
        % Explicitly sort by DTW (should already be sorted, but ensures correctness)
        seg_viz_dtw = sortrows(seg_final, 'dtw_distance', 'ascend');
        colors_seg_dtw = lines(5);
        
        for i = 1:num_viz
            candidate_segment_id = seg_viz_dtw.segment_id{i};
            candidate_seg_vis_data = loadSegmentData(conn, schema, candidate_segment_id, 'position');
            
            if ~isempty(candidate_seg_vis_data)
                candidate_seg_vis_data = candidate_seg_vis_data - candidate_seg_vis_data(1, :);
                
                plot3(candidate_seg_vis_data(:,1), candidate_seg_vis_data(:,2), candidate_seg_vis_data(:,3), ...
                    '-', 'Color', colors_seg_dtw(i,:), 'LineWidth', 1.5, ...
                    'DisplayName', sprintf('#%d: %s (DTW: %.1f)', i, candidate_segment_id, seg_viz_dtw.dtw_distance(i)));
            end
        end
        
        grid on;
        axis equal;
        view(3);
        xlabel('X (mm)', 'FontSize', 11);
        ylabel('Y (mm)', 'FontSize', 11);
        zlabel('Z (mm)', 'FontSize', 11);
        title(sprintf('Segment %d - Top-5 by DTW Distance', seg_idx), 'FontSize', 13, 'FontWeight', 'bold');
        legend('Location', 'bestoutside', 'FontSize', 9);
        
        % Super title
        sgtitle(sprintf('Query Segment: %s | Pipeline: pgvector → LB_Keogh → DTW', query_segment_id), ...
            'FontSize', 14, 'FontWeight', 'bold');
        
        drawnow;
        
        fprintf('    ✓ Segment %d visualization created\n', seg_idx);
    end
end
fprintf('\n✓ All visualizations completed!\n');

% ========================================================================
% 6.6: Export Results (Optional)
% ========================================================================

fprintf('\n--- EXPORTING RESULTS ---\n');

% Create results directory if it doesn't exist
if ~exist('results', 'dir')
    mkdir('results');
end

% Export trajectory-level results
traj_results_file = sprintf('similarity/results/trajectory_results_%s.csv', query_bahn_id);
writetable(final_top_k, traj_results_file);
fprintf('✓ Trajectory results exported: %s\n', traj_results_file);

% Export segment-level results
if ~isempty(segment_results)
    for seg_idx = 1:num_query_segments
        % ✅ FIXED: Verwende query_segment_embeddings statt query_segments_table
        query_segment_id = query_segment_embeddings{seg_idx}.segment_id;
        seg_final = segment_results{seg_idx};
        
        if ~isempty(seg_final)
            seg_results_file = sprintf('results/segment_results_%s.csv', query_segment_id);
            writetable(seg_final, seg_results_file);
            fprintf('✓ Segment %d results exported: %s\n', seg_idx, seg_results_file);
        end
    end
end

fprintf('\n');

% ========================================================================
% 6.7: Close Database Connection
% ========================================================================

fprintf('--- CLEANUP ---\n');
close(conn);
fprintf('✓ Database connection closed\n');

fprintf('\n=== SEARCH COMPLETED ===\n\n');

function result = iif(condition, true_val, false_val)
    % Inline if-function (ternary operator)
    if condition
        result = true_val;
    else
        result = false_val;
    end
end

function data = loadTrajectoryPosition(conn, schema, bahn_id)
    % Load position data (x, y, z) for a trajectory
    query = sprintf(...
        ['SELECT x_soll, y_soll, z_soll ' ...
         'FROM robotervermessung.%s.bahn_position_soll ' ...
         'WHERE bahn_id = ''%s'' ' ...
         'ORDER BY timestamp'], ...
        schema, bahn_id);
    
    result = fetch(conn, query);
    data = table2array(result);
end

function data = loadTrajectoryJointStates(conn, schema, bahn_id)
    % Load joint states data (6 joints) for a trajectory
    query = sprintf(...
        ['SELECT joint_1, joint_2, joint_3, joint_4, joint_5, joint_6 ' ...
         'FROM robotervermessung.%s.bahn_joint_states ' ...
         'WHERE bahn_id = ''%s'' ' ...
         'ORDER BY timestamp'], ...
        schema, bahn_id);
    
    result = fetch(conn, query);
    data = table2array(result);
end

function data = loadSegmentData(conn, schema, segment_id, mode)
    % Helper function to load time-series data for a single segment
    % mode: 'position' or 'joint_states'
    
    if strcmp(mode, 'position')
        query = sprintf(...
            ['SELECT x_soll, y_soll, z_soll ' ...
             'FROM robotervermessung.%s.bahn_position_soll ' ...
             'WHERE segment_id = ''%s'' ' ...
             'ORDER BY timestamp ASC'], ...
            schema, segment_id);
    else
        % joint_states
        query = sprintf(...
            ['SELECT joint_1, joint_2, joint_3, joint_4, joint_5, joint_6 ' ...
             'FROM robotervermessung.%s.bahn_joint_states ' ...
             'WHERE segment_id = ''%s'' ' ...
             'ORDER BY timestamp ASC'], ...
            schema, segment_id);
    end
    
    result = fetch(conn, query);
    
    if isempty(result)
        data = [];
    else
        data = table2array(result);
    end
end

function embedding = parseEmbeddingString(embedding_str)
    % Parse pgvector text representation '[1.2, 3.4, ...]' to MATLAB array
    
    if isempty(embedding_str) || strcmp(embedding_str, '')
        embedding = [];
        return;
    end
    
    % Remove brackets and parse
    embedding_str = strrep(embedding_str, '[', '');
    embedding_str = strrep(embedding_str, ']', '');
    
    % Split by comma and convert to double
    parts = strsplit(embedding_str, ',');
    embedding = cellfun(@str2double, parts);
end
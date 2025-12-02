%  DTW BASELINE - TRAJECTORY SIMILARITY SEARCH (CACHE-OPTIMIZED)
%  ========================================================================
%  Hierarchical similarity search using Dynamic Time Warping (DTW)
%  
%  This version works with pre-loaded data_cache, dtw_cache, embeddings_cache
%  No database loading - all data must be in cache!
%  
%  Author: Gustavo Barros
%  Date: 02.12.2025 (Embeddings cache support added)
%  
%  ========================================================================

%% SECTION 0: CHECK FOR DATA CACHE
%  ========================================================================

if ~exist('data_cache', 'var')
    error(['data_cache not found in workspace!\n', ...
           'This version of dtw_baseline.m requires pre-loaded data.\n', ...
           'Please use loadDataExperiment.m first or pass data_cache via config.']);
end

fprintf('\n✓ Using pre-loaded data cache (no database queries will be made)\n');

% Check for DTW cache
if exist('dtw_cache', 'var') && ~isempty(dtw_cache)
    use_precomputed_dtw = true;
    fprintf('✓ Using pre-computed DTW cache (SKIPPING DTW computation)\n');
else
    use_precomputed_dtw = false;
    fprintf('⚠ No DTW cache found - will compute DTW\n');
end

% Check for Embeddings cache
if exist('embeddings_cache', 'var') && ~isempty(embeddings_cache)
    use_precomputed_embeddings = true;
    fprintf('✓ Using pre-computed Embeddings cache (SKIPPING Embedding computation)\n');
else
    use_precomputed_embeddings = false;
    fprintf('⚠ No Embeddings cache found - will compute embeddings\n');
end

%% SECTION 1: CONFIGURATION & PARAMETERS
%  ========================================================================

% === Query Parameters ===
if ~exist('query_bahn_id', 'var')
    query_bahn_id = '1763740056';
end

% === Database Sampling ===
if ~exist('use_database_sampling', 'var')
    use_database_sampling = true;
end
if ~exist('database_sample_size', 'var')
    database_sample_size = 30;
end
if ~exist('random_seed', 'var')
    random_seed = 42;
end

% === DTW Configuration ===
if ~exist('dtw_mode', 'var')
    dtw_mode = 'position';
end
if ~exist('normalize_dtw', 'var')
    normalize_dtw = false;
end
if ~exist('use_rotation_alignment', 'var')
    use_rotation_alignment = false;
end
if ~exist('cdtw_window', 'var')
    cdtw_window = 0.10;
end

% === Lower Bound Configuration ===
if ~exist('lb_kim_keep_ratio', 'var')
    lb_kim_keep_ratio = 0.40;
end
if ~exist('lb_keogh_keep_ratio', 'var')
    lb_keogh_candidates = 100;
end

% === Embedding Architecture ===
if ~exist('use_multi_scale', 'var')
    use_multi_scale = true;
end
if ~exist('n_coarse', 'var')
    n_coarse = 50;
end
if ~exist('n_medium', 'var')
    n_medium = 0;
end
if ~exist('n_fine', 'var')
    n_fine = 250;
end
if ~exist('norm_strategy', 'var')
    norm_strategy = 'max_extent';
end

% === Embedding Configuration ===
if ~exist('rrf_k', 'var')
    rrf_k = 60;
end

% === Weights (CRITICAL FOR EXPERIMENTS) ===
if ~exist('weights', 'var')
    weights = struct();
    weights.position = 1.0;
    weights.joint = 1.0;
    weights.orientation = 1.0;
    weights.velocity = 1.0;
    weights.metadata = 1.0;
end

% === Output Parameters ===
if ~exist('top_k_trajectories', 'var')
    top_k_trajectories = 50;
end

% === Database Configuration ===
if ~exist('schema', 'var')
    schema = 'bewegungsdaten';
end

% === Embedding Config Name (for cache lookup) ===
if ~exist('embedding_config_name', 'var')
    % Try to construct from parameters
    if use_multi_scale && n_coarse > 0
        embedding_config_name = sprintf('Multi-Dense-%d', n_coarse + n_fine);
    else
        embedding_config_name = sprintf('Single-Fine-%d', n_fine);
    end
end

fprintf('\n=== Configuration Summary ===\n');
fprintf('Query Trajectory: %s\n', query_bahn_id);
fprintf('DTW Mode: %s\n', dtw_mode);

if use_multi_scale
    fprintf('Architecture: Multi-Scale (%d + %d + %d samples)\n', n_coarse, n_medium, n_fine);
else
    fprintf('Architecture: Single-Scale\n');
end

fprintf('Normalization: %s\n', norm_strategy);
fprintf('Embedding Config: %s\n', embedding_config_name);
fprintf('\n');

%% SECTION 2: LOAD QUERY DATA FROM CACHE
%  ========================================================================

fprintf('=== Loading Query Trajectory from Cache: %s ===\n', query_bahn_id);

% Get query field name (convert ID to valid MATLAB field name)
query_field_name = ['q_' strrep(query_bahn_id, '-', '_')];

if ~isfield(data_cache.queries, query_field_name)
    error('Query %s not found in data_cache.queries', query_bahn_id);
end

query_cache = data_cache.queries.(query_field_name);

% === Load trajectory data for DTW ===
if strcmp(dtw_mode, 'position')
    query_data = query_cache.position;
    query_start = query_data(1, :);
    query_data = query_data - query_start;
    fprintf('✓ DTW data loaded from cache: Position (%d points, x/y/z)\n', size(query_data, 1));
    
elseif strcmp(dtw_mode, 'joint_states')
    query_data = query_cache.joint;
    fprintf('✓ DTW data loaded from cache: Joint States (%d points, 6 joints)\n', size(query_data, 1));
    
else
    error('Invalid dtw_mode. Use "position" or "joint_states"');
end

% === Load query metadata from cache ===
query_metadata = query_cache.metadata;

fprintf('  Query length: %.2f mm\n', query_metadata.length);
fprintf('  Query duration: %.2f s\n', query_metadata.duration);
fprintf('  Query segments: %d\n', query_metadata.num_segments);
fprintf('✓ Query metadata loaded from cache\n\n');

% === Load Query Segments from cache ===
query_segments_data = query_cache.segments;
num_query_segments = query_metadata.num_segments;

% Extract segment data based on DTW mode
if strcmp(dtw_mode, 'position')
    query_segments = query_segments_data.position;
elseif strcmp(dtw_mode, 'joint_states')
    query_segments = query_segments_data.joint;
end

% Extract segment metadata
query_segments_metadata = query_segments_data.metadata;

fprintf('✓ Loaded %d query segments from cache\n', num_query_segments);

%% SECTION 3: GET CANDIDATE DATA FROM CACHE
%  ========================================================================

fprintf('\n=== Loading Candidate Data from Cache ===\n');

% Get candidate metadata and trajectories
candidate_metadata = data_cache.candidates.metadata;
num_candidates = height(candidate_metadata);

fprintf('  Candidates in cache: %d trajectories\n', num_candidates);

% Get candidate trajectory data based on DTW mode
if strcmp(dtw_mode, 'position')
    candidate_trajectories = data_cache.candidates.position;
elseif strcmp(dtw_mode, 'joint_states')
    candidate_trajectories = data_cache.candidates.joint;
end

fprintf('✓ Candidate data loaded from cache\n');

% === Load Segment Metadata from cache ===
fprintf('\n=== Loading Segment Data from Cache ===\n');

all_segments_metadata = cell(num_query_segments, 1);

for seg_idx = 1:num_query_segments
    query_segment_id = sprintf('%s_%d', query_bahn_id, seg_idx);
    
    % Filter segments from cache for this query segment
    % Get all segments that are NOT from the query trajectory
    segment_mask = ~strcmp(data_cache.segments.segment_ids, query_segment_id);
    
    % Get filtered segment metadata
    all_segments_metadata{seg_idx} = data_cache.segments.metadata(segment_mask, :);
    
    fprintf('  Segment %d: %d candidates\n', seg_idx, height(all_segments_metadata{seg_idx}));
end

fprintf('\n');

%% SECTION 4: DTW COMPUTATION - TRAJECTORY LEVEL
%  ========================================================================

if use_precomputed_dtw
    % ====================================================================
    % USE PRE-COMPUTED DTW FROM CACHE
    % ====================================================================
    fprintf('\n=== Loading Pre-Computed DTW Results (Trajectory Level) ===\n');
    
    % Get query field name
    query_field_name = ['q_' strrep(query_bahn_id, '-', '_')];
    
    % Check if query exists in cache
    if ~isfield(dtw_cache, query_field_name)
        error('Query %s not found in dtw_cache', query_bahn_id);
    end
    
    % Check if mode exists
    if ~isfield(dtw_cache.(query_field_name), dtw_mode)
        error('DTW mode %s not found in cache for query %s', dtw_mode, query_bahn_id);
    end
    
    % Load from cache
    trajectory_table = dtw_cache.(query_field_name).(dtw_mode).trajectory_ranking;
    dtw_trajectory_time = dtw_cache.(query_field_name).(dtw_mode).dtw_time;
    num_candidates = dtw_cache.(query_field_name).(dtw_mode).num_candidates;
    
    fprintf('  ✓ Loaded trajectory DTW results from cache\n');
    fprintf('  Candidates: %d\n', num_candidates);
    fprintf('  DTW time (pre-computed): %.2fs\n', dtw_trajectory_time);
    fprintf('  Best match: %s (distance: %.4f)\n', ...
        trajectory_table.bahn_id{1}, trajectory_table.dtw_distance(1));
    
else
    % ====================================================================
    % COMPUTE DTW (FALLBACK - should rarely happen)
    % ====================================================================
    fprintf('\n=== Computing DTW Distances (Trajectory Level) ===\n');
    fprintf('  Total candidates: %d\n', num_candidates);
    fprintf('  ✓ Data already in memory (from cache)\n');

    tic;

    % ========================================================================
    % PHASE 1: LB_Kim Pre-filtering
    % ========================================================================
    fprintf('Phase 1: LB_Kim filtering...\n');

    lb_kim_distances = zeros(num_candidates, 1);

    for i = 1:num_candidates
        candidate_data = candidate_trajectories{i};
        lb_kim_distances(i) = LB_Kim(query_data, candidate_data, dtw_mode, use_rotation_alignment, normalize_dtw);
    end

    [~, kim_order] = sort(lb_kim_distances);
    kim_keep_count = round(num_candidates * lb_kim_keep_ratio);

    candidates_after_kim = candidate_metadata(kim_order(1:kim_keep_count), :);
    candidate_trajectories_after_kim = candidate_trajectories(kim_order(1:kim_keep_count));

    fprintf('  ✓ LB_Kim: %d → %d candidates (%.1f%% pruned)\n', ...
        num_candidates, kim_keep_count, 100*(1 - lb_kim_keep_ratio));

    % ========================================================================
    % PHASE 2: LB_Keogh Pre-filtering
    % ========================================================================
    fprintf('Phase 2: LB_Keogh filtering...\n');

    keogh_keep_count = min(lb_keogh_candidates, kim_keep_count);

    fprintf('  Candidates from LB_Kim: %d\n', kim_keep_count);
    fprintf('  Target for LB_Keogh: %d (fixed)\n', keogh_keep_count);

    if kim_keep_count <= keogh_keep_count
        fprintf('  ⚠ LB_Kim already filtered to target size - skipping LB_Keogh\n');
        candidates_after_keogh = candidates_after_kim;
        candidate_trajectories_after_keogh = candidate_trajectories_after_kim;
        fprintf('  ✓ Kept all %d candidates\n', kim_keep_count);
    else
        lb_keogh_distances = zeros(kim_keep_count, 1);
        
        for i = 1:kim_keep_count
            candidate_data = candidate_trajectories_after_kim{i};
            lb_keogh_distances(i) = LB_Keogh(query_data, candidate_data, ...
                cdtw_window, dtw_mode, use_rotation_alignment, normalize_dtw);
        end
        
        [~, keogh_order] = sort(lb_keogh_distances);
        candidates_after_keogh = candidates_after_kim(keogh_order(1:keogh_keep_count), :);
        candidate_trajectories_after_keogh = candidate_trajectories_after_kim(keogh_order(1:keogh_keep_count));
        
        fprintf('  ✓ LB_Keogh: %d → %d candidates (%.1f%% pruned)\n', ...
            kim_keep_count, keogh_keep_count, 100*(1 - keogh_keep_count/kim_keep_count));
    end

    total_pruning = 1 - (keogh_keep_count / num_candidates);
    fprintf('  ✓ Total Pruning: %.1f%% (from %d to %d before DTW)\n', ...
        total_pruning * 100, num_candidates, keogh_keep_count);

    % ========================================================================
    % PHASE 3: DTW Refinement
    % ========================================================================
    fprintf('Phase 3: DTW computation on %d candidates...\n', keogh_keep_count);

    dtw_limit = min(top_k_trajectories, keogh_keep_count);
    dtw_traj_results = struct();
    dtw_traj_results.bahn_id = cell(dtw_limit, 1);
    dtw_traj_results.dtw_distance = inf(dtw_limit, 1);
    dtw_traj_results.length = zeros(dtw_limit, 1);
    dtw_traj_results.duration = zeros(dtw_limit, 1);

    for i = 1:dtw_limit
        candidate_id = candidates_after_keogh.bahn_id{i};
        candidate_data = candidate_trajectories_after_keogh{i};
        
        if i <= top_k_trajectories
            best_so_far = inf;
        else
            sorted_dists = sort(dtw_traj_results.dtw_distance(1:i-1));
            best_so_far = sorted_dists(top_k_trajectories);
        end
        
        dist = cDTW(query_data, candidate_data, dtw_mode, cdtw_window, ...
            best_so_far, use_rotation_alignment, normalize_dtw);
        
        dtw_traj_results.bahn_id{i} = candidate_id;
        dtw_traj_results.dtw_distance(i) = dist;
        dtw_traj_results.length(i) = candidates_after_keogh.length(i);
        dtw_traj_results.duration(i) = candidates_after_keogh.duration(i);
        
        if mod(i, 10) == 0 || i == dtw_limit
            fprintf('  Progress: %d/%d trajectories\n', i, dtw_limit);
        end
    end

    dtw_trajectory_time = toc;

    trajectory_table = struct2table(dtw_traj_results);
    trajectory_table = sortrows(trajectory_table, 'dtw_distance', 'ascend');

    fprintf('✓ DTW computation completed in %.2f seconds\n', dtw_trajectory_time);
end

%% SECTION 5: DTW COMPUTATION - SEGMENT LEVEL
%  ========================================================================

if use_precomputed_dtw
    % ====================================================================
    % USE PRE-COMPUTED SEGMENT DTW FROM CACHE
    % ====================================================================
    fprintf('\n=== Loading Pre-Computed DTW Results (Segment Level) ===\n');
    
    % Get query field name
    query_field_name = ['q_' strrep(query_bahn_id, '-', '_')];
    
    % Load from cache
    segment_results = dtw_cache.(query_field_name).(dtw_mode).segment_rankings;
    segment_dtw_times = dtw_cache.(query_field_name).(dtw_mode).segment_dtw_times;
    segment_dtw_time = mean(segment_dtw_times);
    num_query_segments = length(segment_results);
    
    fprintf('  ✓ Loaded segment DTW results from cache\n');
    fprintf('  Query segments: %d\n', num_query_segments);
    fprintf('  Avg DTW time per segment (pre-computed): %.2fs\n', segment_dtw_time);
    
    % Show best match for each segment
    for seg_idx = 1:num_query_segments
        if ~isempty(segment_results{seg_idx})
            fprintf('  Segment %d best match: %s (distance: %.4f)\n', ...
                seg_idx, segment_results{seg_idx}.segment_id{1}, ...
                segment_results{seg_idx}.dtw_distance(1));
        end
    end
    
else
    % ====================================================================
    % COMPUTE SEGMENT DTW (FALLBACK)
    % ====================================================================
    fprintf('\n=== Computing DTW Distances (Segment Level) ===\n');
    fprintf('  Query segments: %d\n', num_query_segments);

    segment_results = cell(num_query_segments, 1);
    segment_dtw_times = zeros(num_query_segments, 1);

    for seg_idx = 1:num_query_segments
        query_segment_id = sprintf('%s_%d', query_bahn_id, seg_idx);
        query_segment_data = query_segments{seg_idx};
        
        fprintf('\n  Processing Query Segment %d/%d (ID: %s)...\n', seg_idx, num_query_segments, query_segment_id);
        
        candidate_segments = all_segments_metadata{seg_idx};
        num_candidate_segments = height(candidate_segments);
        
        if num_candidate_segments == 0
            fprintf('    No candidate segments found - skipping\n');
            segment_results{seg_idx} = [];
            continue;
        end
        
        fprintf('    Total segment candidates: %d\n', num_candidate_segments);
        fprintf('    ✓ Data already in memory (from cache)\n');
        
        % Get segment IDs
        seg_candidate_ids = candidate_segments.segment_id;
        
        % Find indices in cache
        [~, seg_cache_idx] = ismember(seg_candidate_ids, data_cache.segments.segment_ids);
        
        % Get segment data from cache
        if strcmp(dtw_mode, 'position')
            seg_candidate_segments = data_cache.segments.position(seg_cache_idx);
        elseif strcmp(dtw_mode, 'joint_states')
            seg_candidate_segments = data_cache.segments.joint(seg_cache_idx);
        end
        
        tic;
        
        % LB_Kim Pre-filtering
        fprintf('    Step 1: LB_Kim screening...\n');
        
        segment_lb_kim_distances = zeros(num_candidate_segments, 1);
        
        for cand_seg_idx = 1:num_candidate_segments
            candidate_segment_data = seg_candidate_segments{cand_seg_idx};
            segment_lb_kim_distances(cand_seg_idx) = LB_Kim(query_segment_data, ...
                candidate_segment_data, dtw_mode, use_rotation_alignment, normalize_dtw);
        end
        
        [~, seg_kim_order] = sort(segment_lb_kim_distances);
        seg_kim_keep_count = round(num_candidate_segments * lb_kim_keep_ratio);
        
        candidates_after_kim_seg = candidate_segments(seg_kim_order(1:seg_kim_keep_count), :);
        seg_candidate_segments_after_kim = seg_candidate_segments(seg_kim_order(1:seg_kim_keep_count));
        
        fprintf('      ✓ LB_Kim: %d → %d candidates (%.1f%% pruned)\n', ...
            num_candidate_segments, seg_kim_keep_count, 100*(1 - lb_kim_keep_ratio));
        
        % LB_Keogh Pre-filtering
        fprintf('    Step 2: LB_Keogh refinement...\n');
        
        seg_keogh_keep_count = min(lb_keogh_candidates, seg_kim_keep_count);
        
        fprintf('      Candidates from LB_Kim: %d\n', seg_kim_keep_count);
        fprintf('      Target for LB_Keogh: %d (fixed)\n', seg_keogh_keep_count);
        
        if seg_kim_keep_count <= seg_keogh_keep_count
            fprintf('      ⚠ LB_Kim already filtered to target size - skipping LB_Keogh\n');
            candidates_after_keogh_seg = candidates_after_kim_seg;
            seg_candidate_segments_after_keogh = seg_candidate_segments_after_kim;
            fprintf('      ✓ Kept all %d segment candidates\n', seg_kim_keep_count);
        else
            segment_lb_keogh_distances = zeros(seg_kim_keep_count, 1);
            
            for cand_seg_idx = 1:seg_kim_keep_count
                candidate_segment_data = seg_candidate_segments_after_kim{cand_seg_idx};
                segment_lb_keogh_distances(cand_seg_idx) = LB_Keogh(query_segment_data, ...
                    candidate_segment_data, cdtw_window, dtw_mode, use_rotation_alignment, normalize_dtw);
            end
            
            [~, seg_keogh_order] = sort(segment_lb_keogh_distances);
            candidates_after_keogh_seg = candidates_after_kim_seg(seg_keogh_order(1:seg_keogh_keep_count), :);
            seg_candidate_segments_after_keogh = seg_candidate_segments_after_kim(seg_keogh_order(1:seg_keogh_keep_count));
            
            fprintf('      ✓ LB_Keogh: %d → %d candidates (%.1f%% pruned)\n', ...
                seg_kim_keep_count, seg_keogh_keep_count, ...
                100*(1 - seg_keogh_keep_count/seg_kim_keep_count));
        end
        
        seg_total_pruning = 1 - (seg_keogh_keep_count / num_candidate_segments);
        fprintf('      ✓ Total Pruning: %.1f%% (from %d to %d before DTW)\n', ...
            seg_total_pruning * 100, num_candidate_segments, seg_keogh_keep_count);
        
        % DTW computation
        fprintf('    Step 3: DTW computation...\n');
        
        seg_dtw_limit = min(top_k_trajectories, seg_keogh_keep_count);
        segment_dtw_distances = inf(seg_dtw_limit, 1);
        
        for cand_seg_idx = 1:seg_dtw_limit
            candidate_segment_data = seg_candidate_segments_after_keogh{cand_seg_idx};
            
            if ~isempty(candidate_segment_data)
                candidate_segment_data = candidate_segment_data - candidate_segment_data(1, :);
            end
            
            if cand_seg_idx <= 10
                best_so_far = inf;
            else
                sorted_seg_dists = sort(segment_dtw_distances(1:cand_seg_idx-1));
                best_so_far = sorted_seg_dists(10);
            end
            
            dist = cDTW(query_segment_data, candidate_segment_data, dtw_mode, cdtw_window, ...
                best_so_far, use_rotation_alignment, normalize_dtw);
            segment_dtw_distances(cand_seg_idx) = dist;
            
            if mod(cand_seg_idx, 50) == 0 || cand_seg_idx == seg_dtw_limit
                fprintf('      Progress: %d/%d segments\n', cand_seg_idx, seg_dtw_limit);
            end
        end
        
        segment_results{seg_idx} = table(...
            candidates_after_keogh_seg.segment_id(1:seg_dtw_limit), ...
            candidates_after_keogh_seg.bahn_id(1:seg_dtw_limit), ...
            segment_dtw_distances(1:seg_dtw_limit), ...
            candidates_after_keogh_seg.length(1:seg_dtw_limit), ...
            candidates_after_keogh_seg.duration(1:seg_dtw_limit), ...
            'VariableNames', {'segment_id', 'bahn_id', 'dtw_distance', 'length', 'duration'});
        
        segment_results{seg_idx} = sortrows(segment_results{seg_idx}, 'dtw_distance', 'ascend');
        segment_dtw_times(seg_idx) = toc;
        
        fprintf('    ✓ Segment %d: Compared against %d candidates in %.2f seconds\n', ...
            seg_idx, num_candidate_segments, segment_dtw_times(seg_idx));
    end

    segment_dtw_time = mean(segment_dtw_times);

    fprintf('\n✓ Segment-level DTW computation completed\n');
end  % End of else block for use_precomputed_dtw

%% SECTION 6: COMPUTE EMBEDDINGS (TRAJECTORY LEVEL)
%  ========================================================================

if use_precomputed_embeddings
    % ====================================================================
    % USE PRE-COMPUTED EMBEDDINGS FROM CACHE
    % ====================================================================
    fprintf('\n=== Loading Pre-Computed Embeddings (Trajectory Level) ===\n');
    
    % Get field names
    query_field_name = ['q_' strrep(query_bahn_id, '-', '_')];
    emb_field_name = strrep(embedding_config_name, '-', '_');
    emb_field_name = strrep(emb_field_name, ' ', '_');
    
    % Check if embeddings exist in cache
    if ~isfield(embeddings_cache, query_field_name)
        error('Query %s not found in embeddings_cache', query_bahn_id);
    end
    
    if ~isfield(embeddings_cache.(query_field_name), emb_field_name)
        error('Embedding config %s not found in cache for query %s', ...
            embedding_config_name, query_bahn_id);
    end
    
    % Load from cache
    query_embeddings = embeddings_cache.(query_field_name).(emb_field_name).query_embeddings;
    candidate_embeddings = embeddings_cache.(query_field_name).(emb_field_name).candidate_embeddings;
    
    fprintf('  ✓ Loaded embeddings from cache\n');
    fprintf('  Embedding config: %s\n', embedding_config_name);
    fprintf('  Query embeddings: 5 modalities\n');
    fprintf('  Candidate embeddings: %d × 5 modalities\n', size(candidate_embeddings.position, 1));
    
    % Extract individual embeddings
    query_pos_embedding = query_embeddings.position;
    query_joint_embedding = query_embeddings.joint;
    query_orient_embedding = query_embeddings.orientation;
    query_vel_embedding = query_embeddings.velocity;
    query_meta_embedding = query_embeddings.metadata;
    
    pos_embeddings = candidate_embeddings.position;
    joint_embeddings = candidate_embeddings.joint;
    orient_embeddings = candidate_embeddings.orientation;
    vel_embeddings = candidate_embeddings.velocity;
    meta_embeddings = candidate_embeddings.metadata;
    
else
    % ====================================================================
    % COMPUTE EMBEDDINGS (FALLBACK)
    % ====================================================================
    fprintf('\n=== Computing Embeddings (Trajectory Level) ===\n');

    % Check which modalities are active
    active_modalities.position = weights(1) > 0;
    active_modalities.joint = weights(2) > 0;
    active_modalities.orientation = weights(3) > 0;
    active_modalities.velocity = weights(4) > 0;
    active_modalities.metadata = weights(5) > 0;

    fprintf('  Active modalities: ');
    active_names = {};
    if active_modalities.position, active_names{end+1} = 'Pos'; end
    if active_modalities.joint, active_names{end+1} = 'Joint'; end
    if active_modalities.orientation, active_names{end+1} = 'Orient'; end
    if active_modalities.velocity, active_names{end+1} = 'Vel'; end
    if active_modalities.metadata, active_names{end+1} = 'Meta'; end
    fprintf('%s\n', strjoin(active_names, ', '));

    % ========================================================================
    % Query Embeddings (from cache)
    % ========================================================================
    fprintf('  Computing query embeddings from cache...\n');

    if active_modalities.position
        query_pos_data = query_cache.position;
        query_pos_embedding = computePositionEmbedding(query_pos_data, n_coarse, n_fine);
    else
        query_pos_embedding = [];
    end

    if active_modalities.joint
        query_joint_data = query_cache.joint;
        query_joint_embedding = computeJointEmbedding(query_joint_data, n_coarse, n_fine);
    else
        query_joint_embedding = [];
    end

    if active_modalities.orientation
        query_orient_data = query_cache.orientation;
        query_orient_embedding = computeOrientationEmbedding(query_orient_data, n_coarse, n_fine);
    else
        query_orient_embedding = [];
    end

    if active_modalities.velocity
        if ~active_modalities.position
            query_pos_data = query_cache.position;
        end
        query_vel_embedding = computeVelocityEmbedding(query_pos_data, n_coarse, n_fine);
    else
        query_vel_embedding = [];
    end

    if active_modalities.metadata
        query_meta_embedding = computeMetadataEmbedding(query_metadata);
    else
        query_meta_embedding = [];
    end

    fprintf('    ✓ Query embeddings computed: ');
    dims = {};
    if active_modalities.position, dims{end+1} = sprintf('Pos(%d)', length(query_pos_embedding)); end
    if active_modalities.joint, dims{end+1} = sprintf('Joint(%d)', length(query_joint_embedding)); end
    if active_modalities.orientation, dims{end+1} = sprintf('Orient(%d)', length(query_orient_embedding)); end
    if active_modalities.velocity, dims{end+1} = sprintf('Vel(%d)', length(query_vel_embedding)); end
    if active_modalities.metadata, dims{end+1} = sprintf('Meta(%d)', length(query_meta_embedding)); end
    fprintf('%s\n', strjoin(dims, ', '));

    % ========================================================================
    % Get Candidate Data from Cache
    % ========================================================================
    fprintf('  Getting candidate data from cache...\n');

    if active_modalities.position || active_modalities.velocity
        pos_data = data_cache.candidates.position;
    else
        pos_data = [];
    end

    if active_modalities.joint
        joint_data = data_cache.candidates.joint;
    else
        joint_data = [];
    end

    if active_modalities.orientation
        orient_data = data_cache.candidates.orientation;
    else
        orient_data = [];
    end

    fprintf('    ✓ Candidate data retrieved from cache\n');

    % ========================================================================
    % Compute Candidate Embeddings
    % ========================================================================
    fprintf('  Computing candidate embeddings for %d trajectories...\n', num_candidates);

    if active_modalities.position
        pos_embeddings = zeros(num_candidates, length(query_pos_embedding));
    end
    if active_modalities.joint
        joint_embeddings = zeros(num_candidates, length(query_joint_embedding));
    end
    if active_modalities.orientation
        orient_embeddings = zeros(num_candidates, length(query_orient_embedding));
    end
    if active_modalities.velocity
        vel_embeddings = zeros(num_candidates, length(query_vel_embedding));
    end
    if active_modalities.metadata
        meta_embeddings = zeros(num_candidates, length(query_meta_embedding));
    end

    for i = 1:num_candidates
        if active_modalities.position
            pos_embeddings(i, :) = computePositionEmbedding(pos_data{i}, n_coarse, n_fine);
        end
        
        if active_modalities.joint
            joint_embeddings(i, :) = computeJointEmbedding(joint_data{i}, n_coarse, n_fine);
        end
        
        if active_modalities.orientation
            orient_embeddings(i, :) = computeOrientationEmbedding(orient_data{i}, n_coarse, n_fine);
        end
        
        if active_modalities.velocity
            vel_embeddings(i, :) = computeVelocityEmbedding(pos_data{i}, n_coarse, n_fine);
        end
        
        if active_modalities.metadata
            candidate_meta = struct();
            candidate_meta.movement_type = candidate_metadata.movement_type{i};
            candidate_meta.length = candidate_metadata.length(i);
            candidate_meta.duration = candidate_metadata.duration(i);
            candidate_meta.min_twist = candidate_metadata.min_twist_ist(i);
            candidate_meta.max_twist = candidate_metadata.max_twist_ist(i);
            candidate_meta.mean_twist = candidate_metadata.mean_twist_ist(i);
            candidate_meta.median_twist = candidate_metadata.median_twist_ist(i);
            candidate_meta.std_twist = candidate_metadata.std_twist_ist(i);
            candidate_meta.min_acceleration = candidate_metadata.min_acceleration_ist(i);
            candidate_meta.max_acceleration = candidate_metadata.max_acceleration_ist(i);
            candidate_meta.mean_acceleration = candidate_metadata.mean_acceleration_ist(i);
            candidate_meta.median_acceleration = candidate_metadata.median_acceleration_ist(i);
            candidate_meta.std_acceleration = candidate_metadata.std_acceleration_ist(i);
            
            meta_embeddings(i, :) = computeMetadataEmbedding(candidate_meta);
        end
        
        if mod(i, 100) == 0
            fprintf('    Progress: %d/%d\n', i, num_candidates);
        end
    end

    fprintf('  ✓ All embeddings computed\n');
end  % End of use_precomputed_embeddings check

% ========================================================================
% Compute Cosine Distances (ALWAYS NEEDED - even with pre-computed embeddings)
% ========================================================================
fprintf('\n=== Computing Cosine Distances & Rankings ===\n');

% Check which modalities are active
active_modalities.position = weights(1) > 0;
active_modalities.joint = weights(2) > 0;
active_modalities.orientation = weights(3) > 0;
active_modalities.velocity = weights(4) > 0;
active_modalities.metadata = weights(5) > 0;

if active_modalities.position
    pos_distances = 1 - (pos_embeddings * query_pos_embedding');
    pos_ranking = createRankingTable(candidate_metadata, pos_distances);
else
    pos_ranking = [];
end

if active_modalities.joint
    joint_distances = 1 - (joint_embeddings * query_joint_embedding');
    joint_ranking = createRankingTable(candidate_metadata, joint_distances);
else
    joint_ranking = [];
end

if active_modalities.orientation
    orient_distances = 1 - (orient_embeddings * query_orient_embedding');
    orient_ranking = createRankingTable(candidate_metadata, orient_distances);
else
    orient_ranking = [];
end

if active_modalities.velocity
    vel_distances = 1 - (vel_embeddings * query_vel_embedding');
    vel_ranking = createRankingTable(candidate_metadata, vel_distances);
else
    vel_ranking = [];
end

if active_modalities.metadata
    meta_distances = 1 - (meta_embeddings * query_meta_embedding');
    meta_ranking = createRankingTable(candidate_metadata, meta_distances);
else
    meta_ranking = [];
end

fprintf('  ✓ Rankings created\n');

% ========================================================================
% RRF Fusion
% ========================================================================
fprintf('\n=== RRF Fusion (Active Modalities) ===\n');

rankings = struct();
if active_modalities.position, rankings.position = pos_ranking; end
if active_modalities.joint, rankings.joint = joint_ranking; end
if active_modalities.orientation, rankings.orientation = orient_ranking; end
if active_modalities.velocity, rankings.velocity = vel_ranking; end
if active_modalities.metadata, rankings.metadata = meta_ranking; end

fused_ranking = fuseRankingsRRF(rankings, weights, rrf_k, 'bahn_id');

fprintf('  ✓ Best match: %s (RRF: %.6f)\n', fused_ranking.bahn_id{1}, fused_ranking.rrf_score(1));

%% SECTION 7: COMPUTE EMBEDDINGS (SEGMENT LEVEL)
%  ========================================================================

if num_query_segments > 0
    
    if use_precomputed_embeddings
        % ================================================================
        % USE PRE-COMPUTED SEGMENT EMBEDDINGS
        % ================================================================
        fprintf('\n=== Loading Pre-Computed Embeddings (Segment Level) ===\n');
        
        segment_embeddings_from_cache = embeddings_cache.(query_field_name).(emb_field_name).segment_embeddings;
        
        fprintf('  ✓ Loaded segment embeddings from cache\n');
        fprintf('  Query segments: %d\n', num_query_segments);
        
        segment_embedding_results = cell(num_query_segments, 1);
        
        for seg_idx = 1:num_query_segments
            seg_emb = segment_embeddings_from_cache{seg_idx};
            
            % Get segment DTW ranking to know which candidates to rank
            segment_ranking = segment_results{seg_idx};
            
            if isempty(segment_ranking)
                segment_embedding_results{seg_idx} = [];
                continue;
            end
            
            % Get query embeddings
            query_seg_pos_embedding = seg_emb.query.position;
            query_seg_joint_embedding = seg_emb.query.joint;
            query_seg_orient_embedding = seg_emb.query.orientation;
            query_seg_vel_embedding = seg_emb.query.velocity;
            query_seg_meta_embedding = seg_emb.query.metadata;
            
            % Get candidate embeddings (already in matrix form!)
            seg_pos_embeddings = seg_emb.candidates.position;
            seg_joint_embeddings = seg_emb.candidates.joint;
            seg_orient_embeddings = seg_emb.candidates.orientation;
            seg_vel_embeddings = seg_emb.candidates.velocity;
            seg_meta_embeddings = seg_emb.candidates.metadata;
            
            % Compute distances
            if active_modalities.position
                seg_pos_distances = 1 - (seg_pos_embeddings * query_seg_pos_embedding');
            end
            
            if active_modalities.joint
                seg_joint_distances = 1 - (seg_joint_embeddings * query_seg_joint_embedding');
            end
            
            if active_modalities.orientation
                seg_orient_distances = 1 - (seg_orient_embeddings * query_seg_orient_embedding');
            end
            
            if active_modalities.velocity
                seg_vel_distances = 1 - (seg_vel_embeddings * query_seg_vel_embedding');
            end
            
            if active_modalities.metadata
                seg_meta_distances = 1 - (seg_meta_embeddings * query_seg_meta_embedding');
            end
            
            % Create rankings
            % Get metadata for segments from cache
            candidate_seg_ids = seg_emb.candidates.segment_ids;
            [~, seg_cache_idx] = ismember(candidate_seg_ids, data_cache.segments.segment_ids);
            candidate_seg_metadata = data_cache.segments.metadata(seg_cache_idx, :);
            
            seg_rankings = struct();
            if active_modalities.position
                seg_rankings.position = createRankingTable(candidate_seg_metadata, seg_pos_distances);
            end
            if active_modalities.joint
                seg_rankings.joint = createRankingTable(candidate_seg_metadata, seg_joint_distances);
            end
            if active_modalities.orientation
                seg_rankings.orientation = createRankingTable(candidate_seg_metadata, seg_orient_distances);
            end
            if active_modalities.velocity
                seg_rankings.velocity = createRankingTable(candidate_seg_metadata, seg_vel_distances);
            end
            if active_modalities.metadata
                seg_rankings.metadata = createRankingTable(candidate_seg_metadata, seg_meta_distances);
            end
            
            % RRF Fusion
            seg_fused_ranking = fuseRankingsRRF(seg_rankings, weights, rrf_k, 'segment_id');
            segment_embedding_results{seg_idx} = seg_fused_ranking;
            
            if mod(seg_idx, 5) == 0 || seg_idx == num_query_segments
                fprintf('  Segment %d/%d: Best match %s (RRF: %.6f)\n', ...
                    seg_idx, num_query_segments, ...
                    seg_fused_ranking.segment_id{1}, seg_fused_ranking.rrf_score(1));
            end
        end
        
        fprintf('✓ Segment-level embeddings (from cache) completed\n');
        
    else
        % ================================================================
        % COMPUTE SEGMENT EMBEDDINGS (FALLBACK)
        % ================================================================
        fprintf('\n=== Computing Embeddings (Segment Level) ===\n');
        
        segment_embedding_results = cell(num_query_segments, 1);
        
        for seg_idx = 1:num_query_segments
            query_segment_id = sprintf('%s_%d', query_bahn_id, seg_idx);
            
            fprintf('\n--- Query Segment %d/%d (ID: %s) ---\n', seg_idx, num_query_segments, query_segment_id);
            
            candidate_segments = all_segments_metadata{seg_idx};
            num_candidate_segments = height(candidate_segments);
            
            if num_candidate_segments == 0
                fprintf('  No candidate segments - skipping\n');
                segment_embedding_results{seg_idx} = [];
                continue;
            end
            
            fprintf('  Candidates: %d segments\n', num_candidate_segments);
            
            % Query segment embeddings from cache
            fprintf('  Computing query segment embeddings from cache...\n');
            
            if active_modalities.position
                query_seg_pos_data = query_segments_data.position{seg_idx};
                query_seg_pos_embedding = computePositionEmbedding(query_seg_pos_data, n_coarse, n_fine);
            else
                query_seg_pos_embedding = [];
            end
            
            if active_modalities.joint
                query_seg_joint_data = query_segments_data.joint{seg_idx};
                query_seg_joint_embedding = computeJointEmbedding(query_seg_joint_data, n_coarse, n_fine);
            else
                query_seg_joint_embedding = [];
            end
            
            if active_modalities.orientation
                query_seg_orient_data = query_segments_data.orientation{seg_idx};
                query_seg_orient_embedding = computeOrientationEmbedding(query_seg_orient_data, n_coarse, n_fine);
            else
                query_seg_orient_embedding = [];
            end
            
            if active_modalities.velocity
                if ~active_modalities.position
                    query_seg_pos_data = query_segments_data.position{seg_idx};
                end
                query_seg_vel_embedding = computeVelocityEmbedding(query_seg_pos_data, n_coarse, n_fine);
            else
                query_seg_vel_embedding = [];
            end
            
            if active_modalities.metadata
                query_seg_meta = query_segments_metadata{seg_idx};
                query_seg_meta_embedding = computeMetadataEmbedding(query_seg_meta);
            else
                query_seg_meta_embedding = [];
            end
            
            % Get candidate segment data from cache
            fprintf('  Getting candidate segment data from cache...\n');
            
            seg_candidate_ids = candidate_segments.segment_id;
            [~, seg_cache_idx] = ismember(seg_candidate_ids, data_cache.segments.segment_ids);
            
            if active_modalities.position || active_modalities.velocity
                seg_pos_data = data_cache.segments.position(seg_cache_idx);
            else
                seg_pos_data = [];
            end
            
            if active_modalities.joint
                seg_joint_data = data_cache.segments.joint(seg_cache_idx);
            else
                seg_joint_data = [];
            end
            
            if active_modalities.orientation
                seg_orient_data = data_cache.segments.orientation(seg_cache_idx);
            else
                seg_orient_data = [];
            end
            
            % Compute candidate segment embeddings
            fprintf('  Computing candidate segment embeddings for %d segments...\n', num_candidate_segments);
            
            if active_modalities.position
                seg_pos_embeddings = zeros(num_candidate_segments, length(query_seg_pos_embedding));
            end
            if active_modalities.joint
                seg_joint_embeddings = zeros(num_candidate_segments, length(query_seg_joint_embedding));
            end
            if active_modalities.orientation
                seg_orient_embeddings = zeros(num_candidate_segments, length(query_seg_orient_embedding));
            end
            if active_modalities.velocity
                seg_vel_embeddings = zeros(num_candidate_segments, length(query_seg_vel_embedding));
            end
            if active_modalities.metadata
                seg_meta_embeddings = zeros(num_candidate_segments, length(query_seg_meta_embedding));
            end
            
            for i = 1:num_candidate_segments
                if active_modalities.position
                    seg_pos_embeddings(i, :) = computePositionEmbedding(seg_pos_data{i}, n_coarse, n_fine);
                end
                
                if active_modalities.joint
                    seg_joint_embeddings(i, :) = computeJointEmbedding(seg_joint_data{i}, n_coarse, n_fine);
                end
                
                if active_modalities.orientation
                    seg_orient_embeddings(i, :) = computeOrientationEmbedding(seg_orient_data{i}, n_coarse, n_fine);
                end
                
                if active_modalities.velocity
                    seg_vel_embeddings(i, :) = computeVelocityEmbedding(seg_pos_data{i}, n_coarse, n_fine);
                end
                
                if active_modalities.metadata
                    candidate_seg_meta = struct();
                    candidate_seg_meta.movement_type = candidate_segments.movement_type{i};
                    candidate_seg_meta.length = candidate_segments.length(i);
                    candidate_seg_meta.duration = candidate_segments.duration(i);
                    candidate_seg_meta.min_twist = candidate_segments.min_twist_ist(i);
                    candidate_seg_meta.max_twist = candidate_segments.max_twist_ist(i);
                    candidate_seg_meta.mean_twist = candidate_segments.mean_twist_ist(i);
                    candidate_seg_meta.median_twist = candidate_segments.median_twist_ist(i);
                    candidate_seg_meta.std_twist = candidate_segments.std_twist_ist(i);
                    candidate_seg_meta.min_acceleration = candidate_segments.min_acceleration_ist(i);
                    candidate_seg_meta.max_acceleration = candidate_segments.max_acceleration_ist(i);
                    candidate_seg_meta.mean_acceleration = candidate_segments.mean_acceleration_ist(i);
                    candidate_seg_meta.median_acceleration = candidate_segments.median_acceleration_ist(i);
                    candidate_seg_meta.std_acceleration = candidate_segments.std_acceleration_ist(i);
                
                    seg_meta_embeddings(i, :) = computeMetadataEmbedding(candidate_seg_meta);
                end
                
                if mod(i, 50) == 0
                    fprintf('    Progress: %d/%d\n', i, num_candidate_segments);
                end
            end
            
            % Compute cosine distances
            fprintf('  Computing cosine distances...\n');
            
            if active_modalities.position
                seg_pos_distances = 1 - (seg_pos_embeddings * query_seg_pos_embedding');
                seg_pos_ranking = createRankingTable(candidate_segments, seg_pos_distances);
            else
                seg_pos_ranking = [];
            end
            
            if active_modalities.joint
                seg_joint_distances = 1 - (seg_joint_embeddings * query_seg_joint_embedding');
                seg_joint_ranking = createRankingTable(candidate_segments, seg_joint_distances);
            else
                seg_joint_ranking = [];
            end
            
            if active_modalities.orientation
                seg_orient_distances = 1 - (seg_orient_embeddings * query_seg_orient_embedding');
                seg_orient_ranking = createRankingTable(candidate_segments, seg_orient_distances);
            else
                seg_orient_ranking = [];
            end
            
            if active_modalities.velocity
                seg_vel_distances = 1 - (seg_vel_embeddings * query_seg_vel_embedding');
                seg_vel_ranking = createRankingTable(candidate_segments, seg_vel_distances);
            else
                seg_vel_ranking = [];
            end
            
            if active_modalities.metadata
                seg_meta_distances = 1 - (seg_meta_embeddings * query_seg_meta_embedding');
                seg_meta_ranking = createRankingTable(candidate_segments, seg_meta_distances);
            else
                seg_meta_ranking = [];
            end
            
            % RRF Fusion
            seg_rankings = struct();
            if active_modalities.position, seg_rankings.position = seg_pos_ranking; end
            if active_modalities.joint, seg_rankings.joint = seg_joint_ranking; end
            if active_modalities.orientation, seg_rankings.orientation = seg_orient_ranking; end
            if active_modalities.velocity, seg_rankings.velocity = seg_vel_ranking; end
            if active_modalities.metadata, seg_rankings.metadata = seg_meta_ranking; end
            
            seg_fused_ranking = fuseRankingsRRF(seg_rankings, weights, rrf_k, 'segment_id');
            
            segment_embedding_results{seg_idx} = seg_fused_ranking;
            
            fprintf('  ✓ Best match: %s (RRF: %.6f)\n', ...
                seg_fused_ranking.segment_id{1}, seg_fused_ranking.rrf_score(1));
        end
        
        fprintf('\n✓ Segment-level embeddings completed\n');
    end
    
else
    fprintf('\n=== No segments to process ===\n');
    segment_embedding_results = [];
end

%% SECTION 8: RESULTS & RANKING
%  ========================================================================

fprintf('\n========================================\n');
fprintf('RESULTS SUMMARY\n');
fprintf('========================================\n\n');

fprintf('Query: %s | DTW: %s | Embeddings: 5-way RRF (k=%d)\n', ...
    query_bahn_id, dtw_mode, rrf_k);
fprintf('Candidates: %d trajectories (from cache)\n\n', num_candidates);

% Trajectory-level
fprintf('--- TRAJECTORY LEVEL ---\n');
fprintf('DTW Time: %.2fs\n', dtw_trajectory_time);

embedding_table = fused_ranking;

comparison_table = innerjoin(trajectory_table, embedding_table, ...
    'Keys', 'bahn_id', ...
    'LeftVariables', {'bahn_id', 'dtw_distance'}, ...
    'RightVariables', {'rrf_score'});

dtw_distances = comparison_table.dtw_distance;
emb_scores = -comparison_table.rrf_score;

[rho_spearman, ~] = corr(dtw_distances, emb_scores, 'Type', 'Spearman');
fprintf('Spearman ρ: %.4f\n', rho_spearman);

fprintf('Precision: ');

K = min([top_k_trajectories, height(trajectory_table), height(embedding_table)]);
dtw_top = trajectory_table.bahn_id(1:K);
emb_top = embedding_table.bahn_id(1:K);
prec_k = length(intersect(dtw_top, emb_top)) / K;
fprintf('P@%d=%.3f', K, prec_k);

if height(trajectory_table) >= 10 && height(embedding_table) >= 10
    dtw_top = trajectory_table.bahn_id(1:10);
    emb_top = embedding_table.bahn_id(1:10);
    prec_10 = length(intersect(dtw_top, emb_top)) / 10;
    fprintf(' | P@10=%.3f', prec_10);
end

if height(trajectory_table) >= 5 && height(embedding_table) >= 5
    dtw_top = trajectory_table.bahn_id(1:5);
    emb_top = embedding_table.bahn_id(1:5);
    prec_5 = length(intersect(dtw_top, emb_top)) / 5;
    fprintf(' | P@5=%.3f', prec_5);
end

if height(trajectory_table) >= 3 && height(embedding_table) >= 3
    dtw_top = trajectory_table.bahn_id(1:3);
    emb_top = embedding_table.bahn_id(1:3);
    prec_3 = length(intersect(dtw_top, emb_top)) / 3;
    fprintf(' | P@3=%.3f', prec_3);
end

dtw_top_1 = trajectory_table.bahn_id{1};
emb_top_1 = embedding_table.bahn_id{1};
prec_1 = double(strcmp(dtw_top_1, emb_top_1));
fprintf(' | P@1=%.0f\n', prec_1);

% Segment-level
if num_query_segments > 0 && ~isempty(segment_embedding_results)
    fprintf('\n--- SEGMENT LEVEL (avg over %d segments) ---\n', num_query_segments);
    
    avg_seg_dtw_time = mean(segment_dtw_times);
    fprintf('DTW Time: %.2fs (avg per segment)\n', avg_seg_dtw_time);
    
    seg_rho_all = [];
    seg_prec_k_all = [];
    seg_prec_10_all = [];
    seg_prec_5_all = [];
    seg_prec_3_all = [];
    seg_prec_1_all = [];
    
    for seg_idx = 1:num_query_segments
        seg_dtw_table = segment_results{seg_idx};
        seg_emb_table = segment_embedding_results{seg_idx};
        
        if isempty(seg_dtw_table) || isempty(seg_emb_table)
            continue;
        end
        
        seg_comparison = innerjoin(seg_dtw_table, seg_emb_table, ...
            'Keys', 'segment_id', ...
            'LeftVariables', {'segment_id', 'dtw_distance'}, ...
            'RightVariables', {'rrf_score'});
        
        if height(seg_comparison) >= 3
            seg_emb_scores = -seg_comparison.rrf_score;
            [seg_rho, ~] = corr(seg_comparison.dtw_distance, seg_emb_scores, 'Type', 'Spearman');
            seg_rho_all(end+1) = seg_rho;
        end
        
        K_seg = min([top_k_trajectories, height(seg_dtw_table), height(seg_emb_table)]);
        if K_seg > 0
            seg_dtw_top = seg_dtw_table.segment_id(1:K_seg);
            seg_emb_top = seg_emb_table.segment_id(1:K_seg);
            seg_prec_k_all(end+1) = length(intersect(seg_dtw_top, seg_emb_top)) / K_seg;
        end
        
        if height(seg_dtw_table) >= 10 && height(seg_emb_table) >= 10
            seg_dtw_top = seg_dtw_table.segment_id(1:10);
            seg_emb_top = seg_emb_table.segment_id(1:10);
            seg_prec_10_all(end+1) = length(intersect(seg_dtw_top, seg_emb_top)) / 10;
        end
        
        if height(seg_dtw_table) >= 5 && height(seg_emb_table) >= 5
            seg_dtw_top = seg_dtw_table.segment_id(1:5);
            seg_emb_top = seg_emb_table.segment_id(1:5);
            seg_prec_5_all(end+1) = length(intersect(seg_dtw_top, seg_emb_top)) / 5;
        end
        
        if height(seg_dtw_table) >= 3 && height(seg_emb_table) >= 3
            seg_dtw_top = seg_dtw_table.segment_id(1:3);
            seg_emb_top = seg_emb_table.segment_id(1:3);
            seg_prec_3_all(end+1) = length(intersect(seg_dtw_top, seg_emb_top)) / 3;
        end
        
        if height(seg_dtw_table) >= 1 && height(seg_emb_table) >= 1
            seg_prec_1_all(end+1) = double(strcmp(seg_dtw_table.segment_id{1}, seg_emb_table.segment_id{1}));
        end
    end
    
    if ~isempty(seg_rho_all)
        fprintf('Spearman ρ: %.4f\n', mean(seg_rho_all));
        fprintf('Precision: P@%d=%.3f', top_k_trajectories, mean(seg_prec_k_all));
        
        if ~isempty(seg_prec_10_all)
            fprintf(' | P@10=%.3f', mean(seg_prec_10_all));
        end
        if ~isempty(seg_prec_5_all)
            fprintf(' | P@5=%.3f', mean(seg_prec_5_all));
        end
        if ~isempty(seg_prec_3_all)
            fprintf(' | P@3=%.3f', mean(seg_prec_3_all));
        end
        if ~isempty(seg_prec_1_all)
            fprintf(' | P@1=%.3f\n', mean(seg_prec_1_all));
        end
    end
end

fprintf('\n========================================\n\n');

%% SECTION 9: COVERAGE ANALYSIS
%  ========================================================================

fprintf('\n=== Coverage Analysis ===\n');
fprintf('Analyzing how many embedding top-X are needed for all DTW top-K...\n\n');

% Define K values to analyze
k_values_to_analyze = [10, 50];  % Can be adjusted based on top_k_trajectories

% Filter k_values to only those <= available results
max_k_traj = min([top_k_trajectories, height(trajectory_table), height(embedding_table)]);
k_values_valid = k_values_to_analyze(k_values_to_analyze <= max_k_traj);

if isempty(k_values_valid)
    fprintf('⚠ Not enough candidates for coverage analysis (need at least K=%d)\n', ...
        min(k_values_to_analyze));
    coverage_traj = [];
else
    % Trajectory-level coverage analysis
    fprintf('--- Trajectory Level ---\n');
    coverage_traj = calculateCoverage(trajectory_table, embedding_table, k_values_valid, 'bahn_id');
    
    for i = 1:length(k_values_valid)
        K = k_values_valid(i);
        fprintf('  K=%d: Coverage@%d (%.2fx expansion, Recall@K=%.2f)\n', ...
            K, coverage_traj.coverage_points(i), ...
            coverage_traj.expansion_ratios(i), ...
            coverage_traj.recall_at_k(i));
    end
end

% Segment-level coverage analysis
if num_query_segments > 0 && ~isempty(segment_embedding_results)
    fprintf('\n--- Segment Level ---\n');
    
    segment_coverage_all = cell(num_query_segments, 1);
    valid_segment_count = 0;
    
    for seg_idx = 1:num_query_segments
        seg_dtw_table = segment_results{seg_idx};
        seg_emb_table = segment_embedding_results{seg_idx};
        
        if isempty(seg_dtw_table) || isempty(seg_emb_table)
            continue;
        end
        
        % Get valid k_values for this segment
        max_k_seg = min([top_k_trajectories, height(seg_dtw_table), height(seg_emb_table)]);
        k_values_seg_valid = k_values_to_analyze(k_values_to_analyze <= max_k_seg);
        
        if isempty(k_values_seg_valid)
            continue;
        end
        
        % Calculate coverage for this segment
        seg_coverage = calculateCoverage(seg_dtw_table, seg_emb_table, ...
            k_values_seg_valid, 'segment_id');
        
        segment_coverage_all{seg_idx} = seg_coverage;
        valid_segment_count = valid_segment_count + 1;
    end
    
    % Average segment coverage
    if valid_segment_count > 0
        % Initialize arrays for averaging
        avg_coverage_points = zeros(length(k_values_to_analyze), 1);
        avg_expansion_ratios = zeros(length(k_values_to_analyze), 1);
        avg_recall_at_k = zeros(length(k_values_to_analyze), 1);
        avg_worst_ranks = zeros(length(k_values_to_analyze), 1);
        
        counts = zeros(length(k_values_to_analyze), 1);
        
        % Accumulate across segments
        for seg_idx = 1:num_query_segments
            if isempty(segment_coverage_all{seg_idx})
                continue;
            end
            
            seg_cov = segment_coverage_all{seg_idx};
            
            % Match k_values
            for k_idx = 1:length(k_values_to_analyze)
                K = k_values_to_analyze(k_idx);
                
                % Find this K in segment's k_values
                seg_k_idx = find(seg_cov.k_values == K, 1);
                
                if ~isempty(seg_k_idx) && ~isinf(seg_cov.expansion_ratios(seg_k_idx))
                    avg_coverage_points(k_idx) = avg_coverage_points(k_idx) + seg_cov.coverage_points(seg_k_idx);
                    avg_expansion_ratios(k_idx) = avg_expansion_ratios(k_idx) + seg_cov.expansion_ratios(seg_k_idx);
                    avg_recall_at_k(k_idx) = avg_recall_at_k(k_idx) + seg_cov.recall_at_k(seg_k_idx);
                    avg_worst_ranks(k_idx) = avg_worst_ranks(k_idx) + seg_cov.worst_ranks(seg_k_idx);
                    counts(k_idx) = counts(k_idx) + 1;
                end
            end
        end
        
        % Compute averages
        avg_coverage_points = avg_coverage_points ./ max(counts, 1);
        avg_expansion_ratios = avg_expansion_ratios ./ max(counts, 1);
        avg_recall_at_k = avg_recall_at_k ./ max(counts, 1);
        avg_worst_ranks = avg_worst_ranks ./ max(counts, 1);
        
        % Store in struct
        coverage_seg_avg = struct();
        coverage_seg_avg.k_values = k_values_to_analyze;
        coverage_seg_avg.coverage_points = avg_coverage_points;
        coverage_seg_avg.expansion_ratios = avg_expansion_ratios;
        coverage_seg_avg.recall_at_k = avg_recall_at_k;
        coverage_seg_avg.worst_ranks = avg_worst_ranks;
        
        % Display
        for k_idx = 1:length(k_values_to_analyze)
            if counts(k_idx) > 0
                K = k_values_to_analyze(k_idx);
                fprintf('  K=%d: Avg Coverage@%.1f (%.2fx expansion, Recall@K=%.2f)\n', ...
                    K, avg_coverage_points(k_idx), ...
                    avg_expansion_ratios(k_idx), ...
                    avg_recall_at_k(k_idx));
            end
        end
    else
        coverage_seg_avg = [];
    end
else
    coverage_seg_avg = [];
end

fprintf('\n✓ Coverage analysis completed\n');
fprintf('========================================\n\n');
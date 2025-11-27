%% ========================================================================
%  DTW BASELINE - TRAJECTORY SIMILARITY SEARCH
%  ========================================================================
%  Hierarchical similarity search using Dynamic Time Warping (DTW)
%  
%  Author: Gustavo Barros
%  Date: 26.11.2025
%  
%  ========================================================================

%% SECTION 1: CONFIGURATION & PARAMETERS
%  ========================================================================
%  Set all parameters here
%  
%  NOTE: These can be overridden by external scripts (e.g., experiments)
%  If a variable already exists in workspace, it won't be overwritten.

% === Query Parameters ===
if ~exist('query_bahn_id', 'var')
    query_bahn_id = '1763740056';        % ID der gesuchten Bahn
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

% === Weights (⭐ CRITICAL FOR EXPERIMENTS) ===
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
if ~exist('db_name', 'var')
    db_name = 'robotervermessung';
end

chunk_size = 100;

% ========================================================================
% ⭐ UPDATED: CONFIGURATION SUMMARY
% ========================================================================

fprintf('\n=== Configuration Summary ===\n');
fprintf('Query Trajectory: %s\n', query_bahn_id);
fprintf('DTW Mode: %s\n', dtw_mode);
fprintf('\n');

% Architecture
if use_multi_scale
    fprintf('Architecture: Multi-Scale (%d + %d + %d samples)\n', n_coarse, n_medium, n_fine);
else
    fprintf('Architecture: Single-Scale (%d samples)\n', n_samples_single);
end

fprintf('Normalization: %s\n', norm_strategy);
fprintf('\n');

% SECTION 2: DATABASE CONNECTION
%  ========================================================================

fprintf('\n=== Connecting to Database ===\n');

conn = connectingToPostgres;

if isopen(conn)
    fprintf('✓ Database connection successful\n');
else
    error('✗ Database connection failed');
end

%% SECTION 3: LOAD QUERY TRAJECTORY DATA
%  ========================================================================

fprintf('\n=== Loading Query Trajectory: %s ===\n', query_bahn_id);

% === Load trajectory data for DTW (single mode) ===
if strcmp(dtw_mode, 'position')
    query_data = loadTrajectoryPosition(conn, schema, query_bahn_id);
    query_start = query_data(1, :);
    query_data = query_data - query_start;

    fprintf('✓ DTW data loaded: Position (%d points, x/y/z)\n', size(query_data, 1));
    
elseif strcmp(dtw_mode, 'joint_states')
    query_data = loadTrajectoryJointStates(conn, schema, query_bahn_id);
    
    fprintf('✓ DTW data loaded: Joint States (%d points, 6 joints)\n', size(query_data, 1));
    
else
    error('Invalid dtw_mode. Use "position" or "joint_states"');
end

% === Load trajectory metadata (shared by all modes) ===
query_metadata = loadTrajectoryMetadata(conn, schema, query_bahn_id);

fprintf('  Query length: %.2f mm\n', query_metadata.length);
fprintf('  Query duration: %.2f s\n', query_metadata.duration);
fprintf('  Query segments: %d\n', query_metadata.num_segments);

% ========================================================================
% ⭐ UPDATED: Load metadata features (ALWAYS needed now)
% ========================================================================
% These are needed for:
%   - Hybrid mode: use_metadata_position/joint = true
%   - Separate mode: compute_metadata_embedding = true

fprintf('  Loading metadata features...\n');

query_metadata_query = sprintf(...
    ['SELECT movement_type, min_twist_ist, max_twist_ist, mean_twist_ist, ' ...
     'median_twist_ist, std_twist_ist, min_acceleration_ist, max_acceleration_ist, ' ...
     'mean_acceleration_ist, median_acceleration_ist, std_acceleration_ist ' ...
     'FROM robotervermessung.%s.bahn_metadata ' ...
     'WHERE bahn_id = ''%s'' AND segment_id = ''%s'''], ...
    schema, query_bahn_id, query_bahn_id);

query_meta_result = fetch(conn, query_metadata_query);

% Add to existing query_metadata struct
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

fprintf('\n');

% === Load Query Segments ===
query_segments_data = loadSegments(conn, schema, query_bahn_id, dtw_mode);
num_query_segments = length(query_segments_data);

% Load metadata for each query segment
query_segments_metadata = cell(num_query_segments, 1);
for seg_idx = 1:num_query_segments
    query_segment_id = sprintf('%s_%d', query_bahn_id, seg_idx);
    query_segments_metadata{seg_idx} = loadSegmentMetadata(conn, schema, query_segment_id);
    
    % Load zusätzliche Metadata-Features
    seg_meta_query = sprintf(...
        ['SELECT movement_type, min_twist_ist, max_twist_ist, mean_twist_ist, ' ...
         'median_twist_ist, std_twist_ist, min_acceleration_ist, max_acceleration_ist, ' ...
     'mean_acceleration_ist, median_acceleration_ist, std_acceleration_ist ' ...
         'FROM robotervermessung.%s.bahn_metadata ' ...
         'WHERE segment_id = ''%s'''], ...
        schema, query_segment_id);
    seg_meta_result = fetch(conn, seg_meta_query);
    
    query_segments_metadata{seg_idx}.movement_type = seg_meta_result.movement_type{1};
    query_segments_metadata{seg_idx}.min_twist = seg_meta_result.min_twist_ist;
    query_segments_metadata{seg_idx}.max_twist = seg_meta_result.max_twist_ist;
    query_segments_metadata{seg_idx}.mean_twist = seg_meta_result.mean_twist_ist;
    query_segments_metadata{seg_idx}.median_twist = seg_meta_result.median_twist_ist;
    query_segments_metadata{seg_idx}.std_twist = seg_meta_result.std_twist_ist;
    query_segments_metadata{seg_idx}.min_acceleration = seg_meta_result.min_acceleration_ist;
    query_segments_metadata{seg_idx}.max_acceleration = seg_meta_result.max_acceleration_ist;
    query_segments_metadata{seg_idx}.mean_acceleration = seg_meta_result.mean_acceleration_ist;
    query_segments_metadata{seg_idx}.median_acceleration = seg_meta_result.median_acceleration_ist;
    query_segments_metadata{seg_idx}.std_acceleration = seg_meta_result.std_acceleration_ist;
end

fprintf('✓ Loaded %d query segments with metadata\n', num_query_segments);

%% SECTION 4: DATABASE SAMPLING & CANDIDATE SELECTION
%  ========================================================================

fprintf('\n=== Database Sampling ===\n');

% ========================================================================
% STEP 1: Load full database metadata
% ========================================================================

fprintf('Loading full database metadata...\n');

full_db_query = sprintf(...
    ['SELECT bahn_id, length, duration, movement_type, ' ...
     'min_twist_ist, max_twist_ist, mean_twist_ist, median_twist_ist, std_twist_ist, ' ...
     'min_acceleration_ist, max_acceleration_ist, mean_acceleration_ist, ' ...
     'median_acceleration_ist, std_acceleration_ist ' ...
     'FROM robotervermessung.%s.bahn_metadata ' ...
     'WHERE bahn_id = segment_id ' ...
     'AND bahn_id != ''%s'''], ...
    schema, query_bahn_id);

full_db_metadata = fetch(conn, full_db_query);
num_total_trajectories = height(full_db_metadata);

fprintf('  Total trajectories in database: %d\n', num_total_trajectories);

% ========================================================================
% STEP 2: Random Sampling (if enabled)
% ========================================================================

if use_database_sampling && database_sample_size > 0 && database_sample_size < num_total_trajectories
    fprintf('\nRandom Sampling Enabled\n');
    fprintf('  Target Sample Size: %d (%.1f%% of full database)\n', ...
        database_sample_size, 100 * database_sample_size / num_total_trajectories);
    fprintf('  Random Seed: %d (for reproducibility)\n', random_seed);
    
    % Set random seed
    rng(random_seed);
    
    % Random sample
    sample_indices = randperm(num_total_trajectories, database_sample_size);
    candidate_metadata = full_db_metadata(sample_indices, :);
    
    fprintf('  ✓ Sampled %d trajectories\n', height(candidate_metadata));
    
    % ========================================================================
    % STEP 3: Verify sample distribution (optional)
    % ========================================================================
    
    fprintf('\n  --- Distribution Comparison ---\n');
    fprintf('  Metric          | Full DB | Sample  | Difference\n');
    fprintf('  ----------------|---------|---------|------------\n');
    
    % Length
    full_mean_length = mean(full_db_metadata.length);
    sample_mean_length = mean(candidate_metadata.length);
    fprintf('  Mean Length     | %7.1f | %7.1f | %+6.1f mm\n', ...
        full_mean_length, sample_mean_length, sample_mean_length - full_mean_length);
    
    full_std_length = std(full_db_metadata.length);
    sample_std_length = std(candidate_metadata.length);
    fprintf('  Std Length      | %7.1f | %7.1f | %+6.1f mm\n', ...
        full_std_length, sample_std_length, sample_std_length - full_std_length);
    
    % Duration
    full_mean_duration = mean(full_db_metadata.duration);
    sample_mean_duration = mean(candidate_metadata.duration);
    fprintf('  Mean Duration   | %7.2f | %7.2f | %+6.2f s\n', ...
        full_mean_duration, sample_mean_duration, sample_mean_duration - full_mean_duration);
    
    fprintf('\n');
    
else
    % No sampling - use full database
    candidate_metadata = full_db_metadata;
    fprintf('  Using full database (no sampling)\n');
end

num_candidates = height(candidate_metadata);
fprintf('✓ Final candidate set: %d trajectories\n\n', num_candidates);

% ========================================================================
% STEP 4: Segment-level candidates (from sampled trajectories only)
% ========================================================================

fprintf('=== Loading Segments (from sampled trajectories) ===\n');

all_segments_metadata = cell(num_query_segments, 1);

for seg_idx = 1:num_query_segments
    query_segment_id = sprintf('%s_%d', query_bahn_id, seg_idx);
    
    if use_database_sampling && database_sample_size > 0
        % Only segments from sampled trajectories
        sampled_bahn_ids = candidate_metadata.bahn_id;
        bahn_id_list = sprintf('''%s''', strjoin(sampled_bahn_ids, ''','''));
        
        seg_filter_query = sprintf(...
            ['SELECT segment_id, bahn_id, length, duration, movement_type, ' ...
             'min_twist_ist, max_twist_ist, mean_twist_ist, median_twist_ist, std_twist_ist, ' ...
             'min_acceleration_ist, max_acceleration_ist, mean_acceleration_ist, ' ...
             'median_acceleration_ist, std_acceleration_ist ' ...
             'FROM robotervermessung.%s.bahn_metadata ' ...
             'WHERE segment_id != bahn_id ' ...
             'AND bahn_id IN (%s) ' ...
             'AND segment_id != ''%s'''], ...
            schema, bahn_id_list, query_segment_id);
    else
        % All segments (full database)
        seg_filter_query = sprintf(...
            ['SELECT segment_id, bahn_id, length, duration, movement_type, ' ...
             'min_twist_ist, max_twist_ist, mean_twist_ist, median_twist_ist, std_twist_ist, ' ...
             'min_acceleration_ist, max_acceleration_ist, mean_acceleration_ist, ' ...
             'median_acceleration_ist, std_acceleration_ist ' ...
             'FROM robotervermessung.%s.bahn_metadata ' ...
             'WHERE segment_id != bahn_id ' ...
             'AND bahn_id IN (%s) ' ...
             'AND segment_id != ''%s'''], ...
            schema, bahn_id_list, query_segment_id);
    end
    
    all_segments_metadata{seg_idx} = fetch(conn, seg_filter_query);
    fprintf('  Segment %d: %d candidates\n', seg_idx, height(all_segments_metadata{seg_idx}));
end

fprintf('\n');


%% SECTION 5: DTW COMPUTATION - TRAJECTORY LEVEL

fprintf('\n=== Computing DTW Distances (Trajectory Level) ===\n');
fprintf('  Total candidates: %d\n', num_candidates);

% ========================================================================
% ⭐ BATCH LOAD ALL CANDIDATES AT ONCE (1 query instead of 1000+!)
% ========================================================================

fprintf('  Loading all %d candidates in batch...\n', num_candidates);

candidate_ids = candidate_metadata.bahn_id;
candidate_trajectories = batchLoadTrajectoriesChunked(conn, schema, candidate_ids, dtw_mode, chunk_size);

fprintf('    ✓ Batch loading completed\n');

tic;

% ========================================================================
% PHASE 1: LB_Kim Pre-filtering (Ultra-Fast, O(1))
% ========================================================================
fprintf('Phase 1: LB_Kim filtering...\n');

lb_kim_distances = zeros(num_candidates, 1);

for i = 1:num_candidates
    % ✅ Data already in memory!
    candidate_data = candidate_trajectories{i};
    
    lb_kim_distances(i) = LB_Kim(query_data, candidate_data, dtw_mode, use_rotation_alignment, normalize_dtw);
end

% Sort by LB_Kim distance
[~, kim_order] = sort(lb_kim_distances);
kim_keep_count = round(num_candidates * lb_kim_keep_ratio);

candidates_after_kim = candidate_metadata(kim_order(1:kim_keep_count), :);

% ✅ Keep corresponding trajectory data (no reloading!)
candidate_trajectories_after_kim = candidate_trajectories(kim_order(1:kim_keep_count));

fprintf('  ✓ LB_Kim: %d → %d candidates (%.1f%% pruned)\n', ...
    num_candidates, kim_keep_count, 100*(1 - lb_kim_keep_ratio));

% ========================================================================
% PHASE 2: LB_Keogh Pre-filtering (Precise, O(n))
% ========================================================================
fprintf('Phase 2: LB_Keogh filtering...\n');

% ⭐ Calculate efficient keep count
keogh_keep_count = min(lb_keogh_candidates, kim_keep_count);

fprintf('  Candidates from LB_Kim: %d\n', kim_keep_count);
fprintf('  Target for LB_Keogh: %d (fixed)\n', keogh_keep_count);

% ⭐ Early exit if already at target
if kim_keep_count <= keogh_keep_count
    fprintf('  ⚠ LB_Kim already filtered to target size - skipping LB_Keogh\n');
    candidates_after_keogh = candidates_after_kim;
    candidate_trajectories_after_keogh = candidate_trajectories_after_kim;
    fprintf('  ✓ Kept all %d candidates\n', kim_keep_count);
else
    % Compute LB_Keogh distances
    lb_keogh_distances = zeros(kim_keep_count, 1);
    
    for i = 1:kim_keep_count
        % ✅ Data already in memory!
        candidate_data = candidate_trajectories_after_kim{i};
        
        lb_keogh_distances(i) = LB_Keogh(query_data, candidate_data, ...
            cdtw_window, dtw_mode, use_rotation_alignment, normalize_dtw);
    end
    
    % Sort and keep top-K (⭐ fixed number)
    [~, keogh_order] = sort(lb_keogh_distances);
    candidates_after_keogh = candidates_after_kim(keogh_order(1:keogh_keep_count), :);
    
    % ✅ Keep corresponding trajectory data (no reloading!)
    candidate_trajectories_after_keogh = candidate_trajectories_after_kim(keogh_order(1:keogh_keep_count));
    
    fprintf('  ✓ LB_Keogh: %d → %d candidates (%.1f%% pruned)\n', ...
        kim_keep_count, keogh_keep_count, ...
        100*(1 - keogh_keep_count/kim_keep_count));
end

% Total pruning stats
total_pruning = 1 - (keogh_keep_count / num_candidates);
fprintf('  ✓ Total Pruning: %.1f%% (from %d to %d before DTW)\n', ...
    total_pruning * 100, num_candidates, keogh_keep_count);

% ========================================================================
% PHASE 3: DTW Refinement (Exact, O(n²))
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
    
    % ✅ Data already in memory!
    candidate_data = candidate_trajectories_after_keogh{i};
    
    % Best-so-far for early abandoning
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

%% SECTION 5.5: DTW COMPUTATION - SEGMENT LEVEL

fprintf('\n=== Computing DTW Distances (Segment Level) ===\n');

% Load query segments
query_segments = loadSegments(conn, schema, query_bahn_id, dtw_mode);
num_query_segments = length(query_segments);

fprintf('  Query segments: %d\n', num_query_segments);

% ===== FÜR JEDES QUERY-SEGMENT: UNABHÄNGIGE SUCHE =====
segment_results = cell(num_query_segments, 1);
segment_dtw_times = zeros(num_query_segments, 1);

for seg_idx = 1:num_query_segments
    query_segment_id = sprintf('%s_%d', query_bahn_id, seg_idx);
    query_segment_data = query_segments{seg_idx};
    
    fprintf('\n  Processing Query Segment %d/%d (ID: %s)...\n', seg_idx, num_query_segments, query_segment_id);
    
    % Get candidate segments for this query segment
    candidate_segments = all_segments_metadata{seg_idx};
    num_candidate_segments = height(candidate_segments);
    
    if num_candidate_segments == 0
        fprintf('    No candidate segments found - skipping\n');
        segment_results{seg_idx} = [];
        continue;
    end
    
    fprintf('    Total segment candidates: %d\n', num_candidate_segments);
    
    % ====================================================================
    % ⭐ Batch Load ALL segment candidates at once (chunked)
    % ====================================================================
    
    fprintf('    Loading all %d segment candidates in batch...\n', num_candidate_segments);
    
    seg_candidate_ids = candidate_segments.segment_id;
    seg_candidate_segments = batchLoadSegmentsChunked(conn, schema, seg_candidate_ids, dtw_mode, chunk_size);
    
    fprintf('      ✓ Batch loading completed\n');
    
    % ====================================================================
    % LB_Kim Pre-filtering für Segmente
    % ====================================================================
    
    fprintf('    Step 1: LB_Kim screening...\n');

    tic;
    
    segment_lb_kim_distances = zeros(num_candidate_segments, 1);
    
    for cand_seg_idx = 1:num_candidate_segments
        % ✅ Data already in memory!
        candidate_segment_data = seg_candidate_segments{cand_seg_idx};
              
        segment_lb_kim_distances(cand_seg_idx) = LB_Kim(query_segment_data, ...
            candidate_segment_data, dtw_mode, use_rotation_alignment, normalize_dtw);
    end
    
    % Sort by LB_Kim
    [~, seg_kim_order] = sort(segment_lb_kim_distances);
    seg_kim_keep_count = round(num_candidate_segments * lb_kim_keep_ratio);
    
    candidates_after_kim_seg = candidate_segments(seg_kim_order(1:seg_kim_keep_count), :);
    
    % ✅ Keep corresponding segment data
    seg_candidate_segments_after_kim = seg_candidate_segments(seg_kim_order(1:seg_kim_keep_count));
    
    fprintf('      ✓ LB_Kim: %d → %d candidates (%.1f%% pruned)\n', ...
        num_candidate_segments, seg_kim_keep_count, 100*(1 - lb_kim_keep_ratio));
    
    % ====================================================================
    % LB_Keogh Pre-filtering für Segmente
    % ====================================================================
    fprintf('    Step 2: LB_Keogh refinement...\n');
    
    % ⭐ Calculate efficient keep count
    seg_keogh_keep_count = min(lb_keogh_candidates, seg_kim_keep_count);
    
    fprintf('      Candidates from LB_Kim: %d\n', seg_kim_keep_count);
    fprintf('      Target for LB_Keogh: %d (fixed)\n', seg_keogh_keep_count);
    
    % ⭐ Early exit if already at target
    if seg_kim_keep_count <= seg_keogh_keep_count
        fprintf('      ⚠ LB_Kim already filtered to target size - skipping LB_Keogh\n');
        candidates_after_keogh_seg = candidates_after_kim_seg;
        seg_candidate_segments_after_keogh = seg_candidate_segments_after_kim;
        fprintf('      ✓ Kept all %d segment candidates\n', seg_kim_keep_count);
    else
        % Compute LB_Keogh distances
        segment_lb_keogh_distances = zeros(seg_kim_keep_count, 1);
        
        for cand_seg_idx = 1:seg_kim_keep_count
            % ✅ Data already in memory!
            candidate_segment_data = seg_candidate_segments_after_kim{cand_seg_idx};
            
            segment_lb_keogh_distances(cand_seg_idx) = LB_Keogh(query_segment_data, ...
                candidate_segment_data, cdtw_window, dtw_mode, use_rotation_alignment, normalize_dtw);
        end
        
        % Sort and keep top-K (⭐ fixed number)
        [~, seg_keogh_order] = sort(segment_lb_keogh_distances);
        candidates_after_keogh_seg = candidates_after_kim_seg(seg_keogh_order(1:seg_keogh_keep_count), :);
        
        % ✅ Keep corresponding segment data
        seg_candidate_segments_after_keogh = seg_candidate_segments_after_kim(seg_keogh_order(1:seg_keogh_keep_count));
        
        fprintf('      ✓ LB_Keogh: %d → %d candidates (%.1f%% pruned)\n', ...
            seg_kim_keep_count, seg_keogh_keep_count, ...
            100*(1 - seg_keogh_keep_count/seg_kim_keep_count));
    end
    
    % Total segment pruning
    seg_total_pruning = 1 - (seg_keogh_keep_count / num_candidate_segments);
    fprintf('      ✓ Total Pruning: %.1f%% (from %d to %d before DTW)\n', ...
        seg_total_pruning * 100, num_candidate_segments, seg_keogh_keep_count);
    
    % ====================================================================
    % DTW mit Early Abandoning
    % ====================================================================
    
    fprintf('    Step 3: DTW computation...\n');
    
    seg_dtw_limit = min(top_k_trajectories, seg_keogh_keep_count);
    segment_dtw_distances = inf(seg_dtw_limit, 1);
    
    for cand_seg_idx = 1:seg_dtw_limit
        candidate_segment_id = candidates_after_keogh_seg.segment_id{cand_seg_idx};
        
        % ✅ Data already in memory!
        candidate_segment_data = seg_candidate_segments_after_keogh{cand_seg_idx};
        
        if ~isempty(candidate_segment_data)
            candidate_segment_data = candidate_segment_data - candidate_segment_data(1, :);
        end
        
        % Best-so-far
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
    
    % Store results
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

%% SECTION 6: COMPUTE EMBEDDINGS (TRAJECTORY LEVEL)

fprintf('\n=== Computing Embeddings (Trajectory Level) ===\n');

% === Check which modalities are active (weight > 0) ===
active_modalities.position = weights(1) > 0;       % Position
active_modalities.joint = weights(2) > 0;          % Joint
active_modalities.orientation = weights(3) > 0;    % Orientation
active_modalities.velocity = weights(4) > 0;       % Velocity
active_modalities.metadata = weights(5) > 0;       % Metadata

fprintf('  Active modalities: ');
active_names = {};
if active_modalities.position, active_names{end+1} = 'Pos'; end
if active_modalities.joint, active_names{end+1} = 'Joint'; end
if active_modalities.orientation, active_names{end+1} = 'Orient'; end
if active_modalities.velocity, active_names{end+1} = 'Vel'; end
if active_modalities.metadata, active_names{end+1} = 'Meta'; end
fprintf('%s\n', strjoin(active_names, ', '));

% === Query Embeddings (nur für aktive Modalitäten) ===
fprintf('  Computing query embeddings...\n');

% Position (⭐ nur wenn weight > 0)
if active_modalities.position
    query_pos_data = loadTrajectoryPosition(conn, schema, query_bahn_id);
    query_pos_embedding = computePositionEmbedding(query_pos_data, n_coarse, n_fine);
else
    query_pos_embedding = [];
end

% Joint States (⭐ nur wenn weight > 0)
if active_modalities.joint
    query_joint_data = loadTrajectoryJointStates(conn, schema, query_bahn_id);
    query_joint_embedding = computeJointEmbedding(query_joint_data, n_coarse, n_fine);
else
    query_joint_embedding = [];
end

% Orientation (⭐ nur wenn weight > 0)
if active_modalities.orientation
    query_orient_data = loadTrajectoryOrientation(conn, schema, query_bahn_id);
    query_orient_embedding = computeOrientationEmbedding(query_orient_data, n_coarse, n_fine);
else
    query_orient_embedding = [];
end

% Velocity (⭐ nur wenn weight > 0)
if active_modalities.velocity
    if ~active_modalities.position
        query_pos_data = loadTrajectoryPosition(conn, schema, query_bahn_id);
    end
    query_vel_embedding = computeVelocityEmbedding(query_pos_data, n_coarse, n_fine);
else
    query_vel_embedding = [];
end

% Metadata (⭐ nur wenn weight > 0)
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

% === Batch Load Candidate Data (⭐ nur für aktive Modalitäten) ===
fprintf('  Batch loading candidate data...\n');

candidate_ids = candidate_metadata.bahn_id;

if active_modalities.position || active_modalities.velocity
    pos_data = batchLoadTrajectoriesChunked(conn, schema, candidate_ids, 'position', chunk_size);
else
    pos_data = [];
end

if active_modalities.joint
    joint_data = batchLoadTrajectoriesChunked(conn, schema, candidate_ids, 'joint_states', chunk_size);
else
    joint_data = [];
end

if active_modalities.orientation
    orient_data = batchLoadTrajectoriesChunked(conn, schema, candidate_ids, 'orientation', chunk_size);
else
    orient_data = [];
end

fprintf('    ✓ Batch loading completed\n');

% === Compute Candidate Embeddings (⭐ nur für aktive Modalitäten) ===
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
    % Position (⭐ nur wenn weight > 0)
    if active_modalities.position
        pos_embeddings(i, :) = computePositionEmbedding(pos_data{i}, n_coarse, n_fine);
    end
    
    % Joint (⭐ nur wenn weight > 0)
    if active_modalities.joint
        joint_embeddings(i, :) = computeJointEmbedding(joint_data{i}, n_coarse, n_fine);
    end
    
    % Orientation (⭐ nur wenn weight > 0)
    if active_modalities.orientation
        orient_embeddings(i, :) = computeOrientationEmbedding(orient_data{i}, n_coarse, n_fine);
    end
    
    % Velocity (⭐ nur wenn weight > 0)
    if active_modalities.velocity
        vel_embeddings(i, :) = computeVelocityEmbedding(pos_data{i}, n_coarse, n_fine);
    end
    
    % Metadata (⭐ nur wenn weight > 0)
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

% === Compute Cosine Distances (⭐ nur für aktive Modalitäten) ===
fprintf('  Computing cosine distances...\n');

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

% === RRF Fusion (5-way, ⭐ nur aktive Modalitäten) ===
fprintf('\n=== RRF Fusion (Active Modalities) ===\n');

rankings = struct();
if active_modalities.position, rankings.position = pos_ranking; end
if active_modalities.joint, rankings.joint = joint_ranking; end
if active_modalities.orientation, rankings.orientation = orient_ranking; end
if active_modalities.velocity, rankings.velocity = vel_ranking; end
if active_modalities.metadata, rankings.metadata = meta_ranking; end

fused_ranking = fuseRankingsRRF(rankings, weights, rrf_k, 'bahn_id');

fprintf('  ✓ Best match: %s (RRF: %.6f)\n', fused_ranking.bahn_id{1}, fused_ranking.rrf_score(1));


%% SECTION 6.5: COMPUTE EMBEDDINGS (SEGMENT LEVEL)

if num_query_segments > 0
    fprintf('\n=== Computing Embeddings (Segment Level) ===\n');
    
    segment_embedding_results = cell(num_query_segments, 1);
    
    % ===== FÜR JEDES QUERY-SEGMENT: EMBEDDING-BERECHNUNG =====
    for seg_idx = 1:num_query_segments
        query_segment_id = sprintf('%s_%d', query_bahn_id, seg_idx);
        
        fprintf('\n--- Query Segment %d/%d (ID: %s) ---\n', seg_idx, num_query_segments, query_segment_id);
        
        % Get candidate segments for this query segment
        candidate_segments = all_segments_metadata{seg_idx};
        num_candidate_segments = height(candidate_segments);
        
        if num_candidate_segments == 0
            fprintf('  No candidate segments - skipping\n');
            segment_embedding_results{seg_idx} = [];
            continue;
        end
        
        fprintf('  Candidates: %d segments\n', num_candidate_segments);
        fprintf('  Active modalities: %s\n', strjoin(active_names, ', '));
        
        % ================================================================
        % QUERY SEGMENT EMBEDDINGS (⭐ nur für aktive Modalitäten)
        % ================================================================
        fprintf('  Computing query segment embeddings...\n');
        
        % Position (⭐ nur wenn weight > 0)
        if active_modalities.position
            query_seg_pos_data = loadTrajectoryPosition(conn, schema, query_segment_id);
            query_seg_pos_embedding = computePositionEmbedding(query_seg_pos_data, n_coarse, n_fine);
        else
            query_seg_pos_embedding = [];
        end
        
        % Joint States (⭐ nur wenn weight > 0)
        if active_modalities.joint
            query_seg_joint_data = loadTrajectoryJointStates(conn, schema, query_segment_id);
            query_seg_joint_embedding = computeJointEmbedding(query_seg_joint_data, n_coarse, n_fine);
        else
            query_seg_joint_embedding = [];
        end
        
        % Orientation (⭐ nur wenn weight > 0)
        if active_modalities.orientation
            query_seg_orient_data = loadTrajectoryOrientation(conn, schema, query_segment_id);
            query_seg_orient_embedding = computeOrientationEmbedding(query_seg_orient_data, n_coarse, n_fine);
        else
            query_seg_orient_embedding = [];
        end
        
        % Velocity (⭐ nur wenn weight > 0)
        if active_modalities.velocity
            if ~active_modalities.position
                query_seg_pos_data = loadTrajectoryPosition(conn, schema, query_segment_id);
            end
            query_seg_vel_embedding = computeVelocityEmbedding(query_seg_pos_data, n_coarse, n_fine);
        else
            query_seg_vel_embedding = [];
        end
        
        % Metadata (⭐ nur wenn weight > 0)
        if active_modalities.metadata
            query_seg_meta = query_segments_metadata{seg_idx};
            query_seg_meta_embedding = computeMetadataEmbedding(query_seg_meta);
        else
            query_seg_meta_embedding = [];
        end
        
        fprintf('    ✓ Query embeddings computed: ');
        seg_dims = {};
        if active_modalities.position, seg_dims{end+1} = sprintf('Pos(%d)', length(query_seg_pos_embedding)); end
        if active_modalities.joint, seg_dims{end+1} = sprintf('Joint(%d)', length(query_seg_joint_embedding)); end
        if active_modalities.orientation, seg_dims{end+1} = sprintf('Orient(%d)', length(query_seg_orient_embedding)); end
        if active_modalities.velocity, seg_dims{end+1} = sprintf('Vel(%d)', length(query_seg_vel_embedding)); end
        if active_modalities.metadata, seg_dims{end+1} = sprintf('Meta(%d)', length(query_seg_meta_embedding)); end
        fprintf('%s\n', strjoin(seg_dims, ', '));
        
        % ================================================================
        % BATCH LOAD CANDIDATE SEGMENT DATA (⭐ nur für aktive Modalitäten)
        % ================================================================
        fprintf('  Batch loading candidate segment data...\n');
        
        seg_candidate_ids = candidate_segments.segment_id;
        
        if active_modalities.position || active_modalities.velocity
            seg_pos_data = batchLoadSegmentsChunked(conn, schema, seg_candidate_ids, 'position', chunk_size);
        else
            seg_pos_data = [];
        end
        
        if active_modalities.joint
            seg_joint_data = batchLoadSegmentsChunked(conn, schema, seg_candidate_ids, 'joint_states', chunk_size);
        else
            seg_joint_data = [];
        end
        
        if active_modalities.orientation
            seg_orient_data = batchLoadSegmentsChunked(conn, schema, seg_candidate_ids, 'orientation', chunk_size);
        else
            seg_orient_data = [];
        end
        
        fprintf('    ✓ Batch loading completed\n');
        
        % ================================================================
        % COMPUTE CANDIDATE SEGMENT EMBEDDINGS (⭐ nur für aktive Modalitäten)
        % ================================================================
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
            % Position (⭐ nur wenn weight > 0)
            if active_modalities.position
                seg_pos_embeddings(i, :) = computePositionEmbedding(seg_pos_data{i}, n_coarse, n_fine);
            end
            
            % Joint (⭐ nur wenn weight > 0)
            if active_modalities.joint
                seg_joint_embeddings(i, :) = computeJointEmbedding(seg_joint_data{i}, n_coarse, n_fine);
            end
            
            % Orientation (⭐ nur wenn weight > 0)
            if active_modalities.orientation
                seg_orient_embeddings(i, :) = computeOrientationEmbedding(seg_orient_data{i}, n_coarse, n_fine);
            end
            
            % Velocity (⭐ nur wenn weight > 0)
            if active_modalities.velocity
                seg_vel_embeddings(i, :) = computeVelocityEmbedding(seg_pos_data{i}, n_coarse, n_fine);
            end
            
            % Metadata (⭐ nur wenn weight > 0)
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
        
        fprintf('  ✓ All segment embeddings computed\n');
        
        % ================================================================
        % COMPUTE COSINE DISTANCES (⭐ nur für aktive Modalitäten)
        % ================================================================
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
        
        fprintf('  ✓ Rankings created\n');
        
        % ================================================================
        % RRF FUSION (5-way, ⭐ nur aktive Modalitäten)
        % ================================================================
        fprintf('  RRF Fusion (Active Modalities)...\n');
        
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
else
    fprintf('\n=== No segments to process ===\n');
    segment_embedding_results = [];
end

%% SECTION 7: RESULTS & RANKING
%  ========================================================================

fprintf('\n========================================\n');
fprintf('RESULTS SUMMARY\n');
fprintf('========================================\n\n');

fprintf('Query: %s | DTW: %s | Embeddings: 5-way RRF (k=%d)\n', ...
    query_bahn_id, dtw_mode, rrf_k);
fprintf('Database: %d/%d trajectories (%.1f%%)\n\n', ...
    num_candidates, num_total_trajectories, 100*num_candidates/num_total_trajectories);

% =========================================================================
% TRAJECTORY-LEVEL
% =========================================================================

fprintf('--- TRAJECTORY LEVEL ---\n');
fprintf('DTW Time: %.2fs\n', dtw_trajectory_time);

% Use fused ranking as embedding table
embedding_table = fused_ranking;

% Join tables for correlation
comparison_table = innerjoin(trajectory_table, embedding_table, ...
    'Keys', 'bahn_id', ...
    'LeftVariables', {'bahn_id', 'dtw_distance'}, ...
    'RightVariables', {'rrf_score'});

dtw_distances = comparison_table.dtw_distance;
emb_scores = -comparison_table.rrf_score;  % Invert (higher RRF = better)

% Spearman correlation
[rho_spearman, ~] = corr(dtw_distances, emb_scores, 'Type', 'Spearman');
fprintf('Spearman ρ: %.4f\n', rho_spearman);

% Precision@K for multiple K values
fprintf('Precision: ');

% P@K (K = top_k_trajectories)
K = min([top_k_trajectories, height(trajectory_table), height(embedding_table)]);
dtw_top = trajectory_table.bahn_id(1:K);
emb_top = embedding_table.bahn_id(1:K);
prec_k = length(intersect(dtw_top, emb_top)) / K;
fprintf('P@%d=%.3f', K, prec_k);

% P@10
if height(trajectory_table) >= 10 && height(embedding_table) >= 10
    dtw_top = trajectory_table.bahn_id(1:10);
    emb_top = embedding_table.bahn_id(1:10);
    prec_10 = length(intersect(dtw_top, emb_top)) / 10;
    fprintf(' | P@10=%.3f', prec_10);
end

% P@5
if height(trajectory_table) >= 5 && height(embedding_table) >= 5
    dtw_top = trajectory_table.bahn_id(1:5);
    emb_top = embedding_table.bahn_id(1:5);
    prec_5 = length(intersect(dtw_top, emb_top)) / 5;
    fprintf(' | P@5=%.3f', prec_5);
end

% P@3
if height(trajectory_table) >= 3 && height(embedding_table) >= 3
    dtw_top = trajectory_table.bahn_id(1:3);
    emb_top = embedding_table.bahn_id(1:3);
    prec_3 = length(intersect(dtw_top, emb_top)) / 3;
    fprintf(' | P@3=%.3f', prec_3);
end

% P@1
dtw_top_1 = trajectory_table.bahn_id{1};
emb_top_1 = embedding_table.bahn_id{1};
prec_1 = double(strcmp(dtw_top_1, emb_top_1));
fprintf(' | P@1=%.0f\n', prec_1);

% =========================================================================
% SEGMENT-LEVEL
% =========================================================================

if num_query_segments > 0 && ~isempty(segment_embedding_results)
    fprintf('\n--- SEGMENT LEVEL (avg over %d segments) ---\n', num_query_segments);
    
    % Average DTW time
    avg_seg_dtw_time = mean(segment_dtw_times);
    fprintf('DTW Time: %.2fs (avg per segment)\n', avg_seg_dtw_time);
    
    % Collect metrics across all segments
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
        
        % Correlation
        seg_comparison = innerjoin(seg_dtw_table, seg_emb_table, ...
            'Keys', 'segment_id', ...
            'LeftVariables', {'segment_id', 'dtw_distance'}, ...
            'RightVariables', {'rrf_score'});
        
        if height(seg_comparison) >= 3
            seg_emb_scores = -seg_comparison.rrf_score;
            [seg_rho, ~] = corr(seg_comparison.dtw_distance, seg_emb_scores, 'Type', 'Spearman');
            seg_rho_all(end+1) = seg_rho;
        end
        
        % P@K
        K_seg = min([top_k_trajectories, height(seg_dtw_table), height(seg_emb_table)]);
        if K_seg > 0
            seg_dtw_top = seg_dtw_table.segment_id(1:K_seg);
            seg_emb_top = seg_emb_table.segment_id(1:K_seg);
            seg_prec_k_all(end+1) = length(intersect(seg_dtw_top, seg_emb_top)) / K_seg;
        end
        
        % P@10
        if height(seg_dtw_table) >= 10 && height(seg_emb_table) >= 10
            seg_dtw_top = seg_dtw_table.segment_id(1:10);
            seg_emb_top = seg_emb_table.segment_id(1:10);
            seg_prec_10_all(end+1) = length(intersect(seg_dtw_top, seg_emb_top)) / 10;
        end
        
        % P@5
        if height(seg_dtw_table) >= 5 && height(seg_emb_table) >= 5
            seg_dtw_top = seg_dtw_table.segment_id(1:5);
            seg_emb_top = seg_emb_table.segment_id(1:5);
            seg_prec_5_all(end+1) = length(intersect(seg_dtw_top, seg_emb_top)) / 5;
        end
        
        % P@3
        if height(seg_dtw_table) >= 3 && height(seg_emb_table) >= 3
            seg_dtw_top = seg_dtw_table.segment_id(1:3);
            seg_emb_top = seg_emb_table.segment_id(1:3);
            seg_prec_3_all(end+1) = length(intersect(seg_dtw_top, seg_emb_top)) / 3;
        end
        
        % P@1
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





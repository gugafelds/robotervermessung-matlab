function dtw_cache = precomputeDTW(data_cache, query_ids, config)
% PRECOMPUTEDTW - Pre-compute DTW for all queries in all modes
%
%   This function computes DTW distances ONCE for all query trajectories
%   in both Position and Joint-States modes, storing results for reuse
%   across multiple experiments with different embedding configurations.
%
%   INPUTS:
%       data_cache - Pre-loaded data from loadDataExperiment.m
%       query_ids  - Cell array of query trajectory IDs
%       config     - Configuration struct with DTW parameters
%
%   OUTPUT:
%       dtw_cache - Struct containing all pre-computed DTW results:
%           .(query_field_name).(dtw_mode)
%               .trajectory_ranking - Table with DTW distances
%               .segment_rankings   - Cell array of segment rankings
%               .segment_dtw_times  - Array of segment computation times
%               .dtw_time          - Total DTW time for this query/mode
%               .num_candidates    - Number of candidate trajectories
%
%   PERFORMANCE:
%       For 4 queries × 2 modes × 500 candidates × 4 segments:
%       - Time: ~10-15 minutes (one-time cost)
%       - Saves: ~3-4 hours over 128 experiments (62% reduction!)
%
%   Author: Gustavo Barros
%   Date: 28.11.2025

%% ========================================================================
%  VALIDATION & SETUP
%  ========================================================================

fprintf('\n========================================\n');
fprintf('PRE-COMPUTING DTW FOR ALL QUERIES\n');
fprintf('========================================\n\n');

% Extract config parameters
if isfield(config, 'top_k_trajectories')
    top_k_trajectories = config.top_k_trajectories;
else
    top_k_trajectories = 50;
end

if isfield(config, 'lb_kim_keep_ratio')
    lb_kim_keep_ratio = config.lb_kim_keep_ratio;
else
    lb_kim_keep_ratio = 0.40;
end

if isfield(config, 'lb_keogh_candidates')
    lb_keogh_candidates = config.lb_keogh_candidates;
else
    lb_keogh_candidates = 100;
end

if isfield(config, 'cdtw_window')
    cdtw_window = config.cdtw_window;
else
    cdtw_window = 0.10;
end

if isfield(config, 'normalize_dtw')
    normalize_dtw = config.normalize_dtw;
else
    normalize_dtw = false;
end

if isfield(config, 'use_rotation_alignment')
    use_rotation_alignment = config.use_rotation_alignment;
else
    use_rotation_alignment = false;
end

if isfield(config, 'ground_truth_map')
    ground_truth_map = config.ground_truth_map;
    has_ground_truth = true;
else
    has_ground_truth = false;
end

% DTW modes to compute
dtw_modes = {'position', 'joint_states'};

num_queries = length(query_ids);
num_modes = length(dtw_modes);
total_combinations = num_queries * num_modes;

fprintf('Configuration:\n');
fprintf('  Queries: %d (%s)\n', num_queries, strjoin(query_ids, ', '));
fprintf('  Modes: %s\n', strjoin(dtw_modes, ', '));
fprintf('  Total combinations: %d\n', total_combinations);
fprintf('  Top-K: %d\n', top_k_trajectories);
fprintf('  LB_Kim keep ratio: %.2f\n', lb_kim_keep_ratio);
fprintf('  LB_Keogh candidates: %d\n\n', lb_keogh_candidates);

dtw_cache = struct();
overall_tic = tic;

%% ========================================================================
%  MAIN LOOP: FOR EACH QUERY × EACH MODE
%  ========================================================================

counter = 0;

for q_idx = 1:num_queries
    query_id = query_ids{q_idx};
    query_field_name = ['q_' strrep(query_id, '-', '_')];
    
    % Check if query exists in cache
    if ~isfield(data_cache.queries, query_field_name)
        warning('Query %s not found in data_cache - skipping', query_id);
        continue;
    end
    
    query_cache = data_cache.queries.(query_field_name);
    query_metadata = query_cache.metadata;
    num_query_segments = query_metadata.num_segments;
    
    fprintf('╔════════════════════════════════════════════════════════════════╗\n');
    fprintf('║  QUERY %d/%d: %s                                        ║\n', q_idx, num_queries, query_id);
    fprintf('║  Segments: %d                                                   ║\n', num_query_segments);
    fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');
    
    dtw_cache.(query_field_name) = struct();
    
    % ====================================================================
    % FOR EACH DTW MODE (Position, Joint States)
    % ====================================================================
    
    for mode_idx = 1:num_modes
        dtw_mode = dtw_modes{mode_idx};
        counter = counter + 1;
        
        fprintf('--- [%d/%d] Mode: %s ---\n', counter, total_combinations, dtw_mode);
        
        mode_tic = tic;
        
        % ================================================================
        % STEP 1: Load Query Data from Cache
        % ================================================================
        
        if strcmp(dtw_mode, 'position')
            query_data = query_cache.position;
        elseif strcmp(dtw_mode, 'joint_states')
            query_data = query_cache.joint;
        end
        
        fprintf('  Query data: %d points\n', size(query_data, 1));
       
        % Load query segments
        if strcmp(dtw_mode, 'position')
            query_segments = query_cache.segments.position;
        elseif strcmp(dtw_mode, 'joint_states')
            query_segments = query_cache.segments.joint;
        end
                
        % ================================================================
        % STEP 2: Get Candidate Data from Cache
        % ================================================================
        
        candidate_metadata = data_cache.candidates.metadata;
        num_candidates = height(candidate_metadata);
        
        if strcmp(dtw_mode, 'position')
            candidate_trajectories = data_cache.candidates.position;
        elseif strcmp(dtw_mode, 'joint_states')
            candidate_trajectories = data_cache.candidates.joint;
        end
        
        fprintf('  Candidates: %d trajectories\n', num_candidates);

        % ================================================================
        % 🔍 DEBUG STEP 3: Direct DTW to GT (BEFORE any filtering!)
        % ================================================================
        query_field = ['q_' strrep(query_id, '-', '_')];
        
        if has_ground_truth && isfield(ground_truth_map, query_field)
            gt_ids = ground_truth_map.(query_field).trajectories;
            
            if ~isempty(gt_ids) && length(gt_ids) > 0
                fprintf('\n  🔍 DEBUG: Direct DTW to GT (before filtering)...\n');
                
                gt_id = gt_ids{1};  % First GT
                
                % Find GT in ORIGINAL candidates (before any filtering!)
                gt_candidate_idx = find(strcmp(candidate_metadata.bahn_id, gt_id));
                
                if ~isempty(gt_candidate_idx)
                    gt_data = candidate_trajectories{gt_candidate_idx};
                    
                    % Compute DTW directly
                    tic;
                    direct_dtw = cDTW(query_data, gt_data, dtw_mode, cdtw_window, ...
                        inf, use_rotation_alignment, normalize_dtw);
                    t = toc;
                    
                    fprintf('    GT 1 (%s):\n', gt_id);
                    fprintf('      Direct DTW distance: %.2f (%.3f sec)\n', direct_dtw, t);
                    fprintf('      Query length: %d, GT length: %d\n', ...
                        size(query_data, 1), size(gt_data, 1));
                    
                    % Compare to random
                    random_idx = randi(num_candidates);
                    random_id = candidate_metadata.bahn_id{random_idx};
                    random_data = candidate_trajectories{random_idx};
                    
                    random_dtw = cDTW(query_data, random_data, dtw_mode, cdtw_window, ...
                        inf, use_rotation_alignment, normalize_dtw);
                    
                    fprintf('      Random (%s): DTW = %.2f\n', random_id, random_dtw);
                    
                    if direct_dtw > 0
                        ratio = random_dtw / direct_dtw;
                        fprintf('      Ratio (Random/GT): %.2fx ', ratio);
                        
                        if ratio >= 5.0
                            fprintf('✅\n');
                        elseif ratio >= 2.0
                            fprintf('⚠️\n');
                        else
                            fprintf('❌\n');
                        end
                    end
                    
                    fprintf('\n');
                else
                    fprintf('    ✗ GT not in candidate pool!\n\n');
                end
            end
        end
        
        % ================================================================
        % STEP 3: TRAJECTORY-LEVEL DTW
        % ================================================================
        
        fprintf('  Computing trajectory-level DTW...\n');
        
        traj_dtw_tic = tic;
        
        % LB_Kim filtering
        lb_kim_distances = zeros(num_candidates, 1);
        for i = 1:num_candidates
            candidate_data = candidate_trajectories{i};
            lb_kim_distances(i) = LB_Kim(query_data, candidate_data, dtw_mode, use_rotation_alignment, normalize_dtw);
        end
        
        [~, kim_order] = sort(lb_kim_distances);
        kim_keep_count = round(num_candidates * lb_kim_keep_ratio);
        
        candidates_after_kim = candidate_metadata(kim_order(1:kim_keep_count), :);
        candidate_trajectories_after_kim = candidate_trajectories(kim_order(1:kim_keep_count));
        
        fprintf('    LB_Kim: %d → %d (%.1f%% pruned)\n', ...
            num_candidates, kim_keep_count, 100*(1 - lb_kim_keep_ratio));
        
        % LB_Keogh filtering
        keogh_keep_count = min(lb_keogh_candidates, kim_keep_count);
        
        if kim_keep_count <= keogh_keep_count
            candidates_after_keogh = candidates_after_kim;
            candidate_trajectories_after_keogh = candidate_trajectories_after_kim;
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
        end
        
        fprintf('    LB_Keogh: %d → %d\n', kim_keep_count, keogh_keep_count);
        
        % DTW computation on ALL candidates after Keogh
        dtw_limit = keogh_keep_count;
        
        dtw_traj_results = struct();
        dtw_traj_results.bahn_id = candidates_after_keogh.bahn_id;
        dtw_traj_results.dtw_distance = inf(dtw_limit, 1);
        dtw_traj_results.length = candidates_after_keogh.length;
        dtw_traj_results.duration = candidates_after_keogh.duration;
        
        for i = 1:dtw_limit
            candidate_data = candidate_trajectories_after_keogh{i};
            
            if i <= top_k_trajectories
                best_so_far = inf;
            else
                sorted_dists = sort(dtw_traj_results.dtw_distance(1:i-1));
                if length(sorted_dists) >= top_k_trajectories
                    best_so_far = sorted_dists(top_k_trajectories);
                else
                    best_so_far = inf;
                end
            end
            
            dist = cDTW(query_data, candidate_data, dtw_mode, cdtw_window, ...
                best_so_far, use_rotation_alignment, normalize_dtw);
            
            dtw_traj_results.dtw_distance(i) = dist;
            
            if mod(i, 100) == 0
                fprintf('      Progress: %d/%d (%.1f%%)\n', i, dtw_limit, 100*i/dtw_limit);
            end
        end
        
        trajectory_table = struct2table(dtw_traj_results);
        trajectory_table = sortrows(trajectory_table, 'dtw_distance', 'ascend');
        
        traj_dtw_time = toc(traj_dtw_tic);
        
        fprintf('    ✓ Trajectory DTW: %.2fs (%d candidates)\n', traj_dtw_time, dtw_limit);

        % ================================================================
        % 🔍 DEBUG: Check GT in final DTW results
        % ================================================================
        
        if has_ground_truth && isfield(ground_truth_map, query_field)
            gt_ids = ground_truth_map.(query_field).trajectories;
            
            fprintf('\n  🔍 DEBUG: GT in final DTW results...\n');
            
            for i = 1:min(3, length(gt_ids))
                gt_id = gt_ids{i};
                gt_idx = find(strcmp(trajectory_table.bahn_id, gt_id));
                
                if ~isempty(gt_idx)
                    fprintf('    ✓ GT %d (%s): Rank %d, DTW %.2f\n', ...
                        i, gt_id, gt_idx, trajectory_table.dtw_distance(gt_idx));
                else
                    fprintf('    ✗ GT %d (%s): NOT FOUND\n', i, gt_id);
                end
            end
            
            fprintf('\n');
        end
        
        % ================================================================
        % STEP 4: SEGMENT-LEVEL DTW
        % ================================================================
        
        fprintf('  Computing segment-level DTW (%d segments)...\n', num_query_segments);
        
        segment_results = cell(num_query_segments, 1);
        segment_dtw_times = zeros(num_query_segments, 1);
        
        all_segments_metadata = cell(num_query_segments, 1);
        
        for seg_idx = 1:num_query_segments
            query_segment_id = sprintf('%s_%d', query_id, seg_idx);
            segment_mask = ~strcmp(data_cache.segments.segment_ids, query_segment_id);
            all_segments_metadata{seg_idx} = data_cache.segments.metadata(segment_mask, :);
        end
        
        for seg_idx = 1:num_query_segments
            query_segment_id = sprintf('%s_%d', query_id, seg_idx);
            query_segment_data = query_segments{seg_idx};
            
            candidate_segments = all_segments_metadata{seg_idx};
            num_candidate_segments = height(candidate_segments);
            
            if num_candidate_segments == 0
                segment_results{seg_idx} = [];
                continue;
            end
            
            seg_candidate_ids = candidate_segments.segment_id;
            [~, seg_cache_idx] = ismember(seg_candidate_ids, data_cache.segments.segment_ids);
            
            if strcmp(dtw_mode, 'position')
                seg_candidate_segments = data_cache.segments.position(seg_cache_idx);
            elseif strcmp(dtw_mode, 'joint_states')
                seg_candidate_segments = data_cache.segments.joint(seg_cache_idx);
            end
            
            seg_tic = tic;
            
            % LB_Kim
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
            
            % LB_Keogh
            seg_keogh_keep_count = min(lb_keogh_candidates, seg_kim_keep_count);
            
            if seg_kim_keep_count <= seg_keogh_keep_count
                candidates_after_keogh_seg = candidates_after_kim_seg;
                seg_candidate_segments_after_keogh = seg_candidate_segments_after_kim;
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
            end
            
            % DTW
            seg_dtw_limit = seg_keogh_keep_count;
            segment_dtw_distances = inf(seg_dtw_limit, 1);
            
            for cand_seg_idx = 1:seg_dtw_limit
                candidate_segment_data = seg_candidate_segments_after_keogh{cand_seg_idx};
                
                if ~isempty(candidate_segment_data)
                    candidate_segment_data = candidate_segment_data - candidate_segment_data(1, :);
                end
                
                if cand_seg_idx <= top_k_trajectories
                    best_so_far = inf;
                else
                    sorted_seg_dists = sort(segment_dtw_distances(1:cand_seg_idx-1));
                    if length(sorted_seg_dists) >= top_k_trajectories
                        best_so_far = sorted_seg_dists(top_k_trajectories);
                    else
                        best_so_far = inf;
                    end
                end
                
                dist = cDTW(query_segment_data, candidate_segment_data, dtw_mode, cdtw_window, ...
                    best_so_far, use_rotation_alignment, normalize_dtw);
                segment_dtw_distances(cand_seg_idx) = dist;
            end
            
            segment_results{seg_idx} = table(...
                candidates_after_keogh_seg.segment_id(1:seg_dtw_limit), ...
                candidates_after_keogh_seg.bahn_id(1:seg_dtw_limit), ...
                segment_dtw_distances(1:seg_dtw_limit), ...
                candidates_after_keogh_seg.length(1:seg_dtw_limit), ...
                candidates_after_keogh_seg.duration(1:seg_dtw_limit), ...
                'VariableNames', {'segment_id', 'bahn_id', 'dtw_distance', 'length', 'duration'});
            
            segment_results{seg_idx} = sortrows(segment_results{seg_idx}, 'dtw_distance', 'ascend');
            segment_dtw_times(seg_idx) = toc(seg_tic);
        end
        
        total_seg_dtw_time = sum(segment_dtw_times);
        
        fprintf('    ✓ Segment DTW: %.2fs (%.2fs avg per segment)\n', ...
            total_seg_dtw_time, mean(segment_dtw_times));
        
        % ================================================================
        % STEP 5: Store Results in Cache
        % ================================================================
        
        mode_time = toc(mode_tic);
        
        dtw_cache.(query_field_name).(dtw_mode) = struct();
        dtw_cache.(query_field_name).(dtw_mode).trajectory_ranking = trajectory_table;
        dtw_cache.(query_field_name).(dtw_mode).segment_rankings = segment_results;
        dtw_cache.(query_field_name).(dtw_mode).segment_dtw_times = segment_dtw_times;
        dtw_cache.(query_field_name).(dtw_mode).dtw_time = traj_dtw_time;
        dtw_cache.(query_field_name).(dtw_mode).num_candidates = num_candidates;
        
        fprintf('  ✓ Mode completed in %.2fs\n\n', mode_time);
    end
    
    fprintf('✓ Query %s completed\n\n', query_id);
end

%% ========================================================================
%  SUMMARY
%  ========================================================================

total_time = toc(overall_tic);

fprintf('========================================\n');
fprintf('DTW PRE-COMPUTATION COMPLETED\n');
fprintf('========================================\n\n');

fprintf('Total Time: %.1f minutes (%.2f seconds)\n\n', total_time/60, total_time);

fprintf('--- Summary ---\n');
fprintf('  Queries processed: %d\n', num_queries);
fprintf('  Modes per query: %d\n', num_modes);
fprintf('  Total combinations: %d\n', counter);
fprintf('  Average time per combination: %.2f seconds\n\n', total_time/counter);

cache_info = whos('dtw_cache');
fprintf('--- Cache Info ---\n');
fprintf('  Cache size: %.1f MB\n', cache_info.bytes / 1e6);

fprintf('\n✓ DTW cache ready for experiments!\n');
fprintf('  All experiments can now skip DTW computation\n');
fprintf('  and only compute embeddings (4-5× faster!)\n');
fprintf('========================================\n\n');

end
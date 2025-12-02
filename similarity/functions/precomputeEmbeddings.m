function embeddings_cache = precomputeEmbeddings(data_cache, query_ids, embedding_configs, config)
% PRECOMPUTEEMBEDDINGS - Pre-compute all embeddings for reuse across experiments
%
%   This function computes embeddings ONCE for all modalities, embedding
%   configurations, and queries. Experiments then only need to apply
%   different weights via RRF fusion.
%
%   INPUTS:
%       data_cache        - Pre-loaded data from loadDataExperiment.m
%       query_ids         - Cell array of query trajectory IDs
%       embedding_configs - Cell array of embedding configurations
%                           {name, n_coarse, n_fine, use_multi_scale}
%       config            - Configuration struct with embedding parameters
%
%   OUTPUT:
%       embeddings_cache - Struct containing all pre-computed embeddings:
%           .(query_field_name).(embedding_config_name)
%               .query_embeddings    - Query embeddings for all modalities
%                   .position, .joint, .orientation, .velocity, .metadata
%               .candidate_embeddings - Candidate embeddings
%                   .position, .joint, .orientation, .velocity, .metadata
%               .segment_embeddings  - Segment embeddings (per query segment)
%                   {seg_idx}.query.position, .joint, etc.
%                   {seg_idx}.candidates.position, .joint, etc.
%
%   PERFORMANCE:
%       For 4 queries × 4 embeddings × 5000 candidates:
%       - Time: ~10-15 minutes (one-time cost)
%       - Saves: ~40+ minutes over 128 experiments (75% reduction!)
%
%   Author: Gustavo Barros
%   Date: 02.12.2025

%% ========================================================================
%  VALIDATION & SETUP
%  ========================================================================

fprintf('\n╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  PRE-COMPUTING EMBEDDINGS                                      ║\n');
fprintf('║  This computes embeddings ONCE for reuse in all experiments   ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

% Validate inputs
if ~exist('data_cache', 'var') || isempty(data_cache)
    error('data_cache is required');
end

if ~exist('query_ids', 'var') || isempty(query_ids)
    error('query_ids is required');
end

if ~exist('embedding_configs', 'var') || isempty(embedding_configs)
    error('embedding_configs is required');
end

% Extract config parameters
if ~exist('config', 'var') || isempty(config)
    config = struct();
end

if isfield(config, 'norm_strategy')
    norm_strategy = config.norm_strategy;
else
    norm_strategy = 'max_extent';
end

% Get dimensions
num_queries = length(query_ids);
num_embeddings = size(embedding_configs, 1);
num_candidates = length(data_cache.candidates.bahn_ids);

fprintf('=== Configuration ===\n');
fprintf('  Queries: %d\n', num_queries);
for i = 1:num_queries
    fprintf('    %d. %s\n', i, query_ids{i});
end
fprintf('\n  Embedding Configs: %d\n', num_embeddings);
for i = 1:num_embeddings
    fprintf('    %d. %-25s (coarse=%3d, fine=%3d, multi=%d)\n', i, ...
        embedding_configs{i,1}, embedding_configs{i,2}, ...
        embedding_configs{i,3}, embedding_configs{i,4});
end
fprintf('\n  Candidates: %d trajectories\n', num_candidates);
fprintf('  Normalization: %s\n', norm_strategy);

total_combinations = num_queries * num_embeddings;
fprintf('\n  Total combinations: %d (query × embedding)\n', total_combinations);
fprintf('  Estimated time: %.1f minutes\n\n', total_combinations * 2.5);

embeddings_cache = struct();
overall_tic = tic;

%% ========================================================================
%  MAIN LOOP: FOR EACH QUERY × EACH EMBEDDING CONFIG
%  ========================================================================

counter = 0;

for q_idx = 1:num_queries
    query_id = query_ids{q_idx};
    query_field_name = ['q_' strrep(query_id, '-', '_')];
    
    fprintf('╔════════════════════════════════════════════════════════════════╗\n');
    fprintf('║  QUERY %d/%d: %-50s║\n', q_idx, num_queries, query_id);
    fprintf('╚════════════════════════════════════════════════════════════════╝\n');
    
    % Check if query exists in cache
    if ~isfield(data_cache.queries, query_field_name)
        warning('Query %s not found in data_cache - skipping', query_id);
        continue;
    end
    
    query_cache = data_cache.queries.(query_field_name);
    query_metadata = query_cache.metadata;
    num_query_segments = query_metadata.num_segments;
    
    fprintf('  Segments: %d\n\n', num_query_segments);
    
    embeddings_cache.(query_field_name) = struct();
    
    % ====================================================================
    % FOR EACH EMBEDDING CONFIGURATION
    % ====================================================================
    
    for emb_idx = 1:num_embeddings
        emb_name = embedding_configs{emb_idx, 1};
        n_coarse = embedding_configs{emb_idx, 2};
        n_fine = embedding_configs{emb_idx, 3};
        use_multi_scale = embedding_configs{emb_idx, 4};
        
        counter = counter + 1;
        
        fprintf('  ┌─────────────────────────────────────────────────────────┐\n');
        fprintf('  │ [%3d/%3d] %-48s │\n', counter, total_combinations, emb_name);
        fprintf('  └─────────────────────────────────────────────────────────┘\n');
        
        emb_tic = tic;
        
        % Create valid field name
        emb_field_name = strrep(emb_name, '-', '_');
        emb_field_name = strrep(emb_field_name, ' ', '_');
        
        embeddings_cache.(query_field_name).(emb_field_name) = struct();
        
        % ============================================================
        % STEP 1: QUERY EMBEDDINGS (5 modalities)
        % ============================================================
        
        fprintf('    Step 1/3: Query embeddings (5 modalities)...\n');
        
        query_embeddings = struct();
        
        % Position
        query_position = query_cache.position;
        query_position_norm = query_position - query_position(1, :);
        query_embeddings.position = computePositionEmbedding(...
            query_position_norm, n_coarse, n_fine);
        
        % Joint States
        query_joint = query_cache.joint;
        query_embeddings.joint = computeJointEmbedding(...
            query_joint, n_coarse, n_fine);
        
        % Orientation
        query_orientation = query_cache.orientation;
        query_embeddings.orientation = computeOrientationEmbedding(...
            query_orientation, n_coarse, n_fine);
        
        % Velocity (from position)
        query_velocity = diff(query_position_norm);
        if ~isempty(query_velocity)
            query_embeddings.velocity = computeVelocityEmbedding(...
                query_position_norm, n_coarse, n_fine);
        else
            query_embeddings.velocity = [];
        end
        
        % Metadata
        query_embeddings.metadata = computeMetadataEmbedding(query_metadata);
        
        fprintf('      ✓ Position: [%d], Joint: [%d], Orient: [%d], Vel: [%d], Meta: [%d]\n', ...
            length(query_embeddings.position), length(query_embeddings.joint), ...
            length(query_embeddings.orientation), length(query_embeddings.velocity), ...
            length(query_embeddings.metadata));
        
        % ============================================================
        % STEP 2: CANDIDATE EMBEDDINGS (5 modalities × num_candidates)
        % ============================================================
        
        fprintf('    Step 2/3: Candidate embeddings (%d trajectories)...\n', num_candidates);
        
        % Pre-allocate embedding matrices (much faster than cells!)
        dim_pos = length(query_embeddings.position);
        dim_joint = length(query_embeddings.joint);
        dim_orient = length(query_embeddings.orientation);
        dim_vel = length(query_embeddings.velocity);
        dim_meta = length(query_embeddings.metadata);
        
        candidate_embeddings = struct();
        candidate_embeddings.position = zeros(num_candidates, dim_pos);
        candidate_embeddings.joint = zeros(num_candidates, dim_joint);
        candidate_embeddings.orientation = zeros(num_candidates, dim_orient);
        candidate_embeddings.velocity = zeros(num_candidates, dim_vel);
        candidate_embeddings.metadata = zeros(num_candidates, dim_meta);
        
        for cand_idx = 1:num_candidates
            % Position
            cand_position = data_cache.candidates.position{cand_idx};
            cand_position_norm = cand_position - cand_position(1, :);
            candidate_embeddings.position(cand_idx, :) = computePositionEmbedding(...
                cand_position_norm, n_coarse, n_fine);
            
            % Joint
            cand_joint = data_cache.candidates.joint{cand_idx};
            candidate_embeddings.joint(cand_idx, :) = computeJointEmbedding(...
                cand_joint, n_coarse, n_fine);
            
            % Orientation
            cand_orientation = data_cache.candidates.orientation{cand_idx};
            candidate_embeddings.orientation(cand_idx, :) = computeOrientationEmbedding(...
                cand_orientation, n_coarse, n_fine);
            
            % Velocity
            candidate_embeddings.velocity(cand_idx, :) = computeVelocityEmbedding(...
                cand_position_norm, n_coarse, n_fine);
            
            % Metadata
            cand_metadata_struct = struct();
            cand_metadata_struct.movement_type = data_cache.candidates.metadata.movement_type{cand_idx};
            cand_metadata_struct.length = data_cache.candidates.metadata.length(cand_idx);
            cand_metadata_struct.duration = data_cache.candidates.metadata.duration(cand_idx);
            cand_metadata_struct.min_twist = data_cache.candidates.metadata.min_twist_ist(cand_idx);
            cand_metadata_struct.max_twist = data_cache.candidates.metadata.max_twist_ist(cand_idx);
            cand_metadata_struct.mean_twist = data_cache.candidates.metadata.mean_twist_ist(cand_idx);
            cand_metadata_struct.median_twist = data_cache.candidates.metadata.median_twist_ist(cand_idx);
            cand_metadata_struct.std_twist = data_cache.candidates.metadata.std_twist_ist(cand_idx);
            cand_metadata_struct.min_acceleration = data_cache.candidates.metadata.min_acceleration_ist(cand_idx);
            cand_metadata_struct.max_acceleration = data_cache.candidates.metadata.max_acceleration_ist(cand_idx);
            cand_metadata_struct.mean_acceleration = data_cache.candidates.metadata.mean_acceleration_ist(cand_idx);
            cand_metadata_struct.median_acceleration = data_cache.candidates.metadata.median_acceleration_ist(cand_idx);
            cand_metadata_struct.std_acceleration = data_cache.candidates.metadata.std_acceleration_ist(cand_idx);
            
            candidate_embeddings.metadata(cand_idx, :) = computeMetadataEmbedding(cand_metadata_struct);
            
            if mod(cand_idx, 100) == 0 || cand_idx == num_candidates
                fprintf('      Progress: %d/%d candidates\n', cand_idx, num_candidates);
            end
        end
        
        fprintf('      ✓ All candidate embeddings computed\n');
        
        % ============================================================
        % STEP 3: SEGMENT EMBEDDINGS (per query segment)
        % ============================================================
        
        fprintf('    Step 3/3: Segment embeddings (%d query segments)...\n', num_query_segments);
        
        segment_embeddings = cell(num_query_segments, 1);
        
        for seg_idx = 1:num_query_segments
            query_segment_id = sprintf('%s_%d', query_id, seg_idx);
            
            % Get segment data from query cache
            query_seg_position = query_cache.segments.position{seg_idx};
            query_seg_joint = query_cache.segments.joint{seg_idx};
            query_seg_orientation = query_cache.segments.orientation{seg_idx};
            query_seg_metadata = query_cache.segments.metadata{seg_idx};
            
            % Normalize position
            if ~isempty(query_seg_position)
                query_seg_position_norm = query_seg_position - query_seg_position(1, :);
            else
                query_seg_position_norm = [];
            end
            
            % Query segment embeddings
            seg_emb = struct();
            seg_emb.query = struct();
            
            seg_emb.query.position = computePositionEmbedding(...
                query_seg_position_norm, n_coarse, n_fine);
            seg_emb.query.joint = computeJointEmbedding(...
                query_seg_joint, n_coarse, n_fine);
            seg_emb.query.orientation = computeOrientationEmbedding(...
                query_seg_orientation, n_coarse, n_fine);
            
            if ~isempty(query_seg_position_norm)
                seg_emb.query.velocity = computeVelocityEmbedding(...
                    query_seg_position_norm, n_coarse, n_fine);
            else
                seg_emb.query.velocity = [];
            end
            
            seg_emb.query.metadata = computeMetadataEmbedding(query_seg_metadata);
            
            % Get candidate segments for this query segment
            % Filter segments from cache for this query segment
            segment_mask = ~strcmp(data_cache.segments.segment_ids, query_segment_id);
            filtered_segment_ids = data_cache.segments.segment_ids(segment_mask);
            num_cand_segments = length(filtered_segment_ids);
            
            if num_cand_segments == 0
                segment_embeddings{seg_idx} = seg_emb;
                continue;
            end
            
            % Pre-allocate segment embedding matrices
            seg_emb.candidates = struct();
            seg_emb.candidates.position = zeros(num_cand_segments, dim_pos);
            seg_emb.candidates.joint = zeros(num_cand_segments, dim_joint);
            seg_emb.candidates.orientation = zeros(num_cand_segments, dim_orient);
            seg_emb.candidates.velocity = zeros(num_cand_segments, dim_vel);
            seg_emb.candidates.metadata = zeros(num_cand_segments, dim_meta);
            seg_emb.candidates.segment_ids = filtered_segment_ids;
            
            % Get indices of candidate segments
            [~, seg_cache_idx] = ismember(filtered_segment_ids, data_cache.segments.segment_ids);
            
            for cand_seg_idx = 1:num_cand_segments
                cache_idx = seg_cache_idx(cand_seg_idx);
                
                % Load segment data
                cand_seg_position = data_cache.segments.position{cache_idx};
                cand_seg_joint = data_cache.segments.joint{cache_idx};
                cand_seg_orientation = data_cache.segments.orientation{cache_idx};
                
                % Normalize position
                if ~isempty(cand_seg_position)
                    cand_seg_position_norm = cand_seg_position - cand_seg_position(1, :);
                else
                    cand_seg_position_norm = [];
                end
                
                % Embeddings
                seg_emb.candidates.position(cand_seg_idx, :) = computePositionEmbedding(...
                    cand_seg_position_norm, n_coarse, n_fine);
                seg_emb.candidates.joint(cand_seg_idx, :) = computeJointEmbedding(...
                    cand_seg_joint, n_coarse, n_fine);
                seg_emb.candidates.orientation(cand_seg_idx, :) = computeOrientationEmbedding(...
                    cand_seg_orientation, n_coarse, n_fine);
                
                if ~isempty(cand_seg_position_norm)
                    seg_emb.candidates.velocity(cand_seg_idx, :) = computeVelocityEmbedding(...
                        cand_seg_position_norm, n_coarse, n_fine);
                end
                
                % Metadata
                cand_seg_metadata_struct = struct();
                cand_seg_metadata = data_cache.segments.metadata(cache_idx, :);
                cand_seg_metadata_struct.movement_type = cand_seg_metadata.movement_type{1};
                cand_seg_metadata_struct.length = cand_seg_metadata.length;
                cand_seg_metadata_struct.duration = cand_seg_metadata.duration;
                cand_seg_metadata_struct.min_twist = cand_seg_metadata.min_twist_ist;
                cand_seg_metadata_struct.max_twist = cand_seg_metadata.max_twist_ist;
                cand_seg_metadata_struct.mean_twist = cand_seg_metadata.mean_twist_ist;
                cand_seg_metadata_struct.median_twist = cand_seg_metadata.median_twist_ist;
                cand_seg_metadata_struct.std_twist = cand_seg_metadata.std_twist_ist;
                cand_seg_metadata_struct.min_acceleration = cand_seg_metadata.min_acceleration_ist;
                cand_seg_metadata_struct.max_acceleration = cand_seg_metadata.max_acceleration_ist;
                cand_seg_metadata_struct.mean_acceleration = cand_seg_metadata.mean_acceleration_ist;
                cand_seg_metadata_struct.median_acceleration = cand_seg_metadata.median_acceleration_ist;
                cand_seg_metadata_struct.std_acceleration = cand_seg_metadata.std_acceleration_ist;
                
                seg_emb.candidates.metadata(cand_seg_idx, :) = computeMetadataEmbedding(cand_seg_metadata_struct);
            end
            
            segment_embeddings{seg_idx} = seg_emb;
            
            if mod(seg_idx, 5) == 0 || seg_idx == num_query_segments
                fprintf('      Segment %d/%d: %d candidates\n', ...
                    seg_idx, num_query_segments, num_cand_segments);
            end
        end
        
        fprintf('      ✓ All segment embeddings computed\n');
        
        % ============================================================
        % STORE IN CACHE
        % ============================================================
        
        embeddings_cache.(query_field_name).(emb_field_name).query_embeddings = query_embeddings;
        embeddings_cache.(query_field_name).(emb_field_name).candidate_embeddings = candidate_embeddings;
        embeddings_cache.(query_field_name).(emb_field_name).segment_embeddings = segment_embeddings;
        
        emb_time = toc(emb_tic);
        
        fprintf('    ✓ Completed in %.2f seconds (%.1f min)\n\n', emb_time, emb_time/60);
    end
    
    fprintf('  ✓ Query %s completed\n\n', query_id);
end

%% ========================================================================
%  SUMMARY
%  ========================================================================

total_time = toc(overall_tic);

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  EMBEDDINGS PRE-COMPUTATION COMPLETED                          ║\n');
fprintf('╠════════════════════════════════════════════════════════════════╣\n');
fprintf('║  Total Time: %7.1f minutes (%7.1f seconds)                  ║\n', total_time/60, total_time);
fprintf('╠════════════════════════════════════════════════════════════════╣\n');
fprintf('║  Queries processed:      %3d                                   ║\n', num_queries);
fprintf('║  Embedding configs:      %3d                                   ║\n', num_embeddings);
fprintf('║  Total combinations:     %3d                                   ║\n', counter);
fprintf('║  Avg time/combination:   %7.2f seconds                        ║\n', total_time/counter);
fprintf('╠════════════════════════════════════════════════════════════════╣\n');

% Memory usage
cache_info = whos('embeddings_cache');
fprintf('║  Cache size:             %7.1f MB                             ║\n', cache_info.bytes / 1e6);
fprintf('╠════════════════════════════════════════════════════════════════╣\n');
fprintf('║  ✓ Embeddings cache ready for experiments!                     ║\n');
fprintf('║  All experiments can now skip embedding computation            ║\n');
fprintf('║  and only apply different weights via RRF fusion!              ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

end
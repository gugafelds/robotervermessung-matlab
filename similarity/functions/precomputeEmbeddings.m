function embeddings_cache = precomputeEmbeddings(data_cache, dtw_cache, query_ids, embedding_configs, config)
% PRECOMPUTEEMBEDDINGS - Pre-compute all embeddings for reuse across experiments
%
%   This function computes embeddings ONCE for all modalities, embedding
%   configurations, queries, and DTW modes. Experiments then only need to
%   apply different weights via RRF fusion.
%
%   INPUTS:
%       data_cache        - Pre-loaded data from loadDataExperiment.m
%       dtw_cache         - Pre-computed DTW from precomputeDTW.m
%       query_ids         - Cell array of query trajectory IDs
%       embedding_configs - Cell array of embedding configurations
%                           {name, n_coarse, n_fine, use_multi_scale}
%       config            - Configuration struct with embedding parameters
%
%   OUTPUT:
%       embeddings_cache - Struct containing all pre-computed embeddings:
%           .(query_field_name).(dtw_mode).(embedding_config_name)
%               .query_embeddings    - Query embeddings for all modalities
%                   .position, .joint, .orientation, .velocity, .metadata
%               .candidate_embeddings - Candidate embeddings
%                   .position, .joint, .orientation, .velocity, .metadata
%               .segment_embeddings  - Segment embeddings (per query segment)
%                   {seg_idx}.position, .joint, .orientation, .velocity, .metadata
%
%   PERFORMANCE:
%       For 4 queries × 2 modes × 4 embeddings × 5000 candidates:
%       - Time: ~15-20 minutes (one-time cost)
%       - Saves: ~50 minutes over 128 experiments (80% reduction in embedding time!)
%
%   Author: Gustavo Barros
%   Date: 29.11.2025

%% ========================================================================
%  VALIDATION & SETUP
%  ========================================================================

fprintf('\n========================================\n');
fprintf('PRE-COMPUTING EMBEDDINGS\n');
fprintf('========================================\n\n');

% Extract config parameters
if isfield(config, 'norm_strategy')
    norm_strategy = config.norm_strategy;
else
    norm_strategy = 'max_extent';
end

% DTW modes to compute embeddings for
dtw_modes = {'position', 'joint_states'};

num_queries = length(query_ids);
num_modes = length(dtw_modes);
num_embeddings = size(embedding_configs, 1);
total_combinations = num_queries * num_modes * num_embeddings;

fprintf('Configuration:\n');
fprintf('  Queries: %d (%s)\n', num_queries, strjoin(query_ids, ', '));
fprintf('  DTW Modes: %s\n', strjoin(dtw_modes, ', '));
fprintf('  Embedding Configs: %d\n', num_embeddings);
for i = 1:num_embeddings
    fprintf('    %d. %s (n_coarse=%d, n_fine=%d)\n', i, ...
        embedding_configs{i,1}, embedding_configs{i,2}, embedding_configs{i,3});
end
fprintf('  Total combinations: %d\n', total_combinations);
fprintf('  Normalization: %s\n\n', norm_strategy);

embeddings_cache = struct();
overall_tic = tic;

%% ========================================================================
%  MAIN LOOP: FOR EACH QUERY × EACH MODE × EACH EMBEDDING CONFIG
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
    
    embeddings_cache.(query_field_name) = struct();
    
    % ====================================================================
    % FOR EACH DTW MODE
    % ====================================================================
    
    for mode_idx = 1:num_modes
        dtw_mode = dtw_modes{mode_idx};
        
        fprintf('--- DTW Mode: %s ---\n', dtw_mode);
        
        embeddings_cache.(query_field_name).(dtw_mode) = struct();
        
        % Get DTW rankings from cache
        trajectory_ranking = dtw_cache.(query_field_name).(dtw_mode).trajectory_ranking;
        segment_rankings = dtw_cache.(query_field_name).(dtw_mode).segment_rankings;
        
        num_candidates = height(trajectory_ranking);
        
        % ================================================================
        % FOR EACH EMBEDDING CONFIGURATION
        % ================================================================
        
        for emb_idx = 1:num_embeddings
            emb_name = embedding_configs{emb_idx, 1};
            n_coarse = embedding_configs{emb_idx, 2};
            n_fine = embedding_configs{emb_idx, 3};
            use_multi_scale = embedding_configs{emb_idx, 4};
            
            counter = counter + 1;
            
            fprintf('  [%d/%d] Embedding: %s\n', counter, total_combinations, emb_name);
            
            emb_tic = tic;
            
            % Create valid field name
            emb_field_name = strrep(emb_name, '-', '_');
            
            embeddings_cache.(query_field_name).(dtw_mode).(emb_field_name) = struct();
            
            % ============================================================
            % STEP 1: QUERY EMBEDDINGS (5 modalities)
            % ============================================================
            
            fprintf('    Computing query embeddings...\n');
            
            query_embeddings = struct();
            
            % Position
            query_position = query_cache.position;
            query_position_norm = query_position - query_position(1, :);
            query_embeddings.position = createTrajectoryEmbedding(...
                query_position_norm, n_coarse, 0, n_fine, norm_strategy, use_multi_scale);
            
            % Joint States
            query_joint = query_cache.joint;
            query_embeddings.joint = createTrajectoryEmbedding(...
                query_joint, n_coarse, 0, n_fine, norm_strategy, use_multi_scale);
            
            % Orientation
            query_orientation = query_cache.orientation;
            query_embeddings.orientation = createTrajectoryEmbedding(...
                query_orientation, n_coarse, 0, n_fine, norm_strategy, use_multi_scale);
            
            % Velocity (from position)
            query_velocity = diff(query_position_norm);
            if ~isempty(query_velocity)
                query_embeddings.velocity = createTrajectoryEmbedding(...
                    query_velocity, n_coarse, 0, n_fine, norm_strategy, use_multi_scale);
            else
                query_embeddings.velocity = [];
            end
            
            % Metadata
            query_embeddings.metadata = createMetadataEmbedding(query_metadata);
            
            fprintf('      ✓ Query embeddings: 5 modalities\n');
            
            % ============================================================
            % STEP 2: CANDIDATE EMBEDDINGS (5 modalities × num_candidates)
            % ============================================================
            
            fprintf('    Computing candidate embeddings (%d trajectories)...\n', num_candidates);
            
            candidate_embeddings = struct();
            candidate_embeddings.position = cell(num_candidates, 1);
            candidate_embeddings.joint = cell(num_candidates, 1);
            candidate_embeddings.orientation = cell(num_candidates, 1);
            candidate_embeddings.velocity = cell(num_candidates, 1);
            candidate_embeddings.metadata = cell(num_candidates, 1);
            
            for cand_idx = 1:num_candidates
                candidate_id = trajectory_ranking.bahn_id{cand_idx};
                
                % Find candidate in data_cache
                cache_idx = find(strcmp(data_cache.candidates.bahn_ids, candidate_id), 1);
                
                if isempty(cache_idx)
                    warning('Candidate %s not found in cache', candidate_id);
                    continue;
                end
                
                % Position
                cand_position = data_cache.candidates.position{cache_idx};
                cand_position_norm = cand_position - cand_position(1, :);
                candidate_embeddings.position{cand_idx} = createTrajectoryEmbedding(...
                    cand_position_norm, n_coarse, 0, n_fine, norm_strategy, use_multi_scale);
                
                % Joint
                cand_joint = data_cache.candidates.joint{cache_idx};
                candidate_embeddings.joint{cand_idx} = createTrajectoryEmbedding(...
                    cand_joint, n_coarse, 0, n_fine, norm_strategy, use_multi_scale);
                
                % Orientation
                cand_orientation = data_cache.candidates.orientation{cache_idx};
                candidate_embeddings.orientation{cand_idx} = createTrajectoryEmbedding(...
                    cand_orientation, n_coarse, 0, n_fine, norm_strategy, use_multi_scale);
                
                % Velocity
                cand_velocity = diff(cand_position_norm);
                if ~isempty(cand_velocity)
                    candidate_embeddings.velocity{cand_idx} = createTrajectoryEmbedding(...
                        cand_velocity, n_coarse, 0, n_fine, norm_strategy, use_multi_scale);
                else
                    candidate_embeddings.velocity{cand_idx} = [];
                end
                
                % Metadata
                cand_metadata = data_cache.candidates.metadata(cache_idx, :);
                candidate_embeddings.metadata{cand_idx} = createMetadataEmbedding(cand_metadata);
                
                if mod(cand_idx, 100) == 0
                    fprintf('      Progress: %d/%d candidates\n', cand_idx, num_candidates);
                end
            end
            
            fprintf('      ✓ Candidate embeddings: %d × 5 modalities\n', num_candidates);
            
            % ============================================================
            % STEP 3: SEGMENT EMBEDDINGS (per query segment)
            % ============================================================
            
            fprintf('    Computing segment embeddings (%d query segments)...\n', num_query_segments);
            
            segment_embeddings = cell(num_query_segments, 1);
            
            for seg_idx = 1:num_query_segments
                query_segment_id = sprintf('%s_%d', query_id, seg_idx);
                
                % Get segment data from query cache
                query_seg_position = query_cache.segments.position{seg_idx};
                query_seg_joint = query_cache.segments.joint{seg_idx};
                query_seg_orientation = query_cache.segments.orientation{seg_idx};
                query_seg_metadata = query_cache.segments.metadata{seg_idx};
                
                % Normalize position
                query_seg_position_norm = query_seg_position - query_seg_position(1, :);
                
                % Query segment embeddings
                seg_emb = struct();
                seg_emb.query = struct();
                
                seg_emb.query.position = createTrajectoryEmbedding(...
                    query_seg_position_norm, n_coarse, 0, n_fine, norm_strategy, use_multi_scale);
                seg_emb.query.joint = createTrajectoryEmbedding(...
                    query_seg_joint, n_coarse, 0, n_fine, norm_strategy, use_multi_scale);
                seg_emb.query.orientation = createTrajectoryEmbedding(...
                    query_seg_orientation, n_coarse, 0, n_fine, norm_strategy, use_multi_scale);
                
                query_seg_velocity = diff(query_seg_position_norm);
                if ~isempty(query_seg_velocity)
                    seg_emb.query.velocity = createTrajectoryEmbedding(...
                        query_seg_velocity, n_coarse, 0, n_fine, norm_strategy, use_multi_scale);
                else
                    seg_emb.query.velocity = [];
                end
                
                seg_emb.query.metadata = createMetadataEmbedding(query_seg_metadata);
                
                % Candidate segment embeddings
                segment_ranking = segment_rankings{seg_idx};
                
                if isempty(segment_ranking)
                    segment_embeddings{seg_idx} = seg_emb;
                    continue;
                end
                
                num_cand_segments = height(segment_ranking);
                
                seg_emb.candidates = struct();
                seg_emb.candidates.position = cell(num_cand_segments, 1);
                seg_emb.candidates.joint = cell(num_cand_segments, 1);
                seg_emb.candidates.orientation = cell(num_cand_segments, 1);
                seg_emb.candidates.velocity = cell(num_cand_segments, 1);
                seg_emb.candidates.metadata = cell(num_cand_segments, 1);
                
                for cand_seg_idx = 1:num_cand_segments
                    cand_segment_id = segment_ranking.segment_id{cand_seg_idx};
                    
                    % Find in segment cache
                    seg_cache_idx = find(strcmp(data_cache.segments.segment_ids, cand_segment_id), 1);
                    
                    if isempty(seg_cache_idx)
                        continue;
                    end
                    
                    % Load segment data
                    cand_seg_position = data_cache.segments.position{seg_cache_idx};
                    cand_seg_joint = data_cache.segments.joint{seg_cache_idx};
                    cand_seg_orientation = data_cache.segments.orientation{seg_cache_idx};
                    cand_seg_metadata = data_cache.segments.metadata(seg_cache_idx, :);
                    
                    % Normalize position
                    if ~isempty(cand_seg_position)
                        cand_seg_position_norm = cand_seg_position - cand_seg_position(1, :);
                    else
                        cand_seg_position_norm = [];
                    end
                    
                    % Embeddings
                    seg_emb.candidates.position{cand_seg_idx} = createTrajectoryEmbedding(...
                        cand_seg_position_norm, n_coarse, 0, n_fine, norm_strategy, use_multi_scale);
                    seg_emb.candidates.joint{cand_seg_idx} = createTrajectoryEmbedding(...
                        cand_seg_joint, n_coarse, 0, n_fine, norm_strategy, use_multi_scale);
                    seg_emb.candidates.orientation{cand_seg_idx} = createTrajectoryEmbedding(...
                        cand_seg_orientation, n_coarse, 0, n_fine, norm_strategy, use_multi_scale);
                    
                    cand_seg_velocity = diff(cand_seg_position_norm);
                    if ~isempty(cand_seg_velocity)
                        seg_emb.candidates.velocity{cand_seg_idx} = createTrajectoryEmbedding(...
                            cand_seg_velocity, n_coarse, 0, n_fine, norm_strategy, use_multi_scale);
                    else
                        seg_emb.candidates.velocity{cand_seg_idx} = [];
                    end
                    
                    seg_emb.candidates.metadata{cand_seg_idx} = createMetadataEmbedding(cand_seg_metadata);
                end
                
                segment_embeddings{seg_idx} = seg_emb;
            end
            
            fprintf('      ✓ Segment embeddings: %d segments × 5 modalities\n', num_query_segments);
            
            % ============================================================
            % STORE IN CACHE
            % ============================================================
            
            embeddings_cache.(query_field_name).(dtw_mode).(emb_field_name).query_embeddings = query_embeddings;
            embeddings_cache.(query_field_name).(dtw_mode).(emb_field_name).candidate_embeddings = candidate_embeddings;
            embeddings_cache.(query_field_name).(dtw_mode).(emb_field_name).segment_embeddings = segment_embeddings;
            
            emb_time = toc(emb_tic);
            fprintf('    ✓ Embedding config completed in %.2fs\n\n', emb_time);
        end
    end
    
    fprintf('✓ Query %s completed\n\n', query_id);
end

%% ========================================================================
%  SUMMARY
%  ========================================================================

total_time = toc(overall_tic);

fprintf('========================================\n');
fprintf('EMBEDDINGS PRE-COMPUTATION COMPLETED\n');
fprintf('========================================\n\n');

fprintf('Total Time: %.1f minutes (%.2f seconds)\n\n', total_time/60, total_time);

fprintf('--- Summary ---\n');
fprintf('  Queries processed: %d\n', num_queries);
fprintf('  DTW modes: %d\n', num_modes);
fprintf('  Embedding configs: %d\n', num_embeddings);
fprintf('  Total combinations: %d\n', counter);
fprintf('  Average time per combination: %.2f seconds\n\n', total_time/counter);

% Memory usage
cache_info = whos('embeddings_cache');
fprintf('--- Cache Info ---\n');
fprintf('  Cache size: %.1f MB\n', cache_info.bytes / 1e6);

fprintf('\n✓ Embeddings cache ready for experiments!\n');
fprintf('  All experiments can now skip embedding computation\n');
fprintf('  and only apply different weights via RRF fusion!\n');
fprintf('========================================\n\n');

end
function emb_cache = precomputeEmbeddings_v2(data_cache, query_ids, embedding_configs)
% Berechnet Embeddings für alle Queries × Configs

fprintf('Computing embeddings...\n');
tic;

emb_cache = struct();
num_queries = length(query_ids);
num_configs = size(embedding_configs, 1);

for q_idx = 1:num_queries
    query_id = query_ids{q_idx};
    qf = sprintf('q_%s', strrep(query_id, '-', '_'));
    
    if ~isfield(data_cache.queries, qf), continue; end
    query = data_cache.queries.(qf);
    
    emb_cache.(qf) = struct();
    
    for cfg_idx = 1:num_configs
        cfg_name = embedding_configs{cfg_idx, 1};
        n_coarse = embedding_configs{cfg_idx, 2};
        n_fine = embedding_configs{cfg_idx, 3};
        
        ef = matlab.lang.makeValidName(cfg_name);
        
        % === Query Embeddings ===
        query_data = struct();
        query_data.position = query.position;
        query_data.joint = query.joint;
        query_data.orientation = query.orientation;
        query_data.metadata = query.metadata;
        
        emb_cache.(qf).(ef).query = computeAllEmbeddings(query_data, n_coarse, n_fine);
        
        % === Candidate Embeddings ===
        num_cands = length(data_cache.candidates.bahn_ids);
        first_cand = getCandidateData(data_cache, 1);
        cand_emb = initEmbeddingMatrices(num_cands, first_cand, n_coarse, n_fine);
        
        for c = 1:num_cands
            cand = getCandidateData(data_cache, c);
            cand_emb = fillEmbeddingRow(cand_emb, c, cand, n_coarse, n_fine);
        end
        cand_emb.bahn_ids = data_cache.candidates.bahn_ids;
        emb_cache.(qf).(ef).candidates = cand_emb;
        
        % === Segment Embeddings ===
        num_segs = query.metadata.num_segments;
        seg_emb = cell(num_segs, 1);
        
        for s = 1:num_segs
            seg_emb{s} = struct();
            
            % Query segment
            seg_data = getSegmentData(query, s);
            seg_emb{s}.query = computeAllEmbeddings(seg_data, n_coarse, n_fine);
            
            % Candidate segments (exclude own)
            own_seg_id = sprintf('%s_%d', query_id, s);
            mask = ~strcmp(data_cache.segments.segment_ids, own_seg_id);
            seg_ids = data_cache.segments.segment_ids(mask);
            num_cand_segs = sum(mask);
            
            if num_cand_segs > 0
                idx = find(mask);
                first_cand_seg = getSegmentDataFromCache(data_cache, idx(1));
                cand_seg_emb = initEmbeddingMatrices(num_cand_segs, first_cand_seg, n_coarse, n_fine);
                
                for cs = 1:num_cand_segs
                    cand_seg = getSegmentDataFromCache(data_cache, idx(cs));
                    cand_seg_emb = fillEmbeddingRow(cand_seg_emb, cs, cand_seg, n_coarse, n_fine);
                end
                cand_seg_emb.segment_ids = seg_ids;
                seg_emb{s}.candidates = cand_seg_emb;
            else
                seg_emb{s}.candidates = struct('segment_ids', {});
            end
        end
        emb_cache.(qf).(ef).segments = seg_emb;
    end
    
    fprintf('  Query %d/%d done\n', q_idx, num_queries);
end

fprintf('Embeddings done: %.1f min\n', toc/60);
end

%% === HELPER FUNCTIONS ===

function emb = computeAllEmbeddings(data, n_coarse, n_fine)
    pos = data.position;
    if ~isempty(pos)
        pos = pos - pos(1,:);  % normalize to origin
    end
    
    emb.position = computePositionEmbedding(pos, n_coarse, n_fine);
    emb.joint = computeJointEmbedding(data.joint, n_coarse, n_fine);
    emb.orientation = computeOrientationEmbedding(data.orientation, n_coarse, n_fine);
    emb.velocity = computeVelocityEmbedding(pos, n_coarse, n_fine);
    emb.metadata = computeMetadataEmbedding(data.metadata);
end

function m = initEmbeddingMatrices(n, sample_data, n_coarse, n_fine)
    sample_emb = computeAllEmbeddings(sample_data, n_coarse, n_fine);
    
    m.position = zeros(n, length(sample_emb.position));
    m.joint = zeros(n, length(sample_emb.joint));
    m.orientation = zeros(n, length(sample_emb.orientation));
    m.velocity = zeros(n, length(sample_emb.velocity));
    m.metadata = zeros(n, length(sample_emb.metadata));
end

function m = fillEmbeddingRow(m, idx, data, n_coarse, n_fine)
    emb = computeAllEmbeddings(data, n_coarse, n_fine);
    m.position(idx,:) = emb.position;
    m.joint(idx,:) = emb.joint;
    m.orientation(idx,:) = emb.orientation;
    m.velocity(idx,:) = emb.velocity;
    m.metadata(idx,:) = emb.metadata;
end

function cand = getCandidateData(data_cache, idx)
    cand.position = data_cache.candidates.position{idx};
    cand.joint = data_cache.candidates.joint{idx};
    cand.orientation = data_cache.candidates.orientation{idx};
    cand.metadata = getMetadataStruct(data_cache.candidates.metadata, idx);
end

function seg = getSegmentData(query, seg_idx)
    seg.position = query.segments.position{seg_idx};
    seg.joint = query.segments.joint{seg_idx};
    seg.orientation = query.segments.orientation{seg_idx};
    seg.metadata = query.segments.metadata{seg_idx};
end

function seg = getSegmentDataFromCache(data_cache, idx)
    seg.position = data_cache.segments.position{idx};
    seg.joint = data_cache.segments.joint{idx};
    seg.orientation = data_cache.segments.orientation{idx};
    seg.metadata = getMetadataStructFromTable(data_cache.segments.metadata, idx);
end

function m = getMetadataStruct(tbl, idx)
    m.movement_type = tbl.movement_type{idx};
    m.length = tbl.length(idx);
    m.duration = tbl.duration(idx);
    m.min_twist = tbl.min_twist_ist(idx);
    m.max_twist = tbl.max_twist_ist(idx);
    m.mean_twist = tbl.mean_twist_ist(idx);
    m.median_twist = tbl.median_twist_ist(idx);
    m.std_twist = tbl.std_twist_ist(idx);
    m.min_acceleration = tbl.min_acceleration_ist(idx);
    m.max_acceleration = tbl.max_acceleration_ist(idx);
    m.mean_acceleration = tbl.mean_acceleration_ist(idx);
    m.median_acceleration = tbl.median_acceleration_ist(idx);
    m.std_acceleration = tbl.std_acceleration_ist(idx);
end

function m = getMetadataStructFromTable(tbl, idx)
    m.movement_type = tbl.movement_type{idx};
    m.length = tbl.length(idx);
    m.duration = tbl.duration(idx);
    m.min_twist = tbl.min_twist_ist(idx);
    m.max_twist = tbl.max_twist_ist(idx);
    m.mean_twist = tbl.mean_twist_ist(idx);
    m.median_twist = tbl.median_twist_ist(idx);
    m.std_twist = tbl.std_twist_ist(idx);
    m.min_acceleration = tbl.min_acceleration_ist(idx);
    m.max_acceleration = tbl.max_acceleration_ist(idx);
    m.mean_acceleration = tbl.mean_acceleration_ist(idx);
    m.median_acceleration = tbl.median_acceleration_ist(idx);
    m.std_acceleration = tbl.std_acceleration_ist(idx);
end
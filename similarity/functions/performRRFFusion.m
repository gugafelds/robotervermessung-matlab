function [top_k_ids, top_k_scores, rrf_ranking] = performRRFFusion(query_embeddings, candidate_embeddings, weights, K, candidate_bahn_ids)
    % Extract embeddings and compute similarities for each modality
    % Then call fuseRankingsRRF
    
    rankings = struct();
    
    % Position
    if weights(1) > 0
        pos_similarities = candidate_embeddings.position * query_embeddings.position';
        [~, sort_idx] = sort(pos_similarities, 'descend');
        rankings.position = table(candidate_bahn_ids(sort_idx), ...
            'VariableNames', {'bahn_id'});
    end
    
    % Joint
    if weights(2) > 0
        joint_similarities = candidate_embeddings.joint * query_embeddings.joint';
        [~, sort_idx] = sort(joint_similarities, 'descend');
        rankings.joint = table(candidate_bahn_ids(sort_idx), ...
            'VariableNames', {'bahn_id'});
    end
    
    % Orientation
    if weights(3) > 0
        orient_similarities = candidate_embeddings.orientation * query_embeddings.orientation';
        [~, sort_idx] = sort(orient_similarities, 'descend');
        rankings.orientation = table(candidate_bahn_ids(sort_idx), ...
            'VariableNames', {'bahn_id'});
    end
    
    % Velocity
    if weights(4) > 0
        vel_similarities = candidate_embeddings.velocity * query_embeddings.velocity';
        [~, sort_idx] = sort(vel_similarities, 'descend');
        rankings.velocity = table(candidate_bahn_ids(sort_idx), ...
            'VariableNames', {'bahn_id'});
    end
    
    % Metadata
    if weights(5) > 0
        meta_similarities = candidate_embeddings.metadata * query_embeddings.metadata';
        [~, sort_idx] = sort(meta_similarities, 'descend');
        rankings.metadata = table(candidate_bahn_ids(sort_idx), ...
            'VariableNames', {'bahn_id'});
    end
    
    % RRF fusion
    weights_struct = struct();
    weights_struct.position = weights(1);
    weights_struct.joint = weights(2);
    weights_struct.orientation = weights(3);
    weights_struct.velocity = weights(4);
    weights_struct.metadata = weights(5);
    
    rrf_ranking = fuseRankingsRRF(rankings, weights_struct, 60, 'bahn_id');
    
    % Extract Top-K
    top_k_ids = rrf_ranking.bahn_id(1:min(K, height(rrf_ranking)));
    top_k_scores = rrf_ranking.rrf_score(1:min(K, height(rrf_ranking)));
end
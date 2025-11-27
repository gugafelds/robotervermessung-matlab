function embedding = computePositionEmbedding(pos_data, n_coarse, n_fine)
    % computePositionEmbedding - Position embedding mit Multi-Scale Resampling
    %
    % Syntax: embedding = computePositionEmbedding(pos_data, n_coarse, n_fine)
    %
    % Input:
    %   pos_data - Nx3 matrix [x, y, z]
    %   n_coarse - Coarse resampling points (default: 50)
    %   n_fine - Fine resampling points (default: 250)
    %
    % Output:
    %   embedding - (n_coarse + n_fine) * 3 dimensional vector
    
    if nargin < 2
        n_coarse = 50;
    end
    if nargin < 3
        n_fine = 250;
    end
    
    % ⭐ Handle empty data
    if isempty(pos_data) || size(pos_data, 1) == 0
        total_dims = (n_coarse + n_fine) * 3;
        embedding = zeros(1, total_dims);
        return;
    end
    
    % ⭐ Normalisierung zum Startpunkt (translationsinvariant)
    traj_normalized = pos_data - pos_data(1, :);
    
    % ⭐ Scale-Normalisierung (größeninvariant)
    norms = sqrt(sum(traj_normalized.^2, 2));  % L2-Norm pro Punkt
    max_extent = max(norms);
    if max_extent > 1e-6
        traj_normalized = traj_normalized / max_extent;
    end
    
    % ⭐ Multi-Scale Resampling
    embedding_parts = {};
    
    % Coarse: Grobe Form
    if n_coarse > 0
        coarse = resampleTrajectory(traj_normalized, n_coarse);
        embedding_parts{end+1} = coarse(:)';  % Flatten to row vector
    end
    
    % Fine: Details
    if n_fine > 0
        fine = resampleTrajectory(traj_normalized, n_fine);
        embedding_parts{end+1} = fine(:)';  % Flatten to row vector
    end
    
    % ⭐ Concatenate all parts
    if isempty(embedding_parts)
        error('Both n_coarse and n_fine are 0 or negative - cannot create embedding');
    end
    
    embedding = [embedding_parts{:}];
    
    % ⭐ L2-Normalisierung für Cosine-Similarity
    norm_val = norm(embedding);
    if norm_val > 1e-10
        embedding = embedding / norm_val;
    end
end
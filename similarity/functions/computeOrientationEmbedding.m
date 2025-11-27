function embedding = computeOrientationEmbedding(orient_data, n_coarse, n_fine)
    % computeOrientationEmbedding - Orientation embedding mit Multi-Scale Resampling
    %
    % Syntax: embedding = computeOrientationEmbedding(orient_data, n_coarse, n_fine)
    %
    % Input:
    %   orient_data - Nx4 matrix [qx, qy, qz, qw]
    %   n_coarse - Coarse resampling points (default: 25)
    %   n_fine - Fine resampling points (default: 75)
    %
    % Output:
    %   embedding - (n_coarse + n_fine) * 3 dimensional vector (rotation vectors)
    
    if nargin < 2
        n_coarse = 25;
    end
    if nargin < 3
        n_fine = 75;
    end
    
    % ⭐ Handle empty data
    if isempty(orient_data) || size(orient_data, 1) == 0
        total_dims = (n_coarse + n_fine) * 3;  % 3 axes (rotation vector)
        embedding = zeros(1, total_dims);
        return;
    end
    
    % Extract quaternion components
    qx = orient_data(:, 1);
    qy = orient_data(:, 2);
    qz = orient_data(:, 3);
    qw = orient_data(:, 4);
    
    % Convert to rotation vectors (axis-angle representation)
    angles = 2 * acos(max(-1, min(1, qw)));  % Clamp für numerische Stabilität
    sin_half = sin(angles / 2);
    sin_half(sin_half == 0) = 1;  % Avoid division by zero
    
    rot_vectors = zeros(size(orient_data, 1), 3);
    rot_vectors(:, 1) = qx ./ sin_half;  % rot_x
    rot_vectors(:, 2) = qy ./ sin_half;  % rot_y
    rot_vectors(:, 3) = qz ./ sin_half;  % rot_z
    
    % ⭐ Normalisierung zur ersten Orientation (rotationsinvariant)
    rot_vectors_normalized = rot_vectors - rot_vectors(1, :);
    
    % ⭐ Multi-Scale Resampling
    embedding_parts = {};
    
    % Coarse: Grobe Orientierungsänderungen
    if n_coarse > 0
        coarse = resampleTrajectory(rot_vectors_normalized, n_coarse);
        embedding_parts{end+1} = coarse(:)';  % Flatten to row vector
    end
    
    % Fine: Detaillierte Orientierungsänderungen
    if n_fine > 0
        fine = resampleTrajectory(rot_vectors_normalized, n_fine);
        embedding_parts{end+1} = fine(:)';  % Flatten to row vector
    end
    
    % ⭐ Concatenate all parts
    if isempty(embedding_parts)
        error('Both n_coarse and n_fine are 0 or negative - cannot create embedding');
    end
    
    embedding = [embedding_parts{:}];
    
    % ⭐ L2-Normalisierung
    norm_val = norm(embedding);
    if norm_val > 1e-10
        embedding = embedding / norm_val;
    end
end
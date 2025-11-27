function embedding = computeJointEmbedding(joint_data, n_coarse, n_fine)
    % computeJointEmbedding - Joint embedding mit Multi-Scale Resampling
    %
    % Syntax: embedding = computeJointEmbedding(joint_data, n_coarse, n_fine)
    %
    % Input:
    %   joint_data - Nx6 matrix [j1, j2, j3, j4, j5, j6]
    %   n_coarse - Coarse resampling points (default: 50)
    %   n_fine - Fine resampling points (default: 250)
    %
    % Output:
    %   embedding - (n_coarse + n_fine) * 6 dimensional vector
    
    if nargin < 2
        n_coarse = 50;
    end
    if nargin < 3
        n_fine = 250;
    end
    
    % ⭐ Handle empty data
    if isempty(joint_data) || size(joint_data, 1) == 0
        total_dims = (n_coarse + n_fine) * 6;  % 6 joints
        embedding = zeros(1, total_dims);
        return;
    end
    
    % ⭐ Optional: Wrap angles to [-pi, pi] für zyklische Joints
    % Uncomment wenn deine Joints zyklisch sind:
    % joint_normalized = wrapToPi(joint_normalized);
    
    % ⭐ Multi-Scale Resampling
    embedding_parts = {};
    
    % Coarse: Grobe Gelenkbewegungen
    if n_coarse > 0
        coarse = resampleTrajectory(joint_data, n_coarse);
        embedding_parts{end+1} = coarse(:)';  % Flatten to row vector
    end
    
    % Fine: Detaillierte Gelenkbewegungen
    if n_fine > 0
        fine = resampleTrajectory(joint_data, n_fine);
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
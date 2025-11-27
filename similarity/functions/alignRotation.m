function [seq1_aligned, seq2_aligned] = alignRotation(seq1, seq2)
    % Procrustes alignment: Find optimal rotation between two trajectories
    % Makes DTW rotation-invariant
    %
    % Input:
    %   seq1, seq2: N×3 or M×3 matrices (can have different lengths!)
    % Output:
    %   seq1_aligned, seq2_aligned: Optimally rotated sequences (original lengths preserved)
    
    if isempty(seq1) || isempty(seq2)
        seq1_aligned = seq1;
        seq2_aligned = seq2;
        return;
    end
    
    % Safety: sequences must have at least 3 points for rotation
    if size(seq1, 1) < 3 || size(seq2, 1) < 3
        seq1_aligned = seq1;
        seq2_aligned = seq2;
        return;
    end
    
    n1 = size(seq1, 1);
    n2 = size(seq2, 1);
    
    if n1 ~= n2
        seq2_resampled = zeros(n1, 3);
        for dim = 1:3
            seq2_resampled(:, dim) = interp1(1:n2, seq2(:, dim), ...
                linspace(1, n2, n1), 'pchip');
        end
    else
        seq2_resampled = seq2;
    end
    
    % === 1. CENTER ===
    seq1_centered = seq1 - mean(seq1, 1);
    seq2_centered = seq2_resampled - mean(seq2_resampled, 1);
    
    % === 2. SCALE NORMALIZATION (for numerical stability) ===
    scale1 = sqrt(sum(sum(seq1_centered.^2)));
    scale2 = sqrt(sum(sum(seq2_centered.^2)));
    
    if scale1 > 1e-6
        seq1_normalized = seq1_centered / scale1;
    else
        seq1_normalized = seq1_centered;
    end
    
    if scale2 > 1e-6
        seq2_normalized = seq2_centered / scale2;
    else
        seq2_normalized = seq2_centered;
    end
    
    % === 3. FIND OPTIMAL ROTATION (SVD method) ===
    % Now both sequences have same length: n1 × 3
    H = seq2_normalized' * seq1_normalized;  % ✅ (3 × n1) * (n1 × 3) = (3 × 3)
    
    % Singular Value Decomposition
    [U, ~, V] = svd(H);
    
    % Optimal rotation matrix
    R = V * U';
    
    % Handle reflection case (det(R) = -1)
    if det(R) < 0
        V(:, end) = -V(:, end);
        R = V * U';
    end
    
    % === 4. APPLY ROTATION TO ORIGINAL seq2 (not resampled!) ===
    % Center original seq2
    seq2_original_centered = seq2 - mean(seq2, 1);
    
    % Normalize
    if scale2 > 1e-6
        seq2_original_normalized = seq2_original_centered / scale2;
    else
        seq2_original_normalized = seq2_original_centered;
    end
    
    % Rotate
    seq2_rotated = (R * seq2_original_normalized')';  % Apply rotation

    seq1_aligned = seq1_normalized * scale1;
    seq2_aligned = seq2_rotated * scale2;
end
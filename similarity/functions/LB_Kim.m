function lb_dist = LB_Kim(seq1, seq2, mode, align_rotation, normalize_dtw)
    % LB_Kim: Ultra-fast lower bound using only 4 critical points
    % Complexity: O(1) vs O(n) for LB_Keogh
    %
    % Args:
    %   seq1: Query sequence (N × D)
    %   seq2: Candidate sequence (M × D)
    %   mode: 'position' or 'joint_states'
    %   align_rotation: Align rotation before comparison (default: false)
    %   normalize_dtw: Apply DTW normalization (default: false)
    %
    % Returns:
    %   lb_dist: Lower bound distance
    
    if nargin < 5
        normalize_dtw = false;
    end
    if nargin < 4
        align_rotation = false;
    end
    if nargin < 3
        mode = 'position';
    end
    
    if isempty(seq1) || isempty(seq2)
        lb_dist = inf;
        return;
    end
    
    % Safety check
    if size(seq1, 1) < 2 || size(seq2, 1) < 2
        lb_dist = inf;
        return;
    end
    
    % Apply normalization if requested (same as LB_Keogh)
    if normalize_dtw || strcmp(mode, 'joint_states')
        seq1 = normalizeForDTW(seq1, mode);
        seq2 = normalizeForDTW(seq2, mode);
    end
    
    % Resample seq2 to match seq1 length (for fair comparison)
    n1 = size(seq1, 1);
    n2 = size(seq2, 1);
    d = size(seq1, 2);
    
    if n1 ~= n2
        seq2_resampled = zeros(n1, d);
        for dim = 1:d
            seq2_resampled(:, dim) = interp1(1:n2, seq2(:, dim), ...
                linspace(1, n2, n1), 'pchip');
        end
        seq2 = seq2_resampled;
    end
    
    % Apply rotation alignment if requested (same as LB_Keogh)
    if align_rotation && strcmp(mode, 'position')
        [seq1, seq2] = alignRotation(seq1, seq2);
    end
    
    % ========================================================================
    % Compute distances for 4 critical points
    % ========================================================================
    
    % 1. First point
    d1 = computeDistance(seq1(1, :), seq2(1, :), mode);
    
    % 2. Last point
    d2 = computeDistance(seq1(end, :), seq2(end, :), mode);
    
    % 3. Minimum points (per dimension)
    [~, min_idx1] = min(seq1, [], 1);
    [~, min_idx2] = min(seq2, [], 1);
    d3 = 0;
    for dim = 1:d
        d3 = d3 + computeDistance(seq1(min_idx1(dim), :), seq2(min_idx2(dim), :), mode);
    end
    
    % 4. Maximum points (per dimension)
    [~, max_idx1] = max(seq1, [], 1);
    [~, max_idx2] = max(seq2, [], 1);
    d4 = 0;
    for dim = 1:d
        d4 = d4 + computeDistance(seq1(max_idx1(dim), :), seq2(max_idx2(dim), :), mode);
    end
    
    % Sum all distances
    lb_dist = d1 + d2 + d3 + d4;
end
function lb_dist = LB_Kim(seq1, seq2, mode, align_rotation, normalize_dtw)
    % LB_Kim: Ultra-fast lower bound using only 4 critical points
    % Based on Kim et al. (2001) - uses first, last, min, max of sequences
    % Complexity: O(d) for finding min/max, effectively O(1) compared to O(n) for LB_Keogh
    %
    % Args:
    %   seq1: Query sequence (N × D)
    %   seq2: Candidate sequence (M × D)
    %   mode: 'position' or 'joint_states'
    %   align_rotation: Align rotation before comparison (default: false)
    %   normalize_dtw: Apply DTW normalization (default: false)
    %
    % Returns:
    %   lb_dist: Lower bound distance (square root of max of 4 squared differences)
    
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
    
    % Apply normalization if requested
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
    
    % Apply rotation alignment if requested
    if align_rotation && strcmp(mode, 'position')
        [seq1, seq2] = alignRotation(seq1, seq2);
    end
    
    % ========================================================================
    % Compute LB_Kim: max of 4 feature distances (vectorized)
    % ========================================================================
    
    % Extract 4 features for each sequence (vectorized across all dimensions)
    % 1. First point
    first1 = seq1(1, :);
    first2 = seq2(1, :);
    
    % 2. Last point
    last1 = seq1(end, :);
    last2 = seq2(end, :);
    
    % 3. Minimum values (per dimension)
    min1 = min(seq1, [], 1);
    min2 = min(seq2, [], 1);
    
    % 4. Maximum values (per dimension)
    max1 = max(seq1, [], 1);
    max2 = max(seq2, [], 1);
    
    % Compute distances for each feature
    if strcmp(mode, 'joint_states')
        % Angular distance with wrap-around
        d1 = sum(angularDistance(first1, first2));
        d2 = sum(angularDistance(last1, last2));
        d3 = sum(angularDistance(min1, min2));
        d4 = sum(angularDistance(max1, max2));
    else
        % Standard Euclidean distance (squared)
        d1 = sum((first1 - first2).^2);
        d2 = sum((last1 - last2).^2);
        d3 = sum((min1 - min2).^2);
        d4 = sum((max1 - max2).^2);
    end
    
    % LB_Kim = sqrt(max of the 4 squared distances)
    lb_dist = sqrt(max([d1, d2, d3, d4]));
end

% Helper function for angular distance
function dist_squared = angularDistance(val1, val2)
    % Compute squared angular distance with wrap-around
    diff = abs(val1 - val2);
    diff = min(diff, 2*pi - diff);
    dist_squared = diff.^2;
end
function lb_dist = LB_Keogh(query, candidate, window_percent, mode, align_rotation, normalize_dtw)
    if nargin < 6
        normalize_dtw = false;
    end
    if nargin < 5
        align_rotation = false;
    end
    if nargin < 4
        mode = 'position';
    end
    if nargin < 3
        window_percent = 0.15;
    end

    if normalize_dtw || strcmp(mode, 'joint_states')
        query = normalizeForDTW(query, mode);
        candidate = normalizeForDTW(candidate, mode);
    end
    
    n = size(query, 1);
    m = size(candidate, 1);
    d = size(query, 2);

    if n ~= m
        candidate_resampled = zeros(n, d);
        for dim = 1:d
            candidate_resampled(:, dim) = interp1( ...
                1:m, candidate(:, dim), linspace(1, m, n), 'pchip');
        end
        candidate = candidate_resampled;
    end

    if align_rotation && strcmp(mode, 'position')
        [query, candidate] = alignRotation(query, candidate);
    end
    
    % --- Envelope construction ---
    window = round(n * window_percent);
    lower_env = zeros(n, d);
    upper_env = zeros(n, d);

    for dim = 1:d
        for i = 1:n
            start_idx = max(1, i - window);
            end_idx = min(n, i + window);
            lower_env(i, dim) = min(candidate(start_idx:end_idx, dim));
            upper_env(i, dim) = max(candidate(start_idx:end_idx, dim));
        end
    end
    
    % --- Compute LB_Keogh distance (vectorized) ---
    if strcmp(mode, 'joint_states')
        % Wrap-around angular distance
        % Points above upper envelope
        above_mask = query > upper_env;
        diff_above = abs(query - upper_env);
        diff_above = min(diff_above, 2*pi - diff_above);
        
        % Points below lower envelope
        below_mask = query < lower_env;
        diff_below = abs(query - lower_env);
        diff_below = min(diff_below, 2*pi - diff_below);
        
        % Sum squared distances
        lb_dist_squared = sum(above_mask .* diff_above.^2, 'all') + ...
                          sum(below_mask .* diff_below.^2, 'all');
    else
        % Standard Euclidean
        % Points above upper envelope
        above_diff = query - upper_env;
        above_diff(above_diff < 0) = 0;
        
        % Points below lower envelope
        below_diff = lower_env - query;
        below_diff(below_diff < 0) = 0;
        
        % Sum squared distances
        lb_dist_squared = sum(above_diff.^2, 'all') + sum(below_diff.^2, 'all');
    end
    
    lb_dist = sqrt(lb_dist_squared);
end

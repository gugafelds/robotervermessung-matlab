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
    
    % --- Compute LB_Keogh distance ---
    lb_dist_squared = 0;

    for dim = 1:d
        for i = 1:n
            q_val = query(i, dim);

            if strcmp(mode, 'joint_states')
                % Wrap-around angular distance
                if q_val > upper_env(i, dim)
                    diff = abs(q_val - upper_env(i, dim));
                    diff = min(diff, 2*pi - diff);
                    lb_dist_squared = lb_dist_squared + diff^2;

                elseif q_val < lower_env(i, dim)
                    diff = abs(q_val - lower_env(i, dim));
                    diff = min(diff, 2*pi - diff);
                    lb_dist_squared = lb_dist_squared + diff^2;
                end

            else
                % Standard Euclidean
                if q_val > upper_env(i, dim)
                    lb_dist_squared = lb_dist_squared + (q_val - upper_env(i, dim))^2;
                elseif q_val < lower_env(i, dim)
                    lb_dist_squared = lb_dist_squared + (q_val - lower_env(i, dim))^2;
                end
            end
        end
    end
    
    lb_dist = sqrt(lb_dist_squared);
end

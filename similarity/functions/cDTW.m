function dist = cDTW(seq1, seq2, mode, window_percent, best_so_far, align_rotation, normalize_dtw)
    if nargin < 7
        normalize_dtw = false;
    end
    if nargin < 6
        align_rotation = false;
    end
    if nargin < 5
        best_so_far = inf;
    end
    if nargin < 4
        window_percent = 0.15;
    end
    if nargin < 3
        mode = 'position';
    end

    if strcmp(mode, 'position')
    % Remove translation (set start to origin)
    seq1 = seq1 - seq1(1, :);
    seq2 = seq2 - seq2(1, :);
    end
    
    % Apply normalization AFTER translation removal
    if normalize_dtw && strcmp(mode, 'position')
        seq1 = normalizeForDTW(seq1);
        seq2 = normalizeForDTW(seq2);
    end
    
    % Apply rotation alignment AFTER normalization
    if align_rotation && strcmp(mode, 'position')
        [seq1, seq2] = alignRotation(seq1, seq2);
    end

    % --- DTW core unchanged ---
    n = size(seq1, 1);
    m = size(seq2, 1);
    window = round(max(n, m) * window_percent);
    window = max(window, abs(n - m));
    dtw_matrix = inf(n+1, m+1);
    dtw_matrix(1, 1) = 0;

    for i = 1:n
        j_center = round(i * m / n);
        j_start = max(1, j_center - window);
        j_end = min(m, j_center + window);

        min_cost_in_row = inf;
        for j = j_start:j_end
            if strcmp(mode, 'joint_states')
                angular_diffs = abs(seq1(i, :) - seq2(j, :));
                angular_diffs = min(angular_diffs, 2*pi - angular_diffs);
                cost = norm(angular_diffs);
            else
                cost = norm(seq1(i, :) - seq2(j, :));
            end

            dtw_matrix(i+1, j+1) = cost + min([ ...
                dtw_matrix(i, j+1), ...
                dtw_matrix(i+1, j), ...
                dtw_matrix(i, j) ...
            ]);

            min_cost_in_row = min(min_cost_in_row, dtw_matrix(i+1, j+1));
        end

        if min_cost_in_row > best_so_far
            dist = best_so_far;
            return;
        end
    end

    dist = dtw_matrix(n+1, m+1);
end




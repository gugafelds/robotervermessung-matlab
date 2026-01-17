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
    
    % ================================================================
    % NORMALIZATION
    % ================================================================
    if normalize_dtw
        % For joint_states: ALWAYS normalize (regardless of flag)
        % For position: only if flag is set
        seq1 = normalizeForDTW(seq1, mode);
        seq2 = normalizeForDTW(seq2, mode);
    end
    
    % ================================================================
    % POSITION-SPECIFIC: Translation removal
    % ================================================================
    if strcmp(mode, 'position')
        % Remove translation AFTER normalization
        % (or before, depending on your preference)
        seq1 = seq1 - seq1(1, :);
        seq2 = seq2 - seq2(1, :);
    end
    
    % ================================================================
    % POSITION-SPECIFIC: Rotation alignment
    % ================================================================
    if align_rotation && strcmp(mode, 'position')
        [seq1, seq2] = alignRotation(seq1, seq2);
    end
    
    % ================================================================
    % DTW CORE
    % ================================================================
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
            % After normalization: simple Euclidean distance
            % (no special handling for joint_states needed)
            cost = norm(seq1(i, :) - seq2(j, :));
            
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
    
    % ================================================================
    % LENGTH NORMALIZATION (for both modes)
    % ================================================================
    %avg_length = (n + m) / 2;
    %dist = dist / avg_length;
end
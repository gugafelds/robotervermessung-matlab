function normalized_seq = normalizeForDTW(seq)
    % Z-normalize a trajectory for DTW
    % Makes DTW translation & scale invariant

    %if isempty(seq) || size(seq, 1) < 2
    %    normalized_seq = seq;
    %    return;
    %end

    %mu = mean(seq, 1);
    %sigma = std(seq, 0, 1);

    % Avoid division by zero for flat dimensions
    %sigma(sigma < 1e-9) = 1;

    %normalized_seq = (seq - mu) ./ sigma;
%end

% Max-extent normalization (same as embeddings)
    % Scales to [0, 1] based on range per dimension
    
    if isempty(seq) || size(seq, 1) < 2
        normalized_seq = seq;
        return;
    end
    
    min_vals = min(seq, [], 1);
    max_vals = max(seq, [], 1);
    range_vals = max_vals - min_vals;
    
    % Avoid division by zero
    range_vals(range_vals < 1e-9) = 1;
    
    normalized_seq = (seq - min_vals) ./ range_vals;
end
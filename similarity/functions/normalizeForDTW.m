function normalized_seq = normalizeForDTW(seq)
    % Z-normalize a trajectory for DTW
    % Makes DTW translation & scale invariant

    if isempty(seq) || size(seq, 1) < 2
        normalized_seq = seq;
        return;
    end

    mu = mean(seq, 1);
    sigma = std(seq, 0, 1);

    % Avoid division by zero for flat dimensions
    sigma(sigma < 1e-9) = 1;

    normalized_seq = (seq - mu) ./ sigma;
end
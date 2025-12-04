function normalized_seq = normalizeForDTW(seq, mode)
    % NORMALIZEFORDTW - Normalize trajectory for DTW computation
    %
    % Normalizes sequences to make DTW distance scale-invariant and
    % ensure all dimensions contribute equally to the distance metric.
    %
    % INPUTS:
    %   seq  - Trajectory sequence (N × D matrix)
    %          For position: [x, y, z, ...]
    %          For joint_states: [j1, j2, j3, ...]
    %   mode - 'position' or 'joint_states'
    %
    % OUTPUT:
    %   normalized_seq - Normalized sequence (N × D matrix)
    %                    All values scaled to [0, 1] per dimension
    %
    % NORMALIZATION STRATEGY:
    %   - Position mode: Max-extent normalization (scale + translation invariant)
    %   - Joint_states mode: Min-Max per joint (ensures equal joint contribution)
    %
    % EXAMPLES:
    %   % Position normalization
    %   pos_norm = normalizeForDTW(position_data, 'position');
    %
    %   % Joint angle normalization
    %   joint_norm = normalizeForDTW(joint_data, 'joint_states');
    %
    % LITERATURE REFERENCE:
    %   - Z-normalization is standard for DTW (Keogh & Ratanamahatana, 2005)
    %   - Min-Max normalization for biomechanics (Blazkiewicz et al., 2021)
    %   - Range normalization for multi-dimensional data with different scales
    %
    % Author: Gustavo Barros
    % Date: 04.12.2025
    
    if nargin < 2
        mode = 'position';  % Default to position mode
    end
    
    % Handle empty or single-point sequences
    if isempty(seq) || size(seq, 1) < 2
        normalized_seq = seq;
        return;
    end
    
    % ====================================================================
    % MODE-SPECIFIC NORMALIZATION
    % ====================================================================
    
    switch lower(mode)
        case 'position'
            % ============================================================
            % POSITION MODE: Max-extent normalization
            % ============================================================
            % Scales each dimension to [0, 1] based on range
            % Makes DTW translation & scale invariant
            % Consistent with position embedding normalization
            
            min_vals = min(seq, [], 1);
            max_vals = max(seq, [], 1);
            range_vals = max_vals - min_vals;
            
            % Avoid division by zero for flat dimensions
            range_vals(range_vals < 1e-9) = 1;
            
            % Normalize to [0, 1]
            normalized_seq = (seq - min_vals) ./ range_vals;
            
        case 'joint_states'
            % ============================================================
            % JOINT_STATES MODE: Per-joint Min-Max normalization
            % ============================================================
            % Normalizes each joint dimension to [0, 1]
            % Ensures all joints contribute equally to DTW distance
            % Critical for joints with different ROM (Range of Motion)
            % Consistent with joint embedding normalization
            
            min_vals = min(seq, [], 1);
            max_vals = max(seq, [], 1);
            range_vals = max_vals - min_vals;
            
            % Avoid division by zero for joints with no movement
            % (e.g., locked joints or static configurations)
            range_vals(range_vals < 1e-9) = 1;
            
            % Normalize to [0, 1]
            normalized_seq = (seq - min_vals) ./ range_vals;
            
        otherwise
            error('Unknown mode "%s". Must be "position" or "joint_states"', mode);
    end
    
    % ====================================================================
    % ALTERNATIVE: Z-Score Normalization (commented out)
    % ====================================================================
    % Z-Score is the classic DTW normalization from literature
    % Scales to mean=0, std=1 (translation & scale invariant)
    % Problem: Can be unstable for joints with little movement
    %
    % if false  % Enable if you want Z-score instead
    %     mu = mean(seq, 1);
    %     sigma = std(seq, 0, 1);
    %     sigma(sigma < 1e-9) = 1;  % Avoid division by zero
    %     normalized_seq = (seq - mu) ./ sigma;
    % end
end
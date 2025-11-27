function embedding = computeVelocityEmbedding(pos_data, n_coarse, n_fine)
    % computeVelocityEmbedding - Velocity embedding mit Savitzky-Golay Glättung
    %
    % Syntax: embedding = computeVelocityEmbedding(pos_data, n_coarse, n_fine)
    %
    % Input:
    %   pos_data - Nx3 matrix [x, y, z]
    %   n_coarse - Coarse resampling points (default: 25)
    %   n_fine - Fine resampling points (default: 75)
    %
    % Output:
    %   embedding - (n_coarse + n_fine) * 3 dimensional vector
    
    if nargin < 2
        n_coarse = 25;
    end
    if nargin < 3
        n_fine = 75;
    end
    
    % ⭐ Handle empty data
    if isempty(pos_data) || size(pos_data, 1) < 10
        total_dims = (n_coarse + n_fine) * 3;
        embedding = zeros(1, total_dims);
        return;
    end
    
    n_points = size(pos_data, 1);
    
    % ========================================================================
    % ⭐ Künstliche Zeitstempel (gleichmäßig verteilt von 0 bis 1)
    % ========================================================================
    timestamps = linspace(0, 1, n_points)';
    
    % ========================================================================
    % ⭐ Savitzky-Golay Filter für Position (glättet vor Ableitung!)
    % ========================================================================
    window_length = min(33, floor(n_points / 2) * 2 + 1);  % Muss ungerade sein
    if window_length < 5
        window_length = 5;
    end
    
    positions_smooth = zeros(size(pos_data));
    for dim = 1:3
        positions_smooth(:, dim) = sgolayfilt(pos_data(:, dim), 3, window_length);
    end
    
    % ========================================================================
    % ✅ Velocity aus geglätteter Position berechnen
    % ========================================================================
    delta_pos = diff(positions_smooth, 1, 1);  % (N-1) × 3
    delta_time = diff(timestamps);  % (N-1) × 1
    
    % Verhindere Division durch Null
    delta_time(delta_time == 0) = 1e-9;
    
    % Velocity berechnen (für jede Dimension separat)
    velocity = zeros(size(delta_pos));
    for dim = 1:3
        velocity(:, dim) = delta_pos(:, dim) ./ delta_time;
    end
    
    % ========================================================================
    % ⭐ Velocity nochmal glätten (stabilere Werte)
    % ========================================================================
    vel_window = min(33, floor(size(velocity, 1) / 2) * 2 + 1);
    if vel_window < 5
        vel_window = 5;
    end
    
    velocity_smooth = zeros(size(velocity));
    for dim = 1:3
        velocity_smooth(:, dim) = sgolayfilt(velocity(:, dim), 2, vel_window);
    end
    
    % ========================================================================
    % ⭐ Multi-Scale Resampling (Coarse + Fine)
    % ========================================================================
    embedding_parts = {};
    
    % Coarse: Grobe Geschwindigkeitsänderungen
    if n_coarse > 0
        vel_coarse = resampleTrajectory(velocity_smooth, n_coarse);
        embedding_parts{end+1} = vel_coarse(:)';  % Flatten to row vector
    end
    
    % Fine: Detaillierte Geschwindigkeitsänderungen
    if n_fine > 0
        vel_fine = resampleTrajectory(velocity_smooth, n_fine);
        embedding_parts{end+1} = vel_fine(:)';  % Flatten to row vector
    end
    
    % ⭐ Concatenate all parts
    if isempty(embedding_parts)
        error('Both n_coarse and n_fine are 0 or negative - cannot create embedding');
    end
    
    embedding = [embedding_parts{:}];
    
    % ⭐ L2-Normalisierung
    norm_val = norm(embedding);
    if norm_val > 1e-10
        embedding = embedding / norm_val;
    end
end
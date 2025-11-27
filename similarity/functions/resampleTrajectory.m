function resampled = resampleTrajectory(trajectory, n_samples)
    % resampleTrajectory - Resample trajectory mit Interpolation
    %
    % Input:
    %   trajectory - Nx3 matrix
    %   n_samples - Anzahl gewünschter Samples
    %
    % Output:
    %   resampled - n_samples x 3 matrix
    
    n_points = size(trajectory, 1);
    n_dims = size(trajectory, 2);
    
    % Padding wenn zu kurz
    if n_points <= n_samples
        pad_length = n_samples - n_points;
        % Wiederhole letzten Punkt (mode='edge' in numpy)
        last_point = repmat(trajectory(end, :), pad_length, 1);
        resampled = [trajectory; last_point];
        return;
    end
    
    % ⭐ Interpolation (wie scipy.interpolate.interp1d)
    x_old = linspace(0, 1, n_points);
    x_new = linspace(0, 1, n_samples);
    
    resampled = zeros(n_samples, n_dims);
    for dim = 1:n_dims
        % 'pchip' ist robust, alternativ 'linear' wie in Python
        resampled(:, dim) = interp1(x_old, trajectory(:, dim), x_new, 'linear');
    end
end
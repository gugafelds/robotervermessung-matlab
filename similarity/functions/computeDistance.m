function dist = computeDistance(point1, point2, mode)
    % Helper function for mode-aware distance computation
    
    if strcmp(mode, 'joint_states')
        % Angular distance with wrap-around
        angular_diffs = abs(point1 - point2);
        angular_diffs = min(angular_diffs, 2*pi - angular_diffs);
        dist = norm(angular_diffs);
    else
        % Standard Euclidean distance
        dist = norm(point1 - point2);
    end
end
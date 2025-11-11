function q_clean = removeGimbalLockArtifacts(q)
    n = size(q, 1);
    q_clean = q;
    
    % Berechne Winkeländerungen
    angular_changes = zeros(n-1, 1);
    for i = 1:n-1
        q1 = q(i,:) / norm(q(i,:));
        q2 = q(i+1,:) / norm(q(i+1,:));
        dot_prod = abs(dot(q1, q2));
        angular_changes(i) = 2 * acos(min(dot_prod, 1.0));
    end
    
    % Fester Threshold in Grad
    threshold_deg = 1.5;  % Anpassen nach Bedarf
    threshold_rad = deg2rad(threshold_deg);
    
    large_jumps = find(angular_changes > threshold_rad);
    
    %disp(large_jumps');
    
    % Gruppiere und interpoliere
    if ~isempty(large_jumps)
        groups = {};
        current_group = large_jumps(1);
        
        for i = 2:length(large_jumps)
            if large_jumps(i) == large_jumps(i-1) + 1
                current_group = [current_group, large_jumps(i)];
            else
                groups{end+1} = current_group;
                current_group = large_jumps(i);
            end
        end
        groups{end+1} = current_group;
        
        for g = 1:length(groups)
            artifact_indices = groups{g};
            artifact_start = artifact_indices(1) - 1;
            artifact_end = artifact_indices(end) + 10;
            artifact_end = min(artifact_end, n);
            
            q_start = q(artifact_start,:) / norm(q(artifact_start,:));
            q_end = q(artifact_end,:) / norm(q(artifact_end,:));
            
            num_points = artifact_end - artifact_start - 1;
            for k = 1:num_points
                t = k / (num_points + 1);
                q_clean(artifact_start + k, :) = slerp(q_start, q_end, t);
            end
        end
    end
end

function q_interp = slerp(q1, q2, t)
    dot_prod = dot(q1, q2);
    
    if dot_prod < 0
        q2 = -q2;
        dot_prod = -dot_prod;
    end
    
    dot_prod = min(max(dot_prod, -1), 1);
    theta = acos(dot_prod);
    
    if theta < 1e-6
        q_interp = (1-t)*q1 + t*q2;
    else
        q_interp = (sin((1-t)*theta)/sin(theta))*q1 + (sin(t*theta)/sin(theta))*q2;
    end
    
    q_interp = q_interp / norm(q_interp);
end
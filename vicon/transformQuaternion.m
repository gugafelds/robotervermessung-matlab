function q_transformed = transformQuaternion(data_ist, data_soll, q_transform, trafo_rot)
    % Extract trajectory quaternions in [w x y z] format
    q_traj = [data_ist.qw_ist, data_ist.qx_ist, data_ist.qy_ist, data_ist.qz_ist];
    
    % Normalize input quaternions for numerical stability
    q_traj = q_traj ./ sqrt(sum(q_traj.^2, 2));
    
    % Convert to Euler angles first (using ZYX convention for better robot compatibility)
    euler_angles = rad2deg(quat2eul(q_traj, 'XYZ'));
    
    % Apply coordinate system rotation in Euler space
    euler_transformed = euler_angles * trafo_rot;
    
    % Convert back to quaternions
    q_traj_rotated = eul2quat(deg2rad(euler_transformed), 'XYZ');
    
    % Initialize output array with correct size
    q_transformed = zeros(size(q_traj_rotated));  % Direkt das Ausgabe-Array initialisieren
    
    % Get reference quaternion from first SOLL point and normalize
    q_ref = [data_soll.qw_soll(1), data_soll.qx_soll(1), data_soll.qy_soll(1), data_soll.qz_soll(1)];
    q_ref = q_ref / norm(q_ref);
    
    % Apply calibration transformation
    for i = 1:size(q_traj_rotated,1)
        % Apply transformation
        q_temp = quatmultiply(q_traj_rotated(i,:), q_transform);
        
        % Normalize result
        q_temp = q_temp / norm(q_temp);
        
        % Check orientation consistency
        dot_product = dot(q_temp, q_ref);
        if dot_product < 0
            q_temp = -q_temp;
        end
        
        q_transformed(i,:) = q_temp;  % Direkt ins Ausgabe-Array schreiben
    end
    
    % Keine zusätzliche Zuweisung am Ende nötig
end
function q_transform = calibrateQuaternion(data_ist, data_soll)
% Calculates transformation quaternion from calibration data
    
% Extract quaternions from data tables
q_ist = [data_ist.qw_ist, data_ist.qx_ist, data_ist.qy_ist, data_ist.qz_ist];
q_soll = [data_soll.qw_soll, data_soll.qx_soll, data_soll.qy_soll, data_soll.qz_soll];

% Calculate mean quaternions
q_mean_ist = mean(q_ist, 1);
q_mean_ist = q_mean_ist / norm(q_mean_ist);
q_mean_soll = mean(q_soll, 1);
q_mean_soll = q_mean_soll / norm(q_mean_soll);

% Calculate transformation quaternion
q_transform = quatmultiply(quatconj(q_mean_ist), q_mean_soll);

% Save to workspace
assignin('base', 'q_transform', q_transform);
end
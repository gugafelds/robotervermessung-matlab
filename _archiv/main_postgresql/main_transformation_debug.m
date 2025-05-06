%% Einstellungen
clear;

%bahn_id_ = '1738682877';
bahn_id_ = '1739799300';% Orientierungsänderung ohne Kalibrierungsdatei
%bahn_id_ = '1720784405';
transform_only = true;    % Nur Transformation und Plot, kein Upload
plots = true;

schema = 'bewegungsdaten';

% Verbindung mit PostgreSQL
datasource = "RobotervermessungMATLAB";
username = "felixthomas";
password = "manager";
conn = postgresql(datasource,username,password);

% Überprüfe Verbindung
if isopen(conn)
    disp('Verbindung erfolgreich hergestellt');
else
    error('Datenbankverbindung fehlgeschlagen');
end

clear datasource username password

%% Vorbereitung der zu verarbeitenden Bahn-IDs

% Initialisiere leere Arrays
bahn_ids_to_process = bahn_id_;
existing_bahn_ids = [];

query = ['SELECT DISTINCT bahn_id FROM robotervermessung.' schema '.bahn_pose_trans'];
existing_bahn_ids = str2double(table2array(fetch(conn, query)));

%% Hauptverarbeitungsschleife
tic;


current_bahn_id = bahn_id_;


% Suche nach Kalibrierungslauf
disp(['Verarbeite Bahn-ID: ' current_bahn_id]);
[calibration_id, is_calibration_run] = findCalibrationRun(conn, current_bahn_id, schema);

% Extrahieren der Kalibrierungs-Daten für die Position
tablename_cal = ['robotervermessung.' schema '.bahn_pose_ist'];
opts_cal = databaseImportOptions(conn,tablename_cal);
opts_cal.RowFilter = opts_cal.RowFilter.bahn_id == calibration_id;
data_cal_ist = sqlread(conn,tablename_cal,opts_cal);
data_cal_ist = sortrows(data_cal_ist,'timestamp');

tablename_cal = ['robotervermessung.' schema '.bahn_events'];
opts_cal = databaseImportOptions(conn,tablename_cal);
opts_cal.RowFilter = opts_cal.RowFilter.bahn_id == calibration_id;
data_cal_soll = sqlread(conn,tablename_cal,opts_cal);
data_cal_soll = sortrows(data_cal_soll,'timestamp');

% Positionsdaten für Koordinatentransformation
calibration(data_cal_ist,data_cal_soll, false)

% Extract calibration data for orientation
tablename_cal = ['robotervermessung.' schema '.bahn_orientation_soll'];
opts_cal = databaseImportOptions(conn,tablename_cal);
opts_cal.RowFilter = opts_cal.RowFilter.bahn_id == calibration_id;
data_cal_soll = sqlread(conn,tablename_cal,opts_cal);
data_cal_soll = sortrows(data_cal_soll,'timestamp');

% Read current trajectory data
query = ['SELECT * FROM robotervermessung.' schema '.bahn_pose_ist ' ...
        'WHERE robotervermessung.' schema '.bahn_pose_ist.bahn_id = ''' current_bahn_id ''''];
data_ist = fetch(conn, query);
data_ist = sortrows(data_ist,'timestamp');

% Read target orientation data
query = ['SELECT * FROM robotervermessung.' schema '.bahn_orientation_soll ' ...
         'WHERE robotervermessung.' schema '.bahn_orientation_soll.bahn_id = ''' current_bahn_id ''''];
data_soll = fetch(conn, query);
data_soll = sortrows(data_soll,'timestamp');

% Extract timestamps
ist_times = str2double(data_ist.timestamp);
soll_times = str2double(data_soll.timestamp);

% Normalize timestamps
t0 = min(min(ist_times), min(soll_times));
ist_times_norm = (ist_times - t0) * 1e-9;  % convert to seconds
soll_times_norm = (soll_times - t0) * 1e-9;


findNaNInData(data_cal_soll)
% Perform quaternion transformation
[q_transform, q_traj_transformed] = calibrateAndTransformQuaternions(data_cal_ist, data_cal_soll, data_ist, trafo_rot, data_soll);

% Convert to Euler angles
euler_ist = quat2eul(q_traj_transformed);
q_soll = [data_soll.qw_soll, data_soll.qx_soll, data_soll.qy_soll, data_soll.qz_soll];
euler_soll = quat2eul(q_soll);

% Convert to degrees and fix gimbal lock
euler_ist_deg = rad2deg(euler_ist);
euler_soll_deg = rad2deg(euler_soll);
euler_ist_deg_fixed = fixGimbalLock(euler_ist_deg);
euler_soll_deg_fixed = fixGimbalLock(euler_soll_deg);

% Visualization
if plots
    % Plot Euler angles
    figure('Name', 'Euler Angles Comparison', 'Color', 'white');
    
    % Yaw angle
    subplot(3,1,1);
    plot(ist_times_norm, euler_ist_deg_fixed(:,1), '--', 'LineWidth', 1.5); hold on;
    plot(soll_times_norm, euler_soll_deg_fixed(:,1), '-', 'LineWidth', 1.5);
    title('Yaw Angle');
    legend('Transformed', 'Target');
    ylabel('Angle (degrees)');
    grid on;
    
    % Pitch angle
    subplot(3,1,2);
    plot(ist_times_norm, euler_ist_deg_fixed(:,2), '--', 'LineWidth', 1.5); hold on;
    plot(soll_times_norm, euler_soll_deg_fixed(:,2), '-', 'LineWidth', 1.5);
    title('Pitch Angle');
    legend('Transformed', 'Target');
    ylabel('Angle (degrees)');
    grid on;
    
    % Roll angle
    subplot(3,1,3);
    plot(ist_times_norm, euler_ist_deg_fixed(:,3), '--', 'LineWidth', 1.5); hold on;
    plot(soll_times_norm, euler_soll_deg_fixed(:,3), '-', 'LineWidth', 1.5);
    title('Roll Angle');
    legend('Transformed', 'Target');
    xlabel('Time (s)');
    ylabel('Angle (degrees)');
    grid on;
end

% Find matching time pairs for error analysis
max_time_diff = 0.1;  % maximum acceptable time difference in seconds
matched_pairs = [];   % [ist_index, soll_index]
errors = [];         % [yaw_error, pitch_error, roll_error]

for i = 1:length(ist_times_norm)
    [min_diff, idx] = min(abs(soll_times_norm - ist_times_norm(i)));
    if min_diff <= max_time_diff
        matched_pairs = [matched_pairs; i, idx];
        errors = [errors; euler_ist_deg_fixed(i,:) - euler_soll_deg_fixed(idx,:)];
    end
end

% Calculate error statistics
error_mean = mean(abs(errors));
error_std = std(errors);
error_max = max(abs(errors));

% Calculate quaternion errors
quaternion_errors = calculateQuaternionErrors(q_traj_transformed(matched_pairs(:,1),:), ...
                                           q_soll(matched_pairs(:,2),:));
timestamps = ist_times_norm(matched_pairs(:,1));

% Plot quaternion error analysis
plotQuaternionAnalysis(timestamps, quaternion_errors);

% Display statistics
fprintf('\nError Statistics (degrees):\n');
fprintf('          Yaw     Pitch    Roll\n');
fprintf('Mean:   %6.2f  %6.2f  %6.2f\n', error_mean);
fprintf('Std:    %6.2f  %6.2f  %6.2f\n', error_std);
fprintf('Max:    %6.2f  %6.2f  %6.2f\n', error_max);
fprintf('Number of matched points: %d\n', size(matched_pairs,1));

toc;

%% Local Functions
function euler_fixed = fixGimbalLock(euler_angles)
    euler_fixed = euler_angles;
    
    for i = 1:3
        angle_data = euler_angles(:,i);
        near_180 = abs(abs(angle_data) - 180) < 5;
        
        if any(near_180)
            mask_neg = angle_data < 0;
            angle_data(mask_neg) = angle_data(mask_neg) + 360;
            euler_fixed(:,i) = angle_data;
        end
    end
end

function [q_transform, q_traj_transformed] = calibrateAndTransformQuaternions(data_cal_ist, data_cal_soll, data_ist, trafo_rot, data_soll)
    % Main function to handle both calibration and transformation
    
    % Step 1: Get calibration quaternion
    q_transform = calculateCalibrationQuaternion(data_cal_ist, data_cal_soll);
    
    % Step 2: Transform trajectory
    q_traj_transformed = transformTrajectory(data_ist, data_soll, q_transform, trafo_rot);
end

function q_transform = calculateCalibrationQuaternion(data_cal_ist, data_cal_soll)
    % Calculate the calibration quaternion from calibration data
    
    % Extract quaternions in [w x y z] format
    q_ist = [data_cal_ist.qw_ist, data_cal_ist.qx_ist, data_cal_ist.qy_ist, data_cal_ist.qz_ist];
    q_soll = [data_cal_soll.qw_soll, data_cal_soll.qx_soll, data_cal_soll.qy_soll, data_cal_soll.qz_soll];
    
    % Normalize input quaternions
    q_ist = normalizeQuaternions(q_ist);
    q_soll = normalizeQuaternions(q_soll);
    
    % Calculate mean quaternions using advanced averaging
    q_mean_ist = calculateAverageQuaternion(q_ist);
    q_mean_soll = calculateAverageQuaternion(q_soll);
    
    % Calculate transformation quaternion
    % Note: Using quatconj instead of quaternion inverse for better numerical stability
    q_transform = quatmultiply(quatconj(q_mean_ist), q_mean_soll);
    
    % Ensure consistent handedness
    if q_transform(1) < 0
        q_transform = -q_transform;
    end
end

function q_avg = calculateAverageQuaternion(quaternions)
    % Calculate average quaternion using eigenvalue method for better accuracy
    % Form the accumulator matrix
    M = zeros(4,4);
    for i = 1:size(quaternions, 1)
        q = quaternions(i,:);
        M = M + (q' * q);
    end
    M = M / size(quaternions, 1);
    
    % Find eigenvector corresponding to largest eigenvalue
    [V, D] = eig(M);
    [~, maxIdx] = max(diag(D));
    q_avg = V(:,maxIdx)';
    
    % Ensure unit quaternion
    q_avg = q_avg / norm(q_avg);
end

function q_traj_transformed = transformTrajectory(data_ist, data_soll, q_transform, trafo_rot)
    % Extract trajectory quaternions
    q_traj = [data_ist.qw_ist, data_ist.qx_ist, data_ist.qy_ist, data_ist.qz_ist];
    
    % Convert to Euler angles first
    euler_angles = rad2deg(quat2eul(q_traj, 'XYZ'));
    
    % Apply coordinate system rotation in Euler space
    euler_transformed = euler_angles * trafo_rot;
    
    % Convert back to quaternions
    q_traj_rotated = eul2quat(deg2rad(euler_transformed), 'XYZ');
    
    % Initialize output array
    q_traj_transformed = zeros(size(q_traj));
    
    % Get reference quaternion
    q_ref = [data_soll.qw_soll(1), data_soll.qx_soll(1), data_soll.qy_soll(1), data_soll.qz_soll(1)];
    q_ref = q_ref / norm(q_ref);
    
    % Apply calibration transformation
    for i = 1:size(q_traj_rotated,1)
        q_transformed = quatmultiply(q_traj_rotated(i,:), q_transform);
        q_transformed = q_transformed / norm(q_transformed);
        
        % Check orientation consistency
        dot_product = dot(q_transformed, q_ref);
        if dot_product < 0
            q_transformed = -q_transformed;
        end
        
        q_traj_transformed(i,:) = q_transformed;
    end
end

function q_normalized = normalizeQuaternions(q)
    % Normalize quaternions while handling potential numerical issues
    norms = sqrt(sum(q.^2, 2));
    
    % Handle potential zero quaternions
    valid_idx = norms > eps;
    q_normalized = q;
    q_normalized(valid_idx,:) = q(valid_idx,:) ./ norms(valid_idx);
end

function quaternion_errors = calculateQuaternionErrors(q1, q2)
    % Calculate the angular difference between two sets of quaternions
    n_points = size(q1, 1);
    quaternion_errors = zeros(n_points, 1);
    
    for i = 1:n_points
        % Normalize quaternions
        q1_norm = q1(i,:) / norm(q1(i,:));
        q2_norm = q2(i,:) / norm(q2(i,:));
        
        % Calculate dot product and clamp to [-1,1]
        dot_product = abs(sum(q1_norm .* q2_norm));
        dot_product = min(max(dot_product, -1), 1);
        
        % Convert to angle in degrees
        quaternion_errors(i) = 2 * acosd(dot_product);
    end
end

function plotQuaternionAnalysis(timestamps, quaternion_errors)
    figure('Position', [100, 100, 1200, 800]);
    
    % Error time series
    subplot(2,2,[1,2]);
    plot(timestamps, quaternion_errors, 'b-', 'LineWidth', 1.5);
    hold on;
    plot(timestamps, movmean(quaternion_errors, 50), 'r-', 'LineWidth', 2);
    title('Quaternion Angular Error Over Time');
    xlabel('Time (s)');
    ylabel('Angular Error (degrees)');
    grid on;
    legend('Raw Error', 'Moving Average (50 points)');
    
    % Error distribution
    subplot(2,2,3);
    histogram(quaternion_errors, 30, 'Normalization', 'probability');
    title('Error Distribution');
    xlabel('Angular Error (degrees)');
    ylabel('Probability');
    grid on;
    
    % Statistics
    error_stats = struct(...
        'mean', mean(quaternion_errors), ...
        'median', median(quaternion_errors), ...
        'std', std(quaternion_errors), ...
        'max', max(quaternion_errors), ...
        'min', min(quaternion_errors));
    
    subplot(2,2,4);
    text(0.1, 0.8, sprintf('Statistics (degrees):\nMean: %.2f\nMedian: %.2f\nStd Dev: %.2f\nMin/Max: %.2f / %.2f', ...
        error_stats.mean, error_stats.median, error_stats.std, error_stats.min, error_stats.max), ...
        'Units', 'normalized');
    axis off;
    
    sgtitle('Quaternion Transformation Analysis');
end

function findNaNInData(data_cal_soll)
    fprintf('\nChecking for NaN values in data_cal_soll:\n');
    
    % Get indices of rows with any NaN values
    nan_rows = find(any(isnan([data_cal_soll.qw_soll, data_cal_soll.qx_soll, ...
                              data_cal_soll.qy_soll, data_cal_soll.qz_soll]), 2));
    
    if isempty(nan_rows)
        fprintf('No NaN values found.\n');
    else
        fprintf('Found NaN values in rows: %s\n', num2str(nan_rows'));
        fprintf('Total rows with NaN: %d out of %d rows\n', ...
                length(nan_rows), height(data_cal_soll));
    end
end
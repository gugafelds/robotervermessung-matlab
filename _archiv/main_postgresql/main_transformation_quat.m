%% Einstellungen
clear;

%bahn_id_ = '1738682877';
bahn_id_ = '1739274902';% Orientierungsänderung ohne Kalibrierungsdatei
%bahn_id_ = '1720784405';
plots = false;              % Plotten der Daten 
upload_all = false;        % Upload aller Bahnen
upload_single = false;     % Nur eine einzelne Bahn
transform_only = true;    % Nur Transformation und Plot, kein Upload

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
bahn_ids_to_process = [];
existing_bahn_ids = [];

% Hole existierende Bahn-IDs nur wenn Upload geplant ist
if upload_all || upload_single
    query = ['SELECT DISTINCT bahn_id FROM robotervermessung.' schema '.bahn_pose_trans'];
    existing_bahn_ids = str2double(table2array(fetch(conn, query)));
end

% Bestimme zu verarbeitende Bahn-IDs je nach Modus
if upload_all
    % Hole alle verfügbaren Bahn-IDs
    query = ['SELECT bahn_id FROM robotervermessung.' schema '.bahn_info'];
    all_bahn_ids = str2double(table2array(fetch(conn, query)));
    bahn_ids_to_process = all_bahn_ids;
elseif upload_single || transform_only
    % Einzelne Bahn verarbeiten
    if isempty(bahn_id_)
        error('Bitte geben Sie eine bahn_id an');
    end
    bahn_ids_to_process = str2double(bahn_id_);
end

% Informative Ausgabe
fprintf('Modus: %s\n', getModeString(transform_only, upload_single, upload_all));
if ~isempty(bahn_ids_to_process)
    fprintf('Zu verarbeitende Bahn-IDs: %s\n', ...
        strjoin(string(bahn_ids_to_process), ', '));
end

%% Hauptverarbeitungsschleife
tic;

% Verarbeite jede Bahn entsprechend dem gewählten Modus
for bahn_id = bahn_ids_to_process'
    current_bahn_id = num2str(bahn_id);
    
    % Bei Upload-Modus: Prüfe ob Bahn bereits existiert
    if (upload_all || upload_single) && ismember(bahn_id, existing_bahn_ids)
        disp(['Bahn-ID ' current_bahn_id ' bereits vorhanden - wird übersprungen']);
        continue;
    end

    try
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
        [trafo_rot, trafo_trans, ~] = calibration(data_cal_ist,data_cal_soll, plots);
        
        % Extrahieren der Kalibrierungs-Daten für die Orientierung
        tablename_cal = ['robotervermessung.' schema '.bahn_orientation_soll'];
        opts_cal = databaseImportOptions(conn,tablename_cal);
        opts_cal.RowFilter = opts_cal.RowFilter.bahn_id == calibration_id;
        data_cal_soll = sqlread(conn,tablename_cal,opts_cal);
        data_cal_soll = sortrows(data_cal_soll,'timestamp');

        % Berechnung der relativen Rotationsmatrix für die Orientierung 
        calibrateQuaternion(data_cal_ist, data_cal_soll);
        
        clear data_cal opts_cal tablename_cal
        
        % Auslesen der gesamten Ist-Daten der zu transformierenden Bahn
        query = ['SELECT * FROM robotervermessung.' schema '.bahn_pose_ist ' ...
                'WHERE robotervermessung.' schema '.bahn_pose_ist.bahn_id = ''' current_bahn_id ''''];
        data_ist = fetch(conn, query);
        data_ist = sortrows(data_ist,'timestamp');
        
        % Auslesen der gesamten Soll-Daten der Orientierung der zu transformierenden Bahn
        query = ['SELECT * FROM robotervermessung.' schema '.bahn_orientation_soll ' ...
                'WHERE robotervermessung.' schema '.bahn_orientation_soll.bahn_id = ''' current_bahn_id ''''];
        data_orientation_soll = fetch(conn, query);
        data_orientation_soll = sortrows(data_orientation_soll,'timestamp');
        
        % Auslesen der gesamten Soll-Daten der Position der zu transformierenden Bahn
        query = ['SELECT * FROM robotervermessung.' schema '.bahn_position_soll ' ...
                'WHERE robotervermessung.' schema '.bahn_position_soll.bahn_id = ''' current_bahn_id ''''];
        data_position_soll = fetch(conn, query);
        data_position_soll = sortrows(data_position_soll,'timestamp');
        
        position_ist = table2array(data_ist(:,5:7));
        position_soll = table2array(data_position_soll(:,5:7));
                
        data_ist_trafo = coord_transformation(position_ist,trafo_rot,trafo_trans);
        
        % Quaternion Transformation
        q_transformed = transformQuaternion(data_ist, data_orientation_soll, q_transform, trafo_rot);

        % Zur Visualisierung:
        % Converter quaternions para ângulos de Euler (em radianos)
        euler_trans = quat2eul(q_transformed);
        q_soll = [data_orientation_soll.qw_soll, data_orientation_soll.qx_soll, data_orientation_soll.qy_soll, data_orientation_soll.qz_soll];
        euler_soll = quat2eul(q_soll);
        
        % Converter para graus
        euler_trans_deg = rad2deg(euler_trans);
        euler_soll_deg = rad2deg(euler_soll);
        
        % Gimbal-Lock
        euler_trans_deg_fixed = fixGimbalLock(euler_trans_deg);
        euler_soll_deg_fixed = fixGimbalLock(euler_soll_deg);
             
        % Wenn Upload-Modus: Speichere in Datenbank
        if upload_all || upload_single
            uploadTransformedData(conn, current_bahn_id, calibration_id, data_ist, data_ist_trafo, q_transformed, schema);
        end

        % Wenn Plots aktiviert: Visualisiere Ergebnisse
        if plots
            plotResults(data_ist, data_ist_trafo, data_orientation_soll, position_soll, euler_soll_deg_fixed, euler_trans_deg_fixed)
        end

    catch ME
        warning(['Fehler bei Bahn-ID ' current_bahn_id ': ' ME.message]);
        continue;
    end
end

toc;

%% Hilfsfunktionen

function mode_str = getModeString(transform_only, upload_single, upload_all)
    if transform_only
        mode_str = 'Nur Transformation';
    elseif upload_single
        mode_str = 'Upload einzelne Bahn';
    elseif upload_all
        mode_str = 'Upload alle Bahnen';
    else
        mode_str = 'Kein Modus ausgewählt';
    end
end

function uploadTransformedData(conn, bahn_id, calibration_id, data_ist, data_ist_trafo, q_transformed, schema)
    try
        % Erstelle Tabelle für transformierte Daten
        bahn_pose_trans = table('Size', [height(data_ist), 11], ...
            'VariableTypes', {'string', 'string', 'string', 'double', 'double', ...
                            'double', 'double', 'double', 'double', 'double', 'string'}, ...
            'VariableNames', {'bahn_id', 'segment_id', 'timestamp', 'x_trans', ...
                            'y_trans', 'z_trans', 'qx_trans', 'qy_trans', ...
                            'qz_trans', 'qw_trans', 'calibration_id'});
        
        % bahn_pose_trans.bahn_id = data_ist.bahn_id;
        calibration_ids = repelem(string(calibration_id),height(bahn_pose_trans))';
        bahn_pose_trans{:,:} = [data_ist{:,2:4}, data_ist_trafo, q_transformed(:,2), q_transformed(:,3), q_transformed(:,4), q_transformed(:,1), calibration_ids];
        
        % Upload zur Datenbank
        sqlwrite(conn, ['robotervermessung.' schema '.bahn_pose_trans'], bahn_pose_trans);
        disp(['Bahn-ID ' bahn_id ' erfolgreich in Datenbank geschrieben']);
        
    catch ME
        error('Fehler beim Datenbank-Upload: %s', ME.message);
    end
end

%% Plots
function plotResults(data_ist, data_ist_trafo, data_orientation_soll, position_soll, euler_soll, euler_trans)
    % Colors
    c1 = [0 0.4470 0.7410];    % Blue - SOLL
    c2 = [0.8500 0.3250 0.0980]; % Orange - IST
    c3 = [0.9290 0.6940 0.1250]; % Yellow - Euler Trans

    % Timestamps in seconds
    time_ist = str2double(data_ist.timestamp);
    time_soll = str2double(data_orientation_soll.timestamp);
    timestamps_ist = (time_ist(:,1)- time_soll(1,1))/1e9;
    timestamps_soll = (time_soll(:,1)- time_soll(1,1))/1e9;

    % Get IST Euler angles
    q_ist = table2array(data_ist(:,8:11));
    q_ist = [q_ist(:,4), q_ist(:,3), q_ist(:,2), q_ist(:,1)];
    euler_ist = rad2deg(quat2eul(q_ist));

    % Plot angles
    figure('Color','white','Name','Euler Angles Comparison', 'Position', [100 100 1200 800])
    
    % Roll
    subplot(3,1,1)
    hold on
    plot(timestamps_soll, euler_soll(:,1), '--', 'Color', c1, 'LineWidth', 5, 'DisplayName', 'SOLL')
    plot(timestamps_ist, euler_ist(:,1), ':', 'Color', c2, 'LineWidth', 1, 'DisplayName', 'IST')
    plot(timestamps_ist, euler_trans(:,1), '-', 'Color', c3, 'LineWidth', 1.5, 'DisplayName', 'Euler Trans')
    title('Yaw Angle')
    ylabel('Angle [°]')
    legend('Location', 'best')
    grid on
    hold off

    % Pitch
    subplot(3,1,2)
    hold on
    plot(timestamps_soll, euler_soll(:,2), '--', 'Color', c1, 'LineWidth', 5, 'DisplayName', 'SOLL')
    plot(timestamps_ist, euler_ist(:,2), ':', 'Color', c2, 'LineWidth', 1, 'DisplayName', 'IST')
    plot(timestamps_ist, euler_trans(:,2), '-', 'Color', c3, 'LineWidth', 1.5, 'DisplayName', 'Euler Trans')
    title('Pitch Angle')
    ylabel('Angle [°]')
    legend('Location', 'best')
    grid on
    hold off

    % Yaw
    subplot(3,1,3)
    hold on
    plot(timestamps_soll, euler_soll(:,3), '--', 'Color', c1, 'LineWidth', 5, 'DisplayName', 'SOLL')
    plot(timestamps_ist, euler_ist(:,3), ':', 'Color', c2, 'LineWidth', 1, 'DisplayName', 'IST')
    plot(timestamps_ist, euler_trans(:,3), '-', 'Color', c3, 'LineWidth', 1.5, 'DisplayName', 'Euler Trans')
    title('Roll Angle')
    xlabel('Time [s]')
    ylabel('Angle [°]')
    legend('Location', 'best')
    grid on
    hold off

    % Plot Position
    figure('Color', 'white', 'Name', 'Position Comparison');
    hold on
    plot3(position_soll(:,1), position_soll(:,2), position_soll(:,3), '--', 'Color', c1, 'LineWidth', 1.5, 'DisplayName', 'SOLL')
    %plot3(data_ist(:,5), data_ist(:,6), data_ist(:,7), ':', 'Color', c2, 'LineWidth', 1.5, 'DisplayName', 'IST')
    plot3(data_ist_trafo(:,1), data_ist_trafo(:,2), data_ist_trafo(:,3), '-', 'Color', c3, 'LineWidth', 2, 'DisplayName', 'Transformed')
    grid on
    xlabel('X [mm]')
    ylabel('Y [mm]')
    zlabel('Z [mm]')
    legend('Location', 'best')
    view(3)
    hold off
end

function euler_fixed = fixGimbalLock(euler_angles)
    euler_fixed = euler_angles;
    
    for i = 1:3  % Check each angle component
        angle_data = euler_angles(:,i);
        
        % Check if we have values close to ±180
        near_180 = abs(abs(angle_data) - 180) < 5;
        
        if any(near_180)
            % If we have values near 180, fix sign flips
            mask_neg = angle_data < 0;
            angle_data(mask_neg) = angle_data(mask_neg) + 360;
            euler_fixed(:,i) = angle_data;
        end
    end
end


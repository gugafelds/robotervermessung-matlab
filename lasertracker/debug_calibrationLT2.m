%% Laser Tracker to Robot Coordinate Transformation
% Diese Funktion integriert den Laser Tracker Ansatz in den bestehenden Workflow
% Anstatt eine Kalibrierungsbahn zu erzeugen, werden direkt die Transformationsmatrizen 
% des Lasertrackers verwendet.

% Autor: 
% Datum: April 2025

%% Konfiguration und Datenbankverbindung
clear;
clc;

% Bahn-ID definieren
bahn_id_ = '1748618982';      % ID der zu verarbeitenden Bahn
plots = true;                % Plotten der Daten 
schema = 'bewegungsdaten';    % Datenbankschema

% Transformationsmatrizen vom Lasertracker
% Offset vom gemessenen Punkt zum TCP (Tool Center Point)
offset_matrix = [
    -0.50667304857499 0.86204924111502 -0.01239063116844 20.23474206984281; -0.01130064121523 0.00773011458295 0.99990626602530 -176.90021947940932; 0.86206421881195 0.50676557817354 0.00582506846212 -0.07879235287299; 0.00000000000000 0.00000000000000 0.00000000000000 1.00000000000000 
];

% Transformation vom Instrumenten- zum Roboterkoordinatensystem
robot_to_instrument = [
    0.19955206802495 -0.97967025979272 0.02061926828551 3933.26511788506969; 0.97975560685082 0.19982596445807 0.01218748426779 41.18685528571314; -0.01605998105042 0.01776980602223 0.99971311434961 1280.41540617110309; 0.00000000000000 0.00000000000000 0.00000000000000 1.00000000000000 
];

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

%% Daten aus der Datenbank laden
try
    % Extrahieren der IST-Daten
    query = ['SELECT * FROM robotervermessung.' schema '.bahn_pose_ist ' ...
            'WHERE robotervermessung.' schema '.bahn_pose_ist.bahn_id = ''' bahn_id_ ''''];
    data_ist = fetch(conn, query);
    data_ist = sortrows(data_ist,'timestamp');
    
    % Auslesen der Soll-Daten für Orientierung
    query = ['SELECT * FROM robotervermessung.' schema '.bahn_orientation_soll ' ...
            'WHERE robotervermessung.' schema '.bahn_orientation_soll.bahn_id = ''' bahn_id_ ''''];
    data_orientation_soll = fetch(conn, query);
    data_orientation_soll = sortrows(data_orientation_soll,'timestamp');
    
    % Auslesen der Soll-Daten für Position
    query = ['SELECT * FROM robotervermessung.' schema '.bahn_position_soll ' ...
            'WHERE robotervermessung.' schema '.bahn_position_soll.bahn_id = ''' bahn_id_ ''''];
    data_position_soll = fetch(conn, query);
    data_position_soll = sortrows(data_position_soll,'timestamp');
    
    % Daten extrahieren
    position_ist = table2array(data_ist(:,5:7));     % IST-Position [x, y, z]
    quaternion_ist = table2array(data_ist(:,8:11));  % IST-Quaternion [qx, qy, qz, qw]
    position_soll = table2array(data_position_soll(:,5:7)); % SOLL-Position
    
    disp(['Erfolgreich ' num2str(height(data_ist)) ' Datenpunkte geladen']);
    
catch ME
    error('Fehler beim Laden der Daten: %s', ME.message);
end

%% Transformation der Positionen und Orientierungen
disp('Beginne Transformation der Daten...');

% Initialisierung der Ergebnisarrays
num_points = size(position_ist, 1);
positions_transformed = zeros(num_points, 3);
quaternions_transformed = zeros(num_points, 4);
orientations_transformed = zeros(num_points, 3);

% Verarbeitung jedes Punktes
for i = 1:num_points
    % Umformatierung des Quaternions (von [qx, qy, qz, qw] zu [qw, qx, qy, qz])
    quaternion_sa = LT2SAQuaternionDirect(quaternion_ist(i,:));

    quaternion_reordered = [quaternion_sa(4), quaternion_sa(1), quaternion_sa(2), quaternion_sa(3)];
    
    % 1. Anwenden des Offsets auf den Punkt und die Orientierung
    [point_offset, ~, quaternion_offset] = applyOffset(position_ist(i,:), [], quaternion_reordered, offset_matrix);
    
    % 2. Transformation in Roboterkoordinaten
    positions_transformed(i,:) = transformToRobotCoords(point_offset, robot_to_instrument);
    quaternions_transformed(i,:) = transformQuaternionToRobotCoords(quaternion_offset, robot_to_instrument);
    
    % 3. Berechnung der Euler-Winkel für die Robotersteuerung (XYZ Fixed Angles)
    orientations_transformed(i,:) = quaternionToEulerXYZFixed(quaternions_transformed(i,:));
end

% Quaternionen zurück ins Format [qx, qy, qz, qw] bringen für die Datenbank
quaternions_db_format = [quaternions_transformed(:,2), quaternions_transformed(:,3), quaternions_transformed(:,4), quaternions_transformed(:,1)];

%% Visualisierung der Ergebnisse
if plots
    % Zeitstempel in Sekunden umrechnen
    time_ist = str2double(data_ist.timestamp);
    time_soll = str2double(data_orientation_soll.timestamp);
    timestamps_ist = (time_ist(:,1)- time_ist(1,1))/1e9;
    timestamps_soll = (time_soll(:,1)- time_soll(1,1))/1e9;
    
    % Quaternionen zu Euler-Winkeln konvertieren
    q_soll = [data_orientation_soll.qw_soll, data_orientation_soll.qx_soll, data_orientation_soll.qy_soll, data_orientation_soll.qz_soll];
    euler_ist = rad2deg(quat2eul([quaternion_ist(:,4), quaternion_ist(:,1), quaternion_ist(:,2), quaternion_ist(:,3)]));
    euler_trans = rad2deg(quat2eul(quaternions_transformed));
    euler_soll = rad2deg(quat2eul(q_soll));
    
    % Plot der Orientierungen (Euler-Winkel)
    figure('Color','white','Name','Euler Angles Comparison', 'Position', [100 100 1200 800])
    
    % Farben
    c1 = [0 0.4470 0.7410];     % Blau - SOLL
    c2 = [0.8500 0.3250 0.0980]; % Orange - IST
    c3 = [0.9290 0.6940 0.1250]; % Gelb - Transformiert
    
    % Yaw
    subplot(3,1,1)
    hold on
    plot(timestamps_soll, euler_soll(:,1), '--', 'Color', c1, 'LineWidth', 2, 'DisplayName', 'SOLL')
    plot(timestamps_ist, euler_ist(:,1), ':', 'Color', c2, 'LineWidth', 1, 'DisplayName', 'IST')
    plot(timestamps_ist, euler_trans(:,1), '-', 'Color', c3, 'LineWidth', 1.5, 'DisplayName', 'Transformiert')
    title('Yaw Angle')
    ylabel('Angle [°]')
    legend('Location', 'best')
    grid on
    hold off
    
    % Pitch
    subplot(3,1,2)
    hold on
    plot(timestamps_soll, euler_soll(:,2), '--', 'Color', c1, 'LineWidth', 2, 'DisplayName', 'SOLL')
    plot(timestamps_ist, euler_ist(:,2), ':', 'Color', c2, 'LineWidth', 1, 'DisplayName', 'IST')
    plot(timestamps_ist, euler_trans(:,2), '-', 'Color', c3, 'LineWidth', 1.5, 'DisplayName', 'Transformiert')
    title('Pitch Angle')
    ylabel('Angle [°]')
    legend('Location', 'best')
    grid on
    hold off
    
    % Roll
    subplot(3,1,3)
    hold on
    plot(timestamps_soll, euler_soll(:,3), '--', 'Color', c1, 'LineWidth', 2, 'DisplayName', 'SOLL')
    plot(timestamps_ist, euler_ist(:,3), ':', 'Color', c2, 'LineWidth', 1, 'DisplayName', 'IST')
    plot(timestamps_ist, euler_trans(:,3), '-', 'Color', c3, 'LineWidth', 1.5, 'DisplayName', 'Transformiert')
    title('Roll Angle')
    xlabel('Time [s]')
    ylabel('Angle [°]')
    legend('Location', 'best')
    grid on
    hold off
    
    % Plot der Positionen
    figure('Color', 'white', 'Name', 'Position Comparison', 'Position', [100 100 800 600]);
    hold on
    plot3(position_soll(:,1), position_soll(:,2), position_soll(:,3), '--', 'Color', c1, 'LineWidth', 2, 'DisplayName', 'SOLL')
    plot3(position_ist(:,1), position_ist(:,2), position_ist(:,3), ':', 'Color', c2, 'LineWidth', 1, 'DisplayName', 'IST')
    plot3(positions_transformed(:,1), positions_transformed(:,2), positions_transformed(:,3), '-', 'Color', c3, 'LineWidth', 1.5, 'DisplayName', 'Transformiert')
    grid on
    xlabel('X [mm]')
    ylabel('Y [mm]')
    zlabel('Z [mm]')
    legend('Location', 'best')
    view(3)
    axis equal
    hold off
end

%% Speichern der transformierten Daten in der Datenbank
% Frage den Benutzer, ob die Daten gespeichert werden sollen
save_data = input('Möchten Sie die transformierten Daten in der Datenbank speichern? (j/n): ', 's');

if strcmpi(save_data, 'j')
    try
        % Erstelle leere Tabelle
        bahn_pose_trans = table('Size', [height(data_ist), 11], ...
            'VariableTypes', {'string', 'string', 'string', 'double', 'double', ...
                            'double', 'double', 'double', 'double', 'double', 'string'}, ...
            'VariableNames', {'bahn_id', 'segment_id', 'timestamp', 'x_trans', ...
                            'y_trans', 'z_trans', 'qx_trans', 'qy_trans', ...
                            'qz_trans', 'qw_trans', 'calibration_id'});
        
        % Füllen der Spalten einzeln
        bahn_pose_trans.bahn_id = string(data_ist.bahn_id);
        bahn_pose_trans.segment_id = string(data_ist.segment_id);
        bahn_pose_trans.timestamp = string(data_ist.timestamp);
        
        % Positions- und Orientierungsdaten
        bahn_pose_trans.x_trans = positions_transformed(:,1);
        bahn_pose_trans.y_trans = positions_transformed(:,2);
        bahn_pose_trans.z_trans = positions_transformed(:,3);
        
        % Quaternionen
        bahn_pose_trans.qx_trans = quaternions_db_format(:,1);
        bahn_pose_trans.qy_trans = quaternions_db_format(:,2);
        bahn_pose_trans.qz_trans = quaternions_db_format(:,3);
        bahn_pose_trans.qw_trans = quaternions_db_format(:,4);
        
        % Kalibrier-ID
        bahn_pose_trans.calibration_id = repelem("spatial_analyzer", height(bahn_pose_trans))';
        
        % Überprüfe, ob bereits Daten für diese Bahn existieren
        query = ['SELECT COUNT(*) FROM robotervermessung.' schema '.bahn_pose_trans ' ...
                'WHERE robotervermessung.' schema '.bahn_pose_trans.bahn_id = ''' bahn_id_ ''''];
        count = fetch(conn, query);
        
        if count.count > 0
            % Lösche vorhandene Daten
            delete_query = ['DELETE FROM robotervermessung.' schema '.bahn_pose_trans ' ...
                           'WHERE robotervermessung.' schema '.bahn_pose_trans.bahn_id = ''' bahn_id_ ''''];
            exec(conn, delete_query);
            disp(['Vorhandene Daten für Bahn-ID ' bahn_id_ ' wurden gelöscht']);
        end
        
        % Upload zur Datenbank
        sqlwrite(conn, ['robotervermessung.' schema '.bahn_pose_trans'], bahn_pose_trans);
        disp(['Bahn-ID ' bahn_id_ ' erfolgreich in Datenbank geschrieben']);
        
    catch ME
        error('Fehler beim Datenbank-Upload: %s', ME.message);
    end
end

% Verbindung schließen
close(conn);
disp('Datenbankverbindung geschlossen');


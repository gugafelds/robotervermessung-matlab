%% Einstellungen

clear;
tic;

% bahn_id_ = '172104917';
bahn_id_ = '171992752'; % Bahn mit wenigen Punkten
% bahn_id_ = '171991250';
% bahn_id_ = '172104925'; % ---> mit Orientierungsänderung 

% bahn_id_ = '172079653'; % ---> mit Orientierungsänderung
% bahn_id_ ='172054053'; % Orientierungsänderung ohne Kalibrierungsdatei
% bahn_id_ = '172079403'; % Kalibrierungsdatei selbst

%%% Standard: --> Berechnung der Metriken für Positionsabweichungen %%%

% Berechnung der Metriken für die Geschwindikeitsabweichungen
evaluate_velocity = false;

% Berechnung der Metriken für die Orientierungsabweichungen
evaluate_orientation = false;

% Berechnung der Metriken für bestimmte Bahnabschnitte
evaluate_segmentwise = false;
if evaluate_segmentwise
    segment_first = 1;
    segment_last = 7;
end
% Berechnung der Metriken für die gesamte Messaufnahme
evaluate_all = false;

% Plotten der Daten 
plots = false;

% Upload in die Datenbank
upload = false;

% % Verbindung mit PostgreSQL
% datasource = "RobotervermessungMATLAB";
% username = "postgres";
% password = "200195Beto";
% conn = postgresql(datasource,username,password);
% Verbindung mit PostgreSQL
datasource = "RobotervermessungMATLAB";
username = "felixthomas";
password = "manager";
conn = postgresql(datasource,username,password);

% Überprüfe Verbindung
if isopen(conn)
    disp('Verbindung erfolgreich hergestellt');
else
    disp('Verbindung fehlgeschlagen');
end

clear datasource username password

%% Suche nach zugehörigem "Calibration Run"

query_cal = 'SELECT * FROM robotervermessung.bewegungsdaten.bahn_info WHERE robotervermessung.bewegungsdaten.bahn_info.calibration_run = true';

% Abfrage ausführen und Ergebnisse abrufen
data_cal_info = fetch(conn, query_cal);

% Finden des zugehörigen Calibration Runs anhand der kürzesten vergangen Zeit
check_bahn_id = str2double(data_cal_info.bahn_id);
diff_bahn_id = check_bahn_id - str2double(bahn_id_);

[~,min_diff_idx] = min(abs(diff_bahn_id));

% Wenn eine Kalibierungsdatei vorliegt wird diese für die
% Koordinatentransformation genutzt, ansonsten die wird die gewählte Datei
% selbst verwendet. 
if diff_bahn_id(min_diff_idx) < 0
    calibration_id = char(data_cal_info{min_diff_idx,'bahn_id'});
    disp('Kalibrierungs-Datei vorhanden! ID der Messaufnahme: ' + string(calibration_id))
else
    calibration_id = bahn_id_;
    disp('Zu dem ausgewählten Datensatz liegt keine Kalibirierungsdatei vor!')
end

% Extrahieren der Kalibrierungs-Daten
tablename_cal = 'robotervermessung.bewegungsdaten.bahn_pose_ist';
opts_cal = databaseImportOptions(conn,tablename_cal);
opts_cal.RowFilter = opts_cal.RowFilter.bahn_id == calibration_id;
data_cal_ist= sqlread(conn,tablename_cal,opts_cal);
data_cal_ist = sortrows(data_cal_ist,'timestamp');


tablename_cal = 'robotervermessung.bewegungsdaten.bahn_events';
opts_cal = databaseImportOptions(conn,tablename_cal);
opts_cal.RowFilter = opts_cal.RowFilter.bahn_id == calibration_id;
data_cal_soll = sqlread(conn,tablename_cal,opts_cal);
data_cal_soll = sortrows(data_cal_soll,'timestamp');

% Positionsdaten für Koordinatentransformation
calibration(data_cal_ist,data_cal_soll)

if evaluate_orientation == true
    % Wenn Orientierung wird andere Collection benötigt
    tablename_cal = 'robotervermessung.bewegungsdaten.bahn_orientation_soll';
    opts_cal = databaseImportOptions(conn,tablename_cal);
    opts_cal.RowFilter = opts_cal.RowFilter.bahn_id == calibration_id;
    data_cal_soll = sqlread(conn,tablename_cal,opts_cal);
    data_cal_soll = sortrows(data_cal_soll,'timestamp');


    
    % Transformation der Quaternionen/Eulerwinkel
    euler_transformation(data_cal_ist,data_cal_soll)
end

% if evaluate_orientation == true
%     % Wenn Orientierung wird andere Collection benötigt
%     tablename_cal = 'robotervermessung.bewegungsdaten.bahn_orientation_soll';
%     opts_cal = databaseImportOptions(conn,tablename_cal);
%     opts_cal.RowFilter = opts_cal.RowFilter.bahn_id == calibration_id;
%     data_cal_soll = sqlread(conn,tablename_cal,opts_cal);
%     data_cal_soll = sortrows(data_cal_soll,'timestamp');
% 
%     % Transformation der Quaternionen/Eulerwinkel
%     euler_transformation(data_cal_ist,data_cal_soll)
% else
%     tablename_cal = 'robotervermessung.bewegungsdaten.bahn_events';
%     opts_cal = databaseImportOptions(conn,tablename_cal);
%     opts_cal.RowFilter = opts_cal.RowFilter.bahn_id == calibration_id;
%     data_cal_soll = sqlread(conn,tablename_cal,opts_cal);
%     data_cal_soll = sortrows(data_cal_soll,'timestamp');
% 
%     % Positionsdaten für Koordinatentransformation
%     calibration(data_cal_ist,data_cal_soll)
% end

clear data_cal data_cal_info diff_bahn_id min_diff_bahn_id min_idx opts_cal tablename_cal check_bahn_id
clear query_cal min_diff_idx
clear data_cal_ist data_cal_soll


%% Auslesen der für die entsprechende Auswertung benötigten Daten

% Anzahl der Segmente der gesamten Messaufnahme bestimmen 
query = ['SELECT * FROM robotervermessung.bewegungsdaten.bahn_info ' ...
         'WHERE robotervermessung.bewegungsdaten.bahn_info.bahn_id = ''' bahn_id_ ''''];
data_info = fetch(conn, query);
num_segments = data_info.np_ereignisse;

% Daten auslesen
if evaluate_velocity == false && evaluate_orientation == true

    % Auslesen der gesamten Ist-Daten
    query = ['SELECT * FROM robotervermessung.bewegungsdaten.bahn_pose_ist ' ...
            'WHERE robotervermessung.bewegungsdaten.bahn_pose_ist.bahn_id = ''' bahn_id_ ''''];
    data_ist = fetch(conn, query);
    data_ist = sortrows(data_ist,'timestamp');
    
    % Auslesen der gesamten Soll-Daten
    query = ['SELECT * FROM robotervermessung.bewegungsdaten.bahn_orientation_soll ' ...
            'WHERE robotervermessung.bewegungsdaten.bahn_orientation_soll.bahn_id = ''' bahn_id_ ''''];
    data_soll = fetch(conn, query);
    data_soll = sortrows(data_soll,'timestamp');

    % Transformation der Quarternionen zu Euler-Winkeln
    q_soll = table2array(data_soll(:,5:8));
    q_soll = [q_soll(:,4), q_soll(:,3), q_soll(:,2), q_soll(:,1)];
    euler_soll = quat2eul(q_soll,"ZYX");
    euler_soll = rad2deg(euler_soll);
    
    q_ist = table2array(data_ist(:,8:11));
    q_ist = [q_ist(:,4), q_ist(:,3), q_ist(:,2), q_ist(:,1)];
    euler_ist = quat2eul(q_ist,"ZYX");
    euler_ist = rad2deg(euler_ist);

elseif evaluate_velocity == true && evaluate_orientation == false 

    % Auslesen der gesamten Ist-Daten
    query = ['SELECT * FROM robotervermessung.bewegungsdaten.bahn_twist_ist ' ...
            'WHERE robotervermessung.bewegungsdaten.bahn_twist_ist.bahn_id = ''' bahn_id_ ''''];
    data_ist = fetch(conn, query);
    data_ist = sortrows(data_ist,'timestamp');
    
    % Auslesen der gesamten Soll-Daten

    % query = ['SELECT * FROM robotervermessung.bewegungsdaten.bahn_twist_soll ' ...
    %         'WHERE robotervermessung.bewegungsdaten.bahn_twist_soll.bahn_id = ''' bahn_id_ ''''];
    query = ['SELECT * FROM robotervermessung.bewegungsdaten.bahn_position_soll ' ...
            'WHERE robotervermessung.bewegungsdaten.bahn_position_soll.bahn_id = ''' bahn_id_ ''''];
    data_soll = fetch(conn, query);
    data_soll = sortrows(data_soll,'timestamp');

%%%%%%%%%%%%%%% Interpolation der Sollgeschwindigkeit über Positionen und Timestamp
%%%%%%%%%%%%%%% (funktioniert noch nicht!)
    % vel_soll_time = data_soll(:,4);
    
    vel_soll = data_soll(:,4:7);
    vel_soll.timestamp = str2double(vel_soll.timestamp);
    % vel_soll = table2array(vel_soll);
    
    % Timestamps in Sekunden
    vel_soll.timestamp = (vel_soll.timestamp - vel_soll.timestamp(1,1))/1e9;


    v = [diff(vel_soll.x_soll)./diff(vel_soll.timestamp), diff(vel_soll.y_soll)./diff(vel_soll.timestamp), diff(vel_soll.z_soll)./diff(vel_soll.timestamp)];
    v = sqrt(v(:,1).^2 + v(:,2).^2 + v(:,3).^2);
    v2 = [diff(vel_soll.x_soll) diff(vel_soll.y_soll) diff(vel_soll.z_soll)];
    v2 = sqrt(v2(:,1).^2 + v2(:,2).^2 + v2(:,3).^2);

    idx_zero = find(v2 == 0);

    v2(idx_zero) = [];
    t = diff(vel_soll.timestamp);
    t(idx_zero) = [];

    v3 = [diff(vel_soll.x_soll) , diff(vel_soll.y_soll), diff(vel_soll.z_soll)];
    v3(idx_zero,:) = [];

    v4 = [v3(:,1)./t, v3(:,2)./t, v3(:,3)./t];
    v5 = sqrt(v4(:,1).^2 + v4(:,2).^2 + v4(:,3).^2);

% window = 37; 
% 
% ma_filter = (1/window)* ones (1,window);
% abbx_filtered = filter(ma_filter,1,abb_positions2(:,1));
% abby_filtered = filter(ma_filter,1,abb_positions2(:,2));
% abbz_filtered = filter(ma_filter,1,abb_positions2(:,3));
% 
% abb_filtered = [abbx_filtered abby_filtered abbz_filtered];

%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%

else
    % Auslesen der gesamten Ist-Daten
    query = ['SELECT * FROM robotervermessung.bewegungsdaten.bahn_pose_ist ' ...
            'WHERE robotervermessung.bewegungsdaten.bahn_pose_ist.bahn_id = ''' bahn_id_ ''''];
    data_ist = fetch(conn, query);
    data_ist = sortrows(data_ist,'timestamp');
    
    % Auslesen der gesamten Soll-Daten
    query = ['SELECT * FROM robotervermessung.bewegungsdaten.bahn_position_soll ' ...
            'WHERE robotervermessung.bewegungsdaten.bahn_position_soll.bahn_id = ''' bahn_id_ ''''];
    data_soll = fetch(conn, query);
    data_soll = sortrows(data_soll,'timestamp');

end

clear q_ist q_soll


%% Extraktion und Separation der Segmente der Gesamtaufname

% Alle Segment-ID's 
query = ['SELECT segment_id FROM robotervermessung.bewegungsdaten.bahn_events ' ...
    'WHERE robotervermessung.bewegungsdaten.bahn_events.bahn_id = ''' bahn_id_ ''''];

segment_ids = fetch(conn,query);

% % % IST-DATEN % % %
% Extraktion der Indizes der Segmente 
seg_id = split(data_ist.segment_id, '_');
seg_id = str2double(seg_id(:,2));
idx_new_seg_ist = zeros(num_segments,1);

% Suche nach den Indizes bei denen sich die Segmentnr. ändert
k = 0;
idx = 1;
for i = 1:1:length(seg_id)
    if seg_id(i) == k
        idx = idx + 1;
    else
        k = k +1;
        idx_new_seg_ist(k) = idx;
        idx = idx+1;
    end
end

% % % SOLL-DATEN % % %
seg_id = split(data_soll.segment_id, '_');
seg_id = str2double(seg_id(:,2));
idx_new_seg_soll = zeros(num_segments,1);

k = 0;
idx = 1;
for i = 1:1:length(seg_id)
    if seg_id(i) == k
        idx = idx + 1;
    else
        k = k +1;
        idx_new_seg_soll(k) = idx;
        idx = idx+1;
    end
end


if evaluate_velocity == true && evaluate_orientation == false 

    disp('Es wird die Geschwindigkeit ausgewertet!')

    % Speichern der einzelnen Semgente in Tabelle
    segments_ist = array2table([{data_ist.segment_id(1)} data_ist.tcp_speed_ist(1:idx_new_seg_ist(1)-1)], "VariableNames",{'segment_id','tcp_speed_ist'});
    
    for i = 1:num_segments
    
        if i == length(idx_new_seg_ist)
            segments_ist(i+1,:) = array2table([{segment_ids{i,:}} data_ist.tcp_speed_ist(idx_new_seg_ist(i):end)]);
        else
            segments_ist(i+1,:) = array2table([{segment_ids{i,:}} data_ist.tcp_speed_ist(idx_new_seg_ist(i):idx_new_seg_ist(i+1)-1)]);
        end
    
    end
    
    segments_soll = array2table([{data_soll.segment_id(1)} data_soll.tcp_speed_soll(1:idx_new_seg_soll(1)-1)], "VariableNames",{'segment_id','tcp_speed_soll'});
    for i = 1:num_segments
        if i == length(idx_new_seg_soll)
            segments_soll(i+1,:) = array2table([{segment_ids{i,:}} data_soll.tcp_speed_soll(idx_new_seg_soll(i):end)]);
        else
            segments_soll(i+1,:)= array2table([{segment_ids{i,:}} data_soll.tcp_speed_soll(idx_new_seg_soll(i):idx_new_seg_soll(i+1)-1)]);
        end    
    end
    
elseif evaluate_velocity == false && evaluate_orientation == true

    disp('Es wird die Orientierung ausgewertet!')

    % euler_soll = mod(euler_soll,360);
    % euler_ist = mod(euler_ist,360);

    % Speichern der einzelnen Semgente in Tabelle
    segments_ist = array2table([{data_ist.segment_id(1)} euler_ist(1:idx_new_seg_ist(1)-1,1) euler_ist(1:idx_new_seg_ist(1)-1,2) euler_ist(1:idx_new_seg_ist(1)-1,3)], "VariableNames",{'segment_id','roll_ist','pitch_ist','yaw_ist'});
    
    for i = 1:num_segments
    
        if i == length(idx_new_seg_ist)
            segments_ist(i+1,:) = array2table([{segment_ids{i,:}} euler_ist(idx_new_seg_ist(i):end,1) euler_ist(idx_new_seg_ist(i):end,2) euler_ist(idx_new_seg_ist(i):end,3)]);
        else
            segments_ist(i+1,:) = array2table([{segment_ids{i,:}} euler_ist(idx_new_seg_ist(i):idx_new_seg_ist(i+1)-1,1) euler_ist(idx_new_seg_ist(i):idx_new_seg_ist(i+1)-1,2) euler_ist(idx_new_seg_ist(i):idx_new_seg_ist(i+1)-1,3)]);
        end
    
    end
    
    segments_soll = array2table([{data_soll.segment_id(1)} euler_soll(1:idx_new_seg_soll(1)-1,1) euler_soll(1:idx_new_seg_soll(1)-1,2) euler_soll(1:idx_new_seg_soll(1)-1,3)], "VariableNames",{'segment_id','roll_soll','pitch_soll','yaw_soll'});
    for i = 1:num_segments
        if i == length(idx_new_seg_soll)
            segments_soll(i+1,:) = array2table([{segment_ids{i,:}} euler_soll(idx_new_seg_soll(i):end,1) euler_soll(idx_new_seg_soll(i):end,2) euler_soll(idx_new_seg_soll(i):end,3)]);
        else
            segments_soll(i+1,:)= array2table([{segment_ids{i,:}} euler_soll(idx_new_seg_soll(i):idx_new_seg_soll(i+1)-1,1) euler_soll(idx_new_seg_soll(i):idx_new_seg_soll(i+1)-1,2) euler_soll(idx_new_seg_soll(i):idx_new_seg_soll(i+1)-1,3)]);
        end    
    end
    
    % Transformation der Eulerwinkel für alle Segemente
    segments_trafo = table();
    for i = 1:1:num_segments+1
        euler_transformation(segments_ist(i,:),segments_soll(i,:), trafo_euler, trafo_rot)
        % % Winkel von 0-360° statt -180° bis 180°
        segments_trafo(i,:) = seg_trafo;
    end

%%%%%%% Sonst automatisch Auswertung von Positionsdaten 
else

    disp('Es wird die Position ausgewertet!')
%%% Timestamp Auswertung
    timestamps_ist = str2double(data_ist.timestamp);
    timestamps_soll = str2double(data_soll.timestamp);

    rel_timestamps_soll = (timestamps_soll - timestamps_ist(1))/1e9;
    rel_timestamps_ist = (timestamps_ist - timestamps_ist(1))/1e9;

    % Speichern der einzelnen Semgente in Tabelle
    segments_ist = array2table([{data_ist.segment_id(1)} data_ist.x_ist(1:idx_new_seg_ist(1)-1) data_ist.y_ist(1:idx_new_seg_ist(1)-1) data_ist.z_ist(1:idx_new_seg_ist(1)-1) rel_timestamps_ist(1:idx_new_seg_ist(1)-1)], "VariableNames",{'segment_id','x_ist','y_ist','z_ist','timestamp_ist'});
    
    for i = 1:num_segments
    
        if i == length(idx_new_seg_ist)
            segments_ist(i+1,:) = array2table([{segment_ids{i,:}} data_ist.x_ist(idx_new_seg_ist(i):end) data_ist.y_ist(idx_new_seg_ist(i):end) data_ist.z_ist(idx_new_seg_ist(i):end) rel_timestamps_ist(idx_new_seg_ist(i):end)]);
        else
            segments_ist(i+1,:) = array2table([{segment_ids{i,:}} data_ist.x_ist(idx_new_seg_ist(i):idx_new_seg_ist(i+1)-1) data_ist.y_ist(idx_new_seg_ist(i):idx_new_seg_ist(i+1)-1) data_ist.z_ist(idx_new_seg_ist(i):idx_new_seg_ist(i+1)-1) rel_timestamps_ist(idx_new_seg_ist(i):idx_new_seg_ist(i+1)-1)]);
        end
    
    end
    
    segments_soll = array2table([{data_soll.segment_id(1)} data_soll.x_soll(1:idx_new_seg_soll(1)-1) data_soll.y_soll(1:idx_new_seg_soll(1)-1) data_soll.z_soll(1:idx_new_seg_soll(1)-1) rel_timestamps_soll(1:idx_new_seg_soll(1)-1)], "VariableNames",{'segment_id','x_soll','y_soll','z_soll','timestamp_soll'});
    for i = 1:num_segments
        if i == length(idx_new_seg_soll)
            segments_soll(i+1,:) = array2table([{segment_ids{i,:}} data_soll.x_soll(idx_new_seg_soll(i):end) data_soll.y_soll(idx_new_seg_soll(i):end) data_soll.z_soll(idx_new_seg_soll(i):end) rel_timestamps_soll(idx_new_seg_soll(i):end)]);
        else
            segments_soll(i+1,:)= array2table([{segment_ids{i,:}} data_soll.x_soll(idx_new_seg_soll(i):idx_new_seg_soll(i+1)-1) data_soll.y_soll(idx_new_seg_soll(i):idx_new_seg_soll(i+1)-1) data_soll.z_soll(idx_new_seg_soll(i):idx_new_seg_soll(i+1)-1) rel_timestamps_soll(idx_new_seg_soll(i):idx_new_seg_soll(i+1)-1)]);
        end    
    end
    
    % Koordinatentransformation für alle Segemente
    segments_trafo = table();
    for i = 1:1:num_segments+1
        coord_transformation(segments_ist(i,:),trafo_rot, trafo_trans)
        segments_trafo(i,:) = pos_ist_trafo;
    end

end

clear idx k seg_id query seg_trafo

%% Berechnung der Metriken

% Tabellen initialisieren
table_sidtw_info = table();
table_sidtw_deviation = cell(num_segments+1,1);
table_dtw_info = table();
table_dtw_deviation = cell(num_segments+1,1);
table_dfd_info = table();
table_dfd_deviation = cell(num_segments+1,1);

% Wenn Geschwindigkeit ausgewertet werden soll sind die 1D-Metriken nicht durchfürbar
if size(segments_soll,2) > 2
    table_euclidean_info = table();
    table_euclidean_deviation = cell(num_segments+1,1);
    table_lcss_info = table();
    table_lcss_deviation = cell(num_segments+1,1);
end

% Berechnung der Metriken für alle Segmente
for i = 1:1:num_segments

    if evaluate_velocity == false && evaluate_orientation == false 
        segment_trafo = [segments_trafo.x_ist{i}, segments_trafo.y_ist{i}, segments_trafo.z_ist{i}, segments_trafo.timestamp_ist{i}];
        segment_soll = [segments_soll.x_soll{i}, segments_soll.y_soll{i}, segments_soll.z_soll{i}, segments_soll.timestamp_soll{i}];
    elseif evaluate_velocity == true && evaluate_orientation == false 
        segment_trafo = segments_ist.tcp_speed_ist{i};
        segment_soll = segments_soll.tcp_speed_soll{i};
    elseif evaluate_velocity == false && evaluate_orientation == true 
        segment_trafo = [segments_trafo.roll_ist{i}, segments_trafo.pitch_ist{i},  segments_trafo.yaw_ist{i}];
        segment_soll = [segments_soll.roll_soll{i}, segments_soll.pitch_soll{i}, segments_soll.yaw_soll{i}];
        % Winkel von 0-360° statt -180° bis 180°
        segment_trafo = mod(segment_trafo,360);
        segment_soll = mod(segment_soll,360);
    end

    if size(segment_soll,2) > 1 % Wird nicht betrachtet wenn Geschwindigkeit ausgewertet wird 

    % Berechnung euklidischer Abstand
    [euclidean_ist,euclidean_distances,~] = distance2curve(segment_trafo,segment_soll,'linear');
    % Berechnung LCSS
    % [~, ~, lcss_distances, ~, ~, lcss_soll, lcss_ist, ~, ~] = fkt_lcss(segment_soll,segment_trafo,false);

    end
    % Berechnung SIDTW
    [sidtw_distances, ~, ~, ~, sidtw_soll, sidtw_ist, ~, ~, ~] = fkt_selintdtw3d(segment_soll,segment_trafo,false);
    % Berechnung DTW
    [dtw_distances, ~, ~, ~, dtw_soll, dtw_ist, ~, ~, ~, ~] = fkt_dtw3d(segment_soll,segment_trafo,false);
    % Berechnung diskrete Frechet
    fkt_discreteFrechet(segment_soll,segment_trafo,false)


    if i == 1
        if size(segment_soll,2) > 1
            % Euklidischer Abstand
            metric2postgresql('euclidean',euclidean_distances, segment_soll, euclidean_ist, bahn_id_, data_ist.segment_id(1))
            table_euclidean_info = seg_euclidean_info;
            order_eucl_first = size(seg_euclidean_distances,1);
            seg_euclidean_distances = [seg_euclidean_distances, table((1:1:order_eucl_first)','VariableNames',{'points_order'})];
            table_euclidean_deviation{1} = seg_euclidean_distances;
            % % LCSS
            % metric2postgresql('lcss',lcss_distances, lcss_soll, lcss_ist, bahn_id_, data_ist.segment_id(1))
            % table_lcss_info = seg_lcss_info;
            % order_lcss_first = size(seg_lcss_distances,1);
            % seg_lcss_distances = [seg_lcss_distances, table((1:1:order_lcss_first)','VariableNames',{'points_order'})];
            % table_lcss_deviation{1} = seg_lcss_distances;
        end
        % SIDTW
        metric2postgresql('sidtw', sidtw_distances, sidtw_soll, sidtw_ist, bahn_id_, data_ist.segment_id(1))
        table_sidtw_info = seg_sidtw_info;
        order_sidtw_first = size(seg_sidtw_distances,1);
        seg_sidtw_distances = [seg_sidtw_distances, table((1:1:order_sidtw_first)','VariableNames',{'points_order'})];
        table_sidtw_deviation{1} = seg_sidtw_distances;

        % DTW
        metric2postgresql('dtw',dtw_distances, dtw_soll, dtw_ist, bahn_id_, data_ist.segment_id(1))
        table_dtw_info = seg_dtw_info;
        order_dtw_first = size(seg_dtw_distances,1);
        seg_dtw_distances = [seg_dtw_distances, table((1:1:order_dtw_first)','VariableNames',{'points_order'})];
        table_dtw_deviation{1} = seg_dtw_distances;
        % DFD
        metric2postgresql('dfd',frechet_distances, frechet_soll, frechet_ist, bahn_id_, data_ist.segment_id(1))
        table_dfd_info = seg_dfd_info;
        order_dfd_first = size(seg_dfd_distances,1);
        seg_dfd_distances = [seg_dfd_distances, table((1:1:order_dfd_first)','VariableNames',{'points_order'})];
        table_dfd_deviation{1} = seg_dfd_distances;       


    else
        if size(segment_soll,2) > 1
            % Euklidischer Abstand
            metric2postgresql('euclidean',euclidean_distances, segment_soll, euclidean_ist, bahn_id_, segment_ids{i-1,:})
            table_euclidean_info(i,:) = seg_euclidean_info;
            order_eucl_last = order_eucl_first + size(seg_euclidean_distances,1);
            seg_euclidean_distances = [seg_euclidean_distances, table((order_eucl_first+1:1:order_eucl_last)','VariableNames',{'points_order'})];
            order_eucl_first = order_eucl_last;
            table_euclidean_deviation{i} = seg_euclidean_distances;
            % LCSS
            % metric2postgresql('lcss',lcss_distances, lcss_soll, lcss_ist, bahn_id_, segment_ids{i-1,:})
            % table_lcss_info(i,:) = seg_lcss_info;
            % order_lcss_last = order_lcss_first + size(seg_lcss_distances,1);
            % seg_lcss_distances = [seg_lcss_distances, table((order_lcss_first+1:1:order_lcss_last)','VariableNames',{'points_order'})];
            % order_lcss_first = order_lcss_last;
            % table_lcss_deviation{i} = seg_lcss_distances;
        end
        % SIDTW
        metric2postgresql('sidtw',sidtw_distances, sidtw_soll, sidtw_ist, bahn_id_, segment_ids{i-1,:})
        table_sidtw_info(i,:) = seg_sidtw_info;
        order_sidtw_last = order_sidtw_first + size(seg_sidtw_distances,1);
        seg_sidtw_distances = [seg_sidtw_distances, table((order_sidtw_first+1:1:order_sidtw_last)','VariableNames',{'points_order'})];
        order_sidtw_first = order_sidtw_last;
        table_sidtw_deviation{i} = seg_sidtw_distances;
        % DTW
        metric2postgresql('dtw',dtw_distances, dtw_soll, dtw_ist, bahn_id_, segment_ids{i-1,:})
        table_dtw_info(i,:) = seg_dtw_info;
        order_dtw_last = order_dtw_first + size(seg_dtw_distances,1);
        seg_dtw_distances = [seg_dtw_distances, table((order_dtw_first+1:1:order_dtw_last)','VariableNames',{'points_order'})];
        order_dtw_first = order_dtw_last;
        table_dtw_deviation{i} = seg_dtw_distances;
        % DFD
        metric2postgresql('dfd',frechet_distances, frechet_soll, frechet_ist, bahn_id_, segment_ids{i-1,:})
        table_dfd_info(i,:) = seg_dfd_info;
        order_dfd_last = order_dfd_first + size(seg_dfd_distances,1);
        seg_dfd_distances = [seg_dfd_distances, table((order_dfd_first+1:1:order_dfd_last)','VariableNames',{'points_order'})];
        order_dfd_first = order_dfd_last;
        table_dfd_deviation{i} = seg_dfd_distances;


    end
end

% Berechnung der Kennzahlen für die Gesamtmessung
sidtw = table();
sidtw.bahn_id = {bahn_id_};  
sidtw.calibration_id = {calibration_id};  
sidtw.min_distance = min(table_sidtw_info.sidtw_min_distance); 
sidtw.max_distance = max(table_sidtw_info.sidtw_max_distance);
sidtw.average_distance = mean(table_sidtw_info.sidtw_average_distance);
sidtw.standard_deviation = mean(table_sidtw_info.sidtw_standard_deviation);
sidtw.metrik = "sidtw";

dtw = table();
dtw.bahn_id = {bahn_id_};  
dtw.calibration_id = {calibration_id};  
dtw.min_distance = min(table_dtw_info.dtw_min_distance); 
dtw.max_distance = max(table_dtw_info.dtw_max_distance);
dtw.average_distance = mean(table_dtw_info.dtw_average_distance);
dtw.standard_deviation = mean(table_dtw_info.dtw_standard_deviation);
dtw.metrik = "dtw";

dfd = table();
dfd.bahn_id = {bahn_id_};  
dfd.calibration_id = {calibration_id};  
dfd.min_distance = min(table_dfd_info.dfd_min_distance); 
dfd.max_distance = max(table_dfd_info.dfd_max_distance);
dfd.average_distance = mean(table_dfd_info.dfd_average_distance);
dfd.standard_deviation = mean(table_dfd_info.dfd_standard_deviation);
dfd.metrik = "dfd";

if size(segment_soll,2) > 1
% lcss = table();
% lcss.bahn_id = {bahn_id_};  
% lcss.calibration_id = {calibration_id};  
% lcss.min_distance = min(table_lcss_info.lcss_min_distance); 
% lcss.max_distance = max(table_lcss_info.lcss_max_distance);
% lcss.average_distance = mean(table_lcss_info.lcss_average_distance);
% lcss.standard_deviation = mean(table_lcss_info.lcss_standard_deviation);
% lcss.metrik = "lcss";

table_all_info = table();
table_all_info.bahn_id = {bahn_id_};
table_all_info.calibration_id = {calibration_id};
table_all_info.min_distance = min(table_euclidean_info.euclidean_min_distance);
table_all_info.max_distance = max(table_euclidean_info.euclidean_max_distance);
table_all_info.average_distance = mean(table_euclidean_info.euclidean_average_distance);
table_all_info.standard_deviation = mean(table_euclidean_info.euclidean_standard_deviation);
table_all_info.metrik = "euclidean";

% table_all_info = [table_all_info; sidtw; dtw; dfd; lcss];

else
    table_all_info = [sidtw; dtw; dfd];
end

clear order_eucl_first order_eucl_last order_lcss_first order_lcss_last 
clear order_dfd_first order_dfd_last order_dtw_first order_dtw_last order_sidtw_first order_sidtw_last
clear sidtw sidtw_distances sidtw_ist sidtw_soll seg_sidtw_info seg_sidtw_distances
clear dtw dtw_distances dtw_ist dtw_soll seg_dtw_info seg_dtw_distances
clear dfd frechet_distances frechet_ist frechet_soll frechet_path frechet_matrix frechet_dist frechet_av seg_dfd_info seg_dfd_distances
clear lcss lcss_distances lcss_ist lcss_soll seg_lcss_info seg_lcss_distances
clear euclidean_distances euclidean_ist seg_euclidean_info seg_euclidean_distances
clear pos_ist_trafo segment_ist segment_soll segment_trafo i min_diff 


%% Auswertung der gesamten Messaufnahme 

if evaluate_all == true && evaluate_velocity == false

    % segment_id = bahn_id; 


    if evaluate_orientation == true
        euler_transformation(euler_ist,euler_soll, trafo_euler, trafo_rot)
        % !!!! Zur Sicherheit, damit nicht alle Daten ausgewertet werden !!!!!
        data_all_ist = euler_trans(1:1000,:);
        data_ist_trafo = data_all_ist;
        data_all_soll = euler_soll(1:1000,:);
    else 
        % % !!!! Zur Sicherheit, damit nicht alle Daten ausgewertet werden !!!!!
        % data_all_ist = data_ist(1:600,:);
        % data_all_soll = data_soll(1:500,:);

        % Letztes Segment soll nicht ausgewertet werden
        % Zeilenindex des ersten Auftretens finden
        last_row_ist = find(data_ist.segment_id == segment_ids{end,1}, 1)-1;
        data_all_ist = data_ist(1:last_row_ist,:);
        last_row_soll = find(data_soll.segment_id == segment_ids{end,1}, 1)-1;
        data_all_soll = data_soll(1:last_row_soll,:);

        data_all_ist = table2array(data_all_ist(:,5:7));
        data_all_soll = table2array(data_all_soll(:,5:7));
    
        % Koordinatentrafo für alle Daten 
        coord_transformation(data_all_ist,trafo_rot, trafo_trans)
    end
    
    % Euklidischer Abstand
    tic
    [euclidean_ist,euclidean_distances,~] = distance2curve(data_ist_trafo,data_all_soll,'linear');
    toc
    disp('Euklidischer Abstand berechnet -->')

    % SIDTW
    tic
    [sidtw_distances, ~, ~, ~, sidtw_soll, sidtw_ist, ~, ~, ~] = fkt_selintdtw3d(data_all_soll,data_ist_trafo,false);
    toc
    disp('SIDTW berechnet -->')
    % DTW
    tic
    [dtw_distances, ~, ~, ~, dtw_soll, dtw_ist, ~, ~, ~, ~] = fkt_dtw3d(data_all_soll,data_ist_trafo,false);
    toc
    disp('DTW berechnet -->')
    % Frechet 
    tic
    fkt_discreteFrechet(data_all_soll,data_ist_trafo,false);
    toc
    disp('DFD berechnet -->')
    % LCSS
    tic
    [~, ~, lcss_distances, ~, ~, lcss_soll, lcss_ist, ~, ~] = fkt_lcss(data_all_soll,data_ist_trafo,false);
    toc
    disp('LCSS berechnet -->')

    metric2postgresql('euclidean', euclidean_distances, data_all_soll, euclidean_ist, bahn_id_,bahn_id_)
    metric2postgresql('sidtw', sidtw_distances, sidtw_soll, sidtw_ist, bahn_id_,bahn_id_)
    metric2postgresql('dtw', dtw_distances, dtw_soll, dtw_ist, bahn_id_,bahn_id_)
    metric2postgresql('dfd', frechet_distances, frechet_soll, frechet_ist, bahn_id_,bahn_id_)
    metric2postgresql('lcss', lcss_distances, lcss_soll, lcss_ist, bahn_id_,bahn_id_)

    table_euclidean_info = [seg_euclidean_info; table_euclidean_info];
    table_sidtw_info = [seg_sidtw_info; table_sidtw_info];
    table_dtw_info = [seg_dtw_info; table_dtw_info];
    table_dfd_info = [seg_dfd_info; table_dfd_info];
    table_lcss_info = [seg_lcss_info; table_lcss_info];


%%%%%%%%%% Für die Auswertung in Matlab (für Datenbank irrelevant)
    % Anpassung der Spaltennamen für jede Tabelle
    seg_euclidean_info.Properties.VariableNames = {'bahn_id','segment_id','min_distances', 'max_distance', 'average_distance', 'standard_deviation'};
    seg_sidtw_info.Properties.VariableNames = {'bahn_id','segment_id','min_distances', 'max_distance', 'average_distance', 'standard_deviation'};
    seg_dtw_info.Properties.VariableNames = {'bahn_id','segment_id','min_distances', 'max_distance', 'average_distance', 'standard_deviation'};
    seg_dfd_info.Properties.VariableNames = {'bahn_id','segment_id','min_distances', 'max_distance', 'average_distance', 'standard_deviation'};
    seg_lcss_info.Properties.VariableNames = {'bahn_id','segment_id','min_distances', 'max_distance', 'average_distance', 'standard_deviation'};
    
    table_all_info_2 = [seg_euclidean_info(1,:); seg_sidtw_info(1,:); seg_dtw_info(1,:); seg_dfd_info(1,:); seg_lcss_info(1,:)];

    table_all_info_2.metrik = {'euclidean'; 'sidtw'; 'dtw'; 'dfd'; 'lcss'};

    % Hinzufügen der calibration_id wenn diese existiert
    if exist('calibration_id', 'var') == 1
        calibration_ids = repelem(calibration_id, height(table_all_info_2),1);
        table_all_info_2.calibration_id = calibration_ids;
        table_all_info_2 = table_all_info_2(:,[{'bahn_id'},{'calibration_id'},{'min_distances'},{'max_distance'},{'average_distance'},{'metrik'}]);
        clear calibration_ids
    end
%%%%%%%%%%
    clear sidtw_distances sidtw_ist sidtw_soll seg_sidtw_info seg_sidtw_distances
    clear dtw_distances dtw_ist dtw_soll seg_dtw_info seg_dtw_distances
    clear frechet_distances frechet_ist frechet_soll frechet_path frechet_matrix frechet_dist frechet_av seg_dfd_info seg_dfd_distances
    clear lcss_distances lcss_ist lcss_soll seg_lcss_info seg_lcss_distances
    clear euclidean_distances euclidean_ist seg_euclidean_info seg_euclidean_distances
    clear data_ist_trafo

end

%% Auswertung bestimmter Bahnabschnitte
 
if evaluate_segmentwise == true && evaluate_velocity == false

    n = abs(segment_first - segment_last) + 1;

    if evaluate_orientation == true
        euler_transformation(euler_ist,euler_soll, trafo_euler, trafo_rot)
        data_ist_ = data_ist;
        data_soll_ = data_soll; 
        data_ist = euler_trans;
        data_soll = euler_soll;
    end

    if segment_first == 0 && segment_last == height(segments_ist) - 1
        trajectory_ist = data_ist(1:end,:);
        trajectory_soll = data_soll(1:end,:);

    elseif segment_last == height(segments_ist) - 1 && segment_first > 0
        trajectory_ist = data_ist(idx_new_seg_ist(segment_first):end,:);
        trajectory_soll = data_soll(idx_new_seg_soll(segment_first):end,:);

    elseif segment_first == 0 && segment_last < height(segments_ist) - 1
        trajectory_ist = data_ist(1:idx_new_seg_ist(segment_last+1)-1,:);
        trajectory_soll = data_soll(1:idx_new_seg_soll(segment_last+1)-1,:);

    elseif segment_first > 0 && segment_first < height(segments_ist) - 2 && segment_last < height(segments_ist) - 1

        trajectory_ist = data_ist(idx_new_seg_ist(segment_first):idx_new_seg_ist(segment_last+1)-1,:);
        trajectory_soll = data_soll(idx_new_seg_soll(segment_first):idx_new_seg_soll(segment_last+1)-1,:);
    else
        error('Die ausgewählten Bahnabschnitte sind in den Daten nicht vorhanden!')
    end
    
    if evaluate_orientation == true
        data_ist_trafo = trajectory_ist;
        data_ist = data_ist_; data_soll = data_soll_;
        clear data_ist_ data_soll_ 
    else
        trajectory_ist = table2array(trajectory_ist(:,5:7));
        trajectory_soll = table2array(trajectory_soll(:,5:7));
    
        % Koordinatentransformation
        coord_transformation(trajectory_ist,trafo_rot, trafo_trans)
    end

    % Euklidischer Abstand 
    [euclidean_ist,euclidean_distances,~] = distance2curve(data_ist_trafo,trajectory_soll,'linear');
    % % SIDTW
    [sidtw_distances, ~, ~, ~, sidtw_soll, sidtw_ist, ~, ~, ~] = fkt_selintdtw3d(trajectory_soll,data_ist_trafo,false);
    % DTW
    [dtw_distances, ~, ~, ~, dtw_soll, dtw_ist, ~, ~, ~, ~] = fkt_dtw3d(trajectory_soll,data_ist_trafo,false);
    % Frechet 
    fkt_discreteFrechet(trajectory_soll,data_ist_trafo,false);
    % LCSS
    [~, ~, lcss_distances, ~, ~, lcss_soll, lcss_ist, ~, ~] = fkt_lcss(trajectory_soll,data_ist_trafo,false);

    metric2postgresql('euclidean', euclidean_distances, trajectory_soll, euclidean_ist, bahn_id_)
    metric2postgresql('sidtw', sidtw_distances, sidtw_soll, sidtw_ist, bahn_id_)
    metric2postgresql('dtw', dtw_distances, dtw_soll, dtw_ist, bahn_id_)
    metric2postgresql('dfd', frechet_distances, frechet_soll, frechet_ist, bahn_id_)
    metric2postgresql('lcss', lcss_distances, lcss_soll, lcss_ist, bahn_id_)



    % Anpassung der Spaltennamen für jede Tabelle
    seg_euclidean_info.Properties.VariableNames = {'bahn_id','min_distances', 'max_distance', 'average_distance', 'standard_deviation'};
    seg_sidtw_info.Properties.VariableNames = {'bahn_id','min_distances', 'max_distance', 'average_distance', 'standard_deviation'};
    seg_dtw_info.Properties.VariableNames = {'bahn_id','min_distances', 'max_distance', 'average_distance', 'standard_deviation'};
    seg_dfd_info.Properties.VariableNames = {'bahn_id','min_distances', 'max_distance', 'average_distance', 'standard_deviation'};
    seg_lcss_info.Properties.VariableNames = {'bahn_id','min_distances', 'max_distance', 'average_distance', 'standard_deviation'};
    
    trajectory_info = [seg_euclidean_info(1,:); seg_sidtw_info(1,:); seg_dtw_info(1,:); seg_dfd_info(1,:); seg_lcss_info(1,:)];

    trajectory_info.metrik = {'euclidean'; 'sidtw'; 'dtw'; 'dfd'; 'lcss'};

    clear sidtw_distances sidtw_ist sidtw_soll seg_sidtw_info seg_sidtw_distances
    clear dtw_distances dtw_ist dtw_soll seg_dtw_info seg_dtw_distances
    clear frechet_distances frechet_ist frechet_soll frechet_path frechet_matrix frechet_dist frechet_av seg_dfd_info seg_dfd_distances
    clear lcss_distances lcss_ist lcss_soll seg_lcss_info seg_lcss_distances
    clear euclidean_distances euclidean_ist seg_euclidean_info seg_euclidean_distances
    clear data_ist_trafo

end

%% Plotten
plots = 1;
if plots == true
    % Farben
    c1 = [0 0.4470 0.7410];
    c2 = [0.8500 0.3250 0.0980];
    c3 = [0.9290 0.6940 0.1250];
    c4 = [0.4940 0.1840 0.5560];
    c5 = [0.4660 0.6740 0.1880];
    c6 = [0.3010 0.7450 0.9330];
    c7 = [0.6350 0.0780 0.1840];
end

if plots == true && evaluate_velocity == false && evaluate_orientation == false
    
    ist = table2array(data_ist(:,5:7));
    soll = table2array(data_soll(:,5:7));
    
    % Koordinatentrafo für alle Daten 
    coord_transformation(ist,trafo_rot, trafo_trans)
    ist = data_ist_trafo; 
    clear data_ist_trafo
    
    % Plot der Gesamten Bahn
    f0 = figure('Color','white','Name','Soll und Istbahn (gesamte Messung)');
    f0.Position(3:4) = [1520 840];
    hold on 
    plot3(soll(:,1),soll(:,2),soll(:,3),Color=c1,LineWidth=1.5)
    plot3(ist(:,1),ist(:,2),ist(:,3),Color=c2,LineWidth=1.5)
    xlabel('x','FontWeight','bold');
    ylabel('y','FontWeight','bold');
    zlabel('z','FontWeight','bold','Rotation',0);
    legend('Sollbahn (ABB)','Istbahn (VICON)')
    grid on 
    view(3)

    % Plotten der ausgewählten Segmente und deren Abweichungen 
    if evaluate_segmentwise == true 
    
        f1 = figure('Color','white','Name','Soll- und Istbahnen (Bahnsegmente)');
        f1.Position(3:4) = [1520 840];
        hold on
        if segment_last-segment_first == 0
            plot3(segments_trafo.x_ist{segment_first+1,1},segments_trafo.y_ist{segment_first+1,1},segments_trafo.z_ist{segment_first+1,1},Color=c1,LineWidth=1.5)
            plot3(segments_soll.x_soll{segment_first+1,1},segments_soll.y_soll{segment_first+1,1},segments_soll.z_soll{segment_first+1,1},Color=c2,LineWidth=1.5)
        else
            for i = segment_first:1:segment_last-segment_first
                plot3(segments_trafo.x_ist{i+1,1},segments_trafo.y_ist{i+1,1},segments_trafo.z_ist{i+1,1},Color=c1,LineWidth=1.5)
                plot3(segments_soll.x_soll{i+1,1},segments_soll.y_soll{i+1,1},segments_soll.z_soll{i+1,1},Color=c2,LineWidth=1.5)
            end
        end
        title("Bahnabschnitte " + num2str(segment_first) + " bis " + num2str(segment_last))
        xlabel('x','FontWeight','bold'); ylabel('y','FontWeight','bold'); zlabel('z','FontWeight','bold','Rotation',0);
        legend('Sollbahn (ABB)','Istbahn (VICON)')
        grid on 
        view(3)
        
        f2 = figure('Color','white','Name','Mittlere Abweichungen (Bahnsegmente)');
        f2.Position(3:4) = [1520 840];
        subplot(2,2,1)
        title('Mittlere Abweichungen zwischen Soll- und Istbahn')
        hold on 
        plot(linspace(segment_first,segment_last,segment_last-segment_first+1),table_dtw_info.dtw_average_distance(segment_first+1:segment_last+1,:),LineWidth=2.5,Color=c1)
        plot(linspace(segment_first,segment_last,segment_last-segment_first+1),table_dfd_info.dfd_average_distance(segment_first+1:segment_last+1,:),LineWidth=2.5,Color=c2)
        plot(linspace(segment_first,segment_last,segment_last-segment_first+1),table_lcss_info.lcss_average_distance(segment_first+1:segment_last+1,:),LineWidth=2.5,Color=c3)
        plot(linspace(segment_first,segment_last,segment_last-segment_first+1),table_sidtw_info.sidtw_average_distance(segment_first+1:segment_last+1,:),LineWidth=2.5,Color=c4)
        plot(linspace(segment_first,segment_last,segment_last-segment_first+1),table_euclidean_info.euclidean_average_distance(segment_first+1:segment_last+1,:),LineWidth=2.5,Color=c5)
        % Plot der Werte der gesamten Bahn
        plot(linspace(segment_first,segment_last,segment_last-segment_first+1),repelem(table_all_info.average_distance(3,:),segment_last-segment_first+1,1),LineWidth=1.2,Color=c1)
        plot(linspace(segment_first,segment_last,segment_last-segment_first+1),repelem(table_all_info.average_distance(4,:),segment_last-segment_first+1,1),LineWidth=1.2,Color=c2)
        plot(linspace(segment_first,segment_last,segment_last-segment_first+1),repelem(table_all_info.average_distance(5,:),segment_last-segment_first+1,1),LineWidth=1.2,Color=c3)
        plot(linspace(segment_first,segment_last,segment_last-segment_first+1),repelem(table_all_info.average_distance(2,:),segment_last-segment_first+1,1),LineWidth=1.2,Color=c4)
        plot(linspace(segment_first,segment_last,segment_last-segment_first+1),repelem(table_all_info.average_distance(1,:),segment_last-segment_first+1,1),LineWidth=1.2,Color=c5)
        xlabel('Bahnsegmente');
        ylabel('Abweichung in mm');
        legend('DTW','DFD','LCSS','SIDTW','Eukl. Dist.')
        grid on
        axis padded
        
        subplot(2,2,2)
        title('Maximale Abweichungen zwischen Soll- und Istbahn')
        hold on 
        plot(linspace(segment_first,segment_last,segment_last-segment_first+1),table_dtw_info.dtw_max_distance(segment_first+1:segment_last+1,:),LineWidth=2.5,Color=c1)
        plot(linspace(segment_first,segment_last,segment_last-segment_first+1),table_dfd_info.dfd_max_distance(segment_first+1:segment_last+1,:),LineWidth=2.5,Color=c2)
        plot(linspace(segment_first,segment_last,segment_last-segment_first+1),table_lcss_info.lcss_max_distance(segment_first+1:segment_last+1,:),LineWidth=2.5,Color=c3)
        plot(linspace(segment_first,segment_last,segment_last-segment_first+1),table_sidtw_info.sidtw_max_distance(segment_first+1:segment_last+1,:),LineWidth=2.5,Color=c4)
        plot(linspace(segment_first,segment_last,segment_last-segment_first+1),table_euclidean_info.euclidean_max_distance(segment_first+1:segment_last+1,:),LineWidth=2.5,Color=c5)
        % Plot der Werte der gesamten Bahn
        plot(linspace(segment_first,segment_last,segment_last-segment_first+1),repelem(table_all_info.max_distance(3,:),segment_last-segment_first+1,1),LineWidth=1.2,Color=c1)
        plot(linspace(segment_first,segment_last,segment_last-segment_first+1),repelem(table_all_info.max_distance(4,:),segment_last-segment_first+1,1),LineWidth=1.2,Color=c2)
        plot(linspace(segment_first,segment_last,segment_last-segment_first+1),repelem(table_all_info.max_distance(5,:),segment_last-segment_first+1,1),LineWidth=1.2,Color=c3)
        plot(linspace(segment_first,segment_last,segment_last-segment_first+1),repelem(table_all_info.max_distance(2,:),segment_last-segment_first+1,1),LineWidth=1.2,Color=c4)
        plot(linspace(segment_first,segment_last,segment_last-segment_first+1),repelem(table_all_info.max_distance(1,:),segment_last-segment_first+1,1),LineWidth=1.2,Color=c5)
        xlabel('Bahnsegmente');
        ylabel('Abweichung in mm');
        legend('DTW','DFD','LCSS','SIDTW','Eukl. Dist.')
        grid on
        axis padded
        
        subplot(2,2,3)
        title('Minimale Abweichungen zwischen Soll- und Istbahn')
        hold on 
        plot(linspace(segment_first,segment_last,segment_last-segment_first+1),table_dtw_info.dtw_min_distance(segment_first+1:segment_last+1,:),LineWidth=2.5,Color=c1)
        plot(linspace(segment_first,segment_last,segment_last-segment_first+1),table_dfd_info.dfd_min_distance(segment_first+1:segment_last+1,:),LineWidth=2.5,Color=c2)
        plot(linspace(segment_first,segment_last,segment_last-segment_first+1),table_lcss_info.lcss_min_distance(segment_first+1:segment_last+1,:),LineWidth=2.5,Color=c3)
        plot(linspace(segment_first,segment_last,segment_last-segment_first+1),table_sidtw_info.sidtw_min_distance(segment_first+1:segment_last+1,:),LineWidth=2.5,Color=c4)
        plot(linspace(segment_first,segment_last,segment_last-segment_first+1),table_euclidean_info.euclidean_min_distance(segment_first+1:segment_last+1,:),LineWidth=2.5,Color=c5)
        % Plot der Werte der gesamten Bahn
        plot(linspace(segment_first,segment_last,segment_last-segment_first+1),repelem(table_all_info.min_distance(3,:),segment_last-segment_first+1,1),LineWidth=1.2,Color=c1)
        plot(linspace(segment_first,segment_last,segment_last-segment_first+1),repelem(table_all_info.min_distance(4,:),segment_last-segment_first+1,1),LineWidth=1.2,Color=c2)
        plot(linspace(segment_first,segment_last,segment_last-segment_first+1),repelem(table_all_info.min_distance(5,:),segment_last-segment_first+1,1),LineWidth=1.2,Color=c3)
        plot(linspace(segment_first,segment_last,segment_last-segment_first+1),repelem(table_all_info.min_distance(2,:),segment_last-segment_first+1,1),LineWidth=1.2,Color=c4)
        plot(linspace(segment_first,segment_last,segment_last-segment_first+1),repelem(table_all_info.min_distance(1,:),segment_last-segment_first+1,1),LineWidth=1.2,Color=c5)
        xlabel('Bahnsegmente');
        ylabel('Abweichung in mm');
        legend('DTW','DFD','LCSS','SIDTW','Eukl. Dist.')
        grid on
        axis padded
        
        subplot(2,2,4)
        title('Standardabweichungen zwischen Soll- und Istbahn')
        hold on 
        plot(linspace(segment_first,segment_last,segment_last-segment_first+1),table_dtw_info.dtw_standard_deviation(segment_first+1:segment_last+1,:),LineWidth=2.5,Color=c1)
        plot(linspace(segment_first,segment_last,segment_last-segment_first+1),table_dfd_info.dfd_standard_deviation(segment_first+1:segment_last+1,:),LineWidth=2.5,Color=c2)
        plot(linspace(segment_first,segment_last,segment_last-segment_first+1),table_lcss_info.lcss_standard_deviation(segment_first+1:segment_last+1,:),LineWidth=2.5,Color=c3)
        plot(linspace(segment_first,segment_last,segment_last-segment_first+1),table_sidtw_info.sidtw_standard_deviation(segment_first+1:segment_last+1,:),LineWidth=2.5,Color=c4)
        plot(linspace(segment_first,segment_last,segment_last-segment_first+1),table_euclidean_info.euclidean_standard_deviation(segment_first+1:segment_last+1,:),LineWidth=2.5,Color=c5)
        % Plot der Werte der gesamten Bahn
        plot(linspace(segment_first,segment_last,segment_last-segment_first+1),repelem(table_all_info.standard_deviation(3,:),segment_last-segment_first+1,1),LineWidth=1.2,Color=c1)
        plot(linspace(segment_first,segment_last,segment_last-segment_first+1),repelem(table_all_info.standard_deviation(4,:),segment_last-segment_first+1,1),LineWidth=1.2,Color=c2)
        plot(linspace(segment_first,segment_last,segment_last-segment_first+1),repelem(table_all_info.standard_deviation(5,:),segment_last-segment_first+1,1),LineWidth=1.2,Color=c3)
        plot(linspace(segment_first,segment_last,segment_last-segment_first+1),repelem(table_all_info.standard_deviation(2,:),segment_last-segment_first+1,1),LineWidth=1.2,Color=c4)
        plot(linspace(segment_first,segment_last,segment_last-segment_first+1),repelem(table_all_info.standard_deviation(1,:),segment_last-segment_first+1,1),LineWidth=1.2,Color=c5)
        xlabel('Bahnsegmente');
        ylabel('Abweichung in mm');
        legend('DTW','DFD','LCSS','SIDTW','Eukl. Dist.')
        grid on
        axis padded   
    end
end

% Plotten der Euler-Winkel von Soll- und Ist-Bahn
if plots == true && evaluate_orientation == true && evaluate_velocity == false

    % Timestamps in Sekunden
    time_ist = str2double(data_ist.timestamp);
    time_soll = str2double(data_soll.timestamp);
    timestamps_ist = (time_ist(:,1)- time_soll(1,1))/1e9;
    timestamps_soll = (time_soll(:,1)- time_soll(1,1))/1e9;
    

    % Transformation aller Winkel 
    euler_transformation(euler_ist,euler_soll, trafo_euler, trafo_rot)
    % Winkel zwischen 0 - 360°
    euler_soll = mod(euler_soll,360);
    euler_trans = mod(euler_trans,360);

    % Plot 
    figure('Color','white','Name','Eulerwinkel von 0° bis 360°')
    hold on 
    plot(timestamps_soll,euler_soll(:,1),Color=c1,LineWidth=1.5)
    plot(timestamps_soll,euler_soll(:,2),Color=c2,LineWidth=1.5)
    plot(timestamps_soll,euler_soll(:,3),Color=c4,LineWidth=1.5)
    plot(timestamps_ist,euler_trans(:,1),Color=c1)
    plot(timestamps_ist,euler_trans(:,2),Color=c2)
    plot(timestamps_ist,euler_trans(:,3),Color=c4)
    xlabel('Zeit [s]'); ylabel('Winkel [°]');
    legend("roll","pitch","yaw")
    hold off

    % % Plot 
    % figure('Color','white','Name','Mit Originalzeit')
    % hold on 
    % plot(str2double(data_soll.timestamp),euler_soll(:,1),Color=c1,LineWidth=1.5)
    % plot(str2double(data_soll.timestamp),euler_soll(:,2),Color=c2,LineWidth=1.5)
    % plot(str2double(data_soll.timestamp),euler_soll(:,3),Color=c4,LineWidth=1.5)
    % plot(str2double(data_ist.timestamp),euler_trans(:,1),Color=c1)
    % plot(str2double(data_ist.timestamp),euler_trans(:,2),Color=c2)
    % plot(str2double(data_ist.timestamp),euler_trans(:,3),Color=c4)
    % xlabel('Zeit [s]'); ylabel('Winkel [°]');
    % legend("roll","pitch","yaw")
    % hold off

    % Segmentweise plotten
    k1 = 1;
    k2 = 1; 
    figure('Color','white','Name','Eulerwinkel aller Segmente (-180° bis 180°)')
    hold on
    for i = 1:1:num_segments
        plot(timestamps_soll(k1:k1+length(segments_soll.roll_soll{i})-1),segments_soll.roll_soll{i},Color=c1,LineWidth=1.5)
        plot(timestamps_soll(k1:k1+length(segments_soll.roll_soll{i})-1),segments_soll.pitch_soll{i},Color=c2,LineWidth=1.5)
        plot(timestamps_soll(k1:k1+length(segments_soll.roll_soll{i})-1),segments_soll.yaw_soll{i},Color=c4,LineWidth=1.5)
        plot(timestamps_ist(k2:k2+length(segments_trafo.roll_ist{i})-1),segments_trafo.roll_ist{i},Color=c1)
        plot(timestamps_ist(k2:k2+length(segments_trafo.roll_ist{i})-1),segments_trafo.pitch_ist{i},Color=c2)
        plot(timestamps_ist(k2:k2+length(segments_trafo.roll_ist{i})-1),segments_trafo.yaw_ist{i},Color=c4)
        k1 = k1 + length(segments_soll.roll_soll{i});
        k2 = k2 + length(segments_trafo.roll_ist{i});
    end
    xlabel('Zeit [s]'); ylabel('Winkel [°]');
    legend("roll","pitch","yaw")
    hold off

    % 3D - Plot    
    figure('Color','white','Name','Eulerwinkel 3D')
    hold on
    plot3(euler_soll(:,1),euler_soll(:,2),euler_soll(:,3),Color=c1,LineWidth=1.5);
    plot3(euler_trans(:,1),euler_trans(:,2),euler_trans(:,3),Color=c2,LineWidth=1.5);
    legend('Soll','Ist')
    xlabel('Roll'); ylabel('Pitch'); zlabel('Yaw')
    view(3)
    hold off
    axis equal
end

clear c1 c2 c3 c4 c5 c6 c7 segment_first segment_last n i f0 f1 f2 k1 k2 time

toc;


%% Geschwindigkeit

% euler_transformation(euler_ist,euler_soll, trafo_euler)
% 
% % vel_soll_time = data_soll(:,4);
% 
% vel_soll = data_soll(:,4:7);
% vel_soll.timestamp = str2double(vel_soll.timestamp);
% % vel_soll = table2array(vel_soll);
% 
% % Timestamps in Sekunden
% vel_soll.timestamp = (vel_soll.timestamp - vel_soll.timestamp(1,1))/1e9;

%% Ergebnisse in die Datenbank hochladen
% 
% upload = false; 
% 
% if upload == true
% 
%     upload = input("Eingabe 'upload' wenn die Ergebnisse hochgeladen werden soll. \n",'s');
%     if strcmp(upload,'upload')
%         disp('Upload erfolgt!')
% 
%         tablename = 'robotervermessung.auswertung.euclidean_info';
%         sqlwrite(conn,tablename,table_euclidean_info)
%         a = sqlread(conn,tablename);
% 
%     else
%         disp('Upload fehlgeschlagen!')
%     end
% 
% end

% bahn_info = sqlread(conn,'robotervermessung.bewegungsdaten.bahn_info');
% 
% % 1. Möglichkeit
% search = find(bahn_info.calibration_run == 1);
% bahn_calibration = bahn_info(search,:); clear search 
% 
% % 2. Möglichkeit 
% bahn_calibration_runs = bahn_info(bahn_info.calibration_run == 1, :);
% %%
% tablename = 'robotervermessung.auswertung.euclidean_info'
% a = sqlread(conn,tablename)
% %%
% sqlquery = strcat("DROP TABLE ",tablename);
% execute(conn,sqlquery)
% %%
% % table_euclidean_info.segment_id = convertStringsToChars(table_euclidean_info.segment_id)
% %%
% for i = 1:1:num_segments
% 
% 
% 
% end


% Der Variablen das erste Segment hinzufügen
first_id = table(bahn_id_+"_0",'VariableNames',"segment_id");
segment_ids = [first_id; segment_ids];
if evaluate_all == true
    first_id = table(string(bahn_id_),'VariableNames',"segment_id");
    segment_ids = [first_id; segment_ids];
end
clear first_id

if upload == true
    tic;
    % upload = input("Eingabe 'upload' wenn die Ergebnisse hochgeladen werden soll. \n",'s');
    % if strcmp(upload,'upload')
        disp('Upload erfolgt!')

        if evaluate_velocity == false && evaluate_orientation == false

            type = {'position'};
            % Info Tabellen 
            table_euclidean_info = addvars(table_euclidean_info,repelem(type,height(table_euclidean_info),1),'NewVariableNames','evaluation');
            table_sidtw_info = addvars(table_sidtw_info,repelem(type,height(table_sidtw_info),1),'NewVariableNames','evaluation');
            table_dtw_info = addvars(table_dtw_info,repelem(type,height(table_dtw_info),1),'NewVariableNames','evaluation');
            table_dfd_info = addvars(table_dfd_info,repelem(type,height(table_dfd_info),1),'NewVariableNames','evaluation');
            table_lcss_info = addvars(table_lcss_info,repelem(type,height(table_lcss_info),1),'NewVariableNames','evaluation');
            % Tabellen mit allen Abweichungen 
            for i = 1:1:size(table_euclidean_deviation,1)-1
                table_euclidean_deviation{i} = addvars(table_euclidean_deviation{i},repelem(type,height(table_euclidean_deviation{i}),1),'NewVariableNames','evaluation');
                table_sidtw_deviation{i} = addvars(table_sidtw_deviation{i},repelem(type,height(table_sidtw_deviation{i}),1),'NewVariableNames','evaluation');
                table_dtw_deviation{i} = addvars(table_dtw_deviation{i},repelem(type,height(table_dtw_deviation{i}),1),'NewVariableNames','evaluation');
                table_dfd_deviation{i} = addvars(table_dfd_deviation{i},repelem(type,height(table_dfd_deviation{i}),1),'NewVariableNames','evaluation');
                table_lcss_deviation{i} = addvars(table_lcss_deviation{i},repelem(type,height(table_lcss_deviation{i}),1),'NewVariableNames','evaluation');
            end




        end

        if evaluate_velocity == false && evaluate_orientation == true  

            type = {'orientation'};
            % Info Tabellen 
            table_euclidean_info = addvars(table_euclidean_info,repelem(type,height(table_euclidean_info),1),'NewVariableNames','evaluation');
            table_sidtw_info = addvars(table_sidtw_info,repelem(type,height(table_sidtw_info),1),'NewVariableNames','evaluation');
            table_dtw_info = addvars(table_dtw_info,repelem(type,height(table_dtw_info),1),'NewVariableNames','evaluation');
            table_dfd_info = addvars(table_dfd_info,repelem(type,height(table_dfd_info),1),'NewVariableNames','evaluation');
            table_lcss_info = addvars(table_lcss_info,repelem(type,height(table_lcss_info),1),'NewVariableNames','evaluation');
            % Tabellen mit allen Abweichungen
            for i = 1:1:size(table_euclidean_deviation,1)-1
                table_euclidean_deviation{i} = addvars(table_euclidean_deviation{i},repelem(type,height(table_euclidean_deviation{i}),1),'NewVariableNames','evaluation');
                table_sidtw_deviation{i} = addvars(table_sidtw_deviation{i},repelem(type,height(table_sidtw_deviation{i}),1),'NewVariableNames','evaluation');
                table_dtw_deviation{i} = addvars(table_dtw_deviation{i},repelem(type,height(table_dtw_deviation{i}),1),'NewVariableNames','evaluation');
                table_dfd_deviation{i} = addvars(table_dfd_deviation{i},repelem(type,height(table_dfd_deviation{i}),1),'NewVariableNames','evaluation');
                table_lcss_deviation{i} = addvars(table_lcss_deviation{i},repelem(type,height(table_lcss_deviation{i}),1),'NewVariableNames','evaluation');
            end
        end

        if evaluate_velocity == true && evaluate_orientation == false

            type = {'speed'};
            % Info Tabellen 
            table_euclidean_info = addvars(table_euclidean_info,repelem(type,height(table_euclidean_info),1),'NewVariableNames','evaluation');
            table_sidtw_info = addvars(table_sidtw_info,repelem(type,height(table_sidtw_info),1),'NewVariableNames','evaluation');
            table_dtw_info = addvars(table_dtw_info,repelem(type,height(table_dtw_info),1),'NewVariableNames','evaluation');
            table_dfd_info = addvars(table_dfd_info,repelem(type,height(table_dfd_info),1),'NewVariableNames','evaluation');
            table_lcss_info = addvars(table_lcss_info,repelem(type,height(table_lcss_info),1),'NewVariableNames','evaluation');
            % Tabellen mit allen Abweichungen
            for i = 1:1:size(table_euclidean_deviation,1)-1
                table_euclidean_deviation{i} = addvars(table_euclidean_deviation{i},repelem(type,height(table_euclidean_deviation{i}),1),'NewVariableNames','evaluation');
                table_sidtw_deviation{i} = addvars(table_sidtw_deviation{i},repelem(type,height(table_sidtw_deviation{i}),1),'NewVariableNames','evaluation');
                table_dtw_deviation{i} = addvars(table_dtw_deviation{i},repelem(type,height(table_dtw_deviation{i}),1),'NewVariableNames','evaluation');
                table_dfd_deviation{i} = addvars(table_dfd_deviation{i},repelem(type,height(table_dfd_deviation{i}),1),'NewVariableNames','evaluation');
                table_lcss_deviation{i} = addvars(table_lcss_deviation{i},repelem(type,height(table_lcss_deviation{i}),1),'NewVariableNames','evaluation');
            end
        end
        %%

        % Schreiben in die Datenbank
        upload2postgresql('robotervermessung.auswertung.sidtw_info',table_sidtw_info,segment_ids,type{1},conn)
        upload2postgresql('robotervermessung.auswertung.euclidean_info',table_euclidean_info,segment_ids,type{1},conn)
        upload2postgresql('robotervermessung.auswertung.dtw_info',table_dtw_info,segment_ids,type{1},conn)
        upload2postgresql('robotervermessung.auswertung.dfd_info',table_dfd_info,segment_ids,type{1},conn)
        upload2postgresql('robotervermessung.auswertung.lcss_info',table_lcss_info,segment_ids,type{1},conn)
%%
        % !!!!! (dauert extrem lange) !!!!!
        if evaluate_all == true 
            segment_ids = segment_ids(2:end,:); % segment_id = bahn_id löschen!
        end
        upload2postgresql('robotervermessung.auswertung.euclidean_deviation',table_euclidean_deviation,segment_ids,type{1},conn)
        disp('Euclidean Deviation hochgeladen')
        upload2postgresql('robotervermessung.auswertung.sidtw_deviation',table_sidtw_deviation,segment_ids,type{1},conn)
        disp('SIDTW Deviation hochgeladen')
        upload2postgresql('robotervermessung.auswertung.dtw_deviation',table_dtw_deviation,segment_ids,type{1},conn)
        disp('DTW Deviation hochgeladen')
        upload2postgresql('robotervermessung.auswertung.dfd_deviation',table_dfd_deviation,segment_ids,type{1},conn)
        disp('DFD Deviation hochgeladen')
        upload2postgresql('robotervermessung.auswertung.lcss_deviation',table_lcss_deviation,segment_ids,type{1},conn)

        disp('Der Upload war erfolgreich!')
    else
        disp('Upload fehlgeschlagen!')
    end
    toc;
% end


%%

function upload2postgresql(tablename,table,segment_ids,evaluation,conn)

% Abfrage ob der Eintrag der gesamten Bahn bereits existiert
checkQuery = sprintf("SELECT COUNT(*) FROM %s WHERE bahn_id = '%s' AND evaluation = '%s'", tablename, convertCharsToStrings(table{1,1}{1,1}), evaluation);
duplicates = fetch(conn, checkQuery);
entryExists = duplicates{1,1} > 0;

if entryExists == false
    
    % Daten der Segmente als eine einzige Tabelle schreiben
    if iscell(table)
        table = vertcat(table{:});
    end

    % Schreiben der gesamten Tabelle in die Datenbank
    sqlwrite(conn,tablename,table)
else
    
    % Segmentweise prüfen wenn Eintrag bereits existiert
    for i = 1:1:size(segment_ids,1)-1

        % Abfrage ob der Eintrag eines Segmentes bereits existiert
        checkQuery = sprintf("SELECT COUNT(*) FROM %s WHERE segment_id = '%s' AND evaluation = '%s'", tablename, segment_ids{i,1}, evaluation);
        duplicates = fetch(conn, checkQuery);
        entryExists = duplicates{1,1} > 0;

        % Lösche Daten falls diese bereits existieren
        if entryExists == true 
            deleteQuery = sprintf("DELETE FROM %s WHERE segment_id = '%s' AND evaluation = '%s'", tablename, segment_ids{i,1}, evaluation);
            execute(conn, deleteQuery);
        end

        if iscell(table)
            % Schreiben/Überschreiben der Daten in die Datenbank
            sqlwrite(conn,tablename,table{i,:})
        else          
            sqlwrite(conn,tablename,table(i,:))
        end
    end
end
end





function [table_euclidean_info, table_sidtw_info, table_dtw_info, table_dfd_info, table_lcss_info, ...
          table_euclidean_deviation, table_sidtw_deviation, table_dtw_deviation, ...
          table_dfd_deviation, table_lcss_deviation, segment_ids, calibration_id] = ...
              processTrajectoryLT(conn, bahn_id, evaluate_velocity, ...
                                      evaluate_orientation, plots)
        
    % Collection aus der die Daten extrahiert werden (evtl. übergeben)
    schema = 'bewegungsdaten';
    
    % Prüfe die Datenquelle (Lasertracker AT960 oder Motion Capture Vicon)
    query = ['SELECT source_data_ist FROM robotervermessung.' schema '.bahn_info ' ...
             'WHERE bahn_id = ''' bahn_id ''''];
    data_source_result = fetch(conn, query);
    
    if isempty(data_source_result)
        error(['Keine Quellenangabe für Bahn-ID ' bahn_id ' gefunden.']);
    end
    
    data_source = data_source_result.source_data_ist{1};
    disp(['Datenquelle für Bahn-ID ' bahn_id ': ' data_source]);
    
    % Lasertracker-Transformation oder MoCap-Transformation
    if strcmpi(data_source, 'leica_at960')
        % Lasertracker-Transformation
        disp('Verwende Lasertracker-Transformation für AT960');
        
        % Extrahiere das Aufnahmedatum aus bahn_info für die Matrix-Suche
        query = ['SELECT recording_date FROM robotervermessung.' schema '.bahn_info WHERE bahn_id = ''' bahn_id ''''];
        date_result = fetch(conn, query);
        
        if ~isempty(date_result) && ~isempty(date_result.recording_date)
            raw_date = date_result.recording_date{1};
            date_parts = split(raw_date, ' '); 
            date_only = date_parts{1};
            date_components = split(date_only, '-');
            record_date = [date_components{3} '.' date_components{2} '.' date_components{1}];
            
            disp(['Aufnahmedatum der Bahn: ' record_date]);
        else
            disp('Kein Aufnahmedatum gefunden, verwende aktuelles Datum für Matrizen');
            record_date = datestr(datetime('now'), 'dd.MM.yyyy');
        end
        
        % Lade Transformationsmatrizen
        try
            [offset_matrix, robot_to_instrument, matrix_file] = findTransfMatrices(record_date);
            disp(['Transformationsmatrizen für Datum ' record_date ' geladen']);
            
            % Setze Kalibrierungs-ID für Lasertracker-Transformation
            calibration_id = matrix_file;
        catch ME
            error(['Fehler beim Laden der Transformationsmatrizen: ' ME.message]);
        end
        
        % Auslesen der zu transformierenden Daten
        query = ['SELECT * FROM robotervermessung.' schema '.bahn_pose_ist ' ...
                'WHERE robotervermessung.' schema '.bahn_pose_ist.bahn_id = ''' bahn_id ''''];
        data_ist = fetch(conn, query);
        data_ist = sortrows(data_ist, 'timestamp');
        
        query = ['SELECT * FROM robotervermessung.' schema '.bahn_orientation_soll ' ...
                'WHERE robotervermessung.' schema '.bahn_orientation_soll.bahn_id = ''' bahn_id ''''];
        data_orientation_soll = fetch(conn, query);
        data_orientation_soll = sortrows(data_orientation_soll, 'timestamp');
        
        query = ['SELECT * FROM robotervermessung.' schema '.bahn_position_soll ' ...
                'WHERE robotervermessung.' schema '.bahn_position_soll.bahn_id = ''' bahn_id ''''];
        data_position_soll = fetch(conn, query);
        data_position_soll = sortrows(data_position_soll, 'timestamp');
        
        % Daten extrahieren
        position_ist = table2array(data_ist(:,5:7));     % IST-Position [x, y, z]
        
        % Synchronisiere die Zeitreihen
        quaternion_ist = table2array(data_ist(:,8:11));  % IST-Quaternion [qx, qy, qz, qw]
        
        % Initialisierung der Ergebnisarrays
        num_points = size(position_ist, 1);
        positions_transformed = zeros(num_points, 3);
        quaternions_transformed = zeros(num_points, 4);
        orientations_transformed = zeros(num_points, 3);
        
        % Lasertracker-Transformation
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
        q_transformed = [quaternions_transformed(:,2), quaternions_transformed(:,3), quaternions_transformed(:,4), quaternions_transformed(:,1)];
        
        % Übergebe die bereits transformierten Daten an die getSegmentsLT-Funktion
        if evaluate_orientation
            [num_segments, segment_ids, data_ist, data_soll, segments_ist, segments_soll, segments_trafo, q_transformed] = ...
            getSegmentsLT(conn, bahn_id, schema, evaluate_orientation, evaluate_velocity, [], [], [], positions_transformed, q_transformed);
        else
            [num_segments, segment_ids, data_ist, data_soll, segments_ist, segments_soll, segments_trafo, ~] = ...
            getSegmentsLT(conn, bahn_id, schema, evaluate_orientation, evaluate_velocity, [], [], [], positions_transformed, []);
        end
        
    else
        % MoCap/Vicon-Transformation (bestehender Code)
        disp('Verwende Motion-Capture-Transformation für Vicon');
        
        % Suche nach passender Kalibrierungsdatei
        [calibration_id, is_calibration_run] = findCalibrationRun(conn, bahn_id, schema);

        % Extraktion der Soll- und Ist-Daten der Kalibrierungsdatei
        tablename_cal = ['robotervermessung.' schema '.bahn_pose_ist'];
        query = sprintf("SELECT * FROM %s WHERE bahn_id = '%s'", tablename_cal, calibration_id);
        data_cal_ist = fetch(conn, query);
        data_cal_ist = sortrows(data_cal_ist,'timestamp');

        tablename_cal = ['robotervermessung.' schema '.bahn_events'];
        query = sprintf("SELECT * FROM %s WHERE bahn_id = '%s'", tablename_cal, calibration_id);
        data_cal_soll = fetch(conn, query);
        data_cal_soll = sortrows(data_cal_soll,'timestamp');

        % Positionsdaten für Koordinatentransformation
        [trafo_rot, trafo_trans, error_metrics] = calibration(data_cal_ist,data_cal_soll, plots);

        % Bei Auswertung der Orientierung wird zusätzlich eine andere Collection benötigt
        if evaluate_orientation
            tablename_cal = ['robotervermessung.' schema '.bahn_orientation_soll'];
            query = sprintf("SELECT * FROM %s WHERE bahn_id = '%s'", tablename_cal, calibration_id);
            data_cal_soll = fetch(conn, query);
            data_cal_soll = sortrows(data_cal_soll,'timestamp');
            
            % Transformation der Quaternionen/Eulerwinkel
            q_transform = calibrateQuaternion(data_cal_ist, data_cal_soll);
            
            % Unterteilung der Bahn in ihre Segmente
            [num_segments, segment_ids, data_ist, data_soll, segments_ist, segments_soll, segments_trafo, q_transformed] = ...
            getSegmentsLT(conn, bahn_id, schema, evaluate_orientation, evaluate_velocity, trafo_rot, trafo_trans, q_transform);
        else
            q_transform = 0;
            % Unterteilung der Bahn in ihre Segmente
            [num_segments, segment_ids, data_ist, data_soll, segments_ist, segments_soll, segments_trafo, ~] = ...
            getSegments(conn, bahn_id, schema, evaluate_orientation, evaluate_velocity, trafo_rot, trafo_trans, q_transform);
        end

        clear tablename_cal data_cal_ist data_cal_soll
    end

    % Berechnung der Metriken für die einzelnen Segmente (gemeinsam für beide Transformationstypen)
    [table_euclidean_info, table_euclidean_deviation, ...
     table_lcss_info, table_lcss_deviation, ...
     table_sidtw_info, table_sidtw_deviation, ...
     table_dtw_info, table_dtw_deviation, ...
     table_dfd_info, table_dfd_deviation] = ...
     calculateMetrics(bahn_id, segment_ids, num_segments, evaluate_velocity, evaluate_orientation, segments_soll, segments_ist, segments_trafo);
    
    % Prüfen, ob Kalibrierungsbahn (kommt von Vicon/MoCap)
    is_calibration_run = false;
    if ~strcmpi(data_source, 'leica_at960')
        % Finde heraus, ob es eine Kalibrierungsbahn ist
        query = ['SELECT calibration_run FROM robotervermessung.' schema '.bahn_info ' ...
                'WHERE bahn_id = ''' bahn_id ''''];
        cal_result = fetch(conn, query);
        if ~isempty(cal_result)
            is_calibration_run = cal_result.calibration_run;
        end
    end
    
    if ~is_calibration_run % Calibration run wird nicht ausgewertet, da zu viele Daten
        % Berechnung der Metriken für die Gesamtbahn
        if strcmpi(data_source, 'leica_at960')
            % Für Lasertracker: Verwende direkt die transformierten Daten
            data_all_soll = table2array(data_position_soll(:,5:7));
            data_all_ist_trafo = positions_transformed;
        else
            % Für MoCap: Verwende die evaluateAll-Funktion
            [data_all_soll, data_all_ist_trafo] = evaluateAll(segment_ids, data_ist, data_soll, evaluate_orientation, q_transform, trafo_rot, trafo_trans);
        end
        
        [table_euclidean_all_info, table_euclidean_all_deviation, ...
         table_lcss_all_info, table_lcss_all_deviation, ...
         table_sidtw_all_info, table_sidtw_all_deviation, ...
         table_dtw_all_info, table_dtw_all_deviation, ...
         table_dfd_all_info, table_dfd_all_deviation] = ...
         calculateMetrics(bahn_id, table(string(bahn_id)), 1, evaluate_velocity, evaluate_orientation, data_all_soll, segments_ist, data_all_ist_trafo);

        % Info und Deviation Tabellen der Gesamtauswertung anhängen
        table_euclidean_info = [table_euclidean_all_info; table_euclidean_info];
        table_sidtw_info = [table_sidtw_all_info; table_sidtw_info];
        table_dtw_info = [table_dtw_all_info; table_dtw_info];
        table_dfd_info = [table_dfd_all_info; table_dfd_info];
        table_lcss_info = [table_lcss_all_info; table_lcss_info];
        
        table_sidtw_deviation = [table_sidtw_all_deviation; table_sidtw_deviation];
        table_dtw_deviation = [table_dtw_all_deviation; table_dtw_deviation];
        table_dfd_deviation = [table_dfd_all_deviation; table_dfd_deviation];
        table_euclidean_deviation = [table_euclidean_all_deviation; table_euclidean_deviation];
        table_lcss_deviation = [table_lcss_all_deviation; table_lcss_deviation];

        % Bestimme den Evaluationstyp basierend auf den Funktionsparametern
        if evaluate_velocity == false && evaluate_orientation == false
            evaluation_type = 'position';
        elseif evaluate_velocity == false && evaluate_orientation == true
            evaluation_type = 'orientation';
        elseif evaluate_velocity == true && evaluate_orientation == false
            evaluation_type = 'speed';
        end

        % Füge die evaluation-Spalte zu allen Info-Tabellen hinzu
        if ~isempty(table_euclidean_info)
            table_euclidean_info = addvars(table_euclidean_info, repelem({evaluation_type}, height(table_euclidean_info), 1), 'NewVariableNames', 'evaluation');
        end
        
        if ~isempty(table_sidtw_info)
            table_sidtw_info = addvars(table_sidtw_info, repelem({evaluation_type}, height(table_sidtw_info), 1), 'NewVariableNames', 'evaluation');
        end
        
        if ~isempty(table_dtw_info)
            table_dtw_info = addvars(table_dtw_info, repelem({evaluation_type}, height(table_dtw_info), 1), 'NewVariableNames', 'evaluation');
        end
        
        if ~isempty(table_dfd_info)
            table_dfd_info = addvars(table_dfd_info, repelem({evaluation_type}, height(table_dfd_info), 1), 'NewVariableNames', 'evaluation');
        end
        
        if ~isempty(table_lcss_info)
            table_lcss_info = addvars(table_lcss_info, repelem({evaluation_type}, height(table_lcss_info), 1), 'NewVariableNames', 'evaluation');
        end

        clear table_euclidean_all_info table_euclidean_all_deviation ...
        table_lcss_all_info table_lcss_all_deviation ...
        table_sidtw_all_info table_sidtw_all_deviation ...
        table_dtw_all_info table_dtw_all_deviation ...
        table_dfd_all_info table_dfd_all_deviation
    end
end
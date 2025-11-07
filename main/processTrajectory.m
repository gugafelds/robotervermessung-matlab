function [table_euclidean_info, table_sidtw_info, table_dtw_info, table_dfd_info, table_lcss_info, table_qad_info, table_qdtw_info, ...
          table_euclidean_deviation, table_sidtw_deviation, table_dtw_deviation, ...
          table_dfd_deviation, table_lcss_deviation, table_qad_deviation, table_qdtw_deviation, segment_ids, calibration_id] = ...
              processTrajectory(conn, bahn_id, ...
                                      evaluate_orientation, plots, ...
                                      use_euclidean, use_sidtw, use_dtw, use_dfd, use_lcss, use_qad, use_qdtw)

    schema = 'bewegungsdaten';

    query = ['SELECT source_data_ist, transformation_matrix FROM robotervermessung.' schema '.bahn_info ' ...
             'WHERE bahn_id = ''' bahn_id ''''];
    data_source_result = fetch(conn, query);
    
    if isempty(data_source_result)
        error(['Keine Quellenangabe für Bahn-ID ' bahn_id ' gefunden.']);
    end
    
    data_source = data_source_result.source_data_ist{1};
    check_transf = data_source_result.transformation_matrix{1};
    disp(['Datenquelle für Bahn-ID ' bahn_id ': ' data_source]);
    
    % Lasertracker-Transformation oder MoCap-Transformation
    if strcmpi(data_source, 'leica_at960')
        % Lasertracker-Transformation
        disp('Verwende Lasertracker-Transformation für AT960');
        
        if ~isempty(check_transf)
            % Bereits transformierte Daten sind vorhanden
            disp('Bereits transformierte Daten gefunden - lade direkt aus bahn_pose_trans');
            
            % Lade bereits transformierte Daten direkt
            query = ['SELECT * FROM robotervermessung.' schema '.bahn_pose_trans ' ...
                    'WHERE bahn_id = ''' bahn_id ''''];
            data_trans = fetch(conn, query);
            
            if isempty(data_trans)
                error(['Keine transformierten Daten für Bahn-ID ' bahn_id ' in bahn_pose_trans gefunden.']);
            end
            
            data_trans = sortrows(data_trans, 'timestamp');
            
            % Extrahiere die bereits transformierten Positionen und Quaternionen
            positions_transformed = data_trans(:, [2, 4, 5, 6]);     % Transformierte Position [x, y, z]
            q_transformed = data_trans(:, [2, 7, 8, 9, 10]);            % Transformierte Quaternion [qx, qy, qz, qw]
            
            % Setze Kalibrierungs-ID auf die verwendete transformation_matrix
            calibration_id = check_transf;
            
            disp(['Transformierte Daten geladen mit Kalibrierungs-ID: ' calibration_id]);
        
        else
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
        
        end

        query = ['SELECT segment_id, timestamp, qx_soll, qy_soll, qz_soll, qw_soll FROM robotervermessung.' schema '.bahn_orientation_soll ' ...
                    'WHERE robotervermessung.' schema '.bahn_orientation_soll.bahn_id = ''' bahn_id ''''];
        data_orientation_soll = fetch(conn, query);
        data_orientation_soll = sortrows(data_orientation_soll, 'timestamp');
        
        query = ['SELECT segment_id, timestamp, x_soll, y_soll, z_soll FROM robotervermessung.' schema '.bahn_position_soll ' ...
                'WHERE robotervermessung.' schema '.bahn_position_soll.bahn_id = ''' bahn_id ''''];
        data_position_soll = fetch(conn, query);
        data_position_soll = sortrows(data_position_soll, 'timestamp');

        if evaluate_orientation
            segments_ist = q_transformed;
            segments_soll = data_orientation_soll(:, [1, 3, 4, 5, 6]); 
        else
            segments_ist = positions_transformed;
            segments_soll = data_position_soll(:, [1, 3, 4, 5]);      
        end

        query = ['SELECT segment_id FROM robotervermessung.' schema '.bahn_events ' ...
        'WHERE robotervermessung.' schema '.bahn_events.bahn_id = ''' bahn_id ''''];
        segment_ids = fetch(conn,query);
    end

    [table_euclidean_info, table_euclidean_deviation, ...
     table_lcss_info, table_lcss_deviation, ...
     table_sidtw_info, table_sidtw_deviation, ...
     table_dtw_info, table_dtw_deviation, ...
     table_dfd_info, table_dfd_deviation, ...
     table_qad_info, table_qad_deviation, ...
     table_qdtw_info, table_qdtw_deviation] = ...
     calculateMetrics(bahn_id, segment_ids, 1, evaluate_orientation, segments_soll, segments_ist, ...
                      use_euclidean, use_sidtw, use_dtw, use_dfd, use_lcss, use_qad, use_qdtw);
    
end
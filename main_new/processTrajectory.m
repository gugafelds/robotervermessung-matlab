function [table_euclidean_info, table_sidtw_info, table_dtw_info, table_dfd_info, table_lcss_info, ...
          table_euclidean_deviation, table_sidtw_deviation, table_dtw_deviation, ...
          table_dfd_deviation, table_lcss_deviation, segment_ids, calibration_id] = ...
              processTrajectory(conn, bahn_id, evaluate_velocity, ...
                                      evaluate_orientation, plots)
        
        % Collection aus der die Daten extrahiert werden (evtl. übergeben)
        schema = 'bewegungsdaten';

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
        end

        clear tablename_cal data_cal_ist data_cal_soll
        
        % Unterteilung der Bahn in ihre Segmente
        if evaluate_orientation
            [num_segments, segment_ids, data_ist, data_soll, segments_ist, segments_soll, segments_trafo, q_transformed] = ...
            getSegments(conn, bahn_id, schema, evaluate_orientation, evaluate_velocity, trafo_rot, trafo_trans, q_transform);
        else
            q_transform = 0;
            [num_segments, segment_ids, data_ist, data_soll, segments_ist, segments_soll, segments_trafo, ~] = ...
            getSegments(conn, bahn_id, schema, evaluate_orientation, evaluate_velocity, trafo_rot, trafo_trans, q_transform);
        end

        % Berechnung der Metriken für die einzelnen Segmente
        [table_euclidean_info, table_euclidean_deviation, ...
         table_lcss_info, table_lcss_deviation, ...
         table_sidtw_info, table_sidtw_deviation, ...
         table_dtw_info, table_dtw_deviation, ...
         table_dfd_info, table_dfd_deviation] = ...
         calculateMetrics(bahn_id, segment_ids, num_segments, evaluate_velocity, evaluate_orientation, segments_soll, segments_ist, segments_trafo);
        

        if ~is_calibration_run % Calibration run wird nicht ausgewertet, da zu viele Daten
            % Berechnung der Metriken für die Gesamtbahn
            [data_all_soll, data_all_ist_trafo] = evaluateAll(segment_ids, data_ist, data_soll, evaluate_orientation, q_transform, trafo_rot, trafo_trans);
            
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
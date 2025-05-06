function [calibration_id, is_calibration_run] = findCalibrationRun(conn, bahn_id, schema)
    % FINDCALIBRATIONRUN Finds the appropriate calibration run for a given bahn_id
    %   This function searches for a calibration run from the same day that was
    %   recorded before the current measurement. If no suitable calibration run
    %   is found, it returns the current bahn_id as a fallback.
    %
    % Inputs:
    %   conn - Database connection object
    %   bahn_id - ID of the current measurement (string)
    %   schema - Database schema name (string, e.g., 'bewegungsdaten')
    %
    % Outputs:
    %   calibration_id - ID of the found calibration run or current bahn_id if none found
    %   is_calibration_run - Boolean indicating if the current bahn is itself a calibration run
    
    % Check if the current measurement is itself a calibration run
    query = sprintf(['SELECT calibration_run FROM robotervermessung.%s.bahn_info ' ...
                    'WHERE bahn_id = ''%s'''], schema, bahn_id);
    is_calibration_run = logical(fetch(conn, query).calibration_run);
    
    if is_calibration_run
        disp('Die aktuelle Bahn ist selbst eine Kalibrierungsbahn');
        calibration_id = bahn_id;
        return;
    end
    
    % Get all calibration runs
    query = sprintf(['SELECT bahn_id, recording_date FROM robotervermessung.%s.bahn_info ' ...
                    'WHERE calibration_run = true'], schema);
    cal_data = fetch(conn, query);
    
    % Get recording date of current measurement
    query = sprintf(['SELECT recording_date FROM robotervermessung.%s.bahn_info ' ...
                    'WHERE bahn_id = ''%s'''], schema, bahn_id);
    current_date = fetch(conn, query);
    current_datetime = datetime(current_date.recording_date, 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSSSSS');
    
    % Find best matching calibration
    best_time_diff = Inf;
    calibration_id = bahn_id;  % Default if no suitable calibration is found
    
    for i = 1:height(cal_data)
        cal_datetime = datetime(cal_data.recording_date(i), 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSSSSS');
        
        % Check if calibration is from the same day and before current measurement
        if dateshift(cal_datetime, 'start', 'day') == dateshift(current_datetime, 'start', 'day') && ...
           cal_datetime < current_datetime
            
            time_diff = seconds(current_datetime - cal_datetime);
            if time_diff < best_time_diff
                best_time_diff = time_diff;
                calibration_id = char(cal_data.bahn_id(i));
            end
        end
    end
    
    % Log results
    if best_time_diff < Inf
        disp(['Passende Kalibrierung gefunden: ' calibration_id]);
        disp(['Zeitlicher Abstand: ' num2str(best_time_diff/60, '%.1f') ' Minuten']);
    else
        disp('Keine passende Kalibrierung gefunden - verwende aktuelle Bahn als Referenz');
    end
end
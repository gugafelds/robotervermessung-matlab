%% Manuelle Eingaben
clear; tic

% Mode selection: 'evaluation' oder 'transformation'
mode = 'evaluation'; % Default Modus

% Falls spezifische Bahn-ID ausgewertet werden soll (höchste Priorität)
bahn_id = '';
% Falls Daten aus einem bestimmten Zeitraum ausgewertet werden sollen
record_date = '04.09.2025'; % Format: dd.mm.yyyy

% Falls mehrere Daten ausgewertet werden sollen
loop_record_date = 1; % (Priorität vor loop_all)
loop_all = 0;

% Falls Daten gelöscht und überschrieben werden sollen (Auswertung wird neu berechnet)
overwrite = 0;

% === NEU: Methodenauswahl ===
use_euclidean = 1;    % Euklidische Distanz
use_sidtw = 1;        % Scale Invariant Dynamic Time Warping
use_dtw = 0;          % Dynamic Time Warping
use_dfd = 0;          % Discrete Fréchet Distance
use_lcss = 0;         % Longest Common Subsequence

% === Falls keine Eingabe erfolgt, wird die Position ausgewertet ===
evaluate_velocity = 0;
evaluate_orientation = 0;

% Hochladen der Daten
upload_info = 1;        % Info-Tabellen hochladen
upload_deviations = 1;  % Abweichungs-Tabellen hochladen

batch_size = 3;         % Stapelgröße für Upload

% Plotten -> Muss noch vernünftig eingebunden werden (erstmal keine Funktion)
plots = 0; 

% Verbindung mit der Datenbank
conn = connectingToPostgres;

% Überprüfung ob das Datum korrekt eingegeben wurde, wenn keine Bahn-Id vorliegt
if isempty(bahn_id) && ~isempty(record_date)
    try
        record_date = datetime(record_date, 'InputFormat', 'dd.MM.yyyy');
        record_date = datestr(record_date, 'yyyymmdd'); % Konvertierung in String-Format
    catch
        error('record_date hat nicht das richtige Format.');    
    end
end

%% Modus-Auswahl und Verarbeitung

if strcmpi(mode, 'transformation')
    % Transformationsmodus
    disp('Ausführung im Transformationsmodus - nur Transformation wird durchgeführt');
    
    % Bestimme die zu transformierenden Bahn-IDs
    bahn_ids_to_transform = [];
    
    % Falls eine spezifische Bahn-ID angegeben wurde
    if ~isempty(bahn_id)
        bahn_ids_to_transform = str2double(bahn_id);
        disp(['Transformiere einzelne Bahn-ID: ' bahn_id]);
        
        % Extrahiere das Aufnahmedatum aus bahn_info für die Matrix-Suche
        query = ['SELECT recording_date FROM robotervermessung.bewegungsdaten.bahn_info WHERE bahn_id = ''' bahn_id ''''];
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
            record_date = datetime("now", 'dd.MM.yyyy');
        end
        
    % Falls ein Datum angegeben wurde
    elseif ~isempty(record_date)
        % Hole alle Bahnen des angegebenen Datums
        query = ['SELECT bahn_id FROM robotervermessung.bewegungsdaten.bahn_info WHERE record_filename LIKE ''%record_' record_date '_%'''];
        result = fetch(conn, query);
        
        if isempty(result)
            error(['Keine Bahnen für das Datum ' record_date ' gefunden!']);
        end
        
        bahn_ids_to_transform = str2double(result.bahn_id);
        disp(['Transformiere ' num2str(length(bahn_ids_to_transform)) ' Bahnen vom Datum ' record_date]);
    else
        error('Für den Transformationsmodus ist entweder eine Bahn-ID oder ein Datum erforderlich!');
    end
    
    % Batch-Verarbeitung für Transformation konfigurieren
    num_batches = ceil(length(bahn_ids_to_transform) / batch_size);
    disp(['Transformation in ', num2str(num_batches), ' Batches mit je max. ', num2str(batch_size), ' Bahnen']);
    
    try
        [offset_matrix, robot_to_instrument, matrix_file] = findTransfMatrices(record_date);
        disp(['Transformationsmatrizen für Datum ' record_date ' geladen']);
    catch ME
        warning(['Fehler beim Laden der Transformationsmatrizen: ' ME.message]);
        % Versuche das aktuelle Datum als Fallback
        try
            current_date = datetime('now', 'dd.MM.yyyy');
            [offset_matrix, robot_to_instrument] = findTransfMatrices(current_date);
            disp(['Verwende Transformationsmatrizen vom aktuellen Datum: ' current_date]);
        catch
            error('Konnte keine Transformationsmatrizen laden. Abbruch der Transformation.');
        end
    end
    
    % Durchlaufe alle Batches
    for batch_num = 1:num_batches
        batch_start_time = datetime('now');
        
        % Anfang und Ende des Stapels bestimmen
        batch_start_idx = (batch_num-1) * batch_size + 1;
        batch_end_idx = min(batch_num * batch_size, length(bahn_ids_to_transform));
        
        % Bahn-IDs im aktuellen Stapel
        current_batch_ids = bahn_ids_to_transform(batch_start_idx:batch_end_idx);
        
        disp(['Verarbeite Batch ', num2str(batch_num), ' von ', num2str(num_batches), ...
              ' (Bahnen ', num2str(batch_start_idx), '-', num2str(batch_end_idx), ' von ', ...
              num2str(length(bahn_ids_to_transform)), ')']);
        
        % Sammler für den aktuellen Batch
        batch_ids = cell(length(current_batch_ids), 1);
        batch_calibration_ids = cell(length(current_batch_ids), 1);
        batch_data_ist_list = cell(length(current_batch_ids), 1);
        batch_data_ist_trafo_list = cell(length(current_batch_ids), 1);
        batch_q_transformed_list = cell(length(current_batch_ids), 1);
        batch_valid_count = 0;
        
        % Durchlaufe alle zu transformierenden Bahnen im aktuellen Batch
        for i = 1:length(current_batch_ids)
            current_bahn_id = num2str(current_batch_ids(i));
            disp(['Transformiere Bahn ' num2str(i) ' von ' num2str(length(current_batch_ids)) ' im aktuellen Batch: Bahn-ID ' current_bahn_id]);
            
            try
                % Prüfe die Datenquelle (Lasertracker AT960 oder Motion Capture Vicon)
                query = ['SELECT source_data_ist FROM robotervermessung.bewegungsdaten.bahn_info ' ...
                         'WHERE bahn_id = ''' current_bahn_id ''''];
                data_source_result = fetch(conn, query);
                
                if isempty(data_source_result)
                    warning(['Keine Quellenangabe für Bahn-ID ' current_bahn_id ' gefunden. Überspringe diese Bahn.']);
                    continue;
                end
                
                data_source = data_source_result.source_data_ist{1};
                disp(['Datenquelle für Bahn-ID ' current_bahn_id ': ' data_source]);
                
                % Auslesen der zu transformierenden Daten
                query = ['SELECT * FROM robotervermessung.bewegungsdaten.bahn_pose_ist ' ...
                        'WHERE robotervermessung.bewegungsdaten.bahn_pose_ist.bahn_id = ''' current_bahn_id ''''];
                data_ist = fetch(conn, query);
                data_ist = sortrows(data_ist, 'timestamp');
                
                query = ['SELECT * FROM robotervermessung.bewegungsdaten.bahn_orientation_soll ' ...
                        'WHERE robotervermessung.bewegungsdaten.bahn_orientation_soll.bahn_id = ''' current_bahn_id ''''];
                data_orientation_soll = fetch(conn, query);
                data_orientation_soll = sortrows(data_orientation_soll, 'timestamp');
                
                query = ['SELECT * FROM robotervermessung.bewegungsdaten.bahn_position_soll ' ...
                        'WHERE robotervermessung.bewegungsdaten.bahn_position_soll.bahn_id = ''' current_bahn_id ''''];
                data_position_soll = fetch(conn, query);
                data_position_soll = sortrows(data_position_soll, 'timestamp');
                
                % Transformationslogik basierend auf der Datenquelle
                if strcmpi(data_source, 'leica_at960') % Lasertracker-Transformation
                    % Überprüfen, ob die Matrizen geladen wurden
                    if isempty(offset_matrix) || isempty(robot_to_instrument)
                        error('Transformationsmatrizen für Lasertracker wurden nicht geladen');
                    end
                    
                    % Daten extrahieren
                    position_ist = table2array(data_ist(:,5:7));     % IST-Position [x, y, z]
                    quaternion_ist = table2array(data_ist(:,8:11));  % IST-Quaternion [qx, qy, qz, qw]
                    
                    % Initialisierung der Ergebnisarrays
                    num_points = size(position_ist, 1);
                    positions_transformed = zeros(num_points, 3);
                    quaternions_transformed = zeros(num_points, 4);
                    
                    % Lasertracker-Transformation
                    for j = 1:num_points
                        % Umformatierung des Quaternions (von [qx, qy, qz, qw] zu [qw, qx, qy, qz])
                        quaternion_sa = LT2SAQuaternionDirect(quaternion_ist(j,:));
                        quaternion_reordered = [quaternion_sa(4), quaternion_sa(1), quaternion_sa(2), quaternion_sa(3)];
                        
                        % 1. Anwenden des Offsets auf den Punkt und die Orientierung
                        [point_offset, ~, quaternion_offset] = applyOffset(position_ist(j,:), [], quaternion_reordered, offset_matrix);
                        
                        % 2. Transformation in Roboterkoordinaten
                        positions_transformed(j,:) = transformToRobotCoords(point_offset, robot_to_instrument);
                        quaternions_transformed(j,:) = transformQuaternionToRobotCoords(quaternion_offset, robot_to_instrument);
                        
                        % 3. Berechnung der Euler-Winkel für die Robotersteuerung (XYZ Fixed Angles)
                        orientations_transformed(j,:) = quaternionToEulerXYZFixed(quaternions_transformed(j,:));
                    end
                                                    
                    % Calibration ID für den Lasertracker
                    calibration_id = matrix_file;
                    
                else % Vicon/MoCap-Transformation
                    % Suche nach Kalibrierungslauf
                    [calibration_id, is_calibration_run] = findCalibrationRun(conn, current_bahn_id, schema);
                    
                    % Extrahieren der Kalibrierungs-Daten für die Position
                    tablename_cal = ['robotervermessung.bewegungsdaten.bahn_pose_ist'];
                    opts_cal = databaseImportOptions(conn, tablename_cal);
                    opts_cal.RowFilter = opts_cal.RowFilter.bahn_id == calibration_id;
                    data_cal_ist = sqlread(conn, tablename_cal, opts_cal);
                    data_cal_ist = sortrows(data_cal_ist, 'timestamp');
                    
                    tablename_cal = ['robotervermessung.bewegungsdaten.bahn_events'];
                    opts_cal = databaseImportOptions(conn, tablename_cal);
                    opts_cal.RowFilter = opts_cal.RowFilter.bahn_id == calibration_id;
                    data_cal_soll = sqlread(conn, tablename_cal, opts_cal);
                    data_cal_soll = sortrows(data_cal_soll, 'timestamp');
                    
                    % Positionsdaten für Koordinatentransformation
                    [trafo_rot, trafo_trans, ~] = calibration(data_cal_ist, data_cal_soll, plots);
                    
                    % Extrahieren der Kalibrierungs-Daten für die Orientierung
                    tablename_cal = ['robotervermessung.bewegungsdaten.bahn_orientation_soll'];
                    opts_cal = databaseImportOptions(conn, tablename_cal);
                    opts_cal.RowFilter = opts_cal.RowFilter.bahn_id == calibration_id;
                    data_cal_soll = sqlread(conn, tablename_cal, opts_cal);
                    data_cal_soll = sortrows(data_cal_soll, 'timestamp');
                    
                    % Berechnung der relativen Rotationsmatrix für die Orientierung 
                    q_transform = calibrateQuaternion(data_cal_ist, data_cal_soll);
                    
                    clear data_cal opts_cal tablename_cal
                    
                    % Transformation durchführen
                    position_ist = table2array(data_ist(:,5:7));
                    positions_transformed = coord_transformation(position_ist, trafo_rot, trafo_trans);
                    q_transformed = transformQuaternion(data_ist, data_orientation_soll, q_transform, trafo_rot);
                end
                
                % Daten in die Batch-Arrays einfügen
                batch_valid_count = batch_valid_count + 1;
                batch_ids{batch_valid_count} = current_bahn_id;
                batch_calibration_ids{batch_valid_count} = calibration_id;
                batch_data_ist_list{batch_valid_count} = data_ist;
                batch_data_ist_trafo_list{batch_valid_count} = positions_transformed;
                batch_q_transformed_list{batch_valid_count} = quaternions_transformed;
                
                disp(['Transformation für Bahn-ID ' current_bahn_id ' erfolgreich abgeschlossen.']);
                
            catch ME
                warning(['Fehler bei der Transformation von Bahn-ID ' current_bahn_id ': ' ME.message]);
                continue; % Mit der nächsten Bahn fortfahren
            end
        end
        
        % Entferne leere Einträge, falls Fehler aufgetreten sind
        if batch_valid_count < length(current_batch_ids)
            batch_ids = batch_ids(1:batch_valid_count);
            batch_calibration_ids = batch_calibration_ids(1:batch_valid_count);
            batch_data_ist_list = batch_data_ist_list(1:batch_valid_count);
            batch_data_ist_trafo_list = batch_data_ist_trafo_list(1:batch_valid_count);
            batch_q_transformed_list = batch_q_transformed_list(1:batch_valid_count);
        end
        
        % Batch-Upload durchführen, wenn Daten vorhanden sind
        if ~isempty(batch_ids)
            upload_start_time = datetime('now');
            disp(['Starte Batch-Upload für ' num2str(length(batch_ids)) ' transformierte Bahnen...']);
            
            try
                uploadTransformedData(conn, batch_ids, batch_calibration_ids, batch_data_ist_list, batch_data_ist_trafo_list, batch_q_transformed_list);
                
                upload_end_time = datetime('now');
                upload_duration = seconds(upload_end_time - upload_start_time);
                disp(['Batch-Upload abgeschlossen in ' num2str(upload_duration) ' Sekunden!']);
            catch ME
                warning(['Fehler beim Batch-Upload: ' ME.message]);
            end
        else
            disp('Keine Daten zum Hochladen in diesem Batch.');
        end
        
        batch_end_time = datetime('now');
        batch_duration = seconds(batch_end_time - batch_start_time);
        disp(['Batch ' num2str(batch_num) ' abgeschlossen in ' num2str(batch_duration) ' Sekunden']);
    end
    
    disp(['Alle ' num2str(length(bahn_ids_to_transform)) ' Transformationen abgeschlossen.']);
    
elseif strcmpi(mode, 'evaluation')
    % Originaler Auswertungs-Code bleibt unverändert
    disp('Ausführung im Auswertungsmodus');

    if loop_all || loop_record_date || ~isempty(bahn_id)
        
        % Prüft anhand der mit SIDTW ausgewerteten Daten welche Bahnen bereits ausgewertet wurden
        [bahn_ids_all, bahn_ids_evaluated, evaluation_type] = getBahnIds(conn,evaluate_orientation,evaluate_velocity);
        
        % WENN DIESER BLOCK NICHT AKTIV WIRD WERDEN ALLE DATEN AUSGEWERTET
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Wenn Daten eines bestimmten Datums ausgewertet werden sollen
        if loop_record_date && isempty(bahn_id)
            % Info Tabellen des relevanten Tages extrahieren (ggf. ohne Kalibrierungsdateien)
            query = 'SELECT * FROM robotervermessung.bewegungsdaten.bahn_info WHERE calibration_run = false';
            all_info = fetch(conn,query);
            % Datum aus dem String extrahieren
            all_dates = extractBetween(all_info.record_filename, "record_", "_");
            bahn_info = all_info(all_dates == record_date, :);
            bahn_ids_all = sort(double(bahn_info.bahn_id),'ascend');
            bahn_ids_all_str = string(bahn_ids_all);
        
            if isempty(bahn_ids_all)
                error('An diesem Datum wurden keine Daten aufgezeichnet!')
            end
        
            clear all_dates all_info bahn_info
        
        % Wenn Daten einer bestimmte Bahn ausgewertet werden soll 
        elseif ~isempty(bahn_id) && length(bahn_id) == 10  
            bahn_ids_all = double(string(bahn_id));
        elseif ~isempty(bahn_id) && length(bahn_id) ~= 10  
            error('Ungültige Eingaben! Bahn-Id ist ungültig! ')
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Wenn Daten überschrieben werden sollen, werden die bereits ausgewerteten Daten annuliert
        if overwrite
            bahn_ids_evaluated = zeros(1);
        end
        
        % Bestimmung der zu verarbeitenden Bahnen
        bahn_ids_unprocessed = bahn_ids_all(~ismember(bahn_ids_all, bahn_ids_evaluated));
        bahn_ids_unprocessed_str = string(bahn_ids_unprocessed);
        
        bahn_ids_processed = bahn_ids_all(ismember(bahn_ids_all, bahn_ids_evaluated));
        disp(string(length(bahn_ids_processed)) +' Bahnen bereits ausgewertet.');
        disp('Bahn-Id: ' + string(bahn_ids_processed));
        
        % %%%%%% Später wieder aktivieren! %%%%%%
        % disp('Anzahl der zu verarbeitenden Bahnen: ' + string(length(bahn_ids_unprocessed)) + '. Fahre fort um diese auszuwerten!')
        % % keyboard;
        
        % Anzahl der Batches 
        num_batches = ceil(length(bahn_ids_unprocessed) / batch_size);
        disp(['Verarbeitung in ', num2str(num_batches), ' Batches mit je max. ', num2str(batch_size), ' Bahnen']);
        
        % Tracking der verarbeiteten IDs
        all_processed_ids = [];
        
        % Batch-Verarbeitung mit Upload nach jedem Batch
        for batch_num = 1:num_batches
    
            % Zeiterfassung
            batch_start_time = datetime('now');
            % Anfang und Ende des Stapels bestimmen
            batch_start_idx = (batch_num-1) * batch_size + 1;
            batch_end_idx = min(batch_num * batch_size, length(bahn_ids_unprocessed));
            
            % Bahn-Ids im aktuellen Stapel
            current_batch_ids = bahn_ids_unprocessed(batch_start_idx:batch_end_idx);
        
            disp(['Verarbeite Batch ', num2str(batch_num), ' von ', num2str(num_batches), ...
                  ' (Bahnen ', num2str(batch_start_idx), '-', num2str(batch_end_idx), ' von ', ...
                  num2str(length(bahn_ids_unprocessed)), ')']);
            
            % Verarbeite den Batch (nur Berechnungen, kein Upload)
            [batch_processed_ids, batch_euclidean_info, batch_sidtw_info, batch_dtw_info, batch_dfd_info, batch_lcss_info, ...
             batch_euclidean_deviations, batch_sidtw_deviations, batch_dtw_deviations, batch_dfd_deviations, batch_lcss_deviations] = ...
                processBatch(conn, current_batch_ids, evaluate_velocity, evaluate_orientation, plots, ...
                 use_euclidean, use_sidtw, use_dtw, use_dfd, use_lcss);
            
            all_processed_ids = [all_processed_ids; batch_processed_ids];
            all_processed_ids_str = string(all_processed_ids);
    
            % Structs einzeln aufbauen
            batch_info = struct();
            batch_deviations = struct();
            
            % Nur aktivierte Methoden hinzufügen
            if use_sidtw
                batch_info.sidtw = batch_sidtw_info;
                batch_deviations.sidtw = batch_sidtw_deviations;
            end
            
            if use_dtw
                batch_info.dtw = batch_dtw_info;
                batch_deviations.dtw = batch_dtw_deviations;
            end
            
            if use_dfd
                batch_info.dfd = batch_dfd_info;
                batch_deviations.dfd = batch_dfd_deviations;
            end
            
            if use_euclidean
                batch_info.euclidean = batch_euclidean_info;
                batch_deviations.euclidean = batch_euclidean_deviations;
            end
            
            if use_lcss
                batch_info.lcss = batch_lcss_info;
                batch_deviations.lcss = batch_lcss_deviations;
            end
                        
            % Upload der aktuellen Batch-Daten
            if (upload_deviations || upload_info) && ~isempty(batch_processed_ids)
        
                upload_start_time = datetime('now');
                disp(['Starte Upload für Batch ', num2str(batch_num), ' mit ', num2str(length(batch_processed_ids)), ' Bahnen...']);
        
                % Upload der Daten 
                uploadBatchData(conn, batch_processed_ids, evaluation_type, ...
                evaluate_velocity, evaluate_orientation, ...
                upload_info, upload_deviations, ...
                batch_info, batch_deviations, ...
                use_euclidean, use_sidtw, use_dtw, use_dfd, use_lcss);
        
        
                upload_end_time = datetime('now');
                upload_duration = seconds(upload_end_time - upload_start_time);
                disp(['Batch-Upload abgeschlossen in ', num2str(upload_duration), ' Sekunden!']);
            end
            
            % Speicher freigeben
            clear batch_euclidean_info batch_sidtw_info batch_dtw_info batch_dfd_info batch_lcss_info;
            clear batch_euclidean_deviations batch_sidtw_deviations batch_dtw_deviations batch_dfd_deviations batch_lcss_deviations;
            
            batch_end_time = datetime('now');
            batch_duration = seconds(batch_end_time - batch_start_time);
            disp(['Batch ', num2str(batch_num), ' abgeschlossen in ', num2str(batch_duration), ' Sekunden']);
        end
        
        total_duration = toc;
        disp(['Gesamte Verarbeitung abgeschlossen in ', num2str(total_duration), ' Sekunden']);
        disp(['Durchschnittliche Zeit pro Bahn: ', num2str(total_duration/length(all_processed_ids)), ' Sekunden']);
    
        clear batch_end_time batch_end_idx batch_duration batch_start_idx batch_start_time batch_num batch_processed_ids ...
            num_batches query current_batch_ids 
    
    % Wenn keine Option ausgewählt wurde
    else
        warning("Es wurden keine Daten für die Auswertung ausgewählt !")
    end
    
else
    error(['Ungültiger Modus: ' mode '. Wählen Sie entweder ''evaluation'' oder ''transformation''.']);
end


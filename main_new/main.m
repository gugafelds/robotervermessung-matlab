%% Manuelle Eingaben
clear; tic

% Falls spezifische Bahn-ID ausgewertet werden soll (höchste Priorität)
bahn_id = '';
% Falls Daten aus einem bestimmten Zeitraum ausgewertet werden sollen
record_date = '15.07.2024'; % Format: dd.mm.yyyy

% Falls mehrere Daten ausgewertet werden sollen
loop_record_date = 1; % (Priorität vor loop_all)
loop_all = 0;

% Falls Daten gelöscht und überschrieben werden sollen (Auswertung wird neu berechnet)
overwrite = 1;

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

%% Auswertung

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
    
    % %%%%%%%% Später entfernen ! 
    % num_batches = 2;
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
            processBatch(conn, current_batch_ids, evaluate_velocity, evaluate_orientation, plots);
        
        all_processed_ids = [all_processed_ids; batch_processed_ids];
        all_processed_ids_str = string(all_processed_ids);

        % Strukts anlegen für Übersicht
        batch_info = struct( ...
            'sidtw', batch_sidtw_info, ...
            'dtw', batch_dtw_info, ...
            'dfd', batch_dfd_info, ...
            'euclidean', batch_euclidean_info, ...
            'lcss', batch_lcss_info);
    
        batch_deviations = struct( ...
            'sidtw', batch_sidtw_deviations, ...
            'dtw', batch_dtw_deviations, ...
            'dfd', batch_dfd_deviations, ...
            'euclidean', batch_euclidean_deviations, ...
            'lcss', batch_lcss_deviations);
        
        % Upload der aktuellen Batch-Daten
        if (upload_deviations || upload_info) && ~isempty(batch_processed_ids)
    
            upload_start_time = datetime('now');
            disp(['Starte Upload für Batch ', num2str(batch_num), ' mit ', num2str(length(batch_processed_ids)), ' Bahnen...']);
    
            % Upload der Daten 
            uploadBatchData(conn, batch_processed_ids, evaluation_type, ...
                evaluate_velocity, evaluate_orientation, ...
                upload_info, upload_deviations, ...
                batch_info, batch_deviations);
    
    
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



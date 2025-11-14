%% Manuelle Eingaben
clear; tic

addpath(genpath('lasertracker'))
addpath(genpath('main'))
addpath(genpath('methods'))

% Mode selection
mode = 'automatic'; % 'manual' oder 'automatic'

% Manual mode parameters
bahn_id = '';
record_date = '31.10.2025'; % Format: dd.mm.yyyy
loop_record_date = 1;
loop_all = 0;
overwrite = 0;

% Methodenauswahl
use_euclidean = 1;
use_sidtw = 1;
use_dtw = 0;
use_dfd = 0;
use_lcss = 0;
use_qad = 0;
use_qdtw = 0;

evaluate_orientation = 0;

% Upload-Einstellungen
upload_info = 1;
upload_deviations = 1;
batch_size = 10;

plots = 0;

% Datenbankverbindung
conn = connectingToPostgres;

% Datumsvalidierung
if isempty(bahn_id) && ~isempty(record_date)
    try
        record_date = datetime(record_date, 'InputFormat', 'dd.MM.yyyy');
        record_date = datestr(record_date, 'yyyymmdd');
    catch
        error('record_date hat nicht das richtige Format.');    
    end
end

%% Modus-Auswahl und Bestimmung der zu verarbeitenden Bahnen

bahn_ids_to_process = [];

if strcmpi(mode, 'manual')
    disp('Ausführung im manuellen Modus');
    
    if ~(loop_all || loop_record_date || ~isempty(bahn_id))
        warning("Es wurden keine Daten für die Auswertung ausgewählt!");
        return;
    end
    
    % Bereits ausgewertete Bahnen identifizieren
    [bahn_ids_all, bahn_ids_evaluated] = getBahnIds(conn, evaluate_orientation);
    
    if loop_record_date && isempty(bahn_id)
        % Bahnen eines bestimmten Datums
        query = 'SELECT * FROM robotervermessung.bewegungsdaten.bahn_info WHERE calibration_run = false';
        all_info = fetch(conn, query);
        all_dates = extractBetween(all_info.record_filename, "record_", "_");
        bahn_info = all_info(all_dates == record_date, :);
        bahn_ids_all = sort(double(bahn_info.bahn_id), 'ascend');
        
        if isempty(bahn_ids_all)
            error('An diesem Datum wurden keine Daten aufgezeichnet!');
        end
    elseif ~isempty(bahn_id)
        if length(bahn_id) ~= 10
            error('Ungültige Bahn-Id!');
        end
        bahn_ids_all = double(string(bahn_id));
    end
    
    % Bestimme zu verarbeitende Bahnen
    if overwrite
        bahn_ids_to_process = bahn_ids_all;
    else
        bahn_ids_to_process = bahn_ids_all(~ismember(bahn_ids_all, bahn_ids_evaluated));
        bahn_ids_processed = bahn_ids_all(ismember(bahn_ids_all, bahn_ids_evaluated));
        disp(string(length(bahn_ids_processed)) + ' Bahnen bereits ausgewertet.');
    end
    
elseif strcmpi(mode, 'automatic')
    disp('Ausführung im automatischen Modus');
    
    % Bestimme welche Tabelle für die Prüfung verwendet wird
    check_table = 'info_sidtw';
    if evaluate_orientation == 1
        check_table = 'info_qad';
    end
    
    % Hole alle noch nicht ausgewerteten Bahnen
    query_missing = sprintf(strjoin([
        "SELECT bahn_id"
        "FROM robotervermessung.bewegungsdaten.bahn_info"
        "WHERE source_data_ist = 'leica_at960'"
        "AND bahn_id NOT IN ("
        "  SELECT DISTINCT bahn_id FROM robotervermessung.auswertung.%s"
        ");"
    ], " "), check_table);
    
    bahn_ids_missing = fetch(conn, query_missing);
    bahn_ids_to_process = sort(double(bahn_ids_missing.bahn_id), 'ascend');
    
    fprintf('\n%d Bahnen müssen ausgewertet werden.\n', length(bahn_ids_to_process));
    
    if isempty(bahn_ids_to_process)
        disp('Keine Bahnen mehr zum Auswerten ✅');
        return;
    end
    
    % Sicherheitsabfrage
    user_confirm = input('Mit der Auswertung fortfahren? (j/n): ', 's');
    if ~strcmpi(user_confirm, 'j')
        disp('Vorgang abgebrochen ❌');
        return;
    end
else
    error(['Ungültiger Modus: ' mode]);
end

%% Batch-Verarbeitung (gemeinsam für beide Modi)

if isempty(bahn_ids_to_process)
    disp('Keine Bahnen zu verarbeiten.');
    return;
end

num_batches = ceil(length(bahn_ids_to_process) / batch_size);
disp(['Verarbeitung in ', num2str(num_batches), ' Batches mit je max. ', num2str(batch_size), ' Bahnen']);

all_processed_ids = [];

for batch_num = 1:num_batches
    batch_start_time = datetime('now');
    
    % Batch-Indizes bestimmen
    batch_start_idx = (batch_num-1) * batch_size + 1;
    batch_end_idx = min(batch_num * batch_size, length(bahn_ids_to_process));
    current_batch_ids = bahn_ids_to_process(batch_start_idx:batch_end_idx);
    
    disp(['Verarbeite Batch ', num2str(batch_num), ' von ', num2str(num_batches), ...
          ' (Bahnen ', num2str(batch_start_idx), '-', num2str(batch_end_idx), ')']);
    
    % Batch verarbeiten
    [batch_processed_ids, batch_euclidean_info, batch_sidtw_info, ...
     batch_dtw_info, batch_dfd_info, batch_lcss_info, ...
     batch_qad_info, batch_qdtw_info, ...
     batch_euclidean_deviations, batch_sidtw_deviations, ...
     batch_dtw_deviations, batch_dfd_deviations, batch_lcss_deviations, ...
     batch_qad_deviations, batch_qdtw_deviations] = ...
        processBatch(conn, current_batch_ids, evaluate_orientation, plots, ...
         use_euclidean, use_sidtw, use_dtw, use_dfd, use_lcss, use_qad, use_qdtw);
    
    all_processed_ids = [all_processed_ids; batch_processed_ids];
    
    % Batch-Structs aufbauen
    batch_info = struct();
    batch_deviations = struct();
    
    if use_euclidean
        batch_info.euclidean = batch_euclidean_info;
        batch_deviations.euclidean = batch_euclidean_deviations;
    end
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
    if use_lcss
        batch_info.lcss = batch_lcss_info;
        batch_deviations.lcss = batch_lcss_deviations;
    end
    if use_qad
        batch_info.qad = batch_qad_info;
        batch_deviations.qad = batch_qad_deviations;
    end
    if use_qdtw
        batch_info.qdtw = batch_qdtw_info;
        batch_deviations.qdtw = batch_qdtw_deviations;
    end
    
    % Upload
    if (upload_deviations || upload_info) && ~isempty(batch_processed_ids)
        upload_start_time = datetime('now');
        disp(['Starte Upload für Batch ', num2str(batch_num), '...']);
        
        uploadBatchData(conn, batch_processed_ids, evaluate_orientation, ...
            upload_info, upload_deviations, batch_info, batch_deviations, ...
            use_euclidean, use_sidtw, use_dtw, use_dfd, use_lcss, use_qad, use_qdtw);
        
        upload_duration = seconds(datetime('now') - upload_start_time);
        disp(['Upload abgeschlossen in ', num2str(upload_duration), ' Sekunden']);
    end
    
    % Cleanup
    clear batch_*_info batch_*_deviations
    
    batch_duration = seconds(datetime('now') - batch_start_time);
    disp(['Batch ', num2str(batch_num), ' abgeschlossen in ', num2str(batch_duration), ' Sekunden']);
end

total_duration = toc;
disp(['Gesamtverarbeitung: ', num2str(total_duration), ' Sekunden']);
disp(['Durchschnitt pro Bahn: ', num2str(total_duration/length(all_processed_ids)), ' Sekunden']);
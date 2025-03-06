%% Manuelle Eingaben
clear; tic

% Falls spezifische Bahn-ID ausgewertet werden soll (höchste Priorität)
bahn_id = '1719927316';
bahn_id = '';
% Falls Daten aus einem bestimmten Zeitraum ausgewertet werden sollen
record_date = '12.02.2025'; % Format: dd.mm.yyyy
% record_date = '02.07.2024';

% Falls mehrere Daten ausgewertet werden sollen
loop_record_date = 1; % (Priorität vor loop_all)
loop_all = 1;
% Falls Daten gelöscht und überschrieben werden sollen (Auswertung wird neu berechnet)
overwrite = 1;

% Berechnung der Metriken für die Geschwindikeitsabweichungen
evaluate_velocity = 0;
% Berechnung der Metriken für die Orientierungsabweichungen
evaluate_orientation = 0;

% ! ! ! Falls keine Eingabe erfolgt, wird die Position ausgewertet ! ! !

% Daten hochladen
%%%%
upload = true;
%%%%
upload_info = true;        % Info-Tabellen hochladen
upload_deviations = true;  % Abweichungs-Tabellen hochladen

% Plotten 
plots = false; 

% Verbindung mit der Datenbank
conn = connectingToPostgres;

% Überprüfung ob das Datum korrekt eingeben wurde wenn keine Bahn-Id vorliegt
if isempty(bahn_id) &&  ~isempty(record_date)
    try
        dt = datetime(record_date, 'InputFormat', 'dd.MM.yyyy');
        isValid = true;
    catch
        isValid = false;
    end
    if  isValid
        record_date = datetime(record_date, 'InputFormat', 'dd.MM.yyyy');
        record_date = datestr(record_date, 'yyyymmdd');
    else
        error('record_date hat nicht das richtige Format.');    
    end
end
clear isValid dt

%% Auswertung

if loop_all || loop_record_date || ~isempty(bahn_id)

% Zähgler für Anzahl der Schleifendurchläufe
counter = 0;

% Prüft anhand der mit SIDTW ausgewerteten Daten welche Bahnen bereits ausgewertet wurden
[bahn_ids_all, bahn_ids_evaluated] = getBahnIds(conn,evaluate_orientation,evaluate_velocity);

% WENN DIESER BLOCK NICHT AKTIV WIRD WERDEN ALLE DATEN AUSGEWERTET
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Wenn Daten eines bestimmten Datums ausgewertet werden sollen
if loop_record_date && isempty(bahn_id)
    % Info Tabellen des relevanten Tages extrahieren
    query = sprintf('SELECT * FROM robotervermessung.bewegungsdaten.bahn_info');
    all_info = fetch(conn,query);
    % Datum aus dem String extrahieren
    all_dates = extractBetween(all_info.record_filename, "record_", "_");
    bahn_info = all_info(all_dates == record_date, :);
    bahn_ids_all = double(bahn_info.bahn_id);
    bahn_ids_all_str = bahn_info.bahn_id;

    if isempty(bahn_ids_all)
        error('An diesem Datum wurden keine Daten aufgezeichnet!')
    end

    clear all_dates all_info bahn_info

% Wenn Daten einer bestimmte Bahn ausgewertet werden soll 
elseif ~isempty(bahn_id) && length(bahn_id) == 10  
    bahn_ids_all = double(string(bahn_id));
else
    error('Ungültige Eingaben! Bahn-Id ist ungültig! ') 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Wenn Daten überschrieben werden sollen, werden die bereits ausgewerteten Daten auf 0 gesetzt
if overwrite
    bahn_ids_evaluated = zeros(1);
end

%%%%%% Später wieder aktivieren! %%%%%%
% disp('Es liegen ' + string(length(bahn_ids_all)) + ' Datensätze vor. Fahre fort um diese auszuwerten!')
% keyboard;

% Collection aus der die Daten extrahiert werden
schema = 'bewegungsdaten';

% Schleife
% for j = 1:1:height(bahn_ids_all)
for j = 1:1:1

    bahn_id = convertStringsToChars(string(bahn_ids_all(j)));

    % Wird ausgeführt falls die Bahn noch nicht ausgewertet wurde oder überschrieben werden soll
    if ~ismember(bahn_ids_all(j),bahn_ids_evaluated)

        disp('Datei: '+string(j)+' Auswertung der Datei: ' + string(bahn_id))

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
        calibration(data_cal_ist,data_cal_soll, plots)

        % Bei Auswertung der Orientierung wird zusätzlich eine andere Collection benötigt
        if evaluate_orientation
            tablename_cal = ['robotervermessung.' schema '.bahn_orientation_soll'];
            query = sprintf("SELECT * FROM %s WHERE bahn_id = '%s'", tablename_cal, calibration_id);
            data_cal_soll = fetch(conn, query);
            data_cal_soll = sortrows(data_cal_soll,'timestamp');
            
            % Transformation der Quaternionen/Eulerwinkel
            calibrateQuaternion(data_cal_ist, data_cal_soll);
        end

        clear tablename_cal data_cal_ist data_cal_soll
        
        % Unterteilung der Bahn in ihre Segmente
        if evaluate_orientation 
            getSegments(conn, bahn_id, schema, evaluate_orientation, evaluate_velocity, trafo_rot, trafo_trans, q_transform)
        else
            q_transform = 0;
            getSegments(conn, bahn_id, schema, evaluate_orientation, evaluate_velocity, trafo_rot, trafo_trans, q_transform)
        end

        % Berechnung der Metriken
        [table_euclidean_info, table_euclidean_deviation, ...
         table_lcss_info, table_lcss_deviation, ...
         table_sidtw_info, table_sidtw_deviation, ...
         table_dtw_info, table_dtw_deviation, ...
         table_dfd_info, table_dfd_deviation] = ...
         evaluationMetrics(bahn_id, segment_ids, num_segments, evaluate_velocity, evaluate_orientation, segments_soll, segments_ist, segments_trafo);

    else
        disp('Datei: '+string(j)+" mit der Bahn-ID "+ string(bahn_id+ " lag bereits vor!"))
    end

    % Hochzählen des Counters
    counter = counter +1;
end

% Wenn keine Option ausgewählt wurde
else
    warning("Es wurden keine Daten für die Auswertung ausgewählt !")
end


toc 





%% Plots

% plots = 1;
if plots
    % Farben
    c1 = [0 0.4470 0.7410];
    c2 = [0.8500 0.3250 0.0980];
    c3 = [0.9290 0.6940 0.1250];
    c4 = [0.4940 0.1840 0.5560];
    c5 = [0.4660 0.6740 0.1880];
    c6 = [0.3010 0.7450 0.9330];
    c7 = [0.6350 0.0780 0.1840];
    ist = table2array(data_ist(:,5:7));
    soll = table2array(data_soll(:,5:7));
    
    % Koordinatentrafo für alle Daten 
    coordTransformation(ist,trafo_rot, trafo_trans)
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

    % Segmnteweise plotten
    f0 = figure('Color','white','Name','Soll und Istbahn (Segmentweise)');
    f0.Position(3:4) = [1520 840];
    hold on 
    for i = 1:1:num_segments
      plot3(segments_soll.x_soll{i,1},segments_soll.y_soll{i,1},segments_soll.z_soll{i,1},Color=c1,LineWidth=1.5)
      % plot3(segments_ist.x_ist{i,1},segments_ist.y_ist{i,1},segments_ist.z_ist{i,1},Color=c2,LineWidth=1.5)
      plot3(segments_trafo.x_ist{i,1},segments_trafo.y_ist{i,1},segments_trafo.z_ist{i,1},Color=c2,LineWidth=1.5)
    end
    xlabel('x','FontWeight','bold');
    ylabel('y','FontWeight','bold');
    zlabel('z','FontWeight','bold','Rotation',0);
    legend('Sollbahn (ABB)','Istbahn (VICON)')
    grid on 
    view(3)
end

%%


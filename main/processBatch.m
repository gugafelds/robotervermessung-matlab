function [processed_ids, batch_euclidean_info, batch_sidtw_info, batch_dtw_info, batch_dfd_info, batch_lcss_info, ...
    batch_euclidean_deviations, batch_sidtw_deviations, batch_dtw_deviations, batch_dfd_deviations, batch_lcss_deviations] = ...
    processBatch(conn, batch_ids, evaluate_velocity, evaluate_orientation, plots, ...
                 use_euclidean, use_sidtw, use_dtw, use_dfd, use_lcss)

% Initialisiere Batch-Sammlungen
all_calibration_ids = [];
all_segment_ids = [];
processed_ids = [];

% Tabellen für jede Metrik initialisieren
batch_euclidean_info = table();
batch_sidtw_info = table();
batch_dtw_info = table();
batch_dfd_info = table();
batch_lcss_info = table();

batch_euclidean_deviations = {};
batch_sidtw_deviations = {};
batch_dtw_deviations = {};
batch_dfd_deviations = {};
batch_lcss_deviations = {};

% Schleife über alle Bahnen im Batch
for i = 1:length(batch_ids)
    bahn_id = num2str(batch_ids(i));
    disp(['Verarbeite Bahn ', num2str(i), ' von ', num2str(length(batch_ids)), ': ', bahn_id]);
    
    try
        % NEUE: Übergabe der Methodenauswahl an processTrajectoryLT
        [table_euclidean_info, table_sidtw_info, table_dtw_info, table_dfd_info, table_lcss_info, ...
         table_euclidean_deviation, table_sidtw_deviation, table_dtw_deviation, ...
         table_dfd_deviation, table_lcss_deviation, segment_ids, calibration_id] = ...
            processTrajectoryLT(conn, bahn_id, evaluate_velocity, evaluate_orientation, plots, ...
                               use_euclidean, use_sidtw, use_dtw, use_dfd, use_lcss);
        
        % Sammle die Daten für den Batch-Upload
        processed_ids = [processed_ids; str2double(bahn_id)];
        all_calibration_ids = [all_calibration_ids; calibration_id];
        all_segment_ids = [all_segment_ids; segment_ids];
        
        % Sammle Info-Tabellen - NUR wenn Methode aktiviert
        if use_euclidean && ~isempty(table_euclidean_info)
            batch_euclidean_info = [batch_euclidean_info; table_euclidean_info];
        end
        if use_sidtw && ~isempty(table_sidtw_info)
            batch_sidtw_info = [batch_sidtw_info; table_sidtw_info];
        end
        if use_dtw && ~isempty(table_dtw_info)
            batch_dtw_info = [batch_dtw_info; table_dtw_info];
        end
        if use_dfd && ~isempty(table_dfd_info)
            batch_dfd_info = [batch_dfd_info; table_dfd_info];
        end
        if use_lcss && ~isempty(table_lcss_info)
            batch_lcss_info = [batch_lcss_info; table_lcss_info];
        end
        
        % Sammle Abweichungsdaten - NUR wenn Methode aktiviert
        if use_euclidean && ~isempty(table_euclidean_deviation)
            batch_euclidean_deviations = [batch_euclidean_deviations; table_euclidean_deviation];
        end
        if use_sidtw && ~isempty(table_sidtw_deviation)
            batch_sidtw_deviations = [batch_sidtw_deviations; table_sidtw_deviation];
        end
        if use_dtw && ~isempty(table_dtw_deviation)
            batch_dtw_deviations = [batch_dtw_deviations; table_dtw_deviation];
        end
        if use_dfd && ~isempty(table_dfd_deviation)
            batch_dfd_deviations = [batch_dfd_deviations; table_dfd_deviation];
        end
        if use_lcss && ~isempty(table_lcss_deviation)
            batch_lcss_deviations = [batch_lcss_deviations; table_lcss_deviation];
        end
        
        disp(['Bahn ', bahn_id, ' erfolgreich verarbeitet']);
    catch e
        warning(['Fehler bei der Verarbeitung von Bahn ', bahn_id, ': ', e.message]);
        continue;
    end
end
end
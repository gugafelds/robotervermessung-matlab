function uploadBatchData(conn, batch_processed_ids, evaluation_type, ...
    evaluate_velocity, evaluate_orientation, ...
    upload_info, upload_deviations, ...
    batchInfo, batchDeviations, ...
    use_euclidean, use_sidtw, use_dtw, use_dfd, use_lcss)

% IDs für SQL IN-Query vorbereiten
quotedIds = strcat("'", string(batch_processed_ids), "'");
idString = strjoin(quotedIds, ',');

% === Info-Tabellen ===
if upload_info
    disp('Lade Info-Tabellen im Batch hoch...');
    infoTables = {};
    
    % NUR aktivierte Methoden hinzufügen
    if use_sidtw
        infoTables(end+1,:) = {'robotervermessung.auswertung.info_sidtw', batchInfo.sidtw};
    end
    if use_dtw
        infoTables(end+1,:) = {'robotervermessung.auswertung.info_dtw', batchInfo.dtw};
    end
    if use_dfd
        infoTables(end+1,:) = {'robotervermessung.auswertung.info_dfd', batchInfo.dfd};
    end
    
    if ~evaluate_velocity
        if use_euclidean
            infoTables(end+1,:) = {'robotervermessung.auswertung.info_euclidean', batchInfo.euclidean};
        end
        if use_lcss
            infoTables(end+1,:) = {'robotervermessung.auswertung.info_lcss', batchInfo.lcss};
        end
    end
    
    for i = 1:size(infoTables,1)
        tableName = infoTables{i,1};
        data = infoTables{i,2};
        if ~isempty(data)
            % Löschen der Daten damit diese nicht doppelt vorliegen
            query = sprintf( ...
                "DELETE FROM %s WHERE bahn_id IN (%s) AND evaluation = '%s';", ...
                tableName, idString, evaluation_type);
            execute(conn, query);
            % Hochladen der Daten
            batchUpload2PostgreSQL(tableName, data, conn);
        end
    end
    disp('Info-Tabellen erfolgreich hochgeladen!');
end

% === Deviation-Tabellen ===
if upload_deviations
    disp('Lade Abweichungs-Tabellen im Batch hoch...');
    deviationTables = {};
    
    if evaluate_orientation && ~evaluate_velocity
        if use_euclidean
            deviationTables(end+1,:) = {'robotervermessung.auswertung.orientation_euclidean', batchDeviations.euclidean};
        end
        if use_sidtw
            deviationTables(end+1,:) = {'robotervermessung.auswertung.orientation_sidtw', batchDeviations.sidtw};
        end
        if use_dtw
            deviationTables(end+1,:) = {'robotervermessung.auswertung.orientation_dtw', batchDeviations.dtw};
        end
        if use_dfd
            deviationTables(end+1,:) = {'robotervermessung.auswertung.orientation_dfd', batchDeviations.dfd};
        end
        if use_lcss
            deviationTables(end+1,:) = {'robotervermessung.auswertung.orientation_lcss', batchDeviations.lcss};
        end
    elseif evaluate_velocity && ~evaluate_orientation
        if use_sidtw
            deviationTables(end+1,:) = {'robotervermessung.auswertung.speed_sidtw', batchDeviations.sidtw};
        end
        if use_dtw
            deviationTables(end+1,:) = {'robotervermessung.auswertung.speed_dtw', batchDeviations.dtw};
        end
        if use_dfd
            deviationTables(end+1,:) = {'robotervermessung.auswertung.speed_dfd', batchDeviations.dfd};
        end
    else
        if use_euclidean
            deviationTables(end+1,:) = {'robotervermessung.auswertung.position_euclidean', batchDeviations.euclidean};
        end
        if use_sidtw
            deviationTables(end+1,:) = {'robotervermessung.auswertung.position_sidtw', batchDeviations.sidtw};
        end
        if use_dtw
            deviationTables(end+1,:) = {'robotervermessung.auswertung.position_dtw', batchDeviations.dtw};
        end
        if use_dfd
            deviationTables(end+1,:) = {'robotervermessung.auswertung.position_dfd', batchDeviations.dfd};
        end
        if use_lcss
            deviationTables(end+1,:) = {'robotervermessung.auswertung.position_lcss', batchDeviations.lcss};
        end
    end
    
    for i = 1:size(deviationTables, 1)
        tableName = deviationTables{i, 1};
        data = deviationTables{i, 2};
        if ~isempty(data)
            % Löschen der Daten damit diese nicht doppelt vorliegen
            query = sprintf( ...
                "DELETE FROM %s WHERE bahn_id IN (%s);", ...
                tableName, idString);
            execute(conn, query);
            % Hochladen der Daten
            batchUpload2PostgreSQL(tableName, data, conn);
        end
    end
    disp('Abweichungs-Tabellen erfolgreich hochgeladen!');
end
end
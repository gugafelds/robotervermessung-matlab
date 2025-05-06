function uploadBatchData(conn, batch_processed_ids, evaluation_type, ...
                         evaluate_velocity, evaluate_orientation, ...
                         upload_info, upload_deviations, ...
                         batchInfo, batchDeviations)

    % IDs für SQL IN-Query vorbereiten
    quotedIds = strcat("'", string(batch_processed_ids), "'");
    idString = strjoin(quotedIds, ',');

    % === Info-Tabellen ===
    if upload_info
        disp('Lade Info-Tabellen im Batch hoch...');

        infoTables = {
            'robotervermessung.auswertung.info_sidtw', batchInfo.sidtw;
            'robotervermessung.auswertung.info_dtw', batchInfo.dtw;
            'robotervermessung.auswertung.info_dfd', batchInfo.dfd
        };

        if ~evaluate_velocity
            infoTables(end+1,:) = {'robotervermessung.auswertung.info_euclidean', batchInfo.euclidean};
            infoTables(end+1,:) = {'robotervermessung.auswertung.info_lcss', batchInfo.lcss};
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
            deviationTables = {
                'robotervermessung.auswertung.orientation_euclidean', batchDeviations.euclidean;
                'robotervermessung.auswertung.orientation_sidtw', batchDeviations.sidtw;
                'robotervermessung.auswertung.orientation_dtw', batchDeviations.dtw;
                'robotervermessung.auswertung.orientation_dfd', batchDeviations.dfd;
                'robotervermessung.auswertung.orientation_lcss', batchDeviations.lcss
            };
        elseif evaluate_velocity && ~evaluate_orientation
            deviationTables = {
                'robotervermessung.auswertung.speed_sidtw', batchDeviations.sidtw;
                'robotervermessung.auswertung.speed_dtw', batchDeviations.dtw;
                'robotervermessung.auswertung.speed_dfd', batchDeviations.dfd
            };
        else
            deviationTables = {
                'robotervermessung.auswertung.position_euclidean', batchDeviations.euclidean;
                'robotervermessung.auswertung.position_sidtw', batchDeviations.sidtw;
                'robotervermessung.auswertung.position_dtw', batchDeviations.dtw;
                'robotervermessung.auswertung.position_dfd', batchDeviations.dfd;
                'robotervermessung.auswertung.position_lcss', batchDeviations.lcss
            };
        end

        for i = 1:size(deviationTables, 1)
            tableName = deviationTables{i, 1};
            data = deviationTables(i,2:end)';
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

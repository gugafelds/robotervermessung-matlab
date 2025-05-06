function uploadTransformedData(conn, batch_ids, calibration_ids, data_ist_list, data_ist_trafo_list, q_transformed_list, schema)
    % UPLOADTRANSFORMEDDATA Uploads multiple transformed datasets in batch mode
    %
    % Eingabe:
    %   conn - Datenbankverbindung
    %   batch_ids - Cell-Array mit Bahn-IDs
    %   calibration_ids - Cell-Array mit Kalibrier-IDs (gleiche Länge wie batch_ids)
    %   data_ist_list - Cell-Array mit data_ist Tabellen
    %   data_ist_trafo_list - Cell-Array mit transformierten Positionen
    %   q_transformed_list - Cell-Array mit transformierten Quaternionen
    %   schema - Datenbankschema
    
    try
        % Überprüfen, ob alle Arrays die gleiche Länge haben
        if length(batch_ids) ~= length(calibration_ids) || ...
           length(batch_ids) ~= length(data_ist_list) || ...
           length(batch_ids) ~= length(data_ist_trafo_list) || ...
           length(batch_ids) ~= length(q_transformed_list)
            error('Alle Eingabe-Arrays müssen die gleiche Länge haben');
        end
        
        % IDs für SQL IN-Query vorbereiten
        quotedIds = strcat("'", string(batch_ids), "'");
        idString = strjoin(quotedIds, ',');
        
        % Daten aus existierenden Tabellen löschen
        query = sprintf("DELETE FROM robotervermessung.%s.bahn_pose_trans WHERE bahn_id IN (%s);", ...
            schema, idString);
        execute(conn, query);
        disp(['Bestehende Einträge für ' num2str(length(batch_ids)) ' Bahnen gelöscht']);
        
        % Gesamtgröße für die Batch-Tabelle berechnen
        total_rows = 0;
        for i = 1:length(data_ist_list)
            total_rows = total_rows + height(data_ist_list{i});
        end
        
        % Große Batch-Tabelle erstellen
        bahn_pose_trans_batch = table('Size', [total_rows, 11], ...
            'VariableTypes', {'string', 'string', 'string', 'double', 'double', ...
                             'double', 'double', 'double', 'double', 'double', 'string'}, ...
            'VariableNames', {'bahn_id', 'segment_id', 'timestamp', 'x_trans', ...
                             'y_trans', 'z_trans', 'qx_trans', 'qy_trans', ...
                             'qz_trans', 'qw_trans', 'calibration_id'});
        
        % Tabelle mit Daten füllen
        row_index = 1;
        for i = 1:length(batch_ids)
            data_ist = data_ist_list{i};
            data_ist_trafo = data_ist_trafo_list{i};
            q_transformed = q_transformed_list{i};
            calibration_id = calibration_ids{i};
            bahn_id = batch_ids{i};
            
            rows = height(data_ist);
            calibration_id_array = repelem(string(calibration_id), rows)';
            bahn_id_array = repelem(string(bahn_id), rows)';
            
            % Daten in die Batch-Tabelle einfügen
            bahn_pose_trans_batch(row_index:row_index+rows-1, :) = table(...
                bahn_id_array, ...                    % bahn_id
                data_ist.segment_id, ...              % segment_id
                data_ist.timestamp, ...               % timestamp
                data_ist_trafo(:,1), ...              % x_trans
                data_ist_trafo(:,2), ...              % y_trans
                data_ist_trafo(:,3), ...              % z_trans
                q_transformed(:,2), ...               % qx_trans
                q_transformed(:,3), ...               % qy_trans
                q_transformed(:,4), ...               % qz_trans
                q_transformed(:,1), ...               % qw_trans
                calibration_id_array ...              % calibration_id
            );
            
            row_index = row_index + rows;
        end
        
        % Die optimierte Batch-Upload-Funktion verwenden statt sqlwrite
        tablename = ['robotervermessung.' schema '.bahn_pose_trans'];
        batchUpload2PostgreSQL(tablename, bahn_pose_trans_batch, conn);
        
        disp([num2str(length(batch_ids)) ' Bahnen erfolgreich in einer Batch-Operation in die Datenbank geschrieben']);
        
    catch ME
        error('Fehler beim Batch-Upload der transformierten Daten: %s', ME.message);
    end
end
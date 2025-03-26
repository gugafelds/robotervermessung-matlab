function batchUpload2PostgreSQL(tablename, data, conn)
    try
        % Combine cell data if needed
        if iscell(data)
            try
                data = vertcat(data{:});
            catch e
                disp(['Fehler beim Kombinieren der Zell-Daten: ' e.message]);
                % Process each cell individually
                for i = 1:length(data)
                    if ~isempty(data{i})
                        sqlwrite(conn, tablename, data{i});
                    end
                end
                return;
            end
        end
        
        % Standardize data
        data = standardizeTableData(data);
        
        % Get total rows for progress tracking
        total_rows = height(data);
        disp(['Lade ' num2str(total_rows) ' Zeilen in ' tablename ' hoch...']);
        
        % Optimierte Einstellungen für Bulk-Import
        execute(conn, 'SET synchronous_commit = off');
        
        try
            % Manuell Transaktion starten
            execute(conn, 'BEGIN');
            
            % Optimierter Multi-Row Insert
            batchInsertMultipleRows(conn, tablename, data);
            % ---> Hier Fehler bei Upload Orientierung
            
            % Commit transaction
            execute(conn, 'COMMIT');
            execute(conn, 'SET synchronous_commit = on');
            disp(['Erfolgreich ' num2str(total_rows) ' Zeilen in ' tablename ' hochgeladen']);
            
        catch inner_e
            % Rollback on error
            execute(conn, 'ROLLBACK');
            execute(conn, 'SET synchronous_commit = on');
            disp(['Fehler beim Batch-Upload, versuche Einzelzeilenmodus: ' inner_e.message]);
            
            % Try again with explicit INSERT statements for smaller batches
            batch_size = 250;
            for i = 1:batch_size:total_rows
                try
                    end_idx = min(i + batch_size - 1, total_rows);
                    current_batch = data(i:end_idx, :);
                    insertWithExplicitSQL(conn, tablename, current_batch);
                    
                    if mod(i, 500) == 0 || i == 1
                        disp(['Alternative Methode: ' num2str(end_idx) ' von ' num2str(total_rows)]);
                    end
                catch batch_e
                    disp(['Fehler bei Zeilen ' num2str(i) '-' num2str(end_idx) ': ' batch_e.message]);
                end
            end
        end
        
    catch e
        disp(['Fehler beim Hochladen der Daten: ' e.message]);
        disp('Stack Trace:');
        disp(getReport(e, 'extended'));
    end
end

function batchInsertMultipleRows(conn, tablename, data)
    % Get column names
    colNames = strjoin(data.Properties.VariableNames, ',');
    
    % Use multiple VALUES sets in a single statement
    batchSize = 10000; % Größe anpassen basierend auf Ihren Daten
    total_rows = height(data);
    
    for i = 1:batchSize:total_rows
        % Start building INSERT statement
        insertSql = ['INSERT INTO ' tablename ' (' colNames ') VALUES '];
        
        % Determine end index for this batch
        end_idx = min(i + batchSize - 1, total_rows);
        current_batch = data(i:end_idx, :);
        
        % Build values part with multiple rows
        valueStrings = cell(height(current_batch), 1);
        for j = 1:height(current_batch)
            rowValues = formatRowValues(current_batch(j,:));
            valueStrings{j} = ['(' rowValues ')'];
        end
        
        % Combine into final SQL
        insertSql = [insertSql strjoin(valueStrings, ', ')];
        
        % Execute the multi-row insert
        execute(conn, insertSql);
        
        % Show progress
        if mod(end_idx, batchSize*5) == 0 || end_idx == total_rows
            disp(['Hochgeladen: ' num2str(end_idx) ' von ' num2str(total_rows) ' Zeilen (' num2str(round(100*end_idx/total_rows)) '%)']);
        end
    end
end

% Alternative method using explicit SQL INSERT for small batches
function insertWithExplicitSQL(conn, tablename, data)
    % Get column names
    colNames = strjoin(data.Properties.VariableNames, ', ');
    
    % Process up to 10 rows at a time
    for i = 1:10:height(data)
        end_idx = min(i + 9, height(data));
        current_rows = data(i:end_idx, :);
        
        % Build INSERT statement
        insertSql = ['INSERT INTO ' tablename ' (' colNames ') VALUES '];
        
        valueStrings = cell(height(current_rows), 1);
        for j = 1:height(current_rows)
            rowValues = formatRowValues(current_rows(j,:));
            valueStrings{j} = ['(' rowValues ')'];
        end
        
        insertSql = [insertSql strjoin(valueStrings, ', ')];
        
        % Execute insert
        execute(conn, insertSql);
    end
end

function valueStr = formatRowValues(row)
    values = cell(1, width(row));
    
    for i = 1:width(row)
        val = row{1,i};
        
        if ischar(val) || isstring(val)
            values{i} = ['''' strrep(char(val), '''', '''''') ''''];
        elseif iscell(val)
            if isempty(val)
                values{i} = 'NULL';
            else
                values{i} = ['''' strrep(char(val{1}), '''', '''''') ''''];
            end
        elseif isnumeric(val)
            if isnan(val)
                values{i} = 'NULL';
            else
                values{i} = num2str(val);
            end
        else
            values{i} = 'NULL';
        end
    end
    
    valueStr = strjoin(values, ', ');
end

% Hilfsfunktion zum Standardisieren der Daten
function standardizedData = standardizeTableData(data)
    % Überprüfe jede Spalte auf konsistente Datentypen
    for col = 1:width(data)
        colName = data.Properties.VariableNames{col};
        colData = data.(colName);
        
        % Wenn es sich um eine Zellen-Spalte handelt
        if iscell(colData)
            % Prüfe auf Zellen, die selbst Zellen enthalten
            hasNestedCells = false;
            for i = 1:length(colData)
                if iscell(colData{i})
                    hasNestedCells = true;
                    break;
                end
            end
            
            if hasNestedCells
                % Unneste die Zellen
                newColData = cell(size(colData));
                for i = 1:length(colData)
                    if iscell(colData{i}) && ~isempty(colData{i})
                        newColData{i} = colData{i}{1};
                    else
                        newColData{i} = colData{i};
                    end
                end
                data.(colName) = newColData;
            end
        end
    end
    
    standardizedData = data;
end
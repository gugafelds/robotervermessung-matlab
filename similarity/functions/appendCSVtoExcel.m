function appendCSVtoExcel(csv_pattern, excel_file, sheet_name, timestamp_col)
% APPEND_CSV_TO_EXCEL - Lädt CSVs, dedupliziert nach Timestamp, schreibt zu Excel
%
% Inputs:
%   csv_pattern   - z.B. 'results/embedding_validation_*.csv'
%   excel_file    - z.B. 'results/experiment_data.xlsx'
%   sheet_name    - z.B. 'embedding_validation'
%   timestamp_col - z.B. 'Timestamp' (Spaltenname für Deduplizierung)

    fprintf('\n Smart Excel Export - %s\n', sheet_name);
    fprintf('────────────────────────────────────────\n\n');
    
    % 1. CSV-Dateien finden
    csv_files = dir(csv_pattern);
    
    if isempty(csv_files)
        fprintf('ℹ Keine CSV-Dateien gefunden für: %s\n\n', csv_pattern);
        return;
    end
    
    fprintf('Gefunden: %d CSV-Datei(en)\n', length(csv_files));
    
    % 2. Alle CSVs laden und kombinieren
    all_data = loadAndCombineCSVs(csv_files, csv_pattern);
    fprintf('Total geladen: %d Zeilen\n', height(all_data));
    
    % 4. Mit Excel synchronisieren
    syncToExcel(all_data, excel_file, sheet_name, timestamp_col);
end


function combined = loadAndCombineCSVs(csv_files, csv_pattern)
% Lädt alle CSVs und kombiniert sie
    combined = [];
    
    for i = 1:length(csv_files)
        filepath = fullfile(csv_files(i).folder, csv_files(i).name);
        temp = readtable(filepath, 'VariableNamingRule', 'preserve');
        
        % String-Spalten normalisieren
        temp = normalizeStringCols(temp);
        
        % Timestamp aus Dateinamen extrahieren
        [~, name, ~] = fileparts(csv_files(i).name);
        pattern_base = regexprep(csv_pattern, '\*\.csv$', '');
        pattern_base = regexprep(pattern_base, '^.*/', '');
        timestamp_str = regexprep(name, ['^' pattern_base], '');
        
        if ~ismember('timestamp', temp.Properties.VariableNames)
            temp.Timestamp = repmat({timestamp_str}, height(temp), 1);
        end
        
        combined = [combined; temp];
        fprintf('  Geladen: %s (%d Zeilen)\n', csv_files(i).name, height(temp));
    end
end


function data = normalizeStringCols(data)
% Konvertiert numerische ID-Spalten zu Strings
    string_cols = {'query_id', 'segment_id', 'level'};
    
    for i = 1:length(string_cols)
        col = string_cols{i};
        if ismember(col, data.Properties.VariableNames) && ~iscell(data.(col))
            data.(col) = cellstr(string(data.(col)));
        end
    end
end


function syncToExcel(new_data, excel_file, sheet_name, timestamp_col)
% Synchronisiert Daten mit Excel (nur neue Timestamps hinzufügen)

    if ~exist(excel_file, 'file')
        fprintf('  ⚠ Excel-Datei existiert nicht: %s\n\n', excel_file);
        return;
    end
    
    try
        % Existierende Daten laden
        existing = readtable(excel_file, 'Sheet', sheet_name, 'VariableNamingRule', 'preserve');
        fprintf('\n  Existierende Daten: %d Zeilen\n', height(existing));
        
        % Timestamps vergleichen
        if ~ismember(timestamp_col, existing.Properties.VariableNames) || ...
           ~ismember(timestamp_col, new_data.Properties.VariableNames)
            fprintf('  ⚠ Timestamp-Spalte "%s" nicht gefunden\n\n', timestamp_col);
            return;
        end
        
        existing_ts = unique(existing.(timestamp_col));
        new_ts = unique(new_data.(timestamp_col));
        ts_to_add = setdiff(new_ts, existing_ts);
        
        if isempty(ts_to_add)
            fprintf('  ℹ Alle Daten bereits vorhanden - nichts hinzuzufügen\n\n');
            return;
        end
        
        % Neue Daten filtern und appenden
        mask = ismember(new_data.(timestamp_col), ts_to_add);
        data_to_add = new_data(mask, :);
        
        fprintf('  → %d neue Timestamps, %d Zeilen hinzufügen\n', length(ts_to_add), height(data_to_add));
        
        combined = [existing; data_to_add];
        writetable(combined, excel_file, 'Sheet', sheet_name, 'WriteRowNames', false);
        fprintf('  ✓ Erfolgreich! Jetzt %d Zeilen total\n\n', height(combined));
        
    catch
        % Sheet existiert nicht - neu erstellen
        fprintf('  ℹ Sheet "%s" existiert nicht - erstelle neu\n', sheet_name);
        writetable(new_data, excel_file, 'Sheet', sheet_name, 'WriteRowNames', false);
        fprintf('  ✓ Erstellt mit %d Zeilen\n\n', height(new_data));
    end
end
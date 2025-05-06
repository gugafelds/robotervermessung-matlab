function [offset_matrix, robot_to_instrument, matrix_file] = findTransfMatrices(input_date)
% FINDMATRICESBYDATE Findet die zum angegebenen Datum passenden Transformationsmatrizen
%
% Eingabe:
%   input_date - Datum in verschiedenen Formaten:
%               - 'yyyyMMdd' (z.B. '20250406')
%               - 'dd.MM.yyyy' (z.B. '06.04.2025') 
%               - datetime-Objekt
%
% Ausgabe:
%   offset_matrix - Die zum Datum oder davor passende Offset-Matrix
%   robot_to_instrument - Die zum Datum oder davor passende Robot-to-Instrument-Matrix
%
% Beispiel:
%   [offset, robot2instr] = findMatricesByDate('20250405');
%   [offset, robot2instr] = findMatricesByDate('06.04.2025');
%   [offset, robot2instr] = findMatricesByDate(datetime('2025-04-05'));

    % Relativer Pfad zu den Matrix-Dateien im Repository
    base_dir = fullfile('lasertracker', 'transformation_matrices');
    
    % Zum Repository-Pfad hinzufügen, falls noch nicht im Pfad
    if ~contains(path, base_dir)
        % Aktuelles Verzeichnis des Scripts finden
        current_script_path = fileparts(mfilename('fullpath'));
        % Wenn wir im 'lasertracker'-Ordner sind, gehen wir einen Schritt zurück
        if contains(current_script_path, 'lasertracker')
            repo_root = fileparts(current_script_path);
        else
            % Ansonsten nehmen wir an, dass wir im Root des Repositories sind
            repo_root = current_script_path;
        end
        % Vollständigen Pfad zum Matrix-Ordner erstellen
        matrix_dir = fullfile(repo_root, base_dir);
        
        % Zum MATLAB-Pfad hinzufügen
        addpath(matrix_dir);
        fprintf('Pfad zu den Matrizen hinzugefügt: %s\n', matrix_dir);
    end
    
    % Formatierung des Eingabedatums
    if isa(input_date, 'datetime')
        % Wenn datetime-Objekt übergeben wurde
        date_str = datetime('now', 'yyyymmdd');
    else
        % Wenn String übergeben wurde
        if contains(input_date, '.')
            % Format: 'dd.MM.yyyy'
            parts = split(input_date, '.');
            if length(parts) == 3
                day = parts{1};
                month = parts{2};
                year = parts{3};
                % Stelle sicher, dass Tag und Monat zweistellig sind
                if length(day) == 1, day = ['0', day]; end
                if length(month) == 1, month = ['0', month]; end
                date_str = [year, month, day];
            else
                error('Ungültiges Datumsformat. Erwartet: dd.MM.yyyy oder yyyyMMdd');
            end
        else
            % Annahme: Format ist bereits 'yyyyMMdd'
            date_str = input_date;
        end
    end
    
    % Datum als Zahl für Vergleiche konvertieren
    target_date_num = str2double(date_str);
    
    % Bestimme den vollständigen Pfad zum Matrix-Ordner
    if exist('matrix_dir', 'var')
        full_base_dir = matrix_dir;
    else
        % Wenn der Pfad bereits in MATLAB ist
        full_base_dir = fileparts(which('matrices_20250404.m'));
    end
    
    % Alle Matrix-Dateien im Verzeichnis auflisten
    files = dir(fullfile(full_base_dir, 'matrices_*.m'));
    
    % Verfügbare Datumsangaben extrahieren und als Zahlen konvertieren
    available_dates = zeros(length(files), 1);
    for i = 1:length(files)
        file_name = files(i).name;
        % Extrahiere das Datum aus dem Dateinamen (Format: matrices_YYYYMMDD.m)
        date_part = regexprep(file_name, 'matrices_(\d+)\.m', '$1');
        available_dates(i) = str2double(date_part);
    end
    
    % Finde das nächste Datum vor oder gleich dem Zieldatum
    valid_dates = available_dates(available_dates <= target_date_num);
    
    if isempty(valid_dates)
        error('Keine Matrix-Dateien vor dem angegebenen Datum gefunden');
    end
    
    % Wähle das neueste gültige Datum
    [~, idx] = max(valid_dates);
    best_date = num2str(valid_dates(idx));
    
    % Matrix-Datei als String
    matrix_file = ['matrices_', best_date];
    
    % Lade die Matrizen aus der Datei
    try
        % Führe das Script aus, um die Matrizen in den Workspace zu laden
        run(matrix_file);
        
        % Überprüfe, ob die Matrizen nach dem Ausführen verfügbar sind
        if ~exist('offset_matrix', 'var') || ~exist('robot_to_instrument', 'var')
            error('Die Matrizen wurden nicht korrekt vom Script definiert');
        end
        
        % Formatiere das Datum für die Ausgabe
        formatted_date = [best_date(7:8), '.', best_date(5:6), '.', best_date(1:4)];
        fprintf('Matrizen vom %s erfolgreich geladen.\n', formatted_date);
    catch ME
        error('Fehler beim Laden der Matrizen: %s', ME.message);
    end
end
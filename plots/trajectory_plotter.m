%% MATLAB Script: Plot Trajectories from Database
clear; clc; close all;

% 1. Verbindung zur Datenbank herstellen
% Ersetze 'deine_datenbank.db' durch den tatsächlichen Pfad zu deiner Datei
conn = connectingToPostgres();

% 2. Auswahl der Bahn-IDs (Hier einfach die gewünschten IDs einkommentieren)
selected_ids = {
    '1765991743';
    '1765991747'; 
    '1765991751';
    '1765991754';
    '1765991758';
    '1765991762';
    '1765991766';
    '1765991770';
    '1765991774';
    '1765991778';
};

figure('Name', 'Trajektorien Plot', 'Color', 'w');
hold on; grid on; axis equal;
xlabel('X Soll [mm]'); ylabel('Y Soll [mm]'); zlabel('Z Soll [mm]');

colors = lines(length(selected_ids));

for i = 1:length(selected_ids)
    current_id = selected_ids{i};
    
    % Query - stelle sicher, dass die Spaltennamen exakt so in der DB stehen
    sqlQuery = sprintf('SELECT x_soll, y_soll, z_soll FROM bewegungsdaten.bahn_position_soll WHERE bahn_id = ''%s'' ORDER BY timestamp ASC', current_id);
    
    try
        data = fetch(conn, sqlQuery);
        
        if isempty(data)
            warning('ID %s ist leer.', current_id);
            continue;
        end
        
        % FIX: Zugriff auf Table-Variablen via Punkt-Notation
        % Das funktioniert, wenn die Spalten in der DB x_soll, y_soll, z_soll heißen
        x = data.x_soll;
        y = data.y_soll;
        z = data.z_soll;
        
        % Falls die Spalten in der DB anders heißen, kannst du auch Indexing nutzen:
        % x = data{:, 1}; % Erste Spalte, alle Zeilen
        % y = data{:, 2}; % Zweite Spalte, alle Zeilen
        % z = data{:, 3}; % Dritte Spalte, alle Zeilen

        plot3(x, y, z, 'LineWidth', 1.5, 'Color', colors(i,:), 'DisplayName', ['ID: ' current_id]);
        
    catch ME
        fprintf('Fehler bei ID %s: %s\n', current_id, ME.message);
    end
end

legend('show');
view(3);
close(conn);
%% MATLAB Script: Multi-Record 3D Analysis with SIDTW Animation
clear; clc; close all;

%% 1. Verbindung zur Datenbank
conn = connectingToPostgres();

%% 2. Konfiguration
target_bahn_id = '1762944505';
seg_1 = [target_bahn_id '_1']; % Erstes Segment
seg_2 = [target_bahn_id '_2']; % Zweites Segment
num_points = 500;               % Punkte pro Segment am Übergang

%% 1. Daten laden (Alles für eine Bahn_ID)
% Wir holen alle Daten sortiert nach Segment und Order
sql_query = sprintf('SELECT * FROM robotervermessung.auswertung.position_sidtw WHERE bahn_id = ''%s'' ORDER BY segment_id, points_order', target_bahn_id);
all_data = fetch(conn, sql_query);

%% 2. Übergang finden (Einfacher Lookup)
% Wir wandeln die segment_id in einen cell-array von strings um
ids = cellstr(all_data.segment_id);

% Wir suchen die Stelle, wo ID(i) ungleich ID(i+1) ist
% Das ist der letzte Punkt des ersten Segments
change_idx = find(~strcmp(ids(1:end-1), ids(2:end)), 1);

if isempty(change_idx)
    error('Kein Segmentwechsel gefunden! Besteht die Bahn nur aus einem Segment?');
end

% Bereich: 300 davor und 300 danach
start_idx = max(1, change_idx - 299);
end_idx   = min(height(all_data), change_idx + 300);

plot_data = all_data(start_idx:end_idx, :);

% Koordinaten extrahieren
soll = [plot_data.sidtw_soll_x, plot_data.sidtw_soll_y, plot_data.sidtw_soll_z];
ist  = [plot_data.sidtw_ist_x, plot_data.sidtw_ist_y, plot_data.sidtw_ist_z];

%% 3. Animation Setup
figure('Color', 'w', 'Position', [100, 100, 2500, 2500]);
ax = axes('NextPlot', 'add');
% Entferne alle weißen Ränder um die Grafik herum
set(gcf, 'Units', 'Normalized', 'InnerPosition', [0 0 1 1]);
view(3); grid off; box off; axis off;

c_motion = [0.86, 0.13, 0.15]; % Rot
c_shape  = [0.15, 0.39, 0.91]; % Blau

% Statische Vorab-Pfade (blass im Hintergrund für Kontext)
plot3(soll(:,1), soll(:,2), soll(:,3), 'Color', [0.8 0.8 0.8], 'LineWidth', 0.8, 'HandleVisibility', 'off');

% Dynamische Objekte
hSoll = plot3(nan, nan, nan, 'k-', 'LineWidth', 5, 'Color', c_motion);
hIst  = plot3(nan, nan, nan, 'b-', 'LineWidth', 6, 'Color', c_shape);
hPointSoll = scatter3(nan, nan, nan, 300, c_motion, 'filled', 'DisplayName', 'Aktueller Soll-Punkt');
hPointIst  = scatter3(nan, nan, nan, 300, c_shape, 'filled', 'DisplayName', 'Aktueller Ist-Punkt');

%% 4. Animation ausführen mit Fokus-Zoom
% Ziel-Zoom Bereich definieren (z.B. 100mm Fenster um den Eckpunkt)
zoom_window = 0.85; 

% Wir berechnen den Index des Eckpunktes innerhalb unserer 'plot_data'
% Da plot_data bei start_idx beginnt, ist der lokale Index:
local_change_idx = change_idx - start_idx + 1;

for k = 1:2:height(plot_data)
    % --- Vorhandene Plot-Updates ---
    set(hSoll, 'XData', soll(1:k, 1), 'YData', soll(1:k, 2), 'ZData', soll(1:k, 3));
    set(hIst,  'XData', ist(1:k, 1),  'YData', ist(1:k, 2),  'ZData', ist(1:k, 3));
    set(hPointSoll, 'XData', soll(k,1), 'YData', soll(k,2), 'ZData', soll(k,3));
    set(hPointIst,  'XData', ist(k,1),  'YData', ist(k,2),  'ZData', ist(k,3));
    line([soll(k,1) ist(k,1)], [soll(k,2) ist(k,2)], [soll(k,3) ist(k,3)], ...
         'Color', [0.6, 0.6, 0.6], 'LineWidth', 1);

    % --- NEU: Zoom Logik ---
    % Wenn wir uns dem Eckpunkt nähern (z.B. ab Punkt 100 der Animation)
    if k > 80 && k < (height(plot_data) - 80)
        % Wir berechnen ein dynamisches Limit, das sich um den aktuellen Ist-Punkt zusammenzieht
        % Der Faktor '0.95' sorgt für ein langsames Heranfahren
        current_x = ist(k,1);
        current_y = ist(k,2);
        current_z = ist(k,3);
        
        % Sanftes Anpassen der Achsen (Interp-Effekt)
        xlim(ax, [current_x - zoom_window, current_x + zoom_window]);
        ylim(ax, [current_y - zoom_window, current_y + zoom_window]);
        zlim(ax, [current_z - zoom_window, current_z + zoom_window]);
    end

    camorbit(-0.4, -0.3); 
    drawnow;
    pause(0.01)
end
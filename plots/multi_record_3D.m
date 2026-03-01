%% MATLAB Script: Multi-Record 3D Analysis with Triple Zoom
% Kompaktes Layout mit engen Plot-Abständen
clear; clc; close all;

%% 1. Verbindung zur Datenbank
conn = connectingToPostgres();

%% 2. Die drei zu vergleichenden Records definieren
target_records = {
    'record_20251217_173659';
    'record_20251217_180722';
    'record_20260130_123353';
};

% Zusätzliche Records (mit Transparenz im Hintergrund)
additional_records = {
     %'record_20260130_130116';
     %'record_20260130_125157';
     %'record_20260130_124526';
     %'record_20260130_122226';
};

additional_alpha = 0.3;  % Transparenz für zusätzliche Records

%% 3. Zoom-Boxen automatisch aus bahn_events laden (je Record ein Eckpunkt)
fprintf('Lade Eckpunkte aus bahn_events...\n');

% Würfelgröße um die Eckpunkte (in mm)
box_half_size = [50, 50, 50];

zoom_boxes = cell(3, 1);

for r = 1:length(target_records)
    rec_name = target_records{r};
    
    sqlIdQuery = sprintf(['SELECT bahn_id FROM bewegungsdaten.bahn_info ' ...
        'WHERE record_filename = ''%s'' ORDER BY bahn_id ASC LIMIT 1'], rec_name);
    refIdData = fetch(conn, sqlIdQuery);
    
    if isempty(refIdData)
        warning('Keine bahn_id für Record: %s', rec_name);
        continue;
    end
    
    ref_bahn_id = string(refIdData.bahn_id(1));
    
    sqlEventsQuery = sprintf(['SELECT x_reached, y_reached, z_reached ' ...
        'FROM bewegungsdaten.bahn_events ' ...
        'WHERE bahn_id = ''%s'' ORDER BY timestamp ASC'], ref_bahn_id);
    eventsData = fetch(conn, sqlEventsQuery);
    
    if isempty(eventsData)
        warning('Keine Events für bahn_id: %s', ref_bahn_id);
        continue;
    end
    
    num_events = height(eventsData);
    event_idx = round(num_events / 2);
    
    if r == 1
        cx = eventsData.x_reached(event_idx)-4;
        cy = eventsData.y_reached(event_idx)+8;
        cz = eventsData.z_reached(event_idx)-8;
    elseif r == 2
        cx = eventsData.x_reached(event_idx)-4;
        cy = eventsData.y_reached(event_idx)-4;
        cz = eventsData.z_reached(event_idx)+4;
    else
        cx = eventsData.x_reached(event_idx-1)+2;
        cy = eventsData.y_reached(event_idx-1)+2;
        cz = eventsData.z_reached(event_idx-1)+2;
    end
    
    
    zoom_boxes{r} = [
        cx - box_half_size(1), cx + box_half_size(1);
        cy - box_half_size(2), cy + box_half_size(2);
        cz - box_half_size(3), cz + box_half_size(3)
    ];
    
    fprintf('  Zoom %d: Record "%s" -> Event %d/%d bei [%.1f, %.1f, %.1f]\n', ...
        r, rec_name, event_idx, num_events, cx, cy, cz);
end

zoom_factors = [0.3, 0.3, 0.3];

%% 4. Farben definieren
base_colors = [
    0.000, 0.447, 0.741;  % Blau
    0.850, 0.325, 0.098;  % Orange
    0.466, 0.674, 0.188   % Grün
];

box_colors = [
    0.000, 0.447, 0.741;  % Blau
    0.850, 0.325, 0.098;  % Orange
    0.466, 0.674, 0.188   % Grün
];
box_alpha = 0.15;
box_edge_alpha = 0.8;

%% 5. Figure mit kompaktem Layout erstellen
fig = figure('Name', 'Multi-Record 3D Detail-Vergleich', ...
    'Color', 'w', ...
    'Position', [0, 0, 1100, 700], ...
    'NumberTitle', 'off', ...
    'PaperPositionMode', 'auto', ...
    'InvertHardcopy', 'off');

% =========================================================================
% KOMPAKTES LAYOUT - Manuelle Positionierung
% Position = [left, bottom, width, height] (normalisiert 0-1)
% =========================================================================

% Minimale Abstände für maximale Nutzung des Platzes
margin_left = 0.05;
margin_right = 0.005;
margin_bottom = 0.08;
margin_top = 0.01;
gap_h = 0.00005;  % Horizontaler Abstand zwischen Plots
gap_v = 0.01;  % Vertikaler Abstand zwischen Plots

% Hauptplot (links, volle Höhe)
main_width = 0.67;
main_height = 1 - margin_bottom - margin_top;
axMain = axes('Position', [margin_left, margin_bottom, main_width, main_height]);
hold on; grid on; axis equal;
view(140,28);
xlim(axMain,[900,2000]);
ylim(axMain,[-600,750]);
zlim(axMain, [300, 1400]);
xticks([900, 1200, 1500, 1800]);
yticks([-600, -300, 0, 300, 600]);
zticks([600,900,1200]);
xlabel('X [mm]', 'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
ylabel('Y [mm]', 'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
zlabel('Z [mm]', 'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
set(axMain, 'FontSize', 20, 'FontWeight', 'bold', 'GridAlpha', 0.4, 'FontName', 'Times New Roman');
axMain.XAxis.LineWidth = 1.5;
axMain.YAxis.LineWidth = 1.5;
axMain.ZAxis.LineWidth = 1.5;

% Zoom-Plots (rechts, übereinander)
zoom_left = margin_left + main_width + gap_h;
zoom_width = 1 - zoom_left - margin_right;
zoom_height = (main_height - 2 * gap_v) / 3;

axZoom = gobjects(3, 1);
for z = 1:3
    % Von oben nach unten: z=1 oben, z=3 unten
    zoom_bottom = margin_bottom + (3 - z) * (zoom_height + gap_v);
    
    axZoom(z) = axes('Position', [zoom_left, zoom_bottom, zoom_width, zoom_height]);
    hold on; grid on; axis equal;
    view(120,25);
    set(axZoom(z), 'FontSize', 16, 'FontWeight', 'bold', 'GridAlpha', 0.4, 'FontName', 'Times New Roman');
    axZoom(z).XAxis.LineWidth = 1.5;
    axZoom(z).YAxis.LineWidth = 1.5;
    axZoom(z).ZAxis.LineWidth = 1.5;
    
    % Farbiger Hintergrund entsprechend der Zoom-Box-Farbe
    set(axZoom(z), 'Color', [box_colors(z, :), box_alpha]);
    set(axZoom(z), 'YTickLabel', []);
    set(axZoom(z), 'XTickLabel', []);
    set(axZoom(z), 'ZTickLabel', []);
end

%% 6. Daten laden und plotten (mit Farbverlauf)
fprintf('Lade Bahndaten aus der Datenbank...\n');

% Downsampling-Faktor
downsample_factor = 15;  % Nur jeden 10. Punkt plotten

% Erst die Haupt-Records (highlighted)
for r = 1:length(target_records)
    rec_name = target_records{r};
    fprintf('  Record %d/%d: %s\n', r, length(target_records), rec_name);
    
    sqlIdQuery = sprintf(['SELECT bahn_id FROM bewegungsdaten.bahn_info ' ...
        'WHERE record_filename = ''%s'' ORDER BY bahn_id ASC'], rec_name);
    idData = fetch(conn, sqlIdQuery);
    
    if isempty(idData)
        warning('Keine Daten für Record: %s', rec_name);
        continue;
    end
    
    ids = cellstr(string(idData.bahn_id));
    num_bahnen = length(ids);
    base_color = base_colors(r, :);
    
    for i = 1:num_bahnen
        current_id = ids{i};
        
        % Farbverlauf: dunkel -> hell
        t = (i - 1) / max(num_bahnen - 1, 1);
        bahn_color = base_color * (0.4 + 0.6 * t);
        
        sqlQuery = sprintf(['SELECT x_soll, y_soll, z_soll ' ...
            'FROM bewegungsdaten.bahn_position_soll ' ...
            'WHERE bahn_id = ''%s'' ORDER BY timestamp ASC'], current_id);
        
        try
            data = fetch(conn, sqlQuery);
            if isempty(data), continue; end
            
            x = data.x_soll;
            y = data.y_soll;
            z = data.z_soll;
            
            % Downsampling
            idx = 1:downsample_factor:length(x);
            x = x(idx);
            y = y(idx);
            z = z(idx);
            
            plot3(axMain, x, y, z, ...
                'Color', [bahn_color, 0.6], ...
                'LineWidth', 3.5, ...
                'HandleVisibility', 'off');
            
            for j = 1:3
                plot3(axZoom(j), x, y, z, ...
                    'Color', [bahn_color, 0.9], ...
                    'LineWidth', 2.0);
            end
            
        catch ME
            warning('Fehler bei Bahn %s: %s', current_id, ME.message);
        end
    end
    
    if r == 1
    plot3(axMain, NaN, NaN, NaN, 'Color', bahn_color, 'LineWidth', 3, ...
        'DisplayName', ' 0 mm');
    
    elseif r == 2

    plot3(axMain, NaN, NaN, NaN, 'Color', bahn_color, 'LineWidth', 3, ...
        'DisplayName', ' \pm 2 mm');

    else

    plot3(axMain, NaN, NaN, NaN, 'Color', bahn_color, 'LineWidth', 3, ...
        'DisplayName', ' \pm 5 mm');
    
    end
    
end

% Dann die zusätzlichen Records (transparent, schwarz im Hintergrund)
for r = 1:length(additional_records)
    rec_name = additional_records{r};
    fprintf('  Zusätzlicher Record %d/%d: %s (transparent)\n', r, length(additional_records), rec_name);
    
    sqlIdQuery = sprintf(['SELECT bahn_id FROM bewegungsdaten.bahn_info ' ...
        'WHERE record_filename = ''%s'' ORDER BY bahn_id ASC'], rec_name);
    idData = fetch(conn, sqlIdQuery);
    
    if isempty(idData)
        warning('Keine Daten für Record: %s', rec_name);
        continue;
    end
    
    ids = cellstr(string(idData.bahn_id));
    
    for i = 1:length(ids)
        current_id = ids{i};
        
        sqlQuery = sprintf(['SELECT x_soll, y_soll, z_soll ' ...
            'FROM bewegungsdaten.bahn_position_soll ' ...
            'WHERE bahn_id = ''%s'' ORDER BY timestamp ASC'], current_id);
        
        try
            data = fetch(conn, sqlQuery);
            if isempty(data), continue; end
            
            x = data.x_soll;
            y = data.y_soll;
            z = data.z_soll;
            
            % Downsampling
            idx = 1:downsample_factor:length(x);
            x = x(idx);
            y = y(idx);
            z = z(idx);
            
            % Nur im Hauptplot, schwarz und transparent
            plot3(axMain, x, y, z, ...
                'Color', [0, 0, 0, additional_alpha], ...
                'LineWidth', 0.5, ...
                'HandleVisibility', 'off');
            
        catch ME
            warning('Fehler bei Bahn %s: %s', current_id, ME.message);
        end
    end
end

%% 7. 3D-Würfel (Zoom-Boxen) im Hauptplot zeichnen
fprintf('Zeichne Zoom-Boxen...\n');

for z = 1:3
    if isempty(zoom_boxes{z})
        continue;
    end
    
    box = zoom_boxes{z};
    
    x_min = box(1, 1); x_max = box(1, 2);
    y_min = box(2, 1); y_max = box(2, 2);
    z_min = box(3, 1); z_max = box(3, 2);
    
    vertices = [
        x_min, y_min, z_min;
        x_max, y_min, z_min;
        x_max, y_max, z_min;
        x_min, y_max, z_min;
        x_min, y_min, z_max;
        x_max, y_min, z_max;
        x_max, y_max, z_max;
        x_min, y_max, z_max
    ];
    
    faces = [1,2,3,4; 5,6,7,8; 1,2,6,5; 3,4,8,7; 1,4,8,5; 2,3,7,6];
    
    patch(axMain, 'Vertices', vertices, 'Faces', faces, ...
        'FaceColor', box_colors(z, :), ...
        'FaceAlpha', box_alpha, ...
        'EdgeColor', box_colors(z, :), ...
        'EdgeAlpha', box_edge_alpha, ...
        'LineWidth', 1.5, ...
        'HandleVisibility', 'off');
    
    edges = [1,2; 2,3; 3,4; 4,1; 5,6; 6,7; 7,8; 8,5; 1,5; 2,6; 3,7; 4,8];
    
    for e = 1:size(edges, 1)
        plot3(axMain, ...
            vertices(edges(e, :), 1), ...
            vertices(edges(e, :), 2), ...
            vertices(edges(e, :), 3), ...
            'Color', [box_colors(z, :), box_edge_alpha], ...
            'LineWidth', 1.0, ...
            'HandleVisibility', 'off');
    end
    
    
    % Zoom-Limits
    center_x = (x_min + x_max) / 2;
    center_y = (y_min + y_max) / 2;
    center_z = (z_min + z_max) / 2;
    
    half_x = (x_max - x_min) / 2 * zoom_factors(z);
    half_y = (y_max - y_min) / 2 * zoom_factors(z);
    half_z = (z_max - z_min) / 2 * zoom_factors(z);
    
    xlim(axZoom(z), [center_x - half_x, center_x + half_x]);
    ylim(axZoom(z), [center_y - half_y, center_y + half_y]);
    zlim(axZoom(z), [center_z - half_z, center_z + half_z]);
end

%% 8. Legende und finale Anpassungen
legend(axMain, 'show', ...
    'Location', 'northeast', ...
    'Orientation', 'horizontal', ...
    'FontSize', 20, ...
    'FontName', 'Times New Roman', ...
    'NumColumns', 1, ...
    'Box', 'on', 'LineWidth', 1.5);

for ax = [axMain; axZoom(:)]'
    camlight(ax, 'headlight');
    lighting(ax, 'gouraud');
end

rotate3d on;

%% Minimiere weiße Ränder
% Setze LooseInset auf Minimum für alle Achsen
for ax = [axMain; axZoom(:)]'
    ax.LooseInset = [0, 0, 0, 0];
end

%% 9. Verbindung schließen
close(conn);
folder_path = fullfile('similarity', 'figs');
file_name = 'multi_record_3d.pdf';
full_path = fullfile(folder_path, file_name);

exportgraphics(fig, full_path, ContentType="vector")

fprintf('Fertig!\n');
function [info_table, distances_table] = metric2postgresql(metric, distances, soll, ist, bahn_id, segment_ids_distances)
    % Id's in Cell-Arrays konvertieren und Vektor mit Replica der ID's erstellen
    bahn_ids = repmat(string(bahn_id), length(distances), 1);

    segment_ids = segment_ids_distances;
    
    % Wenn euklidischer Abstand
    if strcmp(metric, 'euclidean')
        euclidean_distances = distances;
        ea_soll_x = soll(:,1);
        ea_soll_y = soll(:,2);
        ea_soll_z = soll(:,3);
        ea_ist_x = ist(:,1);
        ea_ist_y = ist(:,2);
        ea_ist_z = ist(:,3);
        
        % distances_table für alle Punkte
        distances_table = table(bahn_ids, segment_ids, euclidean_distances, ...
            ea_soll_x, ea_soll_y, ea_soll_z, ea_ist_x, ea_ist_y, ea_ist_z);
        distances_table.Properties.VariableNames = {'bahn_id', 'segment_id', 'euclidean_distances', ...
            'ea_soll_x', 'ea_soll_y', 'ea_soll_z', ...
            'ea_ist_x', 'ea_ist_y', 'ea_ist_z'};
        
        % info_table: eine Zeile pro Segment + eine für Gesamtbahn
        unique_segments = unique(segment_ids, 'stable');
        
        % Preallocate cell arrays für die Tabellendaten
        n_segments = length(unique_segments);
        bahn_id_col = cell(n_segments + 1, 1);
        segment_id_col = cell(n_segments + 1, 1);
        min_dist = zeros(n_segments + 1, 1);
        max_dist = zeros(n_segments + 1, 1);
        avg_dist = zeros(n_segments + 1, 1);
        std_dist = zeros(n_segments + 1, 1);
        
        % Daten für jedes Segment sammeln
        for i = 1:n_segments
            seg_mask = strcmp(segment_ids, unique_segments{i});
            seg_distances = distances(seg_mask);
            
            bahn_id_col{i} = bahn_id;
            segment_id_col{i} = unique_segments{i};
            min_dist(i) = min(seg_distances);
            max_dist(i) = max(seg_distances);
            avg_dist(i) = mean(seg_distances);
            std_dist(i) = std(seg_distances);
        end
        
        % Gesamtbahn: bahn_id als segment_id
        bahn_id_col{n_segments + 1} = bahn_id;
        segment_id_col{n_segments + 1} = bahn_id;
        min_dist(n_segments + 1) = min(distances);
        max_dist(n_segments + 1) = max(distances);
        avg_dist(n_segments + 1) = mean(distances);
        std_dist(n_segments + 1) = std(distances);
        
        % Tabelle erstellen
        info_table = table(bahn_id_col, segment_id_col, min_dist, max_dist, avg_dist, std_dist);
        info_table.Properties.VariableNames = {'bahn_id', 'segment_id', ...
            'euclidean_min_distance', 'euclidean_max_distance', ...
            'euclidean_average_distance', 'euclidean_standard_deviation'};

    % Wenn SIDTW
    elseif strcmp(metric, 'sidtw')
        sidtw_distances = distances;
        sidtw_soll_x = soll(:,1);
        sidtw_soll_y = soll(:,2);
        sidtw_soll_z = soll(:,3);
        sidtw_ist_x = ist(:,1);
        sidtw_ist_y = ist(:,2);
        sidtw_ist_z = ist(:,3);
        
        % distances_table für alle Punkte
        distances_table = table(bahn_ids, segment_ids, sidtw_distances, ...
            sidtw_soll_x, sidtw_soll_y, sidtw_soll_z, sidtw_ist_x, sidtw_ist_y, sidtw_ist_z);
        distances_table.Properties.VariableNames = {'bahn_id', 'segment_id', 'sidtw_distances' ...
            'sidtw_soll_x', 'sidtw_soll_y', 'sidtw_soll_z' ...
            'sidtw_ist_x', 'sidtw_ist_y', 'sidtw_ist_z'};
        
        % info_table: eine Zeile pro Segment + eine für Gesamtbahn
        unique_segments = unique(segment_ids, 'stable');
        info_table = table();
        
        for i = 1:length(unique_segments)
            seg_mask = strcmp(segment_ids, unique_segments{i});
            seg_distances = distances(seg_mask);
            
            info_table = [info_table; bahn_id, unique_segments(i), ...
                min(seg_distances), max(seg_distances), ...
                mean(seg_distances), std(seg_distances)];
        end
        
        % Gesamtbahn: bahn_id als segment_id
        info_table = [info_table; {bahn_id}, {bahn_id}, ...
            min(distances), max(distances), ...
            mean(distances), std(distances)];
        
        info_table.Properties.VariableNames = {'bahn_id', 'segment_id', ...
            'sidtw_min_distance', 'sidtw_max_distance', ...
            'sidtw_average_distance', 'sidtw_standard_deviation'};

    % Wenn DTW
    elseif strcmp(metric, 'dtw')
        dtw_min_distance = min(distances);
        dtw_max_distance = max(distances);
        dtw_average_distance = mean(distances);
        dtw_standard_deviation = std(distances);

        dtw_distances = distances;

        if size(soll, 2) < 3
            dtw_soll_speed = soll(:,1);
            dtw_ist_speed = ist(:,1);

            % Tabellen anlegen
            if nargin == 6 
                info_table = table(bahn_id, segment_id, dtw_min_distance, ...
                                  dtw_max_distance, dtw_average_distance, ...
                                  dtw_standard_deviation);
                bahn_id = bahn_ids;
                segment_id = segment_ids;
        
                distances_table = table(bahn_id, segment_id, dtw_distances, dtw_soll_speed, dtw_ist_speed);
            else
                info_table = table(bahn_id, dtw_min_distance, ...
                                  dtw_max_distance, dtw_average_distance, ...
                                  dtw_standard_deviation);
                bahn_id = bahn_ids;
            
                distances_table = table(bahn_id, dtw_distances, dtw_soll_speed, dtw_ist_speed);
            end
        else
            dtw_soll_x = soll(:,1);
            dtw_soll_y = soll(:,2);
            dtw_soll_z = soll(:,3);
            dtw_ist_x = ist(:,1);
            dtw_ist_y = ist(:,2);
            dtw_ist_z = ist(:,3);

            % Tabellen anlegen
            if nargin == 6 
                info_table = table(bahn_id, segment_id, dtw_min_distance, ...
                                  dtw_max_distance, dtw_average_distance, ...
                                  dtw_standard_deviation);
                bahn_id = bahn_ids;
                segment_id = segment_ids;
        
                distances_table = table(bahn_id, segment_id, dtw_distances, ...
                    dtw_soll_x, dtw_soll_y, dtw_soll_z, dtw_ist_x, dtw_ist_y, dtw_ist_z);
            else
                info_table = table(bahn_id, dtw_min_distance, ...
                                  dtw_max_distance, dtw_average_distance, ...
                                  dtw_standard_deviation);
                bahn_id = bahn_ids;
            
                distances_table = table(bahn_id, dtw_distances, ...
                    dtw_soll_x, dtw_soll_y, dtw_soll_z, dtw_ist_x, dtw_ist_y, dtw_ist_z);
            end
        end

    % Wenn Frechet
    elseif strcmp(metric, 'dfd')
        dfd_min_distance = min(distances);
        dfd_max_distance = max(distances);
        dfd_average_distance = mean(distances);
        dfd_standard_deviation = std(distances);

        dfd_distances = distances;

        if size(soll, 2) < 3
            dfd_soll_speed = soll(:,1);
            dfd_ist_speed = ist(:,1);

            % Tabellen anlegen
            if nargin == 6
                info_table = table(bahn_id, segment_id, dfd_min_distance, ...
                                  dfd_max_distance, dfd_average_distance, ...
                                  dfd_standard_deviation);
                bahn_id = bahn_ids;
                segment_id = segment_ids;
        
                distances_table = table(bahn_id, segment_id, dfd_distances, dfd_soll_speed, dfd_ist_speed);
            else
                info_table = table(bahn_id, dfd_min_distance, ...
                                  dfd_max_distance, dfd_average_distance, ...
                                  dfd_standard_deviation);
                bahn_id = bahn_ids;
        
                distances_table = table(bahn_id, dfd_distances, dfd_soll_speed, dfd_ist_speed);
            end
        else
            dfd_soll_x = soll(:,1);
            dfd_soll_y = soll(:,2);
            dfd_soll_z = soll(:,3);
            dfd_ist_x = ist(:,1);
            dfd_ist_y = ist(:,2);
            dfd_ist_z = ist(:,3);
        
            % Tabellen anlegen
            if nargin == 6
                info_table = table(bahn_id, segment_id, dfd_min_distance, ...
                                  dfd_max_distance, dfd_average_distance, ...
                                  dfd_standard_deviation);
                bahn_id = bahn_ids;
                segment_id = segment_ids;
        
                distances_table = table(bahn_id, segment_id, dfd_distances, ...
                    dfd_soll_x, dfd_soll_y, dfd_soll_z, dfd_ist_x, dfd_ist_y, dfd_ist_z);
            else
                info_table = table(bahn_id, dfd_min_distance, ...
                                  dfd_max_distance, dfd_average_distance, ...
                                  dfd_standard_deviation);
                bahn_id = bahn_ids;
        
                distances_table = table(bahn_id, dfd_distances, ...
                    dfd_soll_x, dfd_soll_y, dfd_soll_z, dfd_ist_x, dfd_ist_y, dfd_ist_z);
            end
        end

    % Wenn LCSS
    elseif strcmp(metric, 'lcss')
        lcss_min_distance = min(distances);
        lcss_max_distance = max(distances);
        lcss_average_distance = mean(distances);
        lcss_standard_deviation = std(distances);

        lcss_soll_x = soll(:,1);
        lcss_soll_y = soll(:,2);
        lcss_soll_z = soll(:,3);
        lcss_ist_x = ist(:,1);
        lcss_ist_y = ist(:,2);
        lcss_ist_z = ist(:,3);

        lcss_distances = distances;

        % Tabellen anlegen
        if nargin == 6 
            info_table = table(bahn_id, segment_id, lcss_min_distance, ...
                              lcss_max_distance, lcss_average_distance, ...
                              lcss_standard_deviation);
            bahn_id = bahn_ids;
            segment_id = segment_ids;
        
            distances_table = table(bahn_id, segment_id, lcss_distances, ...
                lcss_soll_x, lcss_soll_y, lcss_soll_z, lcss_ist_x, lcss_ist_y, lcss_ist_z);
        else
            info_table = table(bahn_id, lcss_min_distance, ...
                              lcss_max_distance, lcss_average_distance, ...
                              lcss_standard_deviation);
            bahn_id = bahn_ids;
        
            distances_table = table(bahn_id, lcss_distances, ...
                lcss_soll_x, lcss_soll_y, lcss_soll_z, lcss_ist_x, lcss_ist_y, lcss_ist_z);
        end

    elseif strcmp(metric, 'qad')
        qad_distances = distances;
        qad_soll_x = soll(:,1);
        qad_soll_y = soll(:,2);
        qad_soll_z = soll(:,3);
        qad_soll_w = soll(:,4);
        qad_ist_x = ist(:,1);
        qad_ist_y = ist(:,2);
        qad_ist_z = ist(:,3);
        qad_ist_w = ist(:,4);
        
        % distances_table für alle Punkte
        distances_table = table(bahn_ids, segment_ids, qad_distances, ...
            qad_soll_x, qad_soll_y, qad_soll_z, qad_soll_w, qad_ist_x, qad_ist_y, qad_ist_z, qad_ist_w);
        distances_table.Properties.VariableNames = {'bahn_id', 'segment_id', 'qad_distances', ...
            'qad_soll_x', 'qad_soll_y', 'qad_soll_z', 'qad_soll_w', ...
            'qad_ist_x', 'qad_ist_y', 'qad_ist_z', 'qad_ist_w'};
        
        % info_table: eine Zeile pro Segment + eine für Gesamtbahn
        unique_segments = unique(segment_ids, 'stable');
        
        % Preallocate cell arrays für die Tabellendaten
        n_segments = length(unique_segments);
        bahn_id_col = cell(n_segments + 1, 1);
        segment_id_col = cell(n_segments + 1, 1);
        min_dist = zeros(n_segments + 1, 1);
        max_dist = zeros(n_segments + 1, 1);
        avg_dist = zeros(n_segments + 1, 1);
        std_dist = zeros(n_segments + 1, 1);
        
        % Daten für jedes Segment sammeln
        for i = 1:n_segments
            seg_mask = strcmp(segment_ids, unique_segments{i});
            seg_distances = distances(seg_mask);
            
            bahn_id_col{i} = bahn_id;
            segment_id_col{i} = unique_segments{i};
            min_dist(i) = min(seg_distances);
            max_dist(i) = max(seg_distances);
            avg_dist(i) = mean(seg_distances);
            std_dist(i) = std(seg_distances);
        end
        
        % Gesamtbahn: bahn_id als segment_id
        bahn_id_col{n_segments + 1} = bahn_id;
        segment_id_col{n_segments + 1} = bahn_id;
        min_dist(n_segments + 1) = min(distances);
        max_dist(n_segments + 1) = max(distances);
        avg_dist(n_segments + 1) = mean(distances);
        std_dist(n_segments + 1) = std(distances);
        
        % Tabelle erstellen
        info_table = table(bahn_id_col, segment_id_col, min_dist, max_dist, avg_dist, std_dist);
        info_table.Properties.VariableNames = {'bahn_id', 'segment_id', ...
            'qad_min_distance', 'qad_max_distance', ...
            'qad_average_distance', 'qad_standard_deviation'};
        
    elseif strcmp(metric, 'qdtw')
        qdtw_distances = distances;
        qdtw_soll_x = soll(:,1);
        qdtw_soll_y = soll(:,2);
        qdtw_soll_z = soll(:,3);
        qdtw_soll_w = soll(:,4);
        qdtw_ist_x = ist(:,1);
        qdtw_ist_y = ist(:,2);
        qdtw_ist_z = ist(:,3);
        qdtw_ist_w = ist(:,4);
        
        % distances_table für alle Punkte
        distances_table = table(bahn_ids, segment_ids, qdtw_distances, ...
            qdtw_soll_x, qdtw_soll_y, qdtw_soll_z, qdtw_soll_w, qdtw_ist_x, qdtw_ist_y, qdtw_ist_z, qdtw_ist_w);
        distances_table.Properties.VariableNames = {'bahn_id', 'segment_id', 'qdtw_distances', ...
            'qdtw_soll_x', 'qdtw_soll_y', 'qdtw_soll_z', 'qdtw_soll_w', ...
            'qdtw_ist_x', 'qdtw_ist_y', 'qdtw_ist_z', 'qdtw_ist_w'};
        
        % info_table: eine Zeile pro Segment + eine für Gesamtbahn
        unique_segments = unique(segment_ids, 'stable');
        
        % Preallocate cell arrays für die Tabellendaten
        n_segments = length(unique_segments);
        bahn_id_col = cell(n_segments + 1, 1);
        segment_id_col = cell(n_segments + 1, 1);
        min_dist = zeros(n_segments + 1, 1);
        max_dist = zeros(n_segments + 1, 1);
        avg_dist = zeros(n_segments + 1, 1);
        std_dist = zeros(n_segments + 1, 1);
        
        % Daten für jedes Segment sammeln
        for i = 1:n_segments
            seg_mask = strcmp(segment_ids, unique_segments{i});
            seg_distances = distances(seg_mask);
            
            bahn_id_col{i} = bahn_id;
            segment_id_col{i} = unique_segments{i};
            min_dist(i) = min(seg_distances);
            max_dist(i) = max(seg_distances);
            avg_dist(i) = mean(seg_distances);
            std_dist(i) = std(seg_distances);
        end
        
        % Gesamtbahn: bahn_id als segment_id
        bahn_id_col{n_segments + 1} = bahn_id;
        segment_id_col{n_segments + 1} = bahn_id;
        min_dist(n_segments + 1) = min(distances);
        max_dist(n_segments + 1) = max(distances);
        avg_dist(n_segments + 1) = mean(distances);
        std_dist(n_segments + 1) = std(distances);
        
        % Tabelle erstellen
        info_table = table(bahn_id_col, segment_id_col, min_dist, max_dist, avg_dist, std_dist);
        info_table.Properties.VariableNames = {'bahn_id', 'segment_id', ...
            'qdtw_min_distance', 'qdtw_max_distance', ...
            'qdtw_average_distance', 'qdtw_standard_deviation'};
    end

end
function [info_table, distances_table] = metric2postgresql(metric, distances, soll, ist, bahn_id_, segment_id_)
    % Id's in Cell-Arrays konvertieren und Vektor mit Replica der ID's erstellen
    bahn_id = {bahn_id_};
    bahn_ids = repelem(bahn_id, length(distances), 1);

    % Falls gesamte Bahn 
    if nargin == 6
        segment_id = {convertStringsToChars(segment_id_)};
        segment_ids = repelem(segment_id, length(distances), 1);
    end

    % Wenn euklidischer Abstand
    if strcmp(metric, 'euclidean')
        euclidean_distances = distances;
        euclidean_average_distance = mean(euclidean_distances);
        euclidean_max_distance = max(euclidean_distances);
        euclidean_min_distance = min(euclidean_distances);
        euclidean_standard_deviation = std(euclidean_distances);
        
        ea_soll_x = soll(:,1);
        ea_soll_y = soll(:,2);
        ea_soll_z = soll(:,3);
        ea_ist_x = ist(:,1);
        ea_ist_y = ist(:,2);
        ea_ist_z = ist(:,3);
        
        % Wenn die Auswertung segmentweise erfolgt
        if nargin == 6    
            % Tabellen anlegen
            info_table = table(bahn_id, segment_id, euclidean_min_distance, ...
                              euclidean_max_distance, euclidean_average_distance, ...
                              euclidean_standard_deviation);
            bahn_id = bahn_ids;
            segment_id = segment_ids;
            
            distances_table = table(bahn_id, segment_id, euclidean_distances, ea_soll_x, ea_soll_y, ea_soll_z, ea_ist_x, ea_ist_y, ea_ist_z);
        % Wenn die Auswertung für die gesamte Messaufnahme erfolgt 
        else
            % Tabellen anlegen
            info_table = table(bahn_id, euclidean_min_distance, ...
                              euclidean_max_distance, euclidean_average_distance, ...
                              euclidean_standard_deviation);
            bahn_id = bahn_ids;
            
            distances_table = table(bahn_id, euclidean_distances);
        end

    % Wenn SIDTW
    elseif strcmp(metric, 'sidtw')
        sidtw_min_distance = min(distances);
        sidtw_max_distance = max(distances);
        sidtw_average_distance = mean(distances);
        sidtw_standard_deviation = std(distances);

        sidtw_distances = distances;

        if size(soll, 2) < 3
            sidtw_soll_speed = soll(:,1);
            sidtw_ist_speed = ist(:,1);
            
            % Tabellen anlegen
            if nargin == 6
                info_table = table(bahn_id, segment_id, sidtw_min_distance, ...
                              sidtw_max_distance, sidtw_average_distance, ...
                              sidtw_standard_deviation);
                bahn_id = bahn_ids;
                segment_id = segment_ids; 
        
                distances_table = table(bahn_id, segment_id, sidtw_distances, sidtw_soll_speed, sidtw_ist_speed);
            else
                info_table = table(bahn_id, sidtw_min_distance, ...
                              sidtw_max_distance, sidtw_average_distance, ...
                              sidtw_standard_deviation);
                bahn_id = bahn_ids;
        
                distances_table = table(bahn_id, sidtw_distances, sidtw_soll_speed, sidtw_ist_speed);
            end
        else
            sidtw_soll_x = soll(:,1);
            sidtw_soll_y = soll(:,2);
            sidtw_soll_z = soll(:,3);
            sidtw_ist_x = ist(:,1);
            sidtw_ist_y = ist(:,2);
            sidtw_ist_z = ist(:,3);

            % Tabellen anlegen
            if nargin == 6
                info_table = table(bahn_id, segment_id, sidtw_min_distance, ...
                              sidtw_max_distance, sidtw_average_distance, ...
                              sidtw_standard_deviation);
                bahn_id = bahn_ids;
                segment_id = segment_ids; 
        
                distances_table = table(bahn_id, segment_id, sidtw_distances, ...
                    sidtw_soll_x, sidtw_soll_y, sidtw_soll_z, sidtw_ist_x, sidtw_ist_y, sidtw_ist_z);
            else
                info_table = table(bahn_id, sidtw_min_distance, ...
                              sidtw_max_distance, sidtw_average_distance, ...
                              sidtw_standard_deviation);
                bahn_id = bahn_ids;
        
                distances_table = table(bahn_id, sidtw_distances, ...
                    sidtw_soll_x, sidtw_soll_y, sidtw_soll_z, sidtw_ist_x, sidtw_ist_y, sidtw_ist_z);
            end
        end

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
    end
    
    % Die Tabellen werden jetzt zurückgegeben statt in den Workspace geladen
end
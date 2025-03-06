%% Berechnung der Metriken
function [table_euclidean_info, table_euclidean_deviation, ...
 table_lcss_info, table_lcss_deviation, ...
 table_sidtw_info, table_sidtw_deviation, ...
 table_dtw_info, table_dtw_deviation, ...
 table_dfd_info, table_dfd_deviation] = ...
 evaluationMetrics(bahn_id, segment_ids, num_segments, evaluate_velocity, evaluate_orientation, segments_soll, segments_ist, segments_trafo)

% Tabellen initialisieren
table_sidtw_info = table();
table_sidtw_deviation = cell(num_segments,1);
table_dtw_info = table();
table_dtw_deviation = cell(num_segments,1);
table_dfd_info = table();
table_dfd_deviation = cell(num_segments,1);

% Wenn Geschwindigkeit ausgewertet werden soll sind die 1D-Metriken nicht durchfürbar
if size(segments_soll,2) > 2
    table_euclidean_info = table();
    table_euclidean_deviation = cell(num_segments,1);
    table_lcss_info = table();
    table_lcss_deviation = cell(num_segments,1);
end

% Berechnung der Metriken für alle Segmente
for i = 1:1:num_segments

    if evaluate_velocity == false && evaluate_orientation == false 
        segment_trafo = [segments_trafo.x_ist{i}, segments_trafo.y_ist{i}, segments_trafo.z_ist{i}];
        segment_soll = [segments_soll.x_soll{i}, segments_soll.y_soll{i}, segments_soll.z_soll{i}];
    elseif evaluate_velocity == true && evaluate_orientation == false 
        segment_trafo = segments_ist.tcp_speed_ist{i};
        segment_soll = segments_soll.tcp_speed_soll{i};
    elseif evaluate_velocity == false && evaluate_orientation == true 
        segment_trafo = [segments_trafo.qx_trans{i}, segments_trafo.qy_trans{i},  segments_trafo.qz_trans{i}, segments_trafo.qw_trans{i}];
        segment_soll = [segments_soll.qx_soll{i}, segments_soll.qy_soll{i}, segments_soll.qz_soll{i}, segments_soll.qw_soll{i}];
        
        segment_trafo = fixGimbalLock(rad2deg(quat2eul(segment_trafo)));
        segment_soll = fixGimbalLock(rad2deg(quat2eul(segment_soll)));    
    end

% Anwendung der Methoden
    % Berechnung SIDTW
    [sidtw_distances, ~, ~, ~, sidtw_soll, sidtw_ist, ~, ~, ~] = fkt_selintdtw3d(segment_soll,segment_trafo,false);
    % Berechnung DTW
    [dtw_distances, ~, ~, ~, dtw_soll, dtw_ist, ~, ~, ~, ~] = fkt_dtw3d(segment_soll,segment_trafo,false);
    % Berechnung diskrete Frechet
    [~, ~, frechet_distances, ~, ~, frechet_soll, frechet_ist] = fkt_discreteFrechet(segment_soll,segment_trafo,false);

    % Berechnung nur für mehrdimensionale Größen (z.B. nicht bei Geschwindigkeit)
    if size(segment_soll,2) > 1 
        % Berechnung euklidischer Abstand
        [euclidean_soll,euclidean_distances,~] = distance2curve(segment_soll,segment_trafo, 'linear');
        % Berechnung LCSS
        [~, ~, lcss_distances, ~, ~, lcss_soll, lcss_ist, ~, ~] = fkt_lcss(segment_soll,segment_trafo,false);
    end


% Schreibe die Ergebnisse in die Tabellen 
    if i == 1
        if size(segment_soll,2) > 1
            % Euklidischer Abstand
            [seg_euclidean_info, seg_euclidean_distances] = metric2postgresql('euclidean',euclidean_distances, euclidean_soll, segment_trafo, bahn_id, segment_ids{i,:});
            table_euclidean_info = seg_euclidean_info;
            order_eucl_first = size(seg_euclidean_distances,1);
            seg_euclidean_distances = [seg_euclidean_distances, table((1:1:order_eucl_first)','VariableNames',{'points_order'})];
            table_euclidean_deviation{1} = seg_euclidean_distances;
            % LCSS
            [seg_lcss_info, seg_lcss_distances] = metric2postgresql('lcss',lcss_distances, lcss_soll, lcss_ist, bahn_id, segment_ids{i,:});
            table_lcss_info = seg_lcss_info;
            order_lcss_first = size(seg_lcss_distances,1);
            seg_lcss_distances = [seg_lcss_distances, table((1:1:order_lcss_first)','VariableNames',{'points_order'})];
            table_lcss_deviation{1} = seg_lcss_distances;
        end
        % SIDTW
        [seg_sidtw_info, seg_sidtw_distances] = metric2postgresql('sidtw', sidtw_distances, sidtw_soll, sidtw_ist, bahn_id, segment_ids{i,:});
        table_sidtw_info = seg_sidtw_info;
        order_sidtw_first = size(seg_sidtw_distances,1);
        seg_sidtw_distances = [seg_sidtw_distances, table((1:1:order_sidtw_first)','VariableNames',{'points_order'})];
        table_sidtw_deviation{1} = seg_sidtw_distances;

        % DTW
        [seg_dtw_info, seg_dtw_distances] = metric2postgresql('dtw',dtw_distances, dtw_soll, dtw_ist, bahn_id, segment_ids{i,:});
        table_dtw_info = seg_dtw_info;
        order_dtw_first = size(seg_dtw_distances,1);
        seg_dtw_distances = [seg_dtw_distances, table((1:1:order_dtw_first)','VariableNames',{'points_order'})];
        table_dtw_deviation{1} = seg_dtw_distances;
        % DFD
        [seg_dfd_info, seg_dfd_distances] = metric2postgresql('dfd',frechet_distances, frechet_soll, frechet_ist, bahn_id, segment_ids{i,:});
        table_dfd_info = seg_dfd_info;
        order_dfd_first = size(seg_dfd_distances,1);
        seg_dfd_distances = [seg_dfd_distances, table((1:1:order_dfd_first)','VariableNames',{'points_order'})];
        table_dfd_deviation{1} = seg_dfd_distances;       


    else
        if size(segment_soll,2) > 1
            % Euklidischer Abstand
            [seg_euclidean_info, seg_euclidean_distances] = metric2postgresql('euclidean',euclidean_distances, euclidean_soll, segment_trafo, bahn_id, segment_ids{i,:});
            table_euclidean_info(i,:) = seg_euclidean_info;
            order_eucl_last = order_eucl_first + size(seg_euclidean_distances,1);
            seg_euclidean_distances = [seg_euclidean_distances, table((order_eucl_first+1:1:order_eucl_last)','VariableNames',{'points_order'})];
            order_eucl_first = order_eucl_last;
            table_euclidean_deviation{i} = seg_euclidean_distances;
            % LCSS
            [seg_lcss_info, seg_lcss_distances] = metric2postgresql('lcss',lcss_distances, lcss_soll, lcss_ist, bahn_id, segment_ids{i,:});
            table_lcss_info(i,:) = seg_lcss_info;
            order_lcss_last = order_lcss_first + size(seg_lcss_distances,1);
            seg_lcss_distances = [seg_lcss_distances, table((order_lcss_first+1:1:order_lcss_last)','VariableNames',{'points_order'})];
            order_lcss_first = order_lcss_last;
            table_lcss_deviation{i} = seg_lcss_distances;
        end
        % SIDTW
        [seg_sidtw_info, seg_sidtw_distances] = metric2postgresql('sidtw',sidtw_distances, sidtw_soll, sidtw_ist, bahn_id, segment_ids{i,:});
        table_sidtw_info(i,:) = seg_sidtw_info;
        order_sidtw_last = order_sidtw_first + size(seg_sidtw_distances,1);
        seg_sidtw_distances = [seg_sidtw_distances, table((order_sidtw_first+1:1:order_sidtw_last)','VariableNames',{'points_order'})];
        order_sidtw_first = order_sidtw_last;
        table_sidtw_deviation{i} = seg_sidtw_distances;
        % DTW
        [seg_dtw_info, seg_dtw_distances] = metric2postgresql('dtw',dtw_distances, dtw_soll, dtw_ist, bahn_id, segment_ids{i,:});
        table_dtw_info(i,:) = seg_dtw_info;
        order_dtw_last = order_dtw_first + size(seg_dtw_distances,1);
        seg_dtw_distances = [seg_dtw_distances, table((order_dtw_first+1:1:order_dtw_last)','VariableNames',{'points_order'})];
        order_dtw_first = order_dtw_last;
        table_dtw_deviation{i} = seg_dtw_distances;
        % DFD
        [seg_dfd_info, seg_dfd_distances] = metric2postgresql('dfd',frechet_distances, frechet_soll, frechet_ist, bahn_id, segment_ids{i,:});
        table_dfd_info(i,:) = seg_dfd_info;
        order_dfd_last = order_dfd_first + size(seg_dfd_distances,1);
        seg_dfd_distances = [seg_dfd_distances, table((order_dfd_first+1:1:order_dfd_last)','VariableNames',{'points_order'})];
        order_dfd_first = order_dfd_last;
        table_dfd_deviation{i} = seg_dfd_distances;
    end
end
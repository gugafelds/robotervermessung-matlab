%% Berechnung der Metriken
function [table_euclidean_info, table_euclidean_deviation, ...
    table_lcss_info, table_lcss_deviation, ...
    table_sidtw_info, table_sidtw_deviation, ...
    table_dtw_info, table_dtw_deviation, ...
    table_dfd_info, table_dfd_deviation, table_qad_info, table_qad_deviation, ...
    table_qdtw_info, table_qdtw_deviation] = ...
    calculateMetrics(bahn_id, segment_ids, num_segments, evaluate_orientation, segments_soll, segments_ist, ...
                     use_euclidean, use_sidtw, use_dtw, use_dfd, use_lcss, use_qad, use_qdtw)

% Tabellen initialisieren - nur für aktivierte Methoden
if use_sidtw
    table_sidtw_info = table();
    table_sidtw_deviation = cell(num_segments,1);
else
    table_sidtw_info = table();
    table_sidtw_deviation = {};
end

if use_dtw
    table_dtw_info = table();
    table_dtw_deviation = cell(num_segments,1);
else
    table_dtw_info = table();
    table_dtw_deviation = {};
end

if use_dfd
    table_dfd_info = table();
    table_dfd_deviation = cell(num_segments,1);
else
    table_dfd_info = table();
    table_dfd_deviation = {};
end

if use_euclidean
    table_euclidean_info = table();
    table_euclidean_deviation = cell(num_segments,1);
else
    table_euclidean_info = table();
    table_euclidean_deviation = {};
end

if use_lcss
    table_lcss_info = table();
    table_lcss_deviation = cell(num_segments,1);
else
    table_lcss_info = table();
    table_lcss_deviation = {};
end

if use_qad
    table_qad_info = table();
    table_qad_deviation = cell(num_segments,1);
else
    table_qad_info = table();
    table_qad_deviation = {};
end

if use_qdtw
    table_qdtw_info = table();
    table_qdtw_deviation = cell(num_segments,1);
else
    table_qdtw_info = table();
    table_qdtw_deviation = {};
end


    if evaluate_orientation
        segment_ist = table2array(segments_ist(:, 2:5));
        segment_ids_ist = segments_ist{:,1};
        segment_soll = table2array(segments_soll(:, 2:5));
        segment_ids_soll = segments_soll{:,1};
    else
        segment_ist = table2array(segments_ist(:, 2:4));
        segment_ids_ist = segments_ist{:,1};
        segment_soll = table2array(segments_soll(:, 2:4));
        segment_ids_soll = segments_soll{:,1};
    end

    % Anwendung der Methoden - NUR wenn aktiviert
    if use_sidtw
        [sidtw_distances, ~, ~, ~, sidtw_soll, sidtw_ist, ~, ~, ~, sidtw_segments_order] = fkt_selintdtw3d_opt(segment_soll, segment_ist, segment_ids_soll,segment_ids_ist,false);
    end

    if use_euclidean
        [euclidean_soll, euclidean_distances,~] = distance2curve(segment_soll, segment_ist, 'linear');
        euclidean_segments_order = segment_ids_ist;
    end

    %if use_dtw
    %    [dtw_distances, ~, ~, ~, dtw_soll, dtw_ist, ~, ~, ~, ~] = fkt_dtw3d(segment_soll,segment_ist,false);
    %end
    
    %if use_dfd
    %    [~, ~, frechet_distances, ~, ~, frechet_soll, frechet_ist] = fkt_discreteFrechet(segment_soll,segment_ist,false);
    %end

    %if use_lcss
        %    [~, ~, lcss_distances, ~, ~, lcss_soll, lcss_ist, ~, ~] = fkt_lcss(segment_soll,segment_ist,false);
    %end

    if use_qad
        [qad_soll, qad_distances] = fkt_quaternionAD(segment_soll, segment_ist);
        qad_segments_order = segment_ids_ist;
    end

    if use_qdtw
        [qdtw_distances, qdtw_segments_order, qdtw_soll, qdtw_ist] = fkt_quaternionDTW(segment_soll, segment_ist, segment_ids_ist,false);
    end

    % Schreibe die Ergebnisse in die Tabellen - NUR für aktivierte Methoden

    if use_euclidean
        [seg_euclidean_info, seg_euclidean_distances] = metric2postgresql('euclidean', euclidean_distances, euclidean_soll, segment_ist, bahn_id, euclidean_segments_order);
        table_euclidean_info = seg_euclidean_info;
        order_eucl_first = size(seg_euclidean_distances,1);
        seg_euclidean_distances = [seg_euclidean_distances, table((1:1:order_eucl_first)','VariableNames',{'points_order'})];
        table_euclidean_deviation{1} = seg_euclidean_distances;
    end
    
    %if use_lcss
    %    [seg_lcss_info, seg_lcss_distances] = metric2postgresql('lcss',lcss_distances, lcss_soll, lcss_ist, bahn_id, segment_ids{i,:});
    %    table_lcss_info = seg_lcss_info;
    %    order_lcss_first = size(seg_lcss_distances,1);
    %    seg_lcss_distances = [seg_lcss_distances, table((1:1:order_lcss_first)','VariableNames',{'points_order'})];
    %    table_lcss_deviation{1} = seg_lcss_distances;
    %end
    
    if use_sidtw
        [seg_sidtw_info, seg_sidtw_distances] = metric2postgresql('sidtw', sidtw_distances, sidtw_soll, sidtw_ist, bahn_id, sidtw_segments_order);
        table_sidtw_info = seg_sidtw_info;
        order_sidtw_first = size(seg_sidtw_distances,1);
        seg_sidtw_distances = [seg_sidtw_distances, table((1:1:order_sidtw_first)','VariableNames',{'points_order'})];
        table_sidtw_deviation{1} = seg_sidtw_distances;
    end

    if use_dtw
        [seg_dtw_info, seg_dtw_distances] = metric2postgresql('dtw',dtw_distances, dtw_soll, dtw_ist, bahn_id, segment_ids{i,:});
        table_dtw_info = seg_dtw_info;
        order_dtw_first = size(seg_dtw_distances,1);
        seg_dtw_distances = [seg_dtw_distances, table((1:1:order_dtw_first)','VariableNames',{'points_order'})];
        table_dtw_deviation{1} = seg_dtw_distances;
    end
    
    if use_dfd
        [seg_dfd_info, seg_dfd_distances] = metric2postgresql('dfd',frechet_distances, frechet_soll, frechet_ist, bahn_id, segment_ids{i,:});
        table_dfd_info = seg_dfd_info;
        order_dfd_first = size(seg_dfd_distances,1);
        seg_dfd_distances = [seg_dfd_distances, table((1:1:order_dfd_first)','VariableNames',{'points_order'})];
        table_dfd_deviation{1} = seg_dfd_distances;
    end

    if use_qad
        [seg_qad_info, seg_qad_distances] = metric2postgresql('qad', qad_distances, qad_soll, segment_ist, bahn_id, qad_segments_order);
        table_qad_info = seg_qad_info;
        order_eucl_first = size(seg_qad_distances,1);
        seg_qad_distances = [seg_qad_distances, table((1:1:order_eucl_first)','VariableNames',{'points_order'})];
        table_qad_deviation{1} = seg_qad_distances;
    end

    if use_qdtw
        [seg_qdtw_info, seg_qdtw_distances] = metric2postgresql('qdtw', qdtw_distances, qdtw_soll, qdtw_ist, bahn_id, qdtw_segments_order);
        table_qdtw_info = seg_qdtw_info;
        order_eucl_first = size(seg_qdtw_distances,1);
        seg_qdtw_distances = [seg_qdtw_distances, table((1:1:order_eucl_first)','VariableNames',{'points_order'})];
        table_qdtw_deviation{1} = seg_qdtw_distances;
    end
end
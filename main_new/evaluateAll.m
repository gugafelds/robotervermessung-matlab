function [data_all_soll, data_all_ist_trafo] = evaluateAll(segment_ids, data_ist, data_soll, evaluate_orientation, q_transform, trafo_rot, trafo_trans)
    % Zuerst die segment_ids richtig sortieren (numerisch nach der Zahl nach dem Unterstrich)
    segment_ids_array = table2array(segment_ids);
    segment_ids_numeric = zeros(size(segment_ids_array));
    
    for i = 1:length(segment_ids_array)
        tokens = regexp(segment_ids_array{i}, '_(\d+)$', 'tokens', 'once');
        if ~isempty(tokens)
            segment_ids_numeric(i) = str2double(tokens{1});
        end
    end
    
    [~, sort_idx] = sort(segment_ids_numeric);
    sorted_segment_ids = segment_ids_array(sort_idx);
    
    % Jetzt das erste und letzte Segment basierend auf den sortierten IDs abschneiden
    first_segment = sorted_segment_ids{1};
    last_segment = sorted_segment_ids{end};
    
    % Daten für IST filtern
    first_row_ist = find(data_ist.segment_id == first_segment, 1);
    last_row_ist = find(data_ist.segment_id == last_segment, 1) - 1;
    data_all_ist = data_ist(first_row_ist:last_row_ist,:);
    
    % Daten für SOLL filtern
    first_row_soll = find(data_soll.segment_id == first_segment, 1);
    last_row_soll = find(data_soll.segment_id == last_segment, 1) - 1;
    data_all_soll = data_soll(first_row_soll:last_row_soll,:);

    if evaluate_orientation == true
        q_transformed_all = transformQuaternion(data_all_ist, data_all_soll, q_transform, trafo_rot);
        
        % Quaternion-Transformation für die weitere Verarbeitung verwenden
        data_all_ist = q_transformed_all;
        data_all_ist_trafo = fixGimbalLock(rad2deg(quat2eul(data_all_ist)));
        data_all_soll = [data_all_soll.qw_soll, data_all_soll.qx_soll, data_all_soll.qy_soll, data_all_soll.qz_soll];
        data_all_soll = fixGimbalLock(rad2deg(quat2eul(data_all_soll)));
               
    else 
        % Rest des Codes bleibt gleich...
        data_all_ist = table2array(data_all_ist(:,5:7));
        data_all_soll = table2array(data_all_soll(:,5:7));
    
        % Koordinatentrafo für alle Daten 
        data_all_ist_trafo = coordTransformation(data_all_ist, trafo_rot, trafo_trans);
    end
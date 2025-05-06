function [num_segments, segment_ids, data_ist, data_soll, segments_ist, segments_soll, segments_trafo, q_transformed] = ...
    getSegmentsLT(conn, bahn_id, schema, evaluate_orientation, evaluate_velocity, trafo_rot, trafo_trans, q_transform, positions_transformed, quaternions_transformed)
    
    % Prüfe, ob Lasertracker-Daten vorhanden sind
    is_lasertracker = nargin > 8 && ~isempty(positions_transformed);
    
    % Anzahl der Segmente der gesamten Messaufnahme bestimmen 
    query = ['SELECT * FROM robotervermessung.' schema '.bahn_info ' ...
             'WHERE robotervermessung.' schema '.bahn_info.bahn_id = ''' bahn_id ''''];
    data_info = fetch(conn, query);
    num_segments = data_info.np_ereignisse;

    % Orientierungsdaten
    if evaluate_velocity == false && evaluate_orientation == true
        % Auslesen der gesamten Ist-Daten
        query = ['SELECT * FROM robotervermessung.' schema '.bahn_pose_ist ' ...
                'WHERE robotervermessung.' schema '.bahn_pose_ist.bahn_id = ''' bahn_id ''''];
        data_ist = fetch(conn, query);
        data_ist = sortrows(data_ist,'timestamp');
        
        % Auslesen der gesamten Soll-Daten
        query = ['SELECT * FROM robotervermessung.' schema '.bahn_orientation_soll ' ...
                'WHERE robotervermessung.' schema '.bahn_orientation_soll.bahn_id = ''' bahn_id ''''];
        data_soll = fetch(conn, query);
        data_soll = sortrows(data_soll,'timestamp');
        
        q_ist = table2array(data_ist(:,8:11));
        q_ist = [q_ist(:,4), q_ist(:,1), q_ist(:,2), q_ist(:,3)];
        euler_ist = quat2eul(q_ist,"XYZ");
        euler_ist = rad2deg(euler_ist);

        position_ist = table2array(data_ist(:,5:7));
        
        %%%%% Nur notwendig für plots 
        % Koordinatentransfromation für MoCap-Daten (nicht Lasertracker)
        if ~is_lasertracker
            coordTransformation(position_ist, trafo_rot, trafo_trans);
            % Winkeltransformation
            q_transformed = transformQuaternion(data_ist, data_soll, q_transform, trafo_rot);
        else
            % Verwende bereits transformierte Quaternionen für Lasertracker
            q_transformed = quaternions_transformed;
        end
        %%%%%
    % Geschwindigkeitsdaten aus Positionsdaten
    elseif evaluate_velocity == true && evaluate_orientation == false 
        % Auslesen der gesamten Ist-Daten
        query = ['SELECT * FROM robotervermessung.' schema '.bahn_twist_ist ' ...
                'WHERE robotervermessung.' schema '.bahn_twist_ist.bahn_id = ''' bahn_id ''''];
        data_ist = fetch(conn, query);
        data_ist = sortrows(data_ist,'timestamp');
        
        % Auslesen der gesamten Soll-Daten
        query = ['SELECT * FROM robotervermessung.' schema '.bahn_position_soll ' ...
                'WHERE robotervermessung.' schema '.bahn_position_soll.bahn_id = ''' bahn_id ''''];
        data_soll = fetch(conn, query);
        data_soll = sortrows(data_soll,'timestamp');

        % Geschwindigkeitsdaten präperieren 
        velocityPreparation(data_soll, data_ist);
    % Positionsdaten
    else
        % Auslesen der gesamten Ist-Daten
        query = ['SELECT * FROM robotervermessung.' schema '.bahn_pose_ist ' ...
                'WHERE robotervermessung.' schema '.bahn_pose_ist.bahn_id = ''' bahn_id ''''];
        data_ist = fetch(conn, query);
        data_ist = sortrows(data_ist,'timestamp');
        
        % Auslesen der gesamten Soll-Daten
        query = ['SELECT * FROM robotervermessung.' schema '.bahn_position_soll ' ...
                'WHERE robotervermessung.' schema '.bahn_position_soll.bahn_id = ''' bahn_id ''''];
        data_soll = fetch(conn, query);
        data_soll = sortrows(data_soll,'timestamp');
    end

    %% Extraktion und Separation der Segmente der Gesamtaufname
    % Alle Segment-ID's 
    query = ['SELECT segment_id FROM robotervermessung.' schema '.bahn_events ' ...
        'WHERE robotervermessung.' schema '.bahn_events.bahn_id = ''' bahn_id ''''];
    segment_ids = fetch(conn,query);

    % % % IST-DATEN % % %
    % Extraktion der Indizes der Segmente 
    seg_id = split(data_ist.segment_id, '_');
    seg_id = double(string(seg_id(:,2)));
    idx_new_seg_ist = zeros(num_segments,1);

    % Suche nach den Indizes bei denen sich die Segmentnr. ändert
    k = 0;
    idx = 1;
    for i = 1:1:length(seg_id)
        if seg_id(i) == k
            idx = idx + 1;
        else
            k = k +1;
            idx_new_seg_ist(k) = idx;
            idx = idx+1;
        end
    end

    % % % SOLL-DATEN % % %
    seg_id = split(data_soll.segment_id, '_');
    seg_id = double(string(seg_id(:,2)));
    idx_new_seg_soll = zeros(num_segments,1);

    k = 0;
    idx = 1;
    for i = 1:1:length(seg_id)
        if seg_id(i) == k
            idx = idx + 1;
        else
            k = k +1;
            idx_new_seg_soll(k) = idx;
            idx = idx+1;
        end
    end

    % Speichern der einzelnen Segmente in einer Tabelle 
    %%%%%%%%%%%%%%%%%
    % Geschwindigkeit
    if evaluate_velocity == true && evaluate_orientation == false 
        disp('Es wird die Geschwindigkeit ausgewertet!')
        % [ursprüngliche Geschwindigkeitscode bleibt unverändert...]

    %%%%%%%%%%%%%%
    % Orientierung
    elseif evaluate_velocity == false && evaluate_orientation == true
        disp('Es wird die Orientierung ausgewertet!')
        
        % Erstes Segment Ist (quaternions)
        segments_ist = array2table([{data_ist.segment_id(1)} ...
                                  data_ist.qw_ist(1:idx_new_seg_ist(1)-1) ...
                                  data_ist.qx_ist(1:idx_new_seg_ist(1)-1) ...
                                  data_ist.qy_ist(1:idx_new_seg_ist(1)-1) ...
                                  data_ist.qz_ist(1:idx_new_seg_ist(1)-1)], ...
                                  'VariableNames', {'segment_id', 'qw_ist', 'qx_ist', 'qy_ist', 'qz_ist'});
        
        % Übrige Segmente Ist
        for i = 1:num_segments
            if i == length(idx_new_seg_ist)
                segments_ist(i+1,:) = array2table([{segment_ids{i,:}} ...
                                                 data_ist.qw_ist(idx_new_seg_ist(i):end) ...
                                                 data_ist.qx_ist(idx_new_seg_ist(i):end) ...
                                                 data_ist.qy_ist(idx_new_seg_ist(i):end) ...
                                                 data_ist.qz_ist(idx_new_seg_ist(i):end)]);
            else
                segments_ist(i+1,:) = array2table([{segment_ids{i,:}} ...
                                                 data_ist.qw_ist(idx_new_seg_ist(i):idx_new_seg_ist(i+1)-1) ...
                                                 data_ist.qx_ist(idx_new_seg_ist(i):idx_new_seg_ist(i+1)-1) ...
                                                 data_ist.qy_ist(idx_new_seg_ist(i):idx_new_seg_ist(i+1)-1) ...
                                                 data_ist.qz_ist(idx_new_seg_ist(i):idx_new_seg_ist(i+1)-1)]);
            end
        end
        
        % Erstes Segment Soll
        segments_soll = array2table([{data_soll.segment_id(1)} ...
                                    data_soll.qw_soll(1:idx_new_seg_soll(1)-1) ...
                                    data_soll.qx_soll(1:idx_new_seg_soll(1)-1) ...
                                    data_soll.qy_soll(1:idx_new_seg_soll(1)-1) ...
                                    data_soll.qz_soll(1:idx_new_seg_soll(1)-1)], ...
                                    'VariableNames', {'segment_id', 'qw_soll', 'qx_soll', 'qy_soll', 'qz_soll'});
        
        % Übrige Semgente Soll
        for i = 1:num_segments
            if i == length(idx_new_seg_soll)
                segments_soll(i+1,:) = array2table([{segment_ids{i,:}} ...
                                                  data_soll.qw_soll(idx_new_seg_soll(i):end) ...
                                                  data_soll.qx_soll(idx_new_seg_soll(i):end) ...
                                                  data_soll.qy_soll(idx_new_seg_soll(i):end) ...
                                                  data_soll.qz_soll(idx_new_seg_soll(i):end)]);
            else
                segments_soll(i+1,:) = array2table([{segment_ids{i,:}} ...
                                                  data_soll.qw_soll(idx_new_seg_soll(i):idx_new_seg_soll(i+1)-1) ...
                                                  data_soll.qx_soll(idx_new_seg_soll(i):idx_new_seg_soll(i+1)-1) ...
                                                  data_soll.qy_soll(idx_new_seg_soll(i):idx_new_seg_soll(i+1)-1) ...
                                                  data_soll.qz_soll(idx_new_seg_soll(i):idx_new_seg_soll(i+1)-1)]);
            end
        end
        
        % Initialisierung Quarternion Transformation
        segments_trafo = table();
        q_transformed_all = [];
        
        % Segmentweise Schleife für Transformation
        if ~is_lasertracker
            % MoCap-Transformation (bestehender Code)
            for i = 1:num_segments+1
                % Extrahiere Quarternion für aktuelles Segment
                segment_ist = table2struct(segments_ist(i,:));
                segment_soll = table2struct(segments_soll(i,:));
                
                % Temporäre Tabelle mit Quarternion Daten
                data_ist_seg = table(segment_ist.qw_ist, segment_ist.qx_ist, segment_ist.qy_ist, segment_ist.qz_ist, ...
                                    'VariableNames', {'qw_ist', 'qx_ist', 'qy_ist', 'qz_ist'});
                data_soll_seg = table(segment_soll.qw_soll, segment_soll.qx_soll, segment_soll.qy_soll, segment_soll.qz_soll, ...
                                     'VariableNames', {'qw_soll', 'qx_soll', 'qy_soll', 'qz_soll'});
                
                % Quarternion-Transformation für aktuelles Segment
                q_transformed_seg = transformQuaternion(data_ist_seg, data_soll_seg, q_transform, trafo_rot);

                % Transformationsdaten der Tabelle hinzufügen
                segments_trafo(i,:) = table({segments_ist.segment_id(i)}, ...
                                       {q_transformed_seg(:,1)}, {q_transformed_seg(:,2)}, ...
                                       {q_transformed_seg(:,3)}, {q_transformed_seg(:,4)}, ...
                                       'VariableNames', {'segment_id', 'qw_trans', 'qx_trans', 'qy_trans', 'qz_trans'});
            
                % Hinzufügen der transformierten Quarternion
                q_transformed_all = [q_transformed_all; q_transformed_seg];
            end
        else
            % Lasertracker-Transformation - verwende bereits transformierte Quaternionen
            for i = 1:num_segments+1
                segment_start = idx_new_seg_ist(i);
                if i == length(idx_new_seg_ist)
                    segment_end = length(quaternions_transformed);
                else
                    segment_end = idx_new_seg_ist(i+1)-1;
                end
                
                % Extrahiere transformierte Quaternionen für das aktuelle Segment
                q_segment = quaternions_transformed(segment_start:segment_end, :);
                
                % Transformationsdaten der Tabelle hinzufügen
                segments_trafo(i,:) = table({segments_ist.segment_id(i)}, ...
                                       {q_segment(:,1)}, {q_segment(:,2)}, ...
                                       {q_segment(:,3)}, {q_segment(:,4)}, ...
                                       'VariableNames', {'segment_id', 'qw_trans', 'qx_trans', 'qy_trans', 'qz_trans'});
                
                % Hinzufügen der transformierten Quaternionen
                q_transformed_all = [q_transformed_all; q_segment];
            end
        end

    %%%%%%%%%%%%%%%%    
    % Positionsdaten 
    else
        disp('Es wird die Position ausgewertet!')

        % Erstes Segment Ist
        segments_ist = array2table([{data_ist.segment_id(1)} data_ist.x_ist(1:idx_new_seg_ist(1)-1) data_ist.y_ist(1:idx_new_seg_ist(1)-1) data_ist.z_ist(1:idx_new_seg_ist(1)-1)], "VariableNames",{'segment_id','x_ist','y_ist','z_ist'});
        
        % Übrige Segmente Ist
        for i = 1:num_segments
            if i == length(idx_new_seg_ist)
                segments_ist(i+1,:) = array2table([{segment_ids{i,:}} data_ist.x_ist(idx_new_seg_ist(i):end) data_ist.y_ist(idx_new_seg_ist(i):end) data_ist.z_ist(idx_new_seg_ist(i):end)]);
            else
                segments_ist(i+1,:) = array2table([{segment_ids{i,:}} data_ist.x_ist(idx_new_seg_ist(i):idx_new_seg_ist(i+1)-1) data_ist.y_ist(idx_new_seg_ist(i):idx_new_seg_ist(i+1)-1) data_ist.z_ist(idx_new_seg_ist(i):idx_new_seg_ist(i+1)-1)]);
            end
        end
        
        % Erstes Segment Soll
        if idx_new_seg_soll(1) == 1
            segments_soll = array2table([{data_soll.segment_id(1)} data_soll.x_soll(1:idx_new_seg_soll(1)) data_soll.y_soll(1:idx_new_seg_soll(1)) data_soll.z_soll(1:idx_new_seg_soll(1))], "VariableNames",{'segment_id','x_soll','y_soll','z_soll'});
        else
            segments_soll = array2table([{data_soll.segment_id(1)} data_soll.x_soll(1:idx_new_seg_soll(1)-1) data_soll.y_soll(1:idx_new_seg_soll(1)-1) data_soll.z_soll(1:idx_new_seg_soll(1)-1)], "VariableNames",{'segment_id','x_soll','y_soll','z_soll'});
        end

        % Übrige Segmente Soll
        for i = 1:num_segments
            if i == length(idx_new_seg_soll)
                segments_soll(i+1,:) = array2table([{segment_ids{i,:}} data_soll.x_soll(idx_new_seg_soll(i):end) data_soll.y_soll(idx_new_seg_soll(i):end) data_soll.z_soll(idx_new_seg_soll(i):end)]);
            else
                segments_soll(i+1,:)= array2table([{segment_ids{i,:}} data_soll.x_soll(idx_new_seg_soll(i):idx_new_seg_soll(i+1)-1) data_soll.y_soll(idx_new_seg_soll(i):idx_new_seg_soll(i+1)-1) data_soll.z_soll(idx_new_seg_soll(i):idx_new_seg_soll(i+1)-1)]);
            end    
        end
        
        % Koordinatentransformation für alle Segemente
        segments_trafo = table();
        if ~is_lasertracker
            % MoCap-Transformation
            for i = 1:1:num_segments+1
                segments_trafo(i,:) = coordTransformation(segments_ist(i,:), trafo_rot, trafo_trans);
            end
        else
            % Lasertracker-Transformation - verwende bereits transformierte Positionen
            for i = 1:1:num_segments+1
                if i == 1
                    segment_indices = 1:idx_new_seg_ist(1)-1;
                elseif i == length(idx_new_seg_ist)+1
                    segment_indices = idx_new_seg_ist(i-1):length(positions_transformed);
                else
                    segment_indices = idx_new_seg_ist(i-1):idx_new_seg_ist(i)-1;
                end
                
                % Erstelle Tabelle mit transformierten Koordinaten
                segments_trafo(i,:) = table({segments_ist.segment_id(i)}, ...
                                         {positions_transformed(segment_indices,1)}, ...
                                         {positions_transformed(segment_indices,2)}, ...
                                         {positions_transformed(segment_indices,3)}, ...
                                         'VariableNames', {'segment_id', 'x_ist', 'y_ist', 'z_ist'});
            end
        end
    end

    % Löschen des Segment 0
    segments_soll = segments_soll(2:end,:);
    segments_ist = segments_ist(2:end,:);
    if evaluate_velocity == false
        segments_trafo = segments_trafo(2:end,:);
    end
    num_segments = num_segments - 1;

    if evaluate_orientation == true
        q_transformed = q_transformed_all;
    else
        q_transformed = 0; % Wird dann sowieso nicht ausgegeben!
    end
end
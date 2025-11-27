function segments = batchLoadSegmentsChunked(conn, schema, segment_ids, mode, chunk_size)
    % Batch load segments in chunks
    % 
    % Modes: 'position', 'joint_states', 'orientation'
    
    if nargin < 5
        chunk_size = 200;
    end
    
    num_total = length(segment_ids);
    segments = cell(num_total, 1);
    num_chunks = ceil(num_total / chunk_size);
    
    fprintf('    Loading %d segment candidates (%s) in %d chunks (chunk size: %d)...\n', ...
        num_total, mode, num_chunks, chunk_size);
    
    for chunk_idx = 1:num_chunks
        start_idx = (chunk_idx - 1) * chunk_size + 1;
        end_idx = min(chunk_idx * chunk_size, num_total);
        chunk_ids = segment_ids(start_idx:end_idx);
        
        id_list = sprintf('''%s''', strjoin(chunk_ids, ''','''));
        
        % ⭐ Build query based on mode
        if strcmp(mode, 'position')
            query = sprintf(...
                ['SELECT segment_id, x_soll, y_soll, z_soll, timestamp ' ...
                 'FROM robotervermessung.%s.bahn_position_soll ' ...
                 'WHERE segment_id IN (%s) ' ...
                 'ORDER BY segment_id, timestamp'], ...
                schema, id_list);
                
        elseif strcmp(mode, 'joint_states')
            query = sprintf(...
                ['SELECT segment_id, joint_1, joint_2, joint_3, joint_4, joint_5, joint_6, timestamp ' ...
                 'FROM robotervermessung.%s.bahn_joint_states ' ...
                 'WHERE segment_id IN (%s) ' ...
                 'ORDER BY segment_id, timestamp'], ...
                schema, id_list);
                
        elseif strcmp(mode, 'orientation')
            % ⭐ NEW: Orientation mode
            query = sprintf(...
                ['SELECT segment_id, qx_soll, qy_soll, qz_soll, qw_soll, timestamp ' ...
                 'FROM robotervermessung.%s.bahn_orientation_soll ' ...
                 'WHERE segment_id IN (%s) ' ...
                 'ORDER BY segment_id, timestamp'], ...
                schema, id_list);
        else
            error('Invalid mode: %s. Must be ''position'', ''joint_states'', or ''orientation''', mode);
        end
        
        result = fetch(conn, query);
        
        % ⭐ Extract data based on mode
        for i = 1:length(chunk_ids)
            segment_id = chunk_ids{i};
            rows = strcmp(result.segment_id, segment_id);
            
            if strcmp(mode, 'position')
                seg_data = [result.x_soll(rows), result.y_soll(rows), result.z_soll(rows)];
                
            elseif strcmp(mode, 'joint_states')
                seg_data = [result.joint_1(rows), result.joint_2(rows), result.joint_3(rows), ...
                            result.joint_4(rows), result.joint_5(rows), result.joint_6(rows)];
                            
            elseif strcmp(mode, 'orientation')
                % ⭐ NEW: Orientation data [qx, qy, qz, qw]
                seg_data = [result.qx_soll(rows), result.qy_soll(rows), ...
                            result.qz_soll(rows), result.qw_soll(rows)];
            end
            
            global_idx = start_idx + i - 1;
            segments{global_idx} = seg_data;
        end
        
        fprintf('      Chunk %d/%d: Loaded %d-%d (%.1f%% complete)\n', ...
            chunk_idx, num_chunks, start_idx, end_idx, 100 * end_idx / num_total);
    end
    
    fprintf('    ✓ All %d segment candidates loaded\n', num_total);
end
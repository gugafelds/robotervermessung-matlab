function trajectories = batchLoadTrajectoriesChunked(conn, schema, ids, mode, chunk_size)
    % Batch load trajectories/segments (chunked)
    % 
    % Input: 
    %   ids - cell array von bahn_id ODER segment_id
    %   mode - 'position', 'joint_states', oder 'orientation'
    %   chunk_size - number of IDs per chunk
    
    num_total = length(ids);
    num_chunks = ceil(num_total / chunk_size);
    trajectories = cell(num_total, 1);
    
    for chunk_idx = 1:num_chunks
        start_idx = (chunk_idx - 1) * chunk_size + 1;
        end_idx = min(chunk_idx * chunk_size, num_total);
        chunk_ids = ids(start_idx:end_idx);
        
        id_list = sprintf('''%s''', strjoin(chunk_ids, ''','''));
        
        % ⭐ Build query based on mode
        if strcmp(mode, 'position')
            query = sprintf(['SELECT COALESCE(bahn_id, segment_id) as id, ' ...
                             'x_soll, y_soll, z_soll, timestamp ' ...
                             'FROM %s.bahn_position_soll ' ...
                             'WHERE bahn_id IN (%s) OR segment_id IN (%s) ' ...
                             'ORDER BY id, timestamp'], schema, id_list, id_list);
                             
        elseif strcmp(mode, 'joint_states')
            query = sprintf(['SELECT COALESCE(bahn_id, segment_id) as id, ' ...
                             'joint_1, joint_2, joint_3, ' ...
                             'joint_4, joint_5, joint_6, timestamp ' ...
                             'FROM %s.bahn_joint_states ' ...
                             'WHERE bahn_id IN (%s) OR segment_id IN (%s) ' ...
                             'ORDER BY id, timestamp'], schema, id_list, id_list);
                             
        elseif strcmp(mode, 'orientation')
            % ⭐ NEW: Orientation mode
            query = sprintf(['SELECT COALESCE(bahn_id, segment_id) as id, ' ...
                             'qx_soll, qy_soll, qz_soll, qw_soll, timestamp ' ...
                             'FROM %s.bahn_orientation_soll ' ...
                             'WHERE bahn_id IN (%s) OR segment_id IN (%s) ' ...
                             'ORDER BY id, timestamp'], schema, id_list, id_list);
        else
            error('Invalid mode: %s. Must be ''position'', ''joint_states'', or ''orientation''', mode);
        end
        
        result = fetch(conn, query);
        
        % ⭐ Split by id and extract data
        for i = 1:length(chunk_ids)
            idx_in_result = strcmp(result.id, chunk_ids{i});
            
            if any(idx_in_result)
                if strcmp(mode, 'position')
                    trajectories{start_idx + i - 1} = [result.x_soll(idx_in_result), ...
                                                        result.y_soll(idx_in_result), ...
                                                        result.z_soll(idx_in_result)];
                                                        
                elseif strcmp(mode, 'joint_states')
                    trajectories{start_idx + i - 1} = [result.joint_1(idx_in_result), ...
                                                        result.joint_2(idx_in_result), ...
                                                        result.joint_3(idx_in_result), ...
                                                        result.joint_4(idx_in_result), ...
                                                        result.joint_5(idx_in_result), ...
                                                        result.joint_6(idx_in_result)];
                                                        
                elseif strcmp(mode, 'orientation')
                    % ⭐ NEW: Orientation data [qx, qy, qz, qw]
                    trajectories{start_idx + i - 1} = [result.qx_soll(idx_in_result), ...
                                                        result.qy_soll(idx_in_result), ...
                                                        result.qz_soll(idx_in_result), ...
                                                        result.qw_soll(idx_in_result)];
                end
            end
        end
    end
end
function data = loadTrajectoryJointStates(conn, schema, id)
    % Lädt Joint States (6 Gelenke) für Trajektorie ODER Segment
    % Input: id kann bahn_id ODER segment_id sein
    
    query = sprintf(['SELECT joint_1, joint_2, joint_3, ' ...
                    'joint_4, joint_5, joint_6 ' ...
                    'FROM %s.bahn_joint_states ' ...
                    'WHERE bahn_id = ''%s'' OR segment_id = ''%s'' ' ...
                    'ORDER BY timestamp'], schema, id, id);
    
    result = fetch(conn, query);
    
    if isempty(result)
        data = [];
        return;
    end
    
    data = [result.joint_1, result.joint_2, result.joint_3, ...
            result.joint_4, result.joint_5, result.joint_6];
end
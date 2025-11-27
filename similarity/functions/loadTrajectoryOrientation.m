function data = loadTrajectoryOrientation(conn, schema, id)
    % Lädt Orientation (Quaternion) für Trajektorie ODER Segment
    % Input: id kann bahn_id ODER segment_id sein
    
    query = sprintf(['SELECT qw_soll, qx_soll, qy_soll, qz_soll ' ...
                    'FROM %s.bahn_orientation_soll ' ...
                    'WHERE bahn_id = ''%s'' OR segment_id = ''%s'' ' ...
                    'ORDER BY timestamp'], schema, id, id);
    
    result = fetch(conn, query);
    
    if isempty(result)
        data = [];
        return;
    end
    
    data = [result.qw_soll, result.qx_soll, result.qy_soll, result.qz_soll];
end
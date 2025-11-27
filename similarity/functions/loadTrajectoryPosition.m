function data = loadTrajectoryPosition(conn, schema, id)
    % Lädt Position (x, y, z) für Trajektorie ODER Segment
    % Automatische Erkennung: '_' im id → segment_id, sonst bahn_id
    
    if contains(id, '_')
        % Segment
        query = sprintf(['SELECT x_soll, y_soll, z_soll ' ...
                        'FROM %s.bahn_position_soll ' ...
                        'WHERE segment_id = ''%s'' ' ...
                        'ORDER BY timestamp'], schema, id);
    else
        % Trajektorie
        query = sprintf(['SELECT x_soll, y_soll, z_soll ' ...
                        'FROM %s.bahn_position_soll ' ...
                        'WHERE bahn_id = ''%s'' ' ...
                        'ORDER BY timestamp'], schema, id);
    end
    
    result = fetch(conn, query);
    
    if isempty(result)
        data = [];
        return;
    end
    
    data = [result.x_soll, result.y_soll, result.z_soll];
end
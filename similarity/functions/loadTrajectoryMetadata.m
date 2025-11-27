function metadata = loadTrajectoryMetadata(conn, schema, bahn_id)
    % Load metadata for trajectory (only the row where bahn_id = segment_id)
    query = sprintf(...
        ['SELECT length, duration ' ...
         'FROM robotervermessung.%s.bahn_metadata ' ...
         'WHERE bahn_id = ''%s'' AND segment_id = ''%s'''], ...    % <-- bahn_id = segment_id
        schema, bahn_id, bahn_id);
    
    result = fetch(conn, query);
    
    metadata = struct();
    metadata.length = result.length;
    metadata.duration = result.duration;
    
    % Count segments: number of rows with same bahn_id
    query_segments = sprintf(...
        ['SELECT COUNT(*) as num_segments ' ...
         'FROM robotervermessung.%s.bahn_metadata ' ...
         'WHERE bahn_id = ''%s'' AND segment_id != bahn_id'], ...  % <-- Nur Segment-Zeilen zählen
        schema, bahn_id);
    
    result_segments = fetch(conn, query_segments);
    metadata.num_segments = result_segments.num_segments;
end
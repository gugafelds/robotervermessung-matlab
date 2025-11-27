function metadata = loadSegmentMetadata(conn, schema, segment_id)
    % Load metadata for a specific segment
    query = sprintf(...
        ['SELECT length, duration ' ...
         'FROM robotervermessung.%s.bahn_metadata ' ...
         'WHERE segment_id = ''%s'''], ...
        schema, segment_id);
    
    result = fetch(conn, query);
    
    metadata = struct();
    metadata.length = result.length;
    metadata.duration = result.duration;
end
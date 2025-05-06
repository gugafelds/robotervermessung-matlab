function [bahn_ids, existing_bahn_ids, evaluation_type] = getBahnIds(conn,evaluate_orientation,evaluate_velocity)

    query = 'SELECT bahn_id FROM robotervermessung.bewegungsdaten.bahn_info';
    bahn_ids = fetch(conn, query);
    bahn_ids = double(string(table2array(bahn_ids)));
    
    % Nur nach 10stelligen BahnIds suchen
    log_idx = floor(log10(bahn_ids)) == 9;
    bahn_ids = bahn_ids(log_idx);
    
    % Check ob Orientierung, Geschwindigkeit oder Position ausgewertet wird 
    if evaluate_orientation == true && evaluate_velocity == false 

        evaluation_type = 'orientation';
        query = sprintf("SELECT bahn_id FROM robotervermessung.auswertung.info_sidtw WHERE evaluation = '%s'",evaluation_type);
        % query = "SELECT bahn_id FROM robotervermessung.auswertung.orientation_sidtw";
    end
    if evaluate_velocity == true && evaluate_orientation == false
        evaluation_type = 'speed';
        query = sprintf("SELECT bahn_id FROM robotervermessung.auswertung.info_sidtw WHERE evaluation = '%s'",evaluation_type);
        % query = "SELECT bahn_id FROM robotervermessung.auswertung.speed_sidtw";
    end
    if evaluate_orientation == false && evaluate_velocity == false
        evaluation_type = 'position';
        query = sprintf("SELECT bahn_id FROM robotervermessung.auswertung.info_sidtw WHERE evaluation = '%s'",evaluation_type);
        % query = "SELECT bahn_id FROM robotervermessung.auswertung.position_sidtw"; % ... Dauert extrem lange
        
    end
    
    % Extrahieren aller Bahn-ID's die in der Datenbank vorliegen
    existing_bahn_ids = fetch(conn,query);
    existing_bahn_ids = double(string(table2array(existing_bahn_ids)));
    existing_bahn_ids = unique(existing_bahn_ids);

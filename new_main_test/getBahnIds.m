function [bahn_ids, existing_bahn_ids] = getBahnIds(conn,evaluate_orientation,evaluate_velocity)

    query = 'SELECT bahn_id FROM robotervermessung.bewegungsdaten.bahn_info';
    bahn_ids = fetch(conn, query);
    bahn_ids = double(string(table2array(bahn_ids)));
    
    % Nur nach 10stelligen BahnIds suchen
    log_idx = floor(log10(bahn_ids)) == 9;
    bahn_ids = bahn_ids(log_idx);
    
    % Check ob Orientierung, Geschwindigkeit oder Position ausgewertet wird 
    if evaluate_orientation == true && evaluate_velocity == false 
        query = "SELECT bahn_id FROM robotervermessung.auswertung.orientation_sidtw";
        % query = "SELECT bahn_id FROM robotervermessung.auswertung.info_euclidean WHERE evaluation = 'position'";
    end
    if evaluate_velocity == true && evaluate_orientation == false
        query = "SELECT bahn_id FROM robotervermessung.auswertung.speed_sidtw";
    end
    if evaluate_orientation == false && evaluate_velocity == false
        query = "SELECT bahn_id FROM robotervermessung.auswertung.info_sidtw";
        % query = "SELECT bahn_id FROM robotervermessung.auswertung.position_sidtw"; % Dauert extrem lange
    end
    
    % Extrahieren aller Bahn-ID's die in der Datenbank vorliegen
    existing_bahn_ids = fetch(conn,query);
    existing_bahn_ids = double(string(table2array(existing_bahn_ids)));
    existing_bahn_ids = unique(existing_bahn_ids);

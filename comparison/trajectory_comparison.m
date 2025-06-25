% Löschen des Workspace
clear;
tic;
% Verbindung mit der Datenbank
conn = connecting_to_postgres;

%% Bestimmung der Änlichkeitsparameter anhand der gegebenen Punkte/Ereignisse
% Festlegen nach einer Bahn für die änliche Bahnen gesucht werden sollen

bahn_id = '1719911918'; % Square mässig
bahn_id = '1721048209'; % andere Routine
bahn_id = '1739874821';

query = sprintf('SELECT * FROM robotervermessung.bewegungsdaten.bahn_events WHERE bahn_id = ''%s''',bahn_id);
bahn_events = fetch(conn, query);

% ggf. hier auch selbst Punkte/Trajektorie vorgeben können
% Alle Berechnungen müssen später über die bahn_ident laufen dann!

% Punkte/Events herrausfiltern
events_position = table2array(bahn_events(:,5:7));
events_orientation = table2array(bahn_events(:,8:11));

% Richtungsvektor, norm. Richtungsvektor und Länge der Bahnen berechnen
direction = zeros(height(bahn_events)-1,3);
distance = zeros(height(bahn_events)-1,1);
ndirection = zeros(height(bahn_events)-1,3);
quaternion_diff = zeros(height(bahn_events)-1,4);
for i = 1: height(bahn_events)-1
    direction(i,:) = events_position(i+1,:) - events_position(i,:);
    distance(i) = vecnorm((direction(i,:)));
    ndirection(i,:) = direction(i,:)./distance(i);
    quaternion_diff(i,:) = events_orientation(i+1,:) - events_orientation(i,:);
end

% Identifizierungsvektor über Länge und Richtung der Bahnabschnitte berechnen
bahn_ident = zeros(8, length(distance));
bahn_ident(1:3,:) = ndirection';
bahn_ident(4,:) = distance;
bahn_ident(5:8,:) = quaternion_diff';

clear direction distance ndirection quaternion_diff i query

%% Datenbank durchsuchen um mögliche Bahnen für den Vergleich zu identifizieren
query = sprintf('SELECT * FROM robotervermessung.bewegungsdaten.bahn_events');
all_events = fetch(conn, query);
all_events = sortrows(all_events,"timestamp");

% Tabelle nach Bahn-Ids filtern, die 10 Zeichen enthalten (neue gesplittete Daten)
if exist('all_events', 'var') && istable(all_events) && any(strcmp('bahn_id', all_events.Properties.VariableNames))
    % Filterkriterium: Länge der Strings in der Spalte 'bahn_id'
    isValidLength = cellfun(@(x) length(x) == 10, all_events.bahn_id); % isValidLength = strlength(all_events.bahn_id) == 10;
    
    % Nur Zeilen mit Strings von 10 Zeichen beibehalten
    all_events_filtered = all_events(isValidLength, :);   
else
    error('Die Tabelle "all_events" oder die Spalte "bahn_id" existiert nicht.');
end

% Extrahieren der Spalte 'bahn_id' der gefilteretden Daten
all_bahn_id_array = all_events_filtered.bahn_id;

% Finden der einzigartigen Einträge und deren Häufigkeit
[uniqueEntries, ~, indices] = unique(all_bahn_id_array, 'stable');
% Zählen der Häufigkeit jedes Eintrags
counts = accumarray(indices, 1); 

all_bahn_ids = table(uniqueEntries, counts, 'VariableNames', {'Bahn_ID', 'Häufigkeit'});
num_events = sum(counts);

clear all_bahn_id_array uniqueEntries indices isValidLength query
clear counts all_events

% Initialisieren des Merkmal-Arrays aller Bahnen
all_ident = zeros(8,num_events);

dir = zeros(num_events,3);
dist = zeros(num_events,1);
ndir = zeros(num_events,3);
quat = zeros(num_events,4);

% Alle Positionen und Orientierungen
all_positions = table2array(all_events_filtered(:,5:7));
all_orientations = table2array(all_events_filtered(:,8:11));

% Finden der Bahnindizes, welche min. genau so viele Positionen wie die gesuchte Bahn enthalten
bahn_ids_idx = find(all_bahn_ids.("Häufigkeit") >= size(bahn_ident,2)+1);

% Berechnung der Merkmalsvektoren für alle Ereignisse
k = 1;
for i = 1:size(all_bahn_ids,1)
    for j = 1:all_bahn_ids.("Häufigkeit")(i)-1
        dir(k,:) = all_positions(k+1,:) - all_positions(k,:);
        dist(k) = vecnorm((dir(k,:)));
        ndir(k,:) = dir(k,:)./dist(k);
        quat(k,:) = all_orientations(k+1,:) - all_orientations(k,:);
        k = k + 1;
    end
    k = k + 1;
end

% Merkmals-Array aller Bahnen in den Stützpunkten/Events
all_ident(1:3,:) = ndir';
all_ident(4,:) = dist;
all_ident(5:8,:) = quat';

% Nullvektor zu Beginn anhängen: Alle Bahnen sind durch 0-Vektoren getrennt!
all_ident = [zeros(8,1), all_ident];
zero_columns = find(sum(all_ident) == 0)';

clear dir dist ndir quat i j k clear num_events

%% Finden der ähnlichen Bahnen

% Finden der Indizes im all_ident-Vektor, bei denen eine neue Bahn beginnt. Es werden nur die
% Bahnen gefunden die mindestens genau so viele Bahnsegmente wie die
% gesuchte Bahn enthalten

found_idx = [];
m = size(zero_columns,1);
for i = 1:m
    if i < size(zero_columns,1)
        difference = zero_columns(m) - zero_columns(m-1);
        if difference >= size(bahn_ident,2)+1
            found_idx = [found_idx; zero_columns(m)];
        end
    else 
        difference = zero_columns(m);
        if difference >= size(bahn_ident,2)+1
            found_idx = [found_idx; zero_columns(m)];
        end
    end
    m = m -1;
end

% Index-Vektor mit allen relevanten Startindizes der Bahnen (bezogen auf all_ident)
found_idx = flip(found_idx); 

% Finden des ersten Elements, ab dem die Bahnen geprüft werden!!!
% Überprüfen ob eine Nullspalte vor dem ersten Index liegt
if ~any(all(all_ident(:, 1:found_idx(1)-1) == 0, 1))
    a = 1;
else
    for i = 1:found_idx(1)-1
        % Prüfen, ob Nullspalte vorhanden ist
        if all(all_ident(:, found_idx(1)-i) == 0)
            a = found_idx(1)-i; 
            break;
        end
    end
end

found_idx = [a; found_idx];

clear a difference m i
%% Finden der gleichen Merkmalvektoren in allen Bahnen 

% Vektor der für gleiche Bahnabschnitte kennzeichnet und den Index der jeweiligen Bahn enthält
similarity = [];

threshold_direction = 0.1; % Maximal 10% Abweichung in der Richtung
threshold_length = 0.2; % Maximal 20% Abweichung in der Länge

sz_ident = size(bahn_ident,2);
for i = 1:size(found_idx,1)-1
    % Nullzeile zur Kennzeichnung, dass eine neue Bahn beginnt
    similarity = [similarity; 0 0 0 0]; 
    c = 0;
    last_idx = max(find(all(all_ident(:,found_idx(1) : found_idx(i+1)-1) == 0))); % Anfang des Bahnabschnitts
    comp = all_ident(:,last_idx+1:found_idx(i+1)-1); % Bahnabschnitt der mit dem verglichen wird
    df = 1 + size(comp,2) - sz_ident; % Anzahl der Kombinations-Möglichkeiten zwischen gesuchter und Vergleichsbahn
    
    % Vergleich der Merkmalsvektoren mit allen Möglichkeiten
    for j = 1:df
        for k = 1:sz_ident
            d = c + k;

% Überprüft auf exakte Gleichheit
            % if isequal(bahn_ident(:,k),comp(:,d))
            %     similarity = [similarity; 1 bahn_ids_idx(i)];
            % else
            %     similarity = [similarity;0 bahn_ids_idx(i)];
            % end

% Überprüfung auf Gleichheit entprechend dem Grenzwert
            dir1 = bahn_ident(1:3, k);
            dir2 = comp(1:3, d);
    
            % Prüfe die Winkelabweichung (1 - Cosinus des Winkels zwischen dir1 und dir2)
            direction_deviation = norm(dir1 - dir2);
    
            % Extrahiere die Längen (Zeile 4)
            len1 = bahn_ident(4, k);
            len2 = comp(4, d);
    
            % Berechne die relative Abweichung der Länge
            length_deviation = abs(len1 - len2) / max(len1, len2);
    
            % Überprüfe, ob beide Abweichungen innerhalb der Toleranzgrenzen liegen
            % Schreiben ob Bahnabschnitt gleich ist, Bahn Index und
            % Segmentnummer in den similarity-Vektor
            if direction_deviation <= threshold_direction && length_deviation <= threshold_length
                similarity = [similarity; 1 bahn_ids_idx(i) d k];
            else
                similarity = [similarity; 0 bahn_ids_idx(i) d k];
            end
        end
        c = c +1;        
    end
end

% Finden von Start und Ende der 1er-Sequenzen in der ersten Spalte
differential = diff([0;similarity(:,1);0]);
start_idx = find(differential == 1); 
end_idx = find(differential == -1) - 1; 

% Länge jeder Sequenz berechnen
sequences_all_lengths = end_idx - start_idx + 1;

% Maximale Länge finden
sequences_max_length = max(sequences_all_lengths);
% sequences_max_length = 2;

% Indizes der längsten Sequenzen speichern
sequences_longest = find(sequences_all_lengths == sequences_max_length);

% Ähnliche Bahnen identifizieren
sequences_longest_idx1 =  start_idx(sequences_longest);
sequences_longest_idx2 =  end_idx(sequences_longest);
aehnliche_bahnen_idx = similarity(sequences_longest_idx1,2);
aehnliche_bahnen = table2array(all_bahn_ids(aehnliche_bahnen_idx,1));

% Segmentnummern der ähnlichen Segmente beider Bahnen
seqs = cell(size(aehnliche_bahnen,1),2);
for i = 1:size(aehnliche_bahnen,1)
    seqs{i} = [similarity(sequences_longest_idx1(i),3),similarity(sequences_longest_idx2(i),3)];
    seqs{i,2} = [similarity(sequences_longest_idx1(i),4),similarity(sequences_longest_idx2(i),4)];
end

% Id's und Segmentnummern in einer Tabelle zusammenfassen
aehnliche_seqs = table(aehnliche_bahnen,seqs,'VariableNames', {'bahn_id', 'bahn_abschnitte'});

% Ausgabe der Bahn-Id's und Indizes aus similarity-Vektor
fprintf('Die größte Anzahl gleicher aufeinander folgender Bahnabschnitte beträgt %d.\n', sequences_max_length);
for i = 1:length(sequences_longest)
    fprintf('Bahn_ID %s : Start bei %d, Ende bei %d.\n', ...
        aehnliche_bahnen(i,1), start_idx(sequences_longest(i)), end_idx(sequences_longest(i)));
    
end

% Tabelle zur besseren Übersicht 
similarity_table = table(similarity(:,1),similarity(:,2),similarity(:,3),similarity(:,4),...
    'VariableNames',{'equal?','bahn index','seg_comp','seg_bahn'});

clear c d df sz_ident i j k last_idx start_idx end_idx comp
clear dir1 dir2 len1 len2 direction_deviation length_deviation differential

%% Plotten aller Bahnen 

plot = 1;

if plot

figure;
plot3(events_position(:,1),events_position(:,2),events_position(:,3),'Marker','x','LineWidth',1.5)
xlabel('x');ylabel('y');zlabel('z');
hold on
legend_entries = {"Gesuchte Bahn"};  % Leeres Cell-Array für Legendentexte

for i = 1:1:length(aehnliche_bahnen)
    query = sprintf("SELECT * FROM robotervermessung.bewegungsdaten.bahn_events WHERE bahn_id = '%s'",aehnliche_bahnen(i)); 
    act_bahn = fetch(conn,query);
    act_seqs = act_bahn(aehnliche_seqs.bahn_abschnitte{i}(1):aehnliche_seqs.bahn_abschnitte{i}(end)+1,:);
    pos = table2array(act_seqs(:,5:7));
    plot3(pos(:,1),pos(:,2),pos(:,3),'Marker','o')
    title("Bahn-ID: "+ aehnliche_bahnen(i)+", Bahnabschnitte: "+ string(aehnliche_seqs.bahn_abschnitte{i}(1))+" bis " +string(aehnliche_seqs.bahn_abschnitte{i}(end)) )
    legend_entries{end+1} = aehnliche_bahnen(i)+", Segmente: "+ string(aehnliche_seqs.bahn_abschnitte{i}(1))+" bis " +string(aehnliche_seqs.bahn_abschnitte{i}(end));
    legend(legend_entries);
    % keyboard; % Code anhalten
end

end

toc;

%% Tabelle mit Informationen zu den ähnlichen Bahnen anlegen
tic;

if ~isempty(aehnliche_bahnen)

    % Tabelle initialisieren
    BAHNVERGLEICH  = table('Size',[size(aehnliche_bahnen,1) 12],'VariableTypes',{'string', 'string', 'string', 'string','double','double','double','double','double','double','double','double'}, ...
        'VariableNames',{'bahn_id','robot_model','source_ist','source_soll','num_equal_points','segments','tcp_speed','max_abweichung_soll','avg_abweichung_soll','euclidean','sidtw','equal_orientation'});

    % Alle Soll-Positionsdaten der Hauptbahn auf einmal abrufen
    query = sprintf("SELECT * FROM robotervermessung.bewegungsdaten.bahn_position_soll WHERE bahn_id = '%s'", bahn_id);
    bahn_points = fetch(conn, query);
    % bahn_points = table2array(bahn_points(:,5:7));

    % Eine Abfrage für alle ähnlichen Bahnen (schneller als Schleife für einzelne Bahnen)
    bahn_ids_str = strjoin(string(aehnliche_bahnen), "','");
    query = sprintf("SELECT * FROM robotervermessung.bewegungsdaten.bahn_twist_ist WHERE bahn_id IN ('%s')", bahn_ids_str);
    speed_data = fetch(conn, query);

    query = sprintf("SELECT * FROM robotervermessung.bewegungsdaten.bahn_info WHERE bahn_id IN ('%s')", bahn_ids_str);
    info_data = fetch(conn, query);

    query = sprintf("SELECT * FROM robotervermessung.bewegungsdaten.bahn_position_soll WHERE bahn_id IN ('%s')", bahn_ids_str);
    position_data = fetch(conn, query);
    position_data = sortrows(position_data,"timestamp",'ascend');

    % Werte in Schleife zuweisen
    for i = 1:size(aehnliche_bahnen,1)

        % Filtere die passenden Daten aus den geladenen Tabellen
        act_speed = speed_data(strcmp(speed_data.bahn_id, aehnliche_bahnen(i)), :);
        max_speed = max(act_speed.tcp_speed_ist);

        act_info = info_data(strcmp(info_data.bahn_id, aehnliche_bahnen(i)), :);
        robot_model = act_info.robot_model;
        source_ist = act_info.source_data_ist;
        source_soll = act_info.source_data_soll;

        act_position = position_data(strcmp(position_data.bahn_id, aehnliche_bahnen(i)), :);

% %% Nur die ähnlichen Segmente vergleichen

        b = aehnliche_seqs.bahn_abschnitte{i,2}(1):aehnliche_seqs.bahn_abschnitte{i,2}(end);

        a = aehnliche_seqs.bahn_abschnitte{i,1}(1):aehnliche_seqs.bahn_abschnitte{i,1}(end);
        n = a(end)-a(1)+2;
        ab = strings(n,1);
        ba = strings(n,1);
        for j = 1:n+1
            if j <= n-1
                ab(j) = act_position.bahn_id(1) + '_'+ string(a(j));
                ba(j) = bahn_points.bahn_id(1) + '_' + string(b(j));
            elseif j == n
                ab(j) = act_position.bahn_id(1) + '_'+string(a(j-1)+1);
                ba(j) = bahn_points.bahn_id(1) + '_' + string(b(j-1)+1);
            end
        end
        
        act_position = position_data(ismember(position_data.segment_id, ab), :);
        bahn_points_2 = bahn_points(ismember(bahn_points.segment_id, ba), :);
% %% 

        points = table2array(act_position(:,5:7));
        [~, dists, ~] = distance2curve(points, table2array(bahn_points_2(:,5:7)), 'linear');
        max_dist = max(dists);
        mean_dist = mean(dists);

        % Anzahl gleicher Punkte
        equal_points = sequences_all_lengths(sequences_longest(i)) + 1;

        % Daten in die Tabelle einfügen
        BAHNVERGLEICH{i, :} = [aehnliche_bahnen(i), robot_model, source_ist, source_soll, equal_points, 0, max_speed, max_dist, mean_dist, 0, 0, 0];
        % BAHNVERGLEICH.segments(i) = aehnliche_seqs{i,2};
    end
    % Ähnliche Segmente hinzufügen
    BAHNVERGLEICH.segments = cell(height(BAHNVERGLEICH),1);
    BAHNVERGLEICH.segments = aehnliche_seqs{:,2};
else
    fprintf("Keine ähnlichen Bahnen gefunden! \n")
end

clear source_ist source_soll robot_model max_speed dists max_dist mean_dist act_bahn act_seqs points equal_points
clear i query seqs

clear act_speed act_position act_info info_data  position_data speed_data
 

%% Überprüfen ob die Orientierungen übereinstimmen

if ~isempty(aehnliche_bahnen)

% Tabellenspalte zur Cell-Array umformatieren
BAHNVERGLEICH.equal_orientation = cell(height(BAHNVERGLEICH), 1);

% Orientierungsdaten der ähnlichen Bewegungsbahnen extrahieren
query = sprintf("SELECT * FROM robotervermessung.bewegungsdaten.bahn_events WHERE bahn_id IN ('%s')", bahn_ids_str);
act_bahn_all = fetch(conn,query);
act_bahn_all = sortrows(act_bahn_all,"timestamp");

for i = 1:size(BAHNVERGLEICH,1)

    bahn_quat = bahn_ident(5:8,aehnliche_seqs.bahn_abschnitte{i,2}(1):aehnliche_seqs.bahn_abschnitte{i,2}(end));

    % Differenzen der Orientierungsdaten mit Identifizierungsvektor vergleichen
    act_bahn = act_bahn_all(act_bahn_all.bahn_id == aehnliche_bahnen(i), :);
    act_seqs = act_bahn(aehnliche_seqs.bahn_abschnitte{i,1}(1):aehnliche_seqs.bahn_abschnitte{i,1}(end)+1,:);
    act_quat = table2array(act_seqs(:,8:11));
    act_quat_diff = diff(act_quat)';
    
    quat_compare = [];
    quat_diff = bahn_quat-act_quat_diff;

    % Schreibt Einsen in den Vektor wenn die Orientierung gleich/ähnlich ist
    for j = 1:size(quat_diff,2)
        if norm(quat_diff(:,j)) <= 0.05 % Grenzwert für die Orientierungsabweichung
            quat_compare(end+1) = 1;
        else
            quat_compare(end+1) = 0;
        end
    end

BAHNVERGLEICH.equal_orientation{i} = quat_compare';
end

end 

clear act_bahnn act_bahn_all act_quat act_quat_diff act_seqs i j query
clear bahn_quat quat_compare quat_diff

toc

%% TO-DO

% True oder False ob Geschwindigkeiten übereinstimmen ? Sieht man ja eigentlich

% Positionen manuell eingeben können

% Ähnlichkeit auch für Teilabschnitte anderer Bahnen? z.b. wenn nur 1-2 Bahnabschnitte übereinstimmen ?

% WENN DATEN MIT GLEICHEN TIMESTAMP DOPPELT VORLIEGEN MÜSSEN DIE
% DOPPPELTEN DATEN GELÖSCHT WERDEN !
% ---> GROßES PROBLEM DA NICHT ALLE SEQUENZEN ÜBERLAPPEN 
 


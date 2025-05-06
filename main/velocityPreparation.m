function velocityPreparation(data_soll, data_ist)

v_ist = double(data_ist.tcp_speed_ist);
t_ist = double(string(data_ist.timestamp));

% Berechnung der Abstände zwischen den Punkten der Sollbahn
p_soll = table2array(data_soll(:,5:7));
d_soll = diff(p_soll);
d_soll = sqrt(d_soll(:,1).^2 + d_soll(:,2).^2 + d_soll(:,3).^2);
t_soll = double(string(data_soll.timestamp));

% % Timestamps in Sekunden
t_soll = (t_soll - t_ist(1))/1e9;
t_ist = (t_ist - t_ist(1))/1e9;

% Berechnung der Frequenz der Solldaten Fensterlänge für Mittelwert-Filter
freq_soll = 1/((t_soll(end)-t_soll(1))/length(t_soll));
window = round(freq_soll/100);

% Indizes der Abstänbde > 0 
idx = find(d_soll ~= 0);

data_soll_new = data_soll;
data_soll_new = data_soll_new(idx,:);

% Erneute Berechnung der Abstände und Zeitstempel
p_soll = p_soll(idx,:);
d_soll = abs(diff(p_soll));
d_soll = sqrt(d_soll(:,1).^2 + d_soll(:,2).^2 + d_soll(:,3).^2);
t_soll = t_soll(idx);
t_soll = linspace(t_soll(1),t_soll(end),length(d_soll)+1)';

% Gleitenden Median berechnen (Fenstergröße: 5)
d_soll = movmedian(d_soll, window);
v_soll = d_soll./diff(t_soll);
t_soll = linspace(t_soll(1),t_soll(end),length(v_soll))';

% Soll und Ist-Daten un neue Tabellen schreiben und in Workspace laden
data_ist_new = [data_ist(:,{'bahn_id','segment_id'}) table(t_ist,v_ist,'VariableNames',{'time_ist','tcp_speed_ist'})];
data_soll_new = [data_soll_new(2:end,{'bahn_id','segment_id'}) table(t_soll,v_soll,'VariableNames',{'time_soll','tcp_speed_soll'})];

assignin("caller","data_ist",data_ist_new)
assignin("caller","data_soll",data_soll_new)

end
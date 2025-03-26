function connection = connectingToPostgres

% Verbinden mit Datenbank
datasource = "RobotervermessungMATLAB";
username = "felixthomas";
password = "manager";
connection = postgresql(datasource,username,password);

% Überprüfe Verbindung
if isopen(connection)
    disp('Verbindung erfolgreich hergestellt');
else
    disp('Verbindung fehlgeschlagen');
end

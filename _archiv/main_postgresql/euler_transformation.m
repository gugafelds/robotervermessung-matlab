function euler_transformation(varargin)

data_ist = varargin{1};
data_soll = varargin{2};
    
% Wenn Kalibrierung erfolgt!
if nargin == 2

    % Transformation der Quarternionen zu Euler-Winkeln
    q_soll = table2array(data_soll(:,5:8));
    q_soll = [q_soll(:,4), q_soll(:,3), q_soll(:,2), q_soll(:,1)];
    euler_soll = quat2eul(q_soll,"ZYX");
    euler_soll = rad2deg(euler_soll);

    
    q_ist = table2array(data_ist(:,8:11));
    q_ist = [q_ist(:,4), q_ist(:,3), q_ist(:,2), q_ist(:,1)]; % original
    euler_ist = quat2eul(q_ist,"XYZ");
    euler_ist = rad2deg(euler_ist);
    
    % Berechnung der Rotationsmatrizen aller Orientierungen auf Soll & Istbahn
    R_ist = eul2rotm(deg2rad(euler_ist),"ZYX");
    R_ist_mean = [mean(R_ist(1,1,:)),mean(R_ist(1,2,:)),mean(R_ist(1,3,:)); ...
                   mean(R_ist(2,1,:)),mean(R_ist(2,2,:)),mean(R_ist(2,3,:)); ...
                   mean(R_ist(3,1,:)),mean(R_ist(3,2,:)),mean(R_ist(3,3,:));];
    R_soll = eul2rotm(deg2rad(euler_soll),"ZYX");
    R_soll_mean = [mean(R_soll(1,1,:)),mean(R_soll(1,2,:)),mean(R_soll(1,3,:)); ...
                   mean(R_soll(2,1,:)),mean(R_soll(2,2,:)),mean(R_soll(2,3,:)); ...
                   mean(R_soll(3,1,:)),mean(R_soll(3,2,:)),mean(R_soll(3,3,:));];

    % Gemittelte relative Rotationsmatrix
    R_rel = R_soll_mean * R_ist_mean;

    
    % Laden in Workspace
    assignin("base","trafo_euler", R_rel)

%  Transformation für die Ist-Daten 
else

    % Unterscheide ob ganze Datei oder segmentweise ausgewertet wird
    if istable(data_ist)
        % Konvertierung von Table mit Cell-Daten zu n x 3 double-Array
        euler_ist = cell2mat(data_ist{:,2:4}(:,1:3));
    else
        euler_ist = data_ist;
    end

    % Rotationsmatrix der Koordinatentransformation
    rot_matrix = varargin{4};
    euler_ist = euler_ist * rot_matrix;

    % Berechnung der Rotationsmatrizen aller Orientierungs-Koordinaten
    R_ist = eul2rotm(deg2rad(euler_ist),"ZYX");
%%%%%%%%%%%
    % Transponieren der Ist-Matrizen; 
    % R_ist =permute(R_ist, [2, 1, 3]);
%%%%%%%%%%%
    R_rel = varargin{3};
    
    % Erstellen eines mehrdimensionalen 3x3-Arrays
    R_rel = repmat(R_rel,1,1,length(R_ist));
    
    % Mehrdimensionale Matrixmulitiplikation
    R_trans = pagemtimes(R_rel,R_ist);
    % R_trans = pagemtimes(R_ist,R_rel);
    
    % Transformierte Eulerwinkel der Istbahn
    euler_trans = rad2deg(rotm2eul(R_trans,"ZYX"));
    
    if istable(data_ist)
        % Als Cell-Struktur konvertieren und Istdaten überschreiben
        data_ist(:,2:4) = array2table([{euler_trans(:,1)} {euler_trans(:,2)} {euler_trans(:,3)}]);
        assignin("base","seg_trafo",data_ist)
    else 
        assignin("base","euler_trans",euler_trans)
    end


end
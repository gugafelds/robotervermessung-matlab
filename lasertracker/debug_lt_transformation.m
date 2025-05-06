% LASERTRACKER ZU ROBOTERKOORDINATEN TRANSFORMATION
% Dieses Skript konvertiert Laser Tracker Messungen in Roboterkoordinaten.
% Es verwendet einen direkten Ansatz über Euler-Winkel, um hochpräzise Umwandlungen zu ermöglichen.
%
% Autor: 
% Datum: März 2025

%% Konfiguration - Anpassen für Ihre spezifische Messaufgabe

% Transformationsmatrizen
% Offset vom gemessenen Punkt zum TCP (Tool Center Point)
offset_matrix = [
    -0.49646794506312 0.86797124803586 -0.01206201093802 20.15923993520914;
    -0.01103757423601 0.00758219542268 0.99991033711406 -176.89448180892975;
    0.86798487975297 0.49655656565556 0.00581598010435 -0.20846192267220;
    0.00000000000000 0.00000000000000 0.00000000000000 1.00000000000000
];

% Transformation vom Instrumenten- zum Roboterkoordinatensystem
robot_to_instrument = [
    0.08464662262294 -0.99622482623969 0.01926252481472 3633.87661777230051;
    0.99624481798547 0.08496977929709 0.01662525919475 699.29196877805589;
    -0.01819922843470 0.01778291848691 0.99967622553228 1277.96127304069205;
    0.00000000000000 0.00000000000000 0.00000000000000 1.00000000000000
];

%% Eingabedaten - Hier Ihre Messdaten eintragen

% Format: [x, y, z] Koordinaten in mm
points_lt = [
    % Ersetzen Sie dies mit Ihren gemessenen Laser Tracker Punkten
    % Punkt 1
    -1933.39538, 1935.3579, -414.162973;
    % Punkt 2
    -1343.44, 2632.65, 258.511;
    % Punkt 3
    -1953.309, 2594.8889, 293.27911;
    % Punkt 4
    -1589.25623364662,	2500.76685219856,	-236.520593219670;
    % Punkt 5
    -910.9465, 2174.84966, 378.690; % Null-Punkt Roboter
];

% Format: [qx, qy, qz, qw] Quaternionen für die ersten 3 Punkte
quaternions_lt = [
    % Punkt 1
    0.43971865153854, -0.664874486437972, 0.591112540721293, 0.123188428738208;
    % Punkt 2
    0.495925203026302, -0.506824145345618, 0.687230701632917, 0.15780190566382;
    % Punkt 3
    -0.323481141718943, 0.760407624606544, 0.228725190620096, 0.514545686101458;
    % Fügen Sie Dummy-Werte für die restlichen Punkte ein - diese werden aus Euler-Winkeln berechnet
    0.401228858948662,	-0.505398836986511,	0.757433680054831,	0.0994064315727732;
    0, 0, 0, 1;
];

% Format: [roll, pitch, yaw] Winkel in Grad
orientations_lt = [
    % Ersetzen Sie dies mit Ihren gemessenen Laser Tracker Orientierungen
    % Punkt 1
    -111.779, -43.12744, -143.0161;
    % Punkt 2
    -90.6072, -57.3074457, -148.0362;
    % Punkt 3
    -44.5318, -68.59512, 177.72158;
    % Punkt 4
    %110.2679, -56.0794, 18.32455;
    % Punkt 5
    87.9872, -30.6799, 6.322099; 
];

%% Verarbeitung - Transformation der Messdaten

% Initialisierung der Ergebnisarrays
num_points = size(points_lt, 1);
points_robot = zeros(num_points, 3);
orientations_robot = zeros(num_points, 3);
quaternions_robot = zeros(num_points, 4);

% Verarbeitung jedes Punktes
for i = 1:num_points
    % 1. Konvertiere LT Orientierung zu SA-kompatiblem Quaternion
    quaternion_sa = lt2sa_quaternion_direct(quaternions_lt(i,:));
    
    % 2. Anwenden des Offsets auf den Punkt und die Orientierung
    [point_offset, ~, quaternion_offset] = applyOffset(points_lt(i,:), [], quaternion_sa, offset_matrix);
    
    % 3. Transformation in Roboterkoordinaten
    points_robot(i,:) = transformToRobotCoords(point_offset, robot_to_instrument);
    quaternions_robot(i,:) = transformQuaternionToRobotCoords(quaternion_offset, robot_to_instrument);
    
    % 4. Berechnung der Euler-Winkel für die Robotersteuerung (XYZ Fixed Angles)
    orientations_robot(i,:) = quaternionToEulerXYZFixed(quaternions_robot(i,:));
end

%% Ausgabe der Ergebnisse

disp('==============================================');
disp('TRANSFORMATION VON LASER TRACKER ZU ROBOTERKOORDINATEN');
disp('==============================================');
disp(' ');

for i = 1:num_points
    disp(['Punkt ', num2str(i), ':']);
    disp(['   Position (Roboter) = [', num2str(points_robot(i,:), '%.3f, '), '] mm']);
    disp(['   Orientierung (Roboter XYZ) = [', num2str(orientations_robot(i,:), '%.3f, '), '] Grad']);
    disp(['   Quaternion (Roboter) = [', num2str(quaternions_robot(i,:), '%.6f, '), ']']);
    disp(' ');
end

% Speichern der Ergebnisse in einer CSV-Datei
results_table = table((1:num_points)', points_robot, orientations_robot, quaternions_robot, ...
    'VariableNames', {'PointNumber', 'Position', 'Orientation', 'Quaternion'});
writetable(results_table, 'robot_coordinates.csv');

disp('Die Ergebnisse wurden in "robot_coordinates.csv" gespeichert.');
disp('==============================================');

%% Hilfsfunktionen

function q_sa = lt2sa_quaternion_direct(q_lt_or_euler_lt)
% LT2SA_QUATERNION_DIRECT konvertiert Laser Tracker Quaternionen oder Euler-Winkel
% in Spatial Analyzer Quaternionen
%
% Eingabe:
%   q_lt_or_euler_lt - Entweder ein 1x4 Array mit Laser Tracker Quaternion [x, y, z, w]
%                    ODER ein 1x3 Array mit Euler-Winkeln [roll, pitch, yaw] in Grad
%
% Ausgabe:
%   q_sa - Spatial Analyzer Quaternion im Format [x, y, z, w]

    % Überprüfen, ob ein Quaternion oder Euler-Winkel übergeben wurden
    if length(q_lt_or_euler_lt) == 4
        % Es wurde ein Quaternion übergeben
        q_lt = q_lt_or_euler_lt;
        
        % Umwandlung von LT zu SA Quaternion mit der MATLAB Quaternion-Methode
        if exist('quaternion', 'file') == 2  % Überprüfen, ob MATLAB Aerospace Toolbox verfügbar ist
            % MATLAB quaternion-Klasse verwenden
            q_lt_matlab = quaternion(q_lt(4), q_lt(1), q_lt(2), q_lt(3));
            
            % Rotation um -90 Grad um die X-Achse
            angle = -pi/2;
            rotation_axis = [1, 0, 0];
            q_rotation = quaternion(cos(angle/2), sin(angle/2)*rotation_axis(1), ...
                                    sin(angle/2)*rotation_axis(2), sin(angle/2)*rotation_axis(3));
            
            % Anwenden der Rotation
            q_rotated = q_rotation * q_lt_matlab;
            
            % Extrahieren der Komponenten
            [w, x, y, z] = parts(q_rotated);
            
            % Permutation und Vorzeichen anpassen
            q_sa = [-w, -y, -z, x];
        else
            % Manuelle Quaternion-Multiplikation ohne MATLAB Toolbox
            % Rotation um -90 Grad um die X-Achse entspricht dem Quaternion [sin(-pi/4), 0, 0, cos(-pi/4)]
            cos_half_angle = 0.7071067811865475;  % cos(-pi/4)
            sin_half_angle = -0.7071067811865475; % sin(-pi/4)
            
            % Quaternion-Komponenten extrahieren
            qx = q_lt(1);
            qy = q_lt(2);
            qz = q_lt(3);
            qw = q_lt(4);
            
            % Quaternion-Multiplikation q_rotation * q_lt manuell berechnen
            rotw = cos_half_angle * qw - sin_half_angle * qx;
            rotx = cos_half_angle * qx + sin_half_angle * qw;
            roty = cos_half_angle * qy + sin_half_angle * qz;
            rotz = cos_half_angle * qz - sin_half_angle * qy;
            
            % Permutation und Vorzeichen anpassen
            q_sa = [-rotw, -roty, -rotz, rotx];
        end
    else
        % Es wurden Euler-Winkel übergeben
        euler_lt = q_lt_or_euler_lt;
        
        % Anpassung des Roll-Winkels um +90 Grad
        euler_sa = euler_lt;
        euler_sa(1) = euler_lt(1) + 90; % Roll-Korrektur

        % Konvertiere zu Quaternion mit der Hilfsfunktion
        q_sa = euler2quaternion(euler_sa);
    end

    % Normalisiere zur Sicherheit
    q_sa = q_sa / norm(q_sa);
end

% Hilfsfunktion für die Umwandlung von Euler-Winkeln zu Quaternionen
% Diese Funktion wird nur verwendet, wenn Euler-Winkel übergeben werden
function q = euler2quaternion(euler)
    % Umrechnung in Radian
    roll = euler(1) * pi / 180;
    pitch = euler(2) * pi / 180;
    yaw = euler(3) * pi / 180;

    % Berechne Sinus und Kosinus der halben Winkel
    cr = cos(roll * 0.5);
    sr = sin(roll * 0.5);
    cp = cos(pitch * 0.5);
    sp = sin(pitch * 0.5);
    cy = cos(yaw * 0.5);
    sy = sin(yaw * 0.5);

    % Quaternion berechnen
    q = [
        sr * cp * cy - cr * sp * sy,
        cr * sp * cy + sr * cp * sy,
        cr * cp * sy - sr * sp * cy,
        cr * cp * cy + sr * sp * sy
    ];

    % Normalisieren
    q = q / norm(q);
end

function [new_position, new_orientation, new_quaternion] = applyOffset(position, orientation, quaternion, offset_matrix)
% APPLYOFFSET Wendet eine Offset-Transformationsmatrix auf Position, Orientierung und Quaternion an
%
% Eingabe:
%   position    - 1x3 Array mit [x, y, z] Koordinaten
%   orientation - 1x3 Array mit [roll, pitch, yaw] Winkeln in Grad (optional, kann [] sein)
%   quaternion  - 1x4 Array mit [qx, qy, qz, qw] Quaternion (optional, kann [] sein)
%   offset_matrix - 4x4 homogene Transformationsmatrix
%
% Ausgabe:
%   new_position    - 1x3 Array mit transformierten [x, y, z] Koordinaten
%   new_orientation - 1x3 Array mit transformierten [roll, pitch, yaw] Winkeln in Grad
%   new_quaternion  - 1x4 Array mit transformiertem [qx, qy, qz, qw] Quaternion

    % Überprüfen, welche Orientierungsangabe verwendet werden soll
    if ~isempty(quaternion)
        % Quaternion in Rotationsmatrix umwandeln
        qx = quaternion(1);
        qy = quaternion(2);
        qz = quaternion(3);
        qw = quaternion(4);
        
        R_original = [
            1-2*(qy^2+qz^2), 2*(qx*qy-qz*qw), 2*(qx*qz+qy*qw);
            2*(qx*qy+qz*qw), 1-2*(qx^2+qz^2), 2*(qy*qz-qx*qw);
            2*(qx*qz-qy*qw), 2*(qy*qz+qx*qw), 1-2*(qx^2+qy^2)
        ];
    elseif ~isempty(orientation)
        % Umwandlung der Orientierung von Grad in Radiant
        o_rad = orientation * pi/180;
        
        % Erstellen der Rotationsmatrizen
        Rx = [1 0 0; 0 cos(o_rad(1)) -sin(o_rad(1)); 0 sin(o_rad(1)) cos(o_rad(1))];
        Ry = [cos(o_rad(2)) 0 sin(o_rad(2)); 0 1 0; -sin(o_rad(2)) 0 cos(o_rad(2))];
        Rz = [cos(o_rad(3)) -sin(o_rad(3)) 0; sin(o_rad(3)) cos(o_rad(3)) 0; 0 0 1];
        
        % Gesamtrotation (Reihenfolge: Z-Y-X oder Roll-Pitch-Yaw)
        R_original = Rz * Ry * Rx;
    else
        error('Entweder orientation oder quaternion muss angegeben werden');
    end
    
    % Homogene Transformationsmatrix erstellen
    measured_matrix = [
        R_original, position';
        0 0 0 1
    ];
    
    % Anwenden der Transformation
    result_matrix = measured_matrix * offset_matrix;
    
    % Extrahieren der neuen Position
    new_position = result_matrix(1:3, 4)';
    
    % Extrahieren der neuen Rotation
    new_rotation = result_matrix(1:3, 1:3);
    
    % Umwandlung der Rotationsmatrix zurück in Roll-Pitch-Yaw (in Grad)
    beta = atan2(-new_rotation(3,1), sqrt(new_rotation(1,1)^2 + new_rotation(2,1)^2));
    alpha = atan2(new_rotation(2,1)/cos(beta), new_rotation(1,1)/cos(beta));
    gamma = atan2(new_rotation(3,2)/cos(beta), new_rotation(3,3)/cos(beta));
    
    new_orientation = [alpha, beta, gamma] * 180/pi;
    
    % Umwandlung der Rotationsmatrix in Quaternion
    trace = new_rotation(1,1) + new_rotation(2,2) + new_rotation(3,3);
    
    if trace > 0
        s = 0.5 / sqrt(trace + 1.0);
        qw = 0.25 / s;
        qx = (new_rotation(3,2) - new_rotation(2,3)) * s;
        qy = (new_rotation(1,3) - new_rotation(3,1)) * s;
        qz = (new_rotation(2,1) - new_rotation(1,2)) * s;
    else
        if new_rotation(1,1) > new_rotation(2,2) && new_rotation(1,1) > new_rotation(3,3)
            s = 2.0 * sqrt(1.0 + new_rotation(1,1) - new_rotation(2,2) - new_rotation(3,3));
            qw = (new_rotation(3,2) - new_rotation(2,3)) / s;
            qx = 0.25 * s;
            qy = (new_rotation(1,2) + new_rotation(2,1)) / s;
            qz = (new_rotation(1,3) + new_rotation(3,1)) / s;
        elseif new_rotation(2,2) > new_rotation(3,3)
            s = 2.0 * sqrt(1.0 + new_rotation(2,2) - new_rotation(1,1) - new_rotation(3,3));
            qw = (new_rotation(1,3) - new_rotation(3,1)) / s;
            qx = (new_rotation(1,2) + new_rotation(2,1)) / s;
            qy = 0.25 * s;
            qz = (new_rotation(2,3) + new_rotation(3,2)) / s;
        else
            s = 2.0 * sqrt(1.0 + new_rotation(3,3) - new_rotation(1,1) - new_rotation(2,2));
            qw = (new_rotation(2,1) - new_rotation(1,2)) / s;
            qx = (new_rotation(1,3) + new_rotation(3,1)) / s;
            qy = (new_rotation(2,3) + new_rotation(3,2)) / s;
            qz = 0.25 * s;
        end
    end
    
    new_quaternion = [qx, qy, qz, qw];
end

function p_robot = transformToRobotCoords(p_instr, instr_to_robot_matrix)
% TRANSFORMTOROBOTCOORDS Transformiert einen Punkt vom Instrument- ins Roboterkoordinatensystem
%
% Eingabe:
%   p_instr             - 1x3 Array mit [x, y, z] Koordinaten im Instrumentensystem
%   instr_to_robot_matrix - 4x4 homogene Transformationsmatrix
%
% Ausgabe:
%   p_robot             - 1x3 Array mit [x, y, z] Koordinaten im Robotersystem

    % Punkt in homogene Koordinaten umwandeln
    p_homogeneous = [p_instr, 1]';
    
    % Transformation durchführen
    p_robot_homogeneous = instr_to_robot_matrix * p_homogeneous;
    
    % Zurück zu kartesischen Koordinaten
    p_robot = p_robot_homogeneous(1:3)';
end

function q_robot = transformQuaternionToRobotCoords(quaternion, instr_to_robot_matrix)
% TRANSFORMQUATERNIONTOROBOTCOORDS Transformiert ein Quaternion ins Roboterkoordinatensystem
%
% Eingabe:
%   quaternion           - 1x4 Array mit [qx, qy, qz, qw] Quaternion
%   instr_to_robot_matrix - 4x4 homogene Transformationsmatrix
%
% Ausgabe:
%   q_robot              - 1x4 Array mit transformiertem [qx, qy, qz, qw] Quaternion

    % Extrahieren der Rotationsmatrix aus der homogenen Transformationsmatrix
    R_instrument_to_robot = instr_to_robot_matrix(1:3, 1:3);
    
    % Quaternion in Rotationsmatrix umwandeln
    qx = quaternion(1);
    qy = quaternion(2);
    qz = quaternion(3);
    qw = quaternion(4);
    
    R_original = [
        1-2*(qy^2+qz^2), 2*(qx*qy-qz*qw), 2*(qx*qz+qy*qw);
        2*(qx*qy+qz*qw), 1-2*(qx^2+qz^2), 2*(qy*qz-qx*qw);
        2*(qx*qz-qy*qw), 2*(qy*qz+qx*qw), 1-2*(qx^2+qy^2)
    ];
    
    % Rotationsmatrizen kombinieren
    % Die neue Rotationsmatrix ist die Kombination der ursprünglichen Rotation 
    % und der Rotation des Koordinatensystems
    R_new = R_instrument_to_robot * R_original;
    
    % Umwandlung der kombinierten Rotationsmatrix in Quaternion
    trace = R_new(1,1) + R_new(2,2) + R_new(3,3);
    
    if trace > 0
        s = 0.5 / sqrt(trace + 1.0);
        qw_new = 0.25 / s;
        qx_new = (R_new(3,2) - R_new(2,3)) * s;
        qy_new = (R_new(1,3) - R_new(3,1)) * s;
        qz_new = (R_new(2,1) - R_new(1,2)) * s;
    else
        if R_new(1,1) > R_new(2,2) && R_new(1,1) > R_new(3,3)
            s = 2.0 * sqrt(1.0 + R_new(1,1) - R_new(2,2) - R_new(3,3));
            qw_new = (R_new(3,2) - R_new(2,3)) / s;
            qx_new = 0.25 * s;
            qy_new = (R_new(1,2) + R_new(2,1)) / s;
            qz_new = (R_new(1,3) + R_new(3,1)) / s;
        elseif R_new(2,2) > R_new(3,3)
            s = 2.0 * sqrt(1.0 + R_new(2,2) - R_new(1,1) - R_new(3,3));
            qw_new = (R_new(1,3) - R_new(3,1)) / s;
            qx_new = (R_new(1,2) + R_new(2,1)) / s;
            qy_new = 0.25 * s;
            qz_new = (R_new(2,3) + R_new(3,2)) / s;
        else
            s = 2.0 * sqrt(1.0 + R_new(3,3) - R_new(1,1) - R_new(2,2));
            qw_new = (R_new(2,1) - R_new(1,2)) / s;
            qx_new = (R_new(1,3) + R_new(3,1)) / s;
            qy_new = (R_new(2,3) + R_new(3,2)) / s;
            qz_new = 0.25 * s;
        end
    end
    
    q_robot = [qx_new, qy_new, qz_new, qw_new];
end

function euler_angles = quaternionToEulerXYZFixed(quaternion)
% QUATERNIONTOEULERRXYZFIXED Wandelt ein Quaternion in XYZ Fixed Angles (extrinsisch) um
%
% Eingabe:
%   quaternion  - 1x4 Array mit [qx, qy, qz, qw] Quaternion
%
% Ausgabe:
%   euler_angles - 1x3 Array mit [roll, pitch, yaw] Winkeln in Grad (XYZ Fixed Angles)

    % Quaternion-Komponenten extrahieren
    qx = quaternion(1);
    qy = quaternion(2);
    qz = quaternion(3);
    qw = quaternion(4);
    
    % Quaternion in Rotationsmatrix umwandeln
    R = [
        1-2*(qy^2+qz^2), 2*(qx*qy-qz*qw), 2*(qx*qz+qy*qw);
        2*(qx*qy+qz*qw), 1-2*(qx^2+qz^2), 2*(qy*qz-qx*qw);
        2*(qx*qz-qy*qw), 2*(qy*qz+qx*qw), 1-2*(qx^2+qy^2)
    ];
    
    % Berechnung der XYZ Fixed Angles (extrinsisch)
    % Dies entspricht der Anwendung von Rotationen um die X-, Y- und Z-Achse in dieser Reihenfolge
    % im festen Koordinatensystem (nicht im mitrotierenden Koordinatensystem)
    
    % Bei XYZ Fixed Angles erfolgt die Extraktion aus der Rotationsmatrix wie folgt:
    % R = Rz(yaw) * Ry(pitch) * Rx(roll)
    
    % Pitch (Y-Achse)
    sin_pitch = -R(3,1);
    if abs(sin_pitch) >= 1
        pitch = sign(sin_pitch) * pi/2; % Gimbal Lock
    else
        pitch = asin(sin_pitch);
    end
    
    % Überprüfung auf Gimbal Lock
    if abs(sin_pitch) > 0.99999
        % Bei Gimbal Lock gibt es unendlich viele Lösungen für Roll und Yaw, 
        % wir setzen Roll auf 0 und berechnen Yaw
        roll = 0;
        yaw = atan2(-R(1,2), R(1,1));
    else
        % Roll (X-Achse)
        roll = atan2(R(3,2), R(3,3));
        
        % Yaw (Z-Achse)
        yaw = atan2(R(2,1), R(1,1));
    end
    
    % Umwandlung von Radiant in Grad
    euler_angles = [roll, pitch, yaw] * 180/pi;
end
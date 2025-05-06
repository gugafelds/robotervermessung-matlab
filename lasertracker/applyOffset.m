function [new_position, new_orientation, new_quaternion] = applyOffset(position, orientation, quaternion, offset_matrix)
% APPLYOFFSET Wendet eine Offset-Transformationsmatrix auf Position, Orientierung und Quaternion an
%
% Eingabe:
%   position    - 1x3 Array mit [x, y, z] Koordinaten
%   orientation - 1x3 Array mit [roll, pitch, yaw] Winkeln in Grad (optional, kann [] sein)
%   quaternion  - 1x4 Array mit [qw, qx, qy, qz] Quaternion (optional, kann [] sein)
%   offset_matrix - 4x4 homogene Transformationsmatrix
%
% Ausgabe:
%   new_position    - 1x3 Array mit transformierten [x, y, z] Koordinaten
%   new_orientation - 1x3 Array mit transformierten [roll, pitch, yaw] Winkeln in Grad
%   new_quaternion  - 1x4 Array mit transformiertem [qw, qx, qy, qz] Quaternion

    % Überprüfen, welche Orientierungsangabe verwendet werden soll
    if ~isempty(quaternion)
        % Quaternion in Rotationsmatrix umwandeln
        qw = quaternion(1);
        qx = quaternion(2);
        qy = quaternion(3);
        qz = quaternion(4);
        
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
    
    new_quaternion = [qw, qx, qy, qz];
end
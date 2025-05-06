function q_robot = transformQuaternionToRobotCoords(quaternion, instr_to_robot_matrix)
% TRANSFORMQUATERNIONTOROBOTCOORDS Transformiert ein Quaternion ins Roboterkoordinatensystem
%
% Eingabe:
%   quaternion           - 1x4 Array mit [qw, qx, qy, qz] Quaternion
%   instr_to_robot_matrix - 4x4 homogene Transformationsmatrix
%
% Ausgabe:
%   q_robot              - 1x4 Array mit transformiertem [qw, qx, qy, qz] Quaternion

    % Extrahieren der Rotationsmatrix aus der homogenen Transformationsmatrix
    R_instrument_to_robot = instr_to_robot_matrix(1:3, 1:3);
    
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
    
    q_robot = [qw_new, qx_new, qy_new, qz_new];
end
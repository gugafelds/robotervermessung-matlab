function euler_angles = quaternionToEulerXYZFixed(quaternion)
% QUATERNIONTOEULERRXYZFIXED Wandelt ein Quaternion in XYZ Fixed Angles (extrinsisch) um
%
% Eingabe:
%   quaternion  - 1x4 Array mit [qw, qx, qy, qz] Quaternion
%
% Ausgabe:
%   euler_angles - 1x3 Array mit [roll, pitch, yaw] Winkeln in Grad (XYZ Fixed Angles)

    % Quaternion-Komponenten extrahieren
    qw = quaternion(1);
    qx = quaternion(2);
    qy = quaternion(3);
    qz = quaternion(4);
    
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
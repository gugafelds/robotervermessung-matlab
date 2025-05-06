function q_sa = LT2SAQuaternionDirect(q_lt_or_euler_lt)
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
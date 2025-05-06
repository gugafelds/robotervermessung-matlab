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
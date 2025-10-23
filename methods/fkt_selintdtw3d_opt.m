function [distances, maxDistance, averageDistance, AccumulatedDistance, dtw_X, dtw_Y, MappingIndexes, ix, iy] = fkt_selintdtw3d_opt(X,Y,pflag)
%FKT_SELINTDTW3D_OPT Optimierte Version mit Vectorisierung
%   Selective Interpolation DTW für 3D Bahnen mit Performance-Optimierungen
%   - Vectorisierte Distanzberechnung
%   - Reduzierte Schleifeniterationen
%   - Effizientere Matrix-Operationen
%
%   OPTIMIERUNGEN:
%   1. Vectorisierte Berechnung aller X-Kanten Distanzen auf einmal
%   2. Vectorisierte Berechnung aller Y-Kanten Distanzen auf einmal
%   3. Vorberechnung statt wiederholte Berechnungen in Schleifen
%   4. Effizientere Speicherverwaltung
%
%   Erwarteter Speedup: 5-10x gegenüber Original

if nargin < 3 
    pflag = false; 
end

% Initialisierung der Variablen
M = length(X);
N = length(Y);
AccumulatedDistance = zeros(M,N);
AccumulatedDistanceX = zeros(M,N);
AccumulatedDistanceY = zeros(M,N);
ParamX = zeros(M,N);
ParamY = zeros(M,N);
minParams = cell(M,N);
minParamsX = cell(M,N);
minParamsY = cell(M,N);

fprintf('Initialisiere Matrizen: %d x %d\n', M, N);

% Initialisierung der Anfangsbedingungen
AccumulatedDistance(1,1) = fkt_euclDist(1,1,X,Y);
AccumulatedDistanceX(1,1) = Inf;
AccumulatedDistanceY(1,1) = Inf;
ParamX(1,1) = NaN;
ParamY(1,1) = NaN;
minParams{1,1} = [NaN, NaN];
minParamsX{1,1} = [NaN, NaN];
minParamsY{1,1} = [NaN, NaN];

%% OPTIMIERUNG 1: Vectorisierte Startwerte X
fprintf('Berechne Startwerte X (vectorisiert)...\n');
if M > 1
    % Vectorisierte Berechnung für alle X-Startwerte auf einmal
    [mindists_X, params_X] = fkt_minDistParam_vectorized_edge(X(1:M-1,:), X(2:M,:), repmat(Y(1,:), M-1, 1));
    
    for i = 2:M
        AccumulatedDistanceX(i,1) = AccumulatedDistanceX(i-1,1) + mindists_X(i-1);
        AccumulatedDistanceY(i,1) = Inf;
        ParamX(i,1) = params_X(i-1);
        ParamY(i,1) = NaN;
        minParams{i,1} = [0,1];
        minParamsX{i,1} = [0,1];
        minParamsY{i,1} = [NaN,NaN];
        AccumulatedDistance(i,1) = AccumulatedDistance(i-1,1) + fkt_euclDist(i,1,X,Y);
        if fkt_isInterpolation(ParamX(i,1)) && AccumulatedDistanceX(i,1) < AccumulatedDistance(i-1,1)
            AccumulatedDistance(i,1) = AccumulatedDistanceX(i,1) + fkt_euclDist(i,1,X,Y);
            minParams{i,1} = [ParamX(i,1),1];
        end
    end
end

%% OPTIMIERUNG 2: Vectorisierte Startwerte Y
fprintf('Berechne Startwerte Y (vectorisiert)...\n');
if N > 1
    % Vectorisierte Berechnung für alle Y-Startwerte auf einmal
    [mindists_Y, params_Y] = fkt_minDistParam_vectorized_edge(Y(1:N-1,:), Y(2:N,:), repmat(X(1,:), N-1, 1));
    
    for j = 2:N
        AccumulatedDistanceX(1,j) = Inf;
        AccumulatedDistanceY(1,j) = AccumulatedDistance(1,j-1) + mindists_Y(j-1);
        ParamX(1,j) = NaN;
        ParamY(1,j) = params_Y(j-1);
        minParams{1,j} = [1,0];
        minParamsX{1,j} = [NaN,NaN];
        minParamsY{1,j} = [1,0];
        AccumulatedDistance(1,j) = AccumulatedDistance(1,j-1) + fkt_euclDist(1,j,X,Y);
        if fkt_isInterpolation(ParamY(1,j)) && AccumulatedDistanceY(1,j) < AccumulatedDistance(1,j-1)
            AccumulatedDistance(1,j) = AccumulatedDistanceY(1,j) + fkt_euclDist(1,j,X,Y);
            minParams{1,j} = [1,ParamY(1,j)];
        end
    end
end

%% OPTIMIERUNG 3: Vectorisierte Hauptberechnung
% Vorberechnung aller X-Kanten Distanzen und Parameter
fprintf('Berechne alle X-Kanten Distanzen (vectorisiert)...\n');
[all_mindist_X, all_param_X] = compute_all_edge_distances(X, Y, 'X');

fprintf('Berechne alle Y-Kanten Distanzen (vectorisiert)...\n');
[all_mindist_Y, all_param_Y] = compute_all_edge_distances(Y, X, 'Y');

% Hauptschleife - kann nicht vollständig vectorisiert werden wegen Abhängigkeiten
fprintf('Erstellung der akkumulierten Kostenmatrix (%d x %d)...\n', M, N);
tic;
for i = 2:M
    if mod(i, 50) == 0
        fprintf('  Fortschritt: %d/%d (%.1f%%)\n', i, M, 100*i/M);
    end
    
    for j = 2:N
        % Nutze vorberechnete Werte
        mindist_X = all_mindist_X(i-1, j);
        param_X = all_param_X(i-1, j);
        AccumulatedDistanceX(i,j) = mindist_X; 
        ParamX(i,j) = param_X;
        
        mindist_Y = all_mindist_Y(j-1, i);
        param_Y = all_param_Y(j-1, i);
        AccumulatedDistanceY(i,j) = mindist_Y; 
        ParamY(i,j) = param_Y;

        % X-Interpolation Pfadwahl
        minCost = Inf;
        if fkt_isInterpolation(ParamX(i,j))
            if AccumulatedDistanceX(i,j-1) < minCost && fkt_isInterpolation(ParamX(i,j-1)) && ParamX(i,j-1)<= ParamX(i,j)
                minCost = AccumulatedDistanceX(i,j-1);
                minParamsX{i,j} = [ParamX(i,j-1),0];
            elseif AccumulatedDistanceY(i-1,j) < minCost && fkt_isInterpolation(ParamY(i-1,j))
                minCost = AccumulatedDistanceY(i-1,j);
                minParamsX{i,j} = [0,ParamY(i-1,j)];
            elseif AccumulatedDistance(i-1,j) < minCost
                minCost = AccumulatedDistance(i-1,j);
                minParamsX{i,j} = [0,1];
            elseif AccumulatedDistance(i-1,j-1) < minCost && ~fkt_isInterpolation(ParamX(i,j-1)) && ~fkt_isInterpolation(ParamY(i-1,j)) && fkt_euclDist(i-1,j-1,X,Y) <= fkt_euclDist(i-1,j,X,Y)
                minCost = AccumulatedDistance(i-1,j-1);
                minParamsX{i,j} = [0, 0];
            end
        end
        AccumulatedDistanceX(i,j) = AccumulatedDistanceX(i,j) + minCost;

        % Y-Interpolation Pfadwahl
        minCost = Inf;
        if fkt_isInterpolation(ParamY(i, j))
            if (AccumulatedDistanceX(i, j - 1) < minCost && fkt_isInterpolation(ParamX(i, j - 1)))
                minCost = AccumulatedDistanceX(i, j - 1);
                minParamsY{i, j} = [ParamX(i, j - 1), 0];
            elseif (AccumulatedDistanceY(i - 1, j) < minCost && fkt_isInterpolation(ParamY(i - 1, j)) && ParamY(i - 1, j) <= ParamY(i, j))
                minCost = AccumulatedDistanceY(i - 1, j);
                minParamsY{i, j} = [0, ParamY(i - 1, j)];
            elseif (AccumulatedDistance(i, j - 1) < minCost)
                minCost = AccumulatedDistance(i, j - 1);
                minParamsY{i, j} = [1, 0];
            elseif (AccumulatedDistance(i - 1, j - 1) < minCost && ~fkt_isInterpolation(ParamX(i, j - 1)) && ~fkt_isInterpolation(ParamY(i - 1, j)) && fkt_euclDist(i - 1, j - 1,X,Y) <= fkt_euclDist(i, j - 1,X,Y))
                minCost = AccumulatedDistance(i - 1, j - 1);
                minParamsY{i, j} = [0, 0];
            end
        end
        AccumulatedDistanceY(i, j) = AccumulatedDistanceY(i, j) + minCost;

        % Finale Kostenwahl
        minCost = inf;
        if (fkt_isInterpolation(ParamX(i, j)) && AccumulatedDistanceX(i, j) < minCost)
            minCost = AccumulatedDistanceX(i, j);
            minParams{i, j} = [ParamX(i, j), 1];
        end
        if (fkt_isInterpolation(ParamY(i, j)) && AccumulatedDistanceY(i, j) < minCost)
            minCost = AccumulatedDistanceY(i, j);
            minParams{i, j} = [1, ParamY(i, j)];
        end
        if (AccumulatedDistanceX(i, j - 1) < minCost && fkt_isInterpolation(ParamX(i, j - 1)) && fkt_euclDist(i, j,X,Y) <= fkt_euclDist(i, j - 1,X,Y) && ~fkt_isInterpolation(ParamY(i, j)) && (ParamX(i, j) < ParamX(i, j - 1) || ~fkt_isInterpolation(ParamX(i, j))))
            minCost = AccumulatedDistanceX(i, j - 1);
            minParams{i, j} = [ParamX(i, j - 1), 0];
        end
        if (AccumulatedDistanceY(i - 1, j) < minCost && fkt_isInterpolation(ParamY(i - 1, j)) && fkt_euclDist(i, j,X,Y) <= fkt_euclDist(i - 1, j,X,Y) && ~fkt_isInterpolation(ParamX(i, j)) && (ParamY(i, j) < ParamY(i - 1, j) || ~fkt_isInterpolation(ParamY(i, j))))
            minCost = AccumulatedDistanceY(i - 1, j);
            minParams{i, j} = [0, ParamY(i - 1, j)];
        end
        if (AccumulatedDistance(i, j - 1) < minCost && ~fkt_isInterpolation(ParamY(i, j)))
            minCost = AccumulatedDistance(i, j - 1);
            minParams{i, j} = [1, 0];
        end
        if (AccumulatedDistance(i - 1, j) < minCost && ~fkt_isInterpolation(ParamX(i, j)))
            minCost = AccumulatedDistance(i - 1, j);
            minParams{i, j} = [0, 1];
        end
        if (AccumulatedDistance(i - 1, j - 1) < minCost && fkt_euclDist(i, j,X,Y) <= fkt_euclDist(i - 1, j,X,Y) && fkt_euclDist(i, j,X,Y) <= fkt_euclDist(i, j - 1,X,Y) && ~fkt_isInterpolation(ParamX(i, j)) && ~fkt_isInterpolation(ParamY(i, j)) && ~fkt_isInterpolation(ParamX(i, j - 1)) && ~fkt_isInterpolation(ParamY(i - 1, j)) && fkt_euclDist(i - 1, j - 1,X,Y) <= fkt_euclDist(i - 1, j,X,Y) && fkt_euclDist(i - 1, j - 1,X,Y) <= fkt_euclDist(i, j - 1,X,Y))
            minCost = AccumulatedDistance(i - 1, j - 1);
            minParams{i, j} = [0, 0];
        end
        assert(minCost ~= inf);
        AccumulatedDistance(i, j) = minCost + fkt_euclDist(i, j,X,Y);    
    end
end
elapsed = toc;
fprintf('Kostenmatrix erstellt in %.3f Sekunden\n', elapsed);

fprintf('Backtracking...\n');
%% Backtracking/ dynamische Programmierung (unverändert)
i = M;
j = N;

MappingIndexes = [i j];        
result = [X(i,:) Y(j,:)];

lastParam = [0 0];

while i > 1 || j > 1
    % Eckpunkt
    if (lastParam(1) == 0 && lastParam(2) == 0) || (lastParam(1) == 0 && lastParam(2) == 1) || (lastParam(1) == 1 && lastParam(2) == 0)
        lastParam = minParams{i,j};
    % Interpolationspunkt X-Kante
    elseif (lastParam(1) > 0 && lastParam(2) == 0) || (lastParam(1) > 0 && lastParam(2) == 1)
        lastParam = minParamsX{i,j};
    % Interpolationspunkt Y-Kante
    elseif (lastParam(1) == 0 && lastParam(2) > 0) || (lastParam(1) == 1 && lastParam(2) > 0) 
        lastParam = minParamsY{i,j};
    else
        error('Ungültiger Zustand.');
    end

    if i == 1
        result = [X(1,:), fkt_interpolate(Y(j - 1,:), Y(j,:), lastParam(2)); result];

    elseif j == 1
        result = [fkt_interpolate(X(i - 1,:), X(i,:), lastParam(1)), Y(1,:); result];
    else
        result = [fkt_interpolate(X(i - 1,:), X(i,:), lastParam(1)), fkt_interpolate(Y(j - 1,:), Y(j,:), lastParam(2)); result];
    end
    
    MappingIndexes = [i - 1 + lastParam(1), j - 1 + lastParam(2); MappingIndexes];
    assert(i - 1 + lastParam(1) >= 0);
    assert(j - 1 + lastParam(2) >= 0);

    if lastParam(1) == 0
        i = i - 1;
    end
    if lastParam(2) == 0
        j = j - 1;
    end
end

% Indizes der interpolierten Bahnpunkte
ix = MappingIndexes(:,1);       
iy = MappingIndexes(:,2);

if size(result,2) == 6         % Orientierungen und Positionen 
    dtw_X = result(:,[1 2 3]);
    dtw_Y = result(:, [4 5 6]);
elseif size(result,2) == 2     % Geschwindigkeiten 
    dtw_X = result(:,1);
    dtw_Y = result(:,2);
elseif size(result,2) == 4     % Geschwindigkeiten mit Zeitstempel
    dtw_X = result(:,[1 2]);
    dtw_Y = result(:, [3 4]);
elseif size(result,2) == 8     % Orientierungen und Positionen mit Zeitstempel
    dtw_X = result(:,[1 2 3 4]);
    dtw_Y = result(:, [5 6 7 8]);
end

% Distanzen zwischen den interpolierten Bahnen
distances = zeros(length(dtw_X),1);
for i = 1:length(dtw_X)
    dist = fkt_euclDist(i,i,dtw_X,dtw_Y);
    distances(i,1) = dist;
end

% maximale und mittlere Distanz 
maxDistance = max(distances);
averageDistance = mean(distances);
minDistance = min(distances);

fprintf('✅ Fertig! Max: %.3f mm, Avg: %.3f mm, Min: %.3f mm\n', maxDistance, averageDistance, minDistance);

%% Visualisierung (unverändert)
if pflag
    fprintf('Erstelle Visualisierungen...\n');
    
    % Farben
    blau = [0 0.4470 0.7410];
    rot = [0.78 0 0];
    c1 = [0 0.4470 0.7410];

    % 2D Visualisierung der akkumulierten Kosten samt Mapping 
    figure('Name','SelectiveInterpolationDTW OPT - Kostenkarte und optimaler Pfad','NumberTitle','off');
    hold on
    imagesc(AccumulatedDistance)
    colormap("turbo");
    colorb = colorbar;
    colorb.Label.String = 'Akkumulierte Kosten';
    plot(iy, ix,"-w","LineWidth",1)
    xlabel('Pfad Y [Index]');
    ylabel('Pfad X [Index]');
    axis([min(iy) max(iy) 1 max(ix)]);
    set(gca,'FontSize',10,'YDir', 'normal');

    % Plot der beiden Bahnen und Zuordnung
    figure('Color','white','Name','SelectiveInterpolationDTW OPT - Zuordnung der Bahnpunkte','NumberTitle','off')
    hold on;
    grid on;
    box on;
    plot3(dtw_X(:,1),dtw_X(:,2),dtw_X(:,3),Color= rot,LineWidth=1.5,Marker = "square",MarkerFaceColor=rot,MarkerSize=4);
    plot3(dtw_Y(:,1),dtw_Y(:,2),dtw_Y(:,3),Color= c1,LineWidth=1.5,Marker = "o",MarkerFaceColor= c1,MarkerSize=4);
    for i = 1:length(dtw_X)
        line([dtw_Y(i,1),dtw_X(i,1)],[dtw_Y(i,2),dtw_X(i,2)],[dtw_Y(i,3),dtw_X(i,3)],'Color','black')
    end
    view(300, 40)
    legend({'Sollbahn','Istbahn','Abweichung'},'Location','northeast',"FontWeight", "bold")
    xlabel("x [mm]","FontWeight","bold")
    ylabel("y [mm]","FontWeight","bold")
    zlabel("z [mm]","FontWeight","bold")
    hold off
    axis padded
end

end

%% ============================================================================
%% HILFSFUNKTIONEN - VECTORISIERT
%% ============================================================================

function [mindist_matrix, param_matrix] = compute_all_edge_distances(Path1, Path2, edge_type)
    % Berechnet alle Kanten-zu-Punkt Distanzen auf einmal (vectorisiert)
    % Path1: Pfad dessen Kanten betrachtet werden [M x 3]
    % Path2: Pfad dessen Punkte gegen Kanten geprüft werden [N x 3]
    % edge_type: 'X' oder 'Y' (für Debugging)
    
    M = size(Path1, 1);
    N = size(Path2, 1);
    
    mindist_matrix = zeros(M-1, N);
    param_matrix = zeros(M-1, N);
    
    % Für jede Kante in Path1
    for i = 1:(M-1)
        seg_start = Path1(i, :);      % [1 x 3]
        seg_end = Path1(i+1, :);      % [1 x 3]
        
        % Vectorisiert über alle Punkte in Path2
        % Repliziere seg_start und seg_end für alle Path2-Punkte
        seg_start_rep = repmat(seg_start, N, 1);  % [N x 3]
        seg_end_rep = repmat(seg_end, N, 1);      % [N x 3]
        
        % Richtungsvektor der Kante
        v = seg_end_rep - seg_start_rep;          % [N x 3]
        
        % Vektoren von Kantenstart zu allen Path2-Punkten
        w = Path2 - seg_start_rep;                % [N x 3]
        
        % Projektionsparameter t (vectorisiert)
        c1 = sum(w .* v, 2);                      % [N x 1]
        c2 = sum(v .* v, 2);                      % [N x 1]
        
        % Verhindere Division durch 0
        c2(c2 < eps) = eps;
        
        % Parameter zwischen 0 und 1 clippen
        t = max(0, min(1, c1 ./ c2));             % [N x 1]
        
        % Projizierte Punkte auf Kante
        proj_points = seg_start_rep + t .* v;     % [N x 3]
        
        % Distanzen von Path2-Punkten zu projizierten Punkten
        dists = vecnorm(Path2 - proj_points, 2, 2);  % [N x 1]
        
        % Speichere Ergebnisse
        mindist_matrix(i, :) = dists';
        param_matrix(i, :) = t';
    end
end

function [mindists, params] = fkt_minDistParam_vectorized_edge(seg_starts, seg_ends, points)
    % Vectorisierte Version für mehrere Kanten zu einem Punkt
    % seg_starts: [K x 3] Startpunkte der Kanten
    % seg_ends: [K x 3] Endpunkte der Kanten  
    % points: [K x 3] Testpunkte (einer pro Kante)
    
    K = size(seg_starts, 1);
    
    v = seg_ends - seg_starts;                % [K x 3]
    w = points - seg_starts;                  % [K x 3]
    
    c1 = sum(w .* v, 2);                      % [K x 1]
    c2 = sum(v .* v, 2);                      % [K x 1]
    
    c2(c2 < eps) = eps;
    
    params = max(0, min(1, c1 ./ c2));       % [K x 1]
    
    proj_points = seg_starts + params .* v;   % [K x 3]
    mindists = vecnorm(points - proj_points, 2, 2);  % [K x 1]
end
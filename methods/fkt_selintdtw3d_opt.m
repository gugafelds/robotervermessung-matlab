function [distances, maxDistance, averageDistance, AccumulatedDistance, dtw_X, dtw_Y, MappingIndexes, ix, iy, segment_ids_out] = fkt_selintdtw3d_opt(X, Y, segment_ids_X, segment_ids_Y, pflag)
%FKT_SELINTDTW3D_OPT Optimierte Version mit segment_id Support
%
%   INPUT:
%   X, Y: [N x 3] Arrays mit [x, y, z] Koordinaten (double)
%   segment_ids_X, segment_ids_Y: [N x 1] cell arrays mit segment_id strings
%   pflag: Visualisierung (optional)
%
%   OUTPUT:
%   segment_ids_out: segment_id für jeden interpolierten Punkt

if nargin < 5
    pflag = false;
end

% Check ob segment_ids vorhanden
has_segment_id = nargin >= 4 && ~isempty(segment_ids_X) && ~isempty(segment_ids_Y);

% Initialisierung
M = size(X, 1);
N = size(Y, 1);
AccumulatedDistance = zeros(M,N);
AccumulatedDistanceX = zeros(M,N);
AccumulatedDistanceY = zeros(M,N);
ParamX = zeros(M,N);
ParamY = zeros(M,N);
minParams = cell(M,N);
minParamsX = cell(M,N);
minParamsY = cell(M,N);

if has_segment_id
    fprintf('segment_id wird durch DTW propagiert\n');
end
fprintf('Initialisiere Matrizen: %d x %d\n', M, N);

% Anfangsbedingungen
AccumulatedDistance(1,1) = fkt_euclDist(1,1,X,Y);
AccumulatedDistanceX(1,1) = Inf;
AccumulatedDistanceY(1,1) = Inf;
ParamX(1,1) = NaN;
ParamY(1,1) = NaN;
minParams{1,1} = [NaN, NaN];
minParamsX{1,1} = [NaN, NaN];
minParamsY{1,1} = [NaN, NaN];

%% Vectorisierte Startwerte X
fprintf('Berechne Startwerte X (vectorisiert)...\n');
if M > 1
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

%% Vectorisierte Startwerte Y
fprintf('Berechne Startwerte Y (vectorisiert)...\n');
if N > 1
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

%% Vectorisierte Hauptberechnung
fprintf('Berechne alle X-Kanten Distanzen (vectorisiert)...\n');
[all_mindist_X, all_param_X] = compute_all_edge_distances(X, Y, 'X');

fprintf('Berechne alle Y-Kanten Distanzen (vectorisiert)...\n');
[all_mindist_Y, all_param_Y] = compute_all_edge_distances(Y, X, 'Y');

fprintf('Erstellung der akkumulierten Kostenmatrix (%d x %d)...\n', M, N);
tic;
for i = 2:M
    if mod(i, 50) == 0
        fprintf('  Fortschritt: %d/%d (%.1f%%)\n', i, M, 100*i/M);
    end
    
    for j = 2:N
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

%% Backtracking mit segment_id
fprintf('Backtracking...\n');
i = M;
j = N;

MappingIndexes = [i j];
result = [X(i,:) Y(j,:)];

if has_segment_id
    segment_ids_result = {segment_ids_X{i}};
end

lastParam = [0 0];

while i > 1 || j > 1
    if (lastParam(1) == 0 && lastParam(2) == 0) || (lastParam(1) == 0 && lastParam(2) == 1) || (lastParam(1) == 1 && lastParam(2) == 0)
        lastParam = minParams{i,j};
    elseif (lastParam(1) > 0 && lastParam(2) == 0) || (lastParam(1) > 0 && lastParam(2) == 1)
        lastParam = minParamsX{i,j};
    elseif (lastParam(1) == 0 && lastParam(2) > 0) || (lastParam(1) == 1 && lastParam(2) > 0) 
        lastParam = minParamsY{i,j};
    else
        error('Ungültiger Zustand.');
    end

    % Interpolation
    if i == 1
        new_row = [X(1,:), fkt_interpolate(Y(j-1,:), Y(j,:), lastParam(2))];
        if has_segment_id
            current_seg_id = segment_ids_Y{j};
        end
    elseif j == 1
        new_row = [fkt_interpolate(X(i-1,:), X(i,:), lastParam(1)), Y(1,:)];
        if has_segment_id
            current_seg_id = segment_ids_X{i};
        end
    else
        new_row = [fkt_interpolate(X(i-1,:), X(i,:), lastParam(1)), ...
                   fkt_interpolate(Y(j-1,:), Y(j,:), lastParam(2))];
        if has_segment_id
            current_seg_id = segment_ids_X{i};
        end
    end
    
    result = [new_row; result];
    
    if has_segment_id
        segment_ids_result = [current_seg_id; segment_ids_result];
    end
    
    MappingIndexes = [i - 1 + lastParam(1), j - 1 + lastParam(2); MappingIndexes];
    
    if lastParam(1) == 0
        i = i - 1;
    end
    if lastParam(2) == 0
        j = j - 1;
    end
end

ix = MappingIndexes(:,1);       
iy = MappingIndexes(:,2);

% Output
dtw_X = result(:, 1:3);
dtw_Y = result(:, 4:6);

if has_segment_id
    segment_ids_out = segment_ids_result;
    fprintf('✅ segment_id propagiert für %d Punkte\n', length(segment_ids_out));
else
    segment_ids_out = {};
end

% Distanzen
distances = zeros(length(dtw_X),1);
for i = 1:length(dtw_X)
    distances(i) = fkt_euclDist(i,i,dtw_X,dtw_Y);
end

maxDistance = max(distances);
averageDistance = mean(distances);
minDistance = min(distances);

fprintf('✅ Fertig! Max: %.3f mm, Avg: %.3f mm, Min: %.3f mm\n', maxDistance, averageDistance, minDistance);

%% Visualisierung
if pflag
    fprintf('Erstelle Visualisierungen...\n');
    
    blau = [0 0.4470 0.7410];
    rot = [0.78 0 0];

    figure('Name','SIDTW - Kostenkarte','NumberTitle','off');
    imagesc(AccumulatedDistance)
    colormap("turbo");
    colorbar;
    hold on;
    plot(iy, ix,"-w","LineWidth",1)
    xlabel('Y-Index');
    ylabel('X-Index');
    set(gca,'YDir', 'normal');

    figure('Name','SIDTW - Zuordnung','NumberTitle','off')
    plot3(dtw_X(:,1),dtw_X(:,2),dtw_X(:,3),'Color',rot,'LineWidth',1.5,'Marker',"square",'MarkerSize',4);
    hold on;
    plot3(dtw_Y(:,1),dtw_Y(:,2),dtw_Y(:,3),'Color',blau,'LineWidth',1.5,'Marker',"o",'MarkerSize',4);
    for i = 1:length(dtw_X)
        line([dtw_Y(i,1),dtw_X(i,1)],[dtw_Y(i,2),dtw_X(i,2)],[dtw_Y(i,3),dtw_X(i,3)],'Color','black')
    end
    view(300, 40);
    legend({'Soll','Ist','Abweichung'});
    xlabel("x [mm]");
    ylabel("y [mm]");
    zlabel("z [mm]");
    grid on;
end

end

%% Hilfsfunktionen
function [mindist_matrix, param_matrix] = compute_all_edge_distances(Path1, Path2, ~)
    M = size(Path1, 1);
    N = size(Path2, 1);
    
    mindist_matrix = zeros(M-1, N);
    param_matrix = zeros(M-1, N);
    
    for i = 1:(M-1)
        seg_start_rep = repmat(Path1(i,:), N, 1);
        seg_end_rep = repmat(Path1(i+1,:), N, 1);
        
        v = seg_end_rep - seg_start_rep;
        w = Path2 - seg_start_rep;
        
        c1 = sum(w .* v, 2);
        c2 = sum(v .* v, 2);
        c2(c2 < eps) = eps;
        
        t = max(0, min(1, c1 ./ c2));
        
        proj_points = seg_start_rep + t .* v;
        dists = vecnorm(Path2 - proj_points, 2, 2);
        
        mindist_matrix(i, :) = dists';
        param_matrix(i, :) = t';
    end
end

function [mindists, params] = fkt_minDistParam_vectorized_edge(seg_starts, seg_ends, points)
    v = seg_ends - seg_starts;
    w = points - seg_starts;
    
    c1 = sum(w .* v, 2);
    c2 = sum(v .* v, 2);
    c2(c2 < eps) = eps;
    
    params = max(0, min(1, c1 ./ c2));
    
    proj_points = seg_starts + params .* v;
    mindists = vecnorm(points - proj_points, 2, 2);
end
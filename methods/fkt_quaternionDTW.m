function [distances, segment_ids_out, qdtw_soll, qdtw_ist] = fkt_quaternionDTW(q_soll, q_ist, segment_ids_X, pflag)
    % QDTW mit geodätischer Distanz
    if nargin < 3, pflag = false; end

    q_soll = removeGimbalLockArtifacts(q_soll);
    
    M = size(q_soll, 1);
    N = size(q_ist, 1);
    
    % Normalisierung
    q_soll = q_soll ./ sqrt(sum(q_soll.^2, 2));
    q_ist = q_ist ./ sqrt(sum(q_ist.^2, 2));
    
    % Cost Matrix mit geodätischer Distanz
    D = zeros(M, N);
    for i = 1:M
        for j = 1:N
            dot_prod = abs(dot(q_soll(i,:), q_ist(j,:)));
            D(i,j) = 2 * acos(min(dot_prod, 1.0));
        end
    end
    D = rad2deg(D);
    
    % DTW Algorithmus (wie in fkt_selintdtw3d)
    AccDist = zeros(M, N);
    AccDist(1,1) = D(1,1);
    
    for i = 2:M
        AccDist(i,1) = AccDist(i-1,1) + D(i,1);
    end
    for j = 2:N
        AccDist(1,j) = AccDist(1,j-1) + D(1,j);
    end
    
    for i = 2:M
        for j = 2:N
            AccDist(i,j) = D(i,j) + min([AccDist(i-1,j), AccDist(i,j-1), AccDist(i-1,j-1)]);
        end
    end
    
    % Backtracking für Pfad
    [ix, iy] = backtrack(AccDist);
    
    % Distanzen entlang des optimalen Pfads
    distances = zeros(length(ix), 1);
    for k = 1:length(ix)
        distances(k) = D(ix(k), iy(k));
    end

    qdtw_soll = q_soll(ix,:);
    qdtw_ist = q_ist(iy,:);
    segment_ids_out = segment_ids_X(iy);
end

function [ix, iy] = backtrack(AccDist)
    [M, N] = size(AccDist);
    i = M; j = N;
    ix = i; iy = j;
    
    while i > 1 || j > 1
        if i == 1
            j = j - 1;
        elseif j == 1
            i = i - 1;
        else
            [~, idx] = min([AccDist(i-1,j), AccDist(i,j-1), AccDist(i-1,j-1)]);
            switch idx
                case 1, i = i - 1;
                case 2, j = j - 1;
                case 3, i = i - 1; j = j - 1;
            end
        end
        ix = [i; ix];
        iy = [j; iy];
    end
end
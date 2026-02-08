function dist = cED(seq1, seq2, mode, window_percent, best_so_far, align_rotation, normalize_dtw)
    % Standardwerte setzen (Legacy-Support)
    if nargin < 7, normalize_dtw = false; end
    if nargin < 6, align_rotation = false; end
    % window_percent und best_so_far werden hier nicht für die Rechnung 
    % benötigt, bleiben aber in der Signatur erhalten.

    % ================================================================
    % 1. RESAMPLING (Auf gleiche Länge bringen)
    % ================================================================
    n = size(seq1, 1);
    m = size(seq2, 1);
    
    if n ~= m
        % Wir resamplen seq2 auf die Länge von seq1
        % linspace erstellt neue Indizes, interp1 berechnet die Werte dazwischen
        old_idx = 1:m;
        new_idx = linspace(1, m, n);
        seq2 = interp1(old_idx, seq2, new_idx, 'linear');
    end

    % ================================================================
    % 2. NORMALISIERUNG & TRANSLATION (Analog zu deinem Code)
    % ================================================================
    if normalize_dtw || strcmp(mode, 'joint_states')
        seq1 = normalizeForDTW(seq1, mode);
        seq2 = normalizeForDTW(seq2, mode);
    end
    
    if strcmp(mode, 'position')
        seq1 = seq1 - seq1(1, :);
        seq2 = seq2 - seq2(1, :);
    end
    
    if align_rotation && strcmp(mode, 'position')
        [seq1, seq2] = alignRotation(seq1, seq2);
    end

    % ================================================================
    % 3. EUKLIDISCHE DISTANZ
    % ================================================================
    % Wir berechnen die Differenz für jeden Punkt
    diffs = seq1 - seq2;
    
    % Quadratwurzel der Summe der Quadrate über alle Dimensionen
    % 'fro' (Frobenius Norm) ist effizient für Matrizen
    dist = norm(diffs, 'fro');
end
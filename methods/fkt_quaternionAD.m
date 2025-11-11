function [qad_soll, qad_distances] = fkt_quaternionAD(q_soll, q_ist)
    % Normalisierung
    q_soll = removeGimbalLockArtifacts(q_soll);

    q_soll = q_soll ./ sqrt(sum(q_soll.^2, 2));
    q_ist = q_ist ./ sqrt(sum(q_ist.^2, 2));
    
    n_soll = size(q_soll, 1);
    n_ist = size(q_ist, 1);
    
    % Interpolation mit SLERP
    t_ist = linspace(1, n_soll, n_ist)';
    idx_low = floor(t_ist);
    idx_high = ceil(t_ist);
    frac = t_ist - idx_low;
    idx_high = min(idx_high, n_soll);
    
    q_interp = zeros(n_ist, 4);
    for i = 1:n_ist
        if idx_low(i) == idx_high(i)
            q_interp(i,:) = q_soll(idx_low(i),:);
        else
            q1 = q_soll(idx_low(i),:);
            q2 = q_soll(idx_high(i),:);
            
            dot_prod = sum(q1 .* q2);
            if dot_prod < 0
                q2 = -q2;
                dot_prod = -dot_prod;
            end
            
            if dot_prod > 0.9995
                q_interp(i,:) = q1 + frac(i) * (q2 - q1);
            else
                theta = acos(dot_prod);
                q_interp(i,:) = (sin((1-frac(i))*theta)/sin(theta)) * q1 + ...
                                (sin(frac(i)*theta)/sin(theta)) * q2;
            end
            q_interp(i,:) = q_interp(i,:) / norm(q_interp(i,:));
        end
    end
    
    % Geodätische Distanzen
    qad_distances = zeros(n_ist, 1);
    for i = 1:n_ist
        dot_product = abs(sum(q_interp(i,:) .* q_ist(i,:)));
        qad_distances(i) = 2 * acosd(min(dot_product, 1.0));
    end
    
    qad_soll = q_interp;
    indices_ist = (1:n_ist)';
end
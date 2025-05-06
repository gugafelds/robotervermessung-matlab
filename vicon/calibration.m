function [trafo_rot, trafo_trans, error_metrics] = calibration(data_cal_ist, data_cal_soll, plots)

    ist_time = double(string(data_cal_ist.timestamp));
    ist = [data_cal_ist.x_ist data_cal_ist.y_ist data_cal_ist.z_ist];
    soll_time = double(string(data_cal_soll.timestamp));
    soll_reference = [data_cal_soll.x_reached data_cal_soll.y_reached data_cal_soll.z_reached];
    
    % Berechne Abstände zwischen aufeinanderfolgenden Punkten
    diffs = diff(ist);
    dists = sqrt(sum(diffs.^2,2));
    
    % Plot der Punktabstände nur wenn plots true ist
    if plots
        figure('Name', 'Punktabstände', 'Color', 'white');
        plot(dists, 'LineWidth', 1.5);
        title('Abstände zwischen aufeinanderfolgenden Punkten');
        xlabel('Punkt Index');
        ylabel('Abstand [mm]');
        grid on;
        
        % Dynamischer Threshold basierend auf Statistiken
        dists_mean = mean(dists);
        dists_std = std(dists);
        threshold = min(0.15, dists_mean + 2*dists_std);
        yline(threshold, 'r--', ['Threshold: ' num2str(threshold, '%.3f')], 'LineWidth', 1.5);
    end
    
    dists_mean = mean(dists);
    dists_std = std(dists);
    threshold = min(0.15, dists_mean + 2*dists_std);
    
    % Bestimmung der naheliegensten Timestamps
    ist_base_idx = zeros(length(soll_time),1);
    for i = 1:length(ist_base_idx)
        [~,idx] = min(abs(soll_time(i)-ist_time));
        ist_base_idx(i) = idx;
    end
    
    ist_base_timestamps = ist_time(ist_base_idx);
    ist_base_points = ist(ist_base_idx,:);
    
    % Verbessertes Punktmittelungsverfahren mit Outlier-Erkennung
    ist_reference = zeros(length(ist_base_idx),3);
    valid_points = true(length(ist_base_idx),1);
    
    for i = 1:length(ist_base_idx)
        idx = ist_base_idx(i);
        buffer = [];
        
        while idx <= length(dists) && dists(idx) < threshold
            buffer(end + 1) = idx;
            idx = idx + 1;
        end
        
        if ~isempty(buffer)
            points = ist(buffer,:);
            distances = sqrt(sum((points - mean(points)).^2, 2));
            valid_buffer = ~isoutlier(distances, 'quartiles');
            ist_reference(i,:) = mean(points(valid_buffer,:));
        else
            valid_points(i) = false;
        end
    end
    
    % Entferne ungültige Punkte
    ist_reference = ist_reference(valid_points,:);
    soll_reference = soll_reference(valid_points,:);
    ist_base_points = ist_base_points(valid_points,:);
    
    % Plot der Punkte vor Transformation
    if plots
        figure('Name', 'Punkte vor Transformation', 'Color', 'white');
        scatter3(ist_reference(:,1), ist_reference(:,2), ist_reference(:,3), 50, 'b', 'filled', 'DisplayName', 'IST');
        hold on;
        scatter3(soll_reference(:,1), soll_reference(:,2), soll_reference(:,3), 50, 'r', 'filled', 'DisplayName', 'SOLL');
        xlabel('X [mm]'); ylabel('Y [mm]'); zlabel('Z [mm]');
        title('Punktwolken vor Transformation');
        legend;
        grid on;
        view(3);
    end
    
    % Qualitätscheck vor Transformation
    initial_diffs = ist_reference - ist_base_points;
    initial_dists = sqrt(sum(initial_diffs.^2,2));
    
    if mean(initial_dists) > 1.0 || std(initial_dists) > 0.5
        warning('Hohe Abweichungen in den Initialdaten: Mean=%.2f, Std=%.2f', ...
            mean(initial_dists), std(initial_dists));
    end
    
    % Koordinatentransformation
    soll_mean = mean(soll_reference);
    ist_mean = mean(ist_reference);
    soll_centered = (soll_reference-soll_mean);
    ist_centered = (ist_reference-ist_mean);
    
    % Plot der zentrierten Punkte
    if plots
        figure('Name', 'Zentrierte Punkte', 'Color', 'white');
        scatter3(soll_centered(:,1), soll_centered(:,2), soll_centered(:,3), 50, 'r', 'filled', 'DisplayName', 'SOLL centered');
        hold on;
        scatter3(ist_centered(:,1), ist_centered(:,2), ist_centered(:,3), 50, 'b', 'filled', 'DisplayName', 'IST centered');
        xlabel('X [mm]'); ylabel('Y [mm]'); zlabel('Z [mm]');
        title('Zentrierte Punktwolken');
        legend;
        grid on;
        view(3);
    end
    
    % Kovarianzmatrix und SVD
    H = soll_centered' * ist_centered;
    [U, S, V] = svd(H);
    
    % Überprüfe Konditionszahl der Transformation
    condition_number = max(diag(S)) / min(diag(S));
    if condition_number > 1000
        warning('Hohe Konditionszahl der Transformation: %.2f', condition_number);
    end
    
    % Berechne Transformation
    trafo_rot = V * U';
    trafo_trans = soll_mean - ist_mean * trafo_rot;
    
    % Transformiere Referenzpunkte
    ist_reference_transformed = ist_reference * trafo_rot + trafo_trans;
    
    % Plot der transformierten Punkte
    if plots
        figure('Name', 'Transformationsergebnis', 'Color', 'white');
        scatter3(ist_reference_transformed(:,1), ist_reference_transformed(:,2), ist_reference_transformed(:,3), 50, 'b', 'filled', 'DisplayName', 'IST transformed');
        hold on;
        scatter3(soll_reference(:,1), soll_reference(:,2), soll_reference(:,3), 50, 'r', 'filled', 'DisplayName', 'SOLL');
        xlabel('X [mm]'); ylabel('Y [mm]'); zlabel('Z [mm]');
        title('Ergebnis der Transformation');
        legend;
        grid on;
        view(3);
    end
    
    % Qualitätsmetrik nach Transformation
    final_diffs = ist_reference_transformed - soll_reference;
    final_dists = sqrt(sum(final_diffs.^2,2));
    
    % Plot der Transformationsfehler
    if plots
        figure('Name', 'Transformationsfehler', 'Color', 'white');
        subplot(2,1,1);
        plot(final_dists, 'LineWidth', 1.5);
        title('Abweichungen nach Transformation');
        xlabel('Punkt Index');
        ylabel('Abweichung [mm]');
        grid on;
        
        subplot(2,1,2);
        histogram(final_dists, 20, 'Normalization', 'probability');
        title('Verteilung der Abweichungen');
        xlabel('Abweichung [mm]');
        ylabel('Relative Häufigkeit');
        grid on;
    end
    
    error_metrics = struct();
    error_metrics.mean_error = mean(final_dists);
    error_metrics.max_error = max(final_dists);
    error_metrics.std_error = std(final_dists);
    
    if error_metrics.std_error > 1.0
        warning(['Transformation möglicherweise ungenau:\n' ...
                'Mittlerer Fehler: %.2f mm\n' ...
                'Max Fehler: %.2f mm\n' ...
                'Std Fehler: %.2f mm'], ...
                error_metrics.mean_error, ...
                error_metrics.max_error, ...
                error_metrics.std_error);
    end
    
end

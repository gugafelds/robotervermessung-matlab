function plotIncrementalContribution(csv_file)
    % PLOTINCREMENTALCONTRIBUTION - Plot incremental modality benefits from CSV
    %
    % Usage:
    %   plotIncrementalContribution('similarity/results/multimodality_ablation_2024-11-27T153045.csv')
    %
    % Input:
    %   csv_file - Path to CSV file with experiment results
    %
    % Output:
    %   Saves incremental contribution plot as PNG and PDF
    
    fprintf('=== Creating Incremental Contribution Plot ===\n');
    fprintf('Reading: %s\n', csv_file);
    
    % Read CSV
    results_table = readtable(csv_file, 'ReadRowNames', true);
    
    % Extract output directory
    [output_dir, ~, ~] = fileparts(csv_file);
    
    % Get all row names
    all_names = results_table.Properties.RowNames;
    
    % === Define Color Scheme ===
    % Joint = Blue shades, Position = Orange shades
    % Trajectory = Solid, Segment = Dashed
    colors = struct();
    colors.joint_traj = [0.2, 0.4, 0.8];      % Dark Blue (solid)
    colors.joint_seg = [0.5, 0.6, 0.9];       % Light Blue (dashed)
    colors.pos_traj = [0.8, 0.4, 0.2];        % Dark Orange (solid)
    colors.pos_seg = [0.9, 0.6, 0.4];         % Light Orange (dashed)
    
    % === Create Figure ===
    figure('Position', [100, 100, 1200, 600]);
    
    % === SUBPLOT 1: Joint-based (both levels) ===
    subplot(1, 2, 1);
    plotIncrementalSubplot(results_table, all_names, 'joint_states', colors, ...
        'Joint-Space: Incremental Modality Addition');
    
    % === SUBPLOT 2: Position-based (both levels) ===
    subplot(1, 2, 2);
    plotIncrementalSubplot(results_table, all_names, 'position', colors, ...
        'Cartesian-Space: Incremental Modality Addition');
    
    sgtitle('Incremental Benefit of Additional Modalities', ...
        'FontSize', 16, 'FontWeight', 'bold');
    
    % === Save ===
    output_file_png = fullfile(output_dir, 'incremental_modality_contribution.png');
    output_file_pdf = fullfile(output_dir, 'incremental_modality_contribution.pdf');
    
    saveas(gcf, output_file_png);
    saveas(gcf, output_file_pdf);
    
    fprintf('✓ Saved: %s\n', output_file_png);
    fprintf('✓ Saved: %s\n', output_file_pdf);
end

function plotIncrementalSubplot(results_table, all_names, dtw_mode, colors, title_str)
    % Helper function to create one incremental subplot with both levels
    
    % Filter by DTW mode
    is_mode = strcmp(results_table.DTW_Mode, dtw_mode);
    mode_table = results_table(is_mode, :);
    mode_names = all_names(is_mode);
    
    if height(mode_table) == 0
        text(0.5, 0.5, 'No data', 'HorizontalAlignment', 'center');
        title(title_str, 'FontSize', 12, 'FontWeight', 'bold');
        return;
    end
    
    % Define expected order based on DTW mode
    if strcmp(dtw_mode, 'joint_states')
        % Joint + Orient/Meta combinations
        base_names = {'Joint only', 'Joint + Metadata', ...
                      'Joint + Orientation', 'Joint + Orient + Meta'};
        x_labels = {'Joint', '+Meta', '+Orient', '+Both'};
    else  % position
        % Position + Velocity/Meta combinations
        base_names = {'Position only', 'Pos + Metadata', ...
                      'Pos + Velocity', 'Pos + Vel + Meta'};
        x_labels = {'Position', '+Meta', '+Velocity', '+Both'};
    end
    
    % Separate by level
    is_traj = strcmp(mode_table.Level, 'Trajectory');
    is_seg = strcmp(mode_table.Level, 'Segment');
    
    traj_table = mode_table(is_traj, :);
    traj_names = mode_names(is_traj);
    
    seg_table = mode_table(is_seg, :);
    seg_names = mode_names(is_seg);
    
    % === Extract Trajectory values ===
    traj_values = nan(1, 4);
    for i = 1:length(base_names)
        match_idx = find(contains(traj_names, base_names{i}), 1);
        if ~isempty(match_idx)
            traj_values(i) = traj_table.Spearman(match_idx);
        end
    end
    
    % === Extract Segment values ===
    seg_values = nan(1, 4);
    for i = 1:length(base_names)
        % Handle _seg suffix
        match_idx = find(contains(seg_names, base_names{i}), 1);
        if ~isempty(match_idx)
            seg_values(i) = seg_table.Spearman(match_idx);
        end
    end
    
    % Check if we have any valid data
    has_traj = any(~isnan(traj_values));
    has_seg = any(~isnan(seg_values));
    
    if ~has_traj && ~has_seg
        text(0.5, 0.5, 'Incomplete data', 'HorizontalAlignment', 'center');
        title(title_str, 'FontSize', 12, 'FontWeight', 'bold');
        return;
    end
    
    % === Determine colors based on DTW mode ===
    if strcmp(dtw_mode, 'joint_states')
        traj_color = colors.joint_traj;
        seg_color = colors.joint_seg;
        marker = 'o';
    else
        traj_color = colors.pos_traj;
        seg_color = colors.pos_seg;
        marker = 's';
    end
    
    % === Plot ===
    x = 1:4;
    hold on;
    
    % Plot Trajectory (solid line)
    if has_traj
        valid_traj = ~isnan(traj_values);
        plot(x(valid_traj), traj_values(valid_traj), ...
            ['-' marker], 'LineWidth', 2, 'MarkerSize', 10, ...
            'Color', traj_color, 'MarkerFaceColor', traj_color, ...
            'DisplayName', 'Trajectory');
        
        % Add improvement percentages for trajectory
        baseline_traj = traj_values(1);
        if ~isnan(baseline_traj)
            for i = 2:4
                if valid_traj(i)
                    improvement = (traj_values(i) - baseline_traj) / baseline_traj * 100;
                    text(i, traj_values(i) + 0.01, sprintf('+%.1f%%', improvement), ...
                        'HorizontalAlignment', 'center', 'FontSize', 9, ...
                        'Color', traj_color, 'FontWeight', 'bold');
                end
            end
        end
    end
    
    % Plot Segment (dashed line)
    if has_seg
        valid_seg = ~isnan(seg_values);
        plot(x(valid_seg), seg_values(valid_seg), ...
            ['--' marker], 'LineWidth', 2, 'MarkerSize', 10, ...
            'Color', seg_color, 'MarkerFaceColor', seg_color, ...
            'DisplayName', 'Segment');
        
        % Add improvement percentages for segment
        baseline_seg = seg_values(1);
        if ~isnan(baseline_seg)
            for i = 2:4
                if valid_seg(i)
                    improvement = (seg_values(i) - baseline_seg) / baseline_seg * 100;
                    text(i, seg_values(i) - 0.01, sprintf('+%.1f%%', improvement), ...
                        'HorizontalAlignment', 'center', 'FontSize', 9, ...
                        'Color', seg_color, 'FontWeight', 'bold');
                end
            end
        end
    end
    
    % === Configure plot ===
    xticks(x);
    xticklabels(x_labels);
    ylabel('Spearman \rho', 'FontSize', 12);
    title(title_str, 'FontSize', 13, 'FontWeight', 'bold');
    
    % Set y-limits with padding
    all_values = [traj_values(~isnan(traj_values)), seg_values(~isnan(seg_values))];
    if ~isempty(all_values)
        y_min = min(all_values) * 0.92;
        y_max = max(all_values) * 1.08;
        ylim([y_min, y_max]);
    end
    
    grid on;
    legend('Location', 'southeast', 'FontSize', 10);
    hold off;
end
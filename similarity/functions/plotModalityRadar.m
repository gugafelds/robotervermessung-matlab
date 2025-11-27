function plotModalityRadar(csv_file)
    % PLOTMODALITYRADAR - Create radar chart from experiment results CSV
    %
    % Usage:
    %   plotModalityRadar('similarity/results/multimodality_ablation_2024-11-27T153045.csv')
    %
    % Input:
    %   csv_file - Path to CSV file with experiment results
    %
    % Output:
    %   Saves radar chart as PNG and PDF in same directory as CSV
    
    fprintf('=== Creating Radar Chart ===\n');
    fprintf('Reading: %s\n', csv_file);
    
    % Read CSV
    results_table = readtable(csv_file, 'ReadRowNames', true);
    
    % Extract output directory
    [output_dir, ~, ~] = fileparts(csv_file);
    
    % Get all row names
    all_names = results_table.Properties.RowNames;
    num_rows = height(results_table);
    
    % === Define Color Scheme ===
    % Joint = Blue shades, Position = Orange shades
    % Trajectory = Solid, Segment = Dashed
    colors = struct();
    colors.joint_traj = [0.2, 0.4, 0.8];      % Dark Blue (solid)
    colors.joint_seg = [0.5, 0.6, 0.9];       % Light Blue (dashed)
    colors.pos_traj = [0.8, 0.4, 0.2];        % Dark Orange (solid)
    colors.pos_seg = [0.9, 0.6, 0.4];         % Light Orange (dashed)
    
    % === Metrics to show ===
    metric_cols = {'Spearman', 'P_100', 'P_10', 'P_5', 'P_3', 'P_1'};
    metric_labels = {'\rho', 'P@100', 'P@10', 'P@5', 'P@3', 'P@1'};
    
    % Check which columns exist
    available_metrics = {};
    available_labels = {};
    for i = 1:length(metric_cols)
        if ismember(metric_cols{i}, results_table.Properties.VariableNames)
            available_metrics{end+1} = metric_cols{i};
            available_labels{end+1} = metric_labels{i};
        end
    end
    
    if isempty(available_metrics)
        error('No metrics available in CSV');
    end
    
    angles = linspace(0, 2*pi, length(available_metrics)+1);
    
    % === Create Figure ===
    figure('Position', [100, 100, 1000, 800]);
    ax = polaraxes;
    hold on;
    
    % === Plot all rows ===
    for i = 1:num_rows
        % Get level and mode
        level = results_table.Level{i};
        dtw_mode = results_table.DTW_Mode{i};
        
        % Determine color and line style
        if strcmp(dtw_mode, 'joint_states') && strcmp(level, 'Trajectory')
            color = colors.joint_traj;
            line_style = '-';
            marker = 'o';
        elseif strcmp(dtw_mode, 'joint_states') && strcmp(level, 'Segment')
            color = colors.joint_seg;
            line_style = '--';
            marker = 'o';
        elseif strcmp(dtw_mode, 'position') && strcmp(level, 'Trajectory')
            color = colors.pos_traj;
            line_style = '-';
            marker = 's';
        else  % position + segment
            color = colors.pos_seg;
            line_style = '--';
            marker = 's';
        end
        
        % Extract metrics
        data = zeros(1, length(available_metrics));
        for j = 1:length(available_metrics)
            val = results_table.(available_metrics{j})(i);
            if isnan(val)
                val = 0;  % Replace NaN with 0 for plotting
            end
            data(j) = val;
        end
        
        % Close the loop
        data = [data, data(1)];
        
        % Clean up row name (remove _seg suffix)
        clean_name = strrep(all_names{i}, '_seg', '');
        
        % Create display name with level and mode
        if strcmp(level, 'Segment')
            display_name = sprintf('%s (Seg)', clean_name);
        else
            display_name = clean_name;
        end
        
        % Plot
        polarplot(angles, data, [line_style marker], 'LineWidth', 1.5, ...
            'MarkerSize', 6, 'Color', color, ...
            'DisplayName', display_name);
    end
    
    % === Configure axes ===
    ax.ThetaTick = rad2deg(angles(1:end-1));
    ax.ThetaTickLabel = available_labels;
    ax.RLim = [0, 1];
    
    % === Title and Legend ===
    title('Modality Contribution', ...
        'FontSize', 12, 'FontWeight', 'bold');
    
    % Create custom legend with grouping
    lgd = legend('Location', 'southoutside', 'Orientation', 'vertical', ...
        'FontSize', 9, 'NumColumns', 2);
    lgd.Title.String = 'Configuration (Solid=Traj, Dashed=Seg | Circle=Joint, Square=Pos)';
    
    % === Save ===
    %output_file_png = fullfile(output_dir, 'radar_modality_contribution.png');
    %output_file_pdf = fullfile(output_dir, 'radar_modality_contribution.pdf');
    
    %saveas(gcf, output_file_png);
    %saveas(gcf, output_file_pdf);
    
    %fprintf('✓ Saved: %s\n', output_file_png);
    %fprintf('✓ Saved: %s\n', output_file_pdf);
end
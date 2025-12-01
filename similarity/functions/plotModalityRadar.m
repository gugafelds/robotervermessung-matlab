function plotModalityRadar(csv_file)
    % PLOTMODALITYRADAR - Create radar chart from experiment results CSV
    %
    % Usage:
    %   plotModalityRadar('similarity/results/embedding_validation_2025-11-28.csv')
    %
    % Input:
    %   csv_file - Path to CSV file with experiment results
    %
    % Output:
    %   Saves radar chart as PNG and PDF in same directory as CSV
    
    fprintf('=== Creating Radar Chart ===\n');
    fprintf('Reading: %s\n', csv_file);
    
    % Read CSV (NEW: no row names, standard table format)
    results_table = readtable(csv_file);
    
    % Extract output directory
    [output_dir, ~, ~] = fileparts(csv_file);
    
    % Get number of rows
    num_rows = height(results_table);
    
    fprintf('Found %d rows in CSV\n', num_rows);
    
    % === Define Color Scheme ===
    % Joint = Blue shades, Position = Orange shades
    % Trajectory = Solid, Segment = Dashed
    colors = struct();
    colors.joint_traj = [0.2, 0.4, 0.8];      % Dark Blue (solid)
    colors.joint_seg = [0.5, 0.6, 0.9];       % Light Blue (dashed)
    colors.pos_traj = [0.8, 0.4, 0.2];        % Dark Orange (solid)
    colors.pos_seg = [0.9, 0.6, 0.4];         % Light Orange (dashed)
    
    % === Metrics to show (NEW column names) ===
    metric_cols = {'Spearman', 'P_50', 'P_10', 'P_5', 'P_3', 'P_1'};
    metric_labels = {'\rho', 'P@50', 'P@10', 'P@5', 'P@3', 'P@1'};
    
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
        error('No metrics available in CSV. Available columns: %s', ...
            strjoin(results_table.Properties.VariableNames, ', '));
    end
    
    fprintf('Using metrics: %s\n', strjoin(available_metrics, ', '));
    
    angles = linspace(0, 2*pi, length(available_metrics)+1);
    
    % === Create Figure ===
    figure('Position', [100, 100, 1200, 900]);
    ax = polaraxes;
    hold on;
    
    % === Group rows by configuration (NEW: based on Embedding_Config + Weight_Mode) ===
    % We'll plot one line per unique combination
    
    unique_configs = unique(results_table.Embedding_Config);
    unique_weights = unique(results_table.Weight_Mode);
    
    fprintf('Unique embedding configs: %d\n', length(unique_configs));
    fprintf('Unique weight modes: %d\n', length(unique_weights));
    
    % === Plot strategy: One line per (Embedding_Config × Weight_Mode × Level) ===
    
    plotted_count = 0;
    
    for emb_idx = 1:length(unique_configs)
        emb_config = unique_configs{emb_idx};
        
        for wm_idx = 1:length(unique_weights)
            weight_mode = unique_weights{wm_idx};
            
            % For each level (Trajectory, Segment)
            for level_idx = 1:2
                if level_idx == 1
                    level = 'Trajectory';
                else
                    level = 'Segment';
                end
                
                % Find matching rows (average over all queries for this config)
                mask = strcmp(results_table.Level, level) & ...
                       strcmp(results_table.Embedding_Config, emb_config) & ...
                       strcmp(results_table.Weight_Mode, weight_mode);
                
                if ~any(mask)
                    continue;
                end
                
                % Get DTW mode from first matching row
                matching_rows = results_table(mask, :);
                dtw_mode = matching_rows.DTW_Mode{1};
                
                % Average metrics across all matching rows (different queries)
                data = zeros(1, length(available_metrics));
                for j = 1:length(available_metrics)
                    vals = matching_rows.(available_metrics{j});
                    % Remove NaN values
                    vals = vals(~isnan(vals));
                    if ~isempty(vals)
                        data(j) = mean(vals);
                    else
                        data(j) = 0;
                    end
                end
                
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
                
                % Close the loop
                data = [data, data(1)];
                
                % Create display name
                if strcmp(level, 'Segment')
                    display_name = sprintf('%s | %s (Seg)', emb_config, weight_mode);
                else
                    display_name = sprintf('%s | %s', emb_config, weight_mode);
                end
                
                % Plot
                polarplot(angles, data, [line_style marker], 'LineWidth', 1.5, ...
                    'MarkerSize', 6, 'Color', color, ...
                    'DisplayName', display_name);
                
                plotted_count = plotted_count + 1;
            end
        end
    end
    
    fprintf('Plotted %d lines\n', plotted_count);
    
    % === Configure axes ===
    ax.ThetaTick = rad2deg(angles(1:end-1));
    ax.ThetaTickLabel = available_labels;
    ax.RLim = [0, 1];
    ax.FontSize = 10;
    
    % === Title and Legend ===
    title('Embedding Validation: Modality Performance', ...
        'FontSize', 14, 'FontWeight', 'bold');
    
    % Create custom legend with grouping
    lgd = legend('Location', 'southoutside', 'Orientation', 'vertical', ...
        'FontSize', 8, 'NumColumns', 3);
    lgd.Title.String = 'Config | Weights (Solid=Traj, Dashed=Seg | Circle=Joint, Square=Pos)';
    
end
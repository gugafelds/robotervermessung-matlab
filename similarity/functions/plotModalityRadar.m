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
    
    % === Metrics to show (DYNAMIC: detect P@k automatically) ===
    
    % Always include these base metrics
    metric_cols = {'Spearman', 'P_10', 'P_5', 'P_3', 'P_1'};
    metric_labels = {'\rho', 'P@10', 'P@5', 'P@3', 'P@1'};
    
    % Find P@k column dynamically (e.g., P_50, P_100, P_250, etc.)
    all_cols = results_table.Properties.VariableNames;
    p_k_pattern = '^P_\d+$';  % Matches P_50, P_100, P_250, etc.
    
    for i = 1:length(all_cols)
        col_name = all_cols{i};
        if ~isempty(regexp(col_name, p_k_pattern, 'once'))
            % Extract the k value (e.g., "P_50" -> 50)
            k_value = str2double(strrep(col_name, 'P_', ''));
            
            % Only add if it's NOT already in the list (not P_10, P_5, etc.)
            if k_value > 10 && ~ismember(col_name, metric_cols)
                % Insert at position 2 (after Spearman, before P@10)
                metric_cols = [metric_cols(1), {col_name}, metric_cols(2:end)];
                metric_labels = [metric_labels(1), {sprintf('P@%d', k_value)}, metric_labels(2:end)];
                fprintf('Detected dynamic P@k metric: %s\n', col_name);
                break;  % Only use the first one found
            end
        end
    end
    
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
    
    % === Get unique embedding configs ===
    unique_configs = unique(results_table.Embedding_Config);
    num_configs = length(unique_configs);
    
    fprintf('Found %d embedding configs\n', num_configs);
    
    % === Create Figure with Subplots ===
    figure('Position', [100, 100, 1400, 1000]);
    
    % Determine subplot layout
    if num_configs <= 2
        rows = 1;
        cols = num_configs;
    elseif num_configs <= 4
        rows = 2;
        cols = 2;
    elseif num_configs <= 6
        rows = 2;
        cols = 3;
    else
        rows = ceil(num_configs / 3);
        cols = 3;
    end
    
    % === Plot each embedding config in its own subplot ===
    for emb_idx = 1:num_configs
        emb_config = unique_configs{emb_idx};
        
        % Create subplot
        subplot(rows, cols, emb_idx, polaraxes);
        ax = gca;
        hold on;
        
        fprintf('\nProcessing subplot %d: %s\n', emb_idx, emb_config);
        
        % === Plot 4 lines: Joint-Traj, Joint-Seg, Pos-Traj, Pos-Seg ===
        % (Averaged over all weight modes)
        
        combinations = {
            'joint_states', 'Trajectory', colors.joint_traj, '-', 'o', 'Joint (Traj)';
            'joint_states', 'Segment',    colors.joint_seg,  '--', 'o', 'Joint (Seg)';
            'position',     'Trajectory', colors.pos_traj,   '-', 's', 'Position (Traj)';
            'position',     'Segment',    colors.pos_seg,    '--', 's', 'Position (Seg)'
        };
        
        for comb_idx = 1:size(combinations, 1)
            dtw_mode = combinations{comb_idx, 1};
            level = combinations{comb_idx, 2};
            color = combinations{comb_idx, 3};
            line_style = combinations{comb_idx, 4};
            marker = combinations{comb_idx, 5};
            display_name = combinations{comb_idx, 6};
            
            % Find all matching rows (this emb_config + this mode + this level)
            % Average over ALL weight modes
            mask = strcmp(results_table.Embedding_Config, emb_config) & ...
                   strcmp(results_table.DTW_Mode, dtw_mode) & ...
                   strcmp(results_table.Level, level);
            
            if ~any(mask)
                fprintf('  No data for %s\n', display_name);
                continue;
            end
            
            matching_rows = results_table(mask, :);
            fprintf('  %s: %d rows (avg over weight modes & queries)\n', ...
                display_name, height(matching_rows));
            
            % Average metrics across all matching rows
            data = zeros(1, length(available_metrics));
            for j = 1:length(available_metrics)
                vals = matching_rows.(available_metrics{j});
                vals = vals(~isnan(vals));
                if ~isempty(vals)
                    data(j) = mean(vals);
                else
                    data(j) = 0;
                end
            end
            
            % Close the loop
            data = [data, data(1)];
            
            % Plot
            polarplot(angles, data, [line_style marker], 'LineWidth', 2, ...
                'MarkerSize', 8, 'Color', color, ...
                'DisplayName', display_name);
        end
        
        % === Configure subplot axes ===
        ax.ThetaTick = rad2deg(angles(1:end-1));
        ax.ThetaTickLabel = available_labels;
        ax.RLim = [0, 1];
        ax.FontSize = 9;
        
        % === Subplot Title ===
        title(emb_config, 'FontSize', 11, 'FontWeight', 'bold');
        
        % === Legend for this subplot ===
        legend('Location', 'southoutside', 'Orientation', 'horizontal', ...
            'FontSize', 8);
    end
    
    % === Main Title ===
    sgtitle('Embedding Validation: Performance by Architecture', ...
        'FontSize', 14, 'FontWeight', 'bold');
    
    % === Save ===
    %output_file_png = fullfile(output_dir, 'radar_embedding_validation.png');
    %output_file_pdf = fullfile(output_dir, 'radar_embedding_validation.pdf');
    
    %saveas(gcf, output_file_png);
    %saveas(gcf, output_file_pdf);
    
    %fprintf('✓ Saved: %s\n', output_file_png);
    %fprintf('✓ Saved: %s\n', output_file_pdf);
    
    fprintf('=== Radar Chart Complete ===\n');
end
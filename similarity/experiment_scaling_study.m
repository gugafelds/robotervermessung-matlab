%% ========================================================================
%  EXPERIMENT 3A: DATABASE SIZE SCALING STUDY
%  ========================================================================
%  Goal: Test scalability and stability of DTW vs Embeddings
%  
%  Tests: N = [100, 500, 1000, 2500, 5000, 10000]
%  
%  Metrics:
%    - Spearman ρ (should stay stable)
%    - DTW Time (should scale linearly O(N))
%    - P@K, P@10 (should stay stable)
%
%  Runtime: ~1-1.5 hours (6 experiments × 10-15 min)
%  Output: results/scaling_study_YYYY-MM-DD_HH-MM-SS.csv + plot
%  ========================================================================

clear; clc;

fprintf('========================================\n');
fprintf('EXPERIMENT 3A: DATABASE SIZE SCALING\n');
fprintf('========================================\n');
fprintf('Runtime: ~1-1.5 hours\n');
fprintf('========================================\n\n');

% === Base Configuration ===
base_config = struct();
base_config.query_bahn_id = '1763740056';
base_config.random_seed = 42;
base_config.weights = [1, 1, 1, 1, 1];  % All modalities (baseline)
base_config.n_coarse = 50;
base_config.n_fine = 250;
base_config.top_k_trajectories = 50;
base_config.dtw_mode = 'joint_states';

% === Experiment Matrix ===
sample_sizes = [200, 500, 1000, 2500, 5000, 10000, 20000];
num_experiments = length(sample_sizes);

% === Storage ===
all_results = cell(num_experiments, 1);

% === Run Experiments ===
experiment_start = tic;

for exp_idx = 1:num_experiments
    sample_size = sample_sizes(exp_idx);
    
    fprintf('\n╔════════════════════════════════════════╗\n');
    fprintf('║  EXPERIMENT %d/%d: N = %-7d         ║\n', exp_idx, num_experiments, sample_size);
    fprintf('╚════════════════════════════════════════╝\n');
    
    % Configure
    config = base_config;
    config.database_sample_size = sample_size;
    config.exp_name = sprintf('N=%d', sample_size);
    
    % Run
    exp_tic = tic;
    results = runExperiment(config);
    exp_time = toc(exp_tic);
    
    results.exp_runtime = exp_time;
    all_results{exp_idx} = results;
    
    fprintf('✓ Experiment %d/%d completed in %.1f minutes\n', ...
        exp_idx, num_experiments, exp_time/60);
    
    % Estimate remaining time
    avg_time_per_exp = toc(experiment_start) / exp_idx;
    remaining_time = avg_time_per_exp * (num_experiments - exp_idx);
    fprintf('  Estimated time remaining: %.1f minutes\n', remaining_time/60);
end

total_time = toc(experiment_start);

% === Create Results Table ===
fprintf('\n=== Creating Results Table ===\n');

results_matrix = zeros(num_experiments, 8);  % ⭐ +1 column for P@K

for i = 1:num_experiments
    r = all_results{i};
    results_matrix(i, :) = [
        r.database_size, ...
        r.spearman, ...
        r.p_at_k, ...        % ⭐ P@K added
        r.p_at_10, ...
        r.p_at_1, ...
        r.dtw_time, ...
        r.exp_runtime, ...
        r.dtw_time / r.database_size  % Time per trajectory
    ];
end

results_table = array2table(results_matrix, ...
    'VariableNames', {'N', 'Spearman', sprintf('P@%d', base_config.top_k_trajectories), 'P@10', 'P@1', 'DTW_Time', 'Total_Time', 'Time_per_Traj'});

% === Save Results with Timestamp ===
output_dir = 'similarity/results';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% ⭐ Create timestamp with date and time
timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd''T''HHmmss'));

output_file = fullfile(output_dir, sprintf('scaling_study_%s.csv', timestamp));
writetable(results_table, output_file);

fprintf('✓ Results saved to: %s\n\n', output_file);

% === Create Plots ===
fprintf('=== Creating Plots ===\n');

plots_dir = fullfile(output_dir, 'plots');
if ~exist(plots_dir, 'dir')
    mkdir(plots_dir);
end

figure('Position', [100, 100, 1200, 400]);

% Plot 1: DTW Time vs. N (Scalability)
subplot(1, 3, 1);
plot(results_matrix(:, 1), results_matrix(:, 6), '-o', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Database Size (N)', 'FontSize', 12);
ylabel('DTW Time (s)', 'FontSize', 12);
title('Scalability: DTW Time vs. N', 'FontSize', 13, 'FontWeight', 'bold');
grid on;

% Add linear fit
p = polyfit(results_matrix(:, 1), results_matrix(:, 6), 1);
hold on;
plot(results_matrix(:, 1), polyval(p, results_matrix(:, 1)), '--r', 'LineWidth', 1.5);
legend('Data', sprintf('Linear Fit: y = %.4fx + %.2f', p(1), p(2)), 'Location', 'northwest');

% Plot 2: Spearman ρ vs. N (Stability)
subplot(1, 3, 2);
plot(results_matrix(:, 1), results_matrix(:, 2), '-o', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Database Size (N)', 'FontSize', 12);
ylabel('Spearman ρ', 'FontSize', 12);
title('Quality: Spearman ρ vs. N', 'FontSize', 13, 'FontWeight', 'bold');
ylim([0, 1]);
grid on;

% Add mean line
mean_spearman = mean(results_matrix(:, 2));
hold on;
yline(mean_spearman, '--r', sprintf('Mean: %.4f', mean_spearman), 'LineWidth', 1.5);

% Plot 3: P@K and P@10 vs. N (Stability)
subplot(1, 3, 3);
plot(results_matrix(:, 1), results_matrix(:, 3), '-o', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', sprintf('P@%d', base_config.top_k_trajectories));
hold on;
plot(results_matrix(:, 1), results_matrix(:, 4), '-s', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'P@10');
xlabel('Database Size (N)', 'FontSize', 12);
ylabel('Precision', 'FontSize', 12);
title('Quality: Precision vs. N', 'FontSize', 13, 'FontWeight', 'bold');
ylim([0, 1]);
grid on;
legend('Location', 'best');

% Add mean lines
mean_p_k = mean(results_matrix(:, 3));
mean_p10 = mean(results_matrix(:, 4));
yline(mean_p_k, '--r', sprintf('Mean P@%d: %.3f', base_config.top_k_trajectories, mean_p_k), 'LineWidth', 1);
yline(mean_p10, '--b', sprintf('Mean P@10: %.3f', mean_p10), 'LineWidth', 1);

% Save plot with timestamp
plot_file = fullfile(plots_dir, sprintf('scaling_study_%s.png', timestamp));
saveas(gcf, plot_file);
fprintf('✓ Plot saved to: %s\n\n', plot_file);

% === Summary ===
fprintf('========================================\n');
fprintf('EXPERIMENT COMPLETED\n');
fprintf('========================================\n');
fprintf('Total Runtime: %.1f minutes\n', total_time/60);
fprintf('Results: %s\n', output_file);
fprintf('Plot:    %s\n\n', plot_file);

fprintf('--- KEY METRICS ---\n');
fprintf('Spearman ρ:  Mean=%.4f, Std=%.4f, Range=[%.4f, %.4f]\n', ...
    mean(results_matrix(:, 2)), std(results_matrix(:, 2)), ...
    min(results_matrix(:, 2)), max(results_matrix(:, 2)));
fprintf('P@%d:       Mean=%.3f, Std=%.3f, Range=[%.3f, %.3f]\n', ...
    base_config.top_k_trajectories, ...
    mean(results_matrix(:, 3)), std(results_matrix(:, 3)), ...
    min(results_matrix(:, 3)), max(results_matrix(:, 3)));
fprintf('P@10:        Mean=%.3f, Std=%.3f, Range=[%.3f, %.3f]\n', ...
    mean(results_matrix(:, 4)), std(results_matrix(:, 4)), ...
    min(results_matrix(:, 4)), max(results_matrix(:, 4)));
fprintf('Time/Traj:   Mean=%.4fs\n', mean(results_matrix(:, 8)));

fprintf('\n--- LINEARITY CHECK ---\n');
fprintf('Linear Fit: DTW_Time = %.4f × N + %.2f\n', p(1), p(2));
r_squared = 1 - sum((results_matrix(:, 6) - polyval(p, results_matrix(:, 1))).^2) / ...
                sum((results_matrix(:, 6) - mean(results_matrix(:, 6))).^2);
fprintf('R² = %.4f (1.0 = perfect linear scaling)\n', r_squared);

fprintf('\n========================================\n\n');
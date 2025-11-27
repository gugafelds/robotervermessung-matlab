%% ========================================================================
%  EXPERIMENT 1B: MULTI-MODALITY ABLATION STUDY
%  ========================================================================
%  Goal: Measure contribution of each modality to ranking quality
%  
%  Tests both Position-based and Joint-based DTW with relevant embeddings
%
%  Runtime: ~30-40 minutes (8 experiments × 4-5 min)
%  Output: results/multimodality_ablation_YYYY-MM-DD.csv
%  ========================================================================

clear; clc;

fprintf('========================================\n');
fprintf('EXPERIMENT: MULTI-MODALITY ABLATION\n');
fprintf('========================================\n');

% === Base Configuration ===
base_config = struct();
base_config.query_bahn_id = '1763740056';
base_config.database_sample_size = 1100;  % Fixed for fair comparison
base_config.random_seed = 42;
base_config.n_coarse = 50;
base_config.n_fine = 150;
base_config.top_k_trajectories = 100;

% === Experiment Matrix ===
% Format: {Name, DTW_Mode, [Pos, Joint, Orient, Vel, Meta]}
experiments = {
    % Joint-based DTW experiments
    'Joint only',              'joint_states',  [0,   1,     0,      0,   0];
    'Joint + Orientation',     'joint_states',  [0,   1,     1,      0,   0];
    'Joint + Metadata',        'joint_states',  [0,   1,     0,      0,   1];
    'Joint + Orient + Meta',   'joint_states',  [0,   1,     1,      0,   1];
    
    % Position-based DTW experiments
    'Position only',           'position',      [1,   0,     0,      0,   0];
    'Pos + Velocity',          'position',      [1,   0,     0,      1,   0];
    'Pos + Metadata',          'position',      [1,   0,     0,      0,   1];
    'Pos + Vel + Meta',        'position',      [1,   0,     0,      1,   1];
};

num_experiments = size(experiments, 1);

% === Storage ===
all_results = cell(num_experiments, 1);

% === Run Experiments ===
experiment_start = tic;

for exp_idx = 1:num_experiments
    exp_name = experiments{exp_idx, 1};
    dtw_mode = experiments{exp_idx, 2};
    weight_array = experiments{exp_idx, 3};
    
    fprintf('\n╔════════════════════════════════════════╗\n');
    fprintf('║  EXPERIMENT %d/%d: %-22s ║\n', exp_idx, num_experiments, exp_name);
    fprintf('╚════════════════════════════════════════╝\n');
    fprintf('DTW Mode: %s\n', dtw_mode);
    fprintf('Weights: Pos=%.1f, Joint=%.1f, Orient=%.1f, Vel=%.1f, Meta=%.1f\n', ...
        weight_array(1), weight_array(2), weight_array(3), weight_array(4), weight_array(5));
    
    % Configure
    config = base_config;
    config.dtw_mode = dtw_mode;
    config.exp_name = exp_name;
    config.weights = weight_array;  % ⭐ Als Array übergeben
    
    % Run
    exp_tic = tic;
    results = runExperiment(config);
    exp_time = toc(exp_tic);
    
    results.exp_runtime = exp_time;
    all_results{exp_idx} = results;
    
    fprintf('✓ Experiment %d/%d completed in %.1f minutes\n', ...
        exp_idx, num_experiments, exp_time/60);
end

total_time = toc(experiment_start);

%% ========================================================================
%  CREATE RESULTS TABLE (TRAJECTORY + SEGMENT LEVEL)
%  ========================================================================

fprintf('\n=== Creating Results Table ===\n');

% Trajectory-level results
traj_results_matrix = zeros(num_experiments, 12);
% Segment-level results
seg_results_matrix = zeros(num_experiments, 12);
exp_names = cell(num_experiments, 1);
dtw_modes = cell(num_experiments, 1);

for i = 1:num_experiments
    r = all_results{i};
    exp_names{i} = r.exp_name;
    dtw_modes{i} = r.dtw_mode;
    
    % Trajectory-level metrics
    traj_results_matrix(i, :) = [
        r.weight_pos, r.weight_joint, r.weight_orient, r.weight_vel, r.weight_meta, ...
        r.spearman, r.p_at_k, r.p_at_10, r.p_at_5, r.p_at_3, r.p_at_1, ...
        r.exp_runtime/60  % Runtime in minutes
    ];
    
    % Segment-level metrics (averaged across segments)
    seg_results_matrix(i, :) = [
        r.weight_pos, r.weight_joint, r.weight_orient, r.weight_vel, r.weight_meta, ...
        r.seg_spearman, r.seg_p_at_k, r.seg_p_at_10, r.seg_p_at_5, r.seg_p_at_3, r.seg_p_at_1, ...
        r.exp_runtime/60
    ];
end

query_ids_traj = repmat({base_config.query_bahn_id}, num_experiments, 1);
query_ids_seg = repmat({base_config.query_bahn_id}, num_experiments, 1);

% Create trajectory-level table
traj_table = array2table(traj_results_matrix, ...
    'VariableNames', {'W_Pos', 'W_Joint', 'W_Orient', 'W_Vel', 'W_Meta', ...
                      'Spearman', sprintf('P@%d', base_config.top_k_trajectories), ...
                      'P@10', 'P@5', 'P@3', 'P@1', 'Runtime_min'}, ...
    'RowNames', exp_names);

% Add columns
traj_table = addvars(traj_table, dtw_modes, 'Before', 'W_Pos', ...
    'NewVariableNames', 'DTW_Mode');
level_traj = repmat({'Trajectory'}, num_experiments, 1);
traj_table = addvars(traj_table, level_traj, 'Before', 'DTW_Mode', ...
    'NewVariableNames', 'Level');
% ⭐ NEU: Query Bahn ID hinzufügen
traj_table = addvars(traj_table, query_ids_traj, 'Before', 'Level', ...
    'NewVariableNames', 'Query_Bahn_ID');

% Create segment-level table
seg_table = array2table(seg_results_matrix, ...
    'VariableNames', {'W_Pos', 'W_Joint', 'W_Orient', 'W_Vel', 'W_Meta', ...
                      'Spearman', sprintf('P@%d', base_config.top_k_trajectories), ...
                      'P@10', 'P@5', 'P@3', 'P@1', 'Runtime_min'}, ...
    'RowNames', strcat(exp_names, '_seg'));

% Add columns
seg_table = addvars(seg_table, dtw_modes, 'Before', 'W_Pos', ...
    'NewVariableNames', 'DTW_Mode');
level_seg = repmat({'Segment'}, num_experiments, 1);
seg_table = addvars(seg_table, level_seg, 'Before', 'DTW_Mode', ...
    'NewVariableNames', 'Level');
% ⭐ NEU: Query Bahn ID hinzufügen
seg_table = addvars(seg_table, query_ids_seg, 'Before', 'Level', ...
    'NewVariableNames', 'Query_Bahn_ID');

% Combine both tables (trajectory rows + segment rows)
combined_table = [traj_table; seg_table];

%% ========================================================================
%  SAVE RESULTS
%  ========================================================================

output_dir = 'similarity/results';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd''T''HHmmss'));
output_file = fullfile(output_dir, sprintf('multimodality_ablation_%s.csv', timestamp));
writetable(combined_table, output_file, 'WriteRowNames', true);

fprintf('✓ Results saved to: %s\n\n', output_file);

%% ========================================================================
%  SUMMARY
%  ========================================================================

fprintf('========================================\n');
fprintf('EXPERIMENT COMPLETED\n');
fprintf('========================================\n');
fprintf('Total Runtime: %.1f minutes\n', total_time/60);
fprintf('Results saved to: %s\n\n', output_file);

fprintf('--- RANKING BY DTW MODE (TRAJECTORY LEVEL) ---\n\n');

% Joint-based ranking
joint_idx = find(strcmp(dtw_modes, 'joint_states'));
if ~isempty(joint_idx)
    fprintf('JOINT-BASED DTW:\n');
    [~, sort_idx] = sort(traj_results_matrix(joint_idx, 6), 'descend');
    for i = 1:length(joint_idx)
        idx = joint_idx(sort_idx(i));
        fprintf('%d. %-25s: ρ=%.4f, P@%d=%.3f, P@10=%.3f\n', ...
            i, exp_names{idx}, traj_results_matrix(idx, 6), ...
            base_config.top_k_trajectories, traj_results_matrix(idx, 7), ...
            traj_results_matrix(idx, 8));
    end
    fprintf('\n');
end

% Position-based ranking
pos_idx = find(strcmp(dtw_modes, 'position'));
if ~isempty(pos_idx)
    fprintf('POSITION-BASED DTW:\n');
    [~, sort_idx] = sort(traj_results_matrix(pos_idx, 6), 'descend');
    for i = 1:length(pos_idx)
        idx = pos_idx(sort_idx(i));
        fprintf('%d. %-25s: ρ=%.4f, P@%d=%.3f, P@10=%.3f\n', ...
            i, exp_names{idx}, traj_results_matrix(idx, 6), ...
            base_config.top_k_trajectories, traj_results_matrix(idx, 7), ...
            traj_results_matrix(idx, 8));
    end
    fprintf('\n');
end

fprintf('--- RANKING BY DTW MODE (SEGMENT LEVEL) ---\n\n');

% Joint-based ranking (segments)
if ~isempty(joint_idx)
    fprintf('JOINT-BASED DTW (SEGMENTS):\n');
    [~, sort_idx] = sort(seg_results_matrix(joint_idx, 6), 'descend');
    for i = 1:length(joint_idx)
        idx = joint_idx(sort_idx(i));
        fprintf('%d. %-25s: ρ=%.4f, P@10=%.3f, P@1=%.3f\n', ...
            i, exp_names{idx}, seg_results_matrix(idx, 6), ...
            seg_results_matrix(idx, 8), seg_results_matrix(idx, 11));
    end
    fprintf('\n');
end

% Position-based ranking (segments)
if ~isempty(pos_idx)
    fprintf('POSITION-BASED DTW (SEGMENTS):\n');
    [~, sort_idx] = sort(seg_results_matrix(pos_idx, 6), 'descend');
    for i = 1:length(pos_idx)
        idx = pos_idx(sort_idx(i));
        fprintf('%d. %-25s: ρ=%.4f, P@10=%.3f, P@1=%.3f\n', ...
            i, exp_names{idx}, seg_results_matrix(idx, 6), ...
            seg_results_matrix(idx, 8), seg_results_matrix(idx, 11));
    end
    fprintf('\n');
end

fprintf('--- KEY INSIGHTS ---\n');
[best_spearman, best_idx] = max(traj_results_matrix(:, 6));
fprintf('Best Overall (Trajectory): %s (%s) - ρ=%.4f\n', ...
    exp_names{best_idx}, dtw_modes{best_idx}, best_spearman);

[best_seg_spearman, best_seg_idx] = max(seg_results_matrix(:, 6));
fprintf('Best Overall (Segment):    %s (%s) - ρ=%.4f\n', ...
    exp_names{best_seg_idx}, dtw_modes{best_seg_idx}, best_seg_spearman);

% Best per mode
if ~isempty(joint_idx)
    [best_joint_rho, best_joint_i] = max(traj_results_matrix(joint_idx, 6));
    best_joint_idx = joint_idx(best_joint_i);
    fprintf('Best Joint (Trajectory):   %s - ρ=%.4f\n', exp_names{best_joint_idx}, best_joint_rho);
    
    [best_joint_seg_rho, best_joint_seg_i] = max(seg_results_matrix(joint_idx, 6));
    best_joint_seg_idx = joint_idx(best_joint_seg_i);
    fprintf('Best Joint (Segment):      %s - ρ=%.4f\n', exp_names{best_joint_seg_idx}, best_joint_seg_rho);
end

if ~isempty(pos_idx)
    [best_pos_rho, best_pos_i] = max(traj_results_matrix(pos_idx, 6));
    best_pos_idx = pos_idx(best_pos_i);
    fprintf('Best Position (Trajectory): %s - ρ=%.4f\n', exp_names{best_pos_idx}, best_pos_rho);
    
    [best_pos_seg_rho, best_pos_seg_i] = max(seg_results_matrix(pos_idx, 6));
    best_pos_seg_idx = pos_idx(best_pos_seg_i);
    fprintf('Best Position (Segment):    %s - ρ=%.4f\n', exp_names{best_pos_seg_idx}, best_pos_seg_rho);
end

fprintf('\n========================================\n\n');

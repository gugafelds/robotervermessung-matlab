%% MAIN ANALYSIS: Load Data and Create Plots
%  ========================================================================
%  This script:
%  1. Loads all experiment CSV files
%  2. Creates 3 main analysis plots
%  3. Prints key findings
%  ========================================================================

clear; clc;

%% === CONFIGURATION ===
% Filter settings
cfg.db_size = 7500;
cfg.top_k = 50;
cfg.dtw_norm = 0;
cfg.weight_modes_tab3 = {'Joint + All', 'Position + All'};

% Composite weights (summe = 1.0)
cfg.w.spearman = 1.0;
cfg.w.ndcg50_dtw = 0.0;
cfg.w.ndcg50_gt = 0.0;
cfg.w.r50_dtw = 0.0;
cfg.w.r50_gt = 0.0;

% Paths
cfg.output_folder = 'results';
cfg.figure_folder = 'figs';
cfg.save_figures = false;

%% EXCEL EXPORT
export_filename = 'results/experiment_data.xlsx';

% Embedding Validation
appendCSVtoExcel('results/embedding_validation_*.csv', ...
                    export_filename, 'embedding_validation', 'timestamp');

% Similarity Search
appendCSVtoExcel('results/similarity_search_*.csv', ...
                    export_filename, 'similarity_search', 'Timestamp');

%% === LOAD & FILTER ===
data_table = readtable('experiment_data.xlsx', 'Sheet', 'embedding_validation', 'VariableNamingRule', 'preserve');
fprintf('Loaded: %d rows\n', height(data_table));

filtered = data_table( data_table.database_size == cfg.db_size & ...
                      data_table.top_k == cfg.top_k & ...
                      data_table.dtw_normalize == cfg.dtw_norm, :);
fprintf('Filtered: %d rows\n\n', height(filtered));

%% FIGURE 1: DIMENSIONALITY (BASELINES)
dtw_modes = {'joint_states', 'position'};

filtered_baselines = [];
for m = 1:length(dtw_modes)
    mask = strcmp(data_table.dtw_mode, dtw_modes{m}) & ...
           data_table.top_k == cfg.top_k;
    filtered_baselines = [filtered_baselines; data_table(mask, :)];
end
fprintf('✓ Baselines: %d experiments\n\n', height(filtered_baselines));

exclude_configs = {'Multi-15', 'Multi-25'};
% Falls deine Spalte anders heißt (z.B. 'emb_name'), hier anpassen:
keep_mask = ~ismember(filtered_baselines.emb_config, exclude_configs);
filtered_baselines = filtered_baselines(keep_mask, :);

fprintf('Creating Figure 1...\n');

fig1 = figure('Position', [200, 100, 750, 500], 'Color', 'w');
hold on;

% Farben
c_motion = [0.86, 0.13, 0.15];
c_shape  = [0.15, 0.39, 0.91];

% Trajectory - Motion
d = filtered_baselines(strcmpi(filtered_baselines.level, 'bahn') & strcmp(filtered_baselines.dtw_mode, 'joint_states'), :);
[dims, comp] = calcDimComposite(d, 6, cfg.w);
h1 = plot(dims, comp, 'o-', 'Color', c_motion, 'MarkerFaceColor', c_motion, 'LineWidth', 3, 'MarkerSize', 11);

% Trajectory - Shape
d = filtered_baselines(strcmpi(filtered_baselines.level, 'bahn') & strcmp(filtered_baselines.dtw_mode, 'position'), :);
[dims, comp] = calcDimComposite(d, 3, cfg.w);
h2 = plot(dims, comp, 's-', 'Color', c_shape, 'MarkerFaceColor', c_shape, 'LineWidth', 3, 'MarkerSize', 11);

% Segment - Motion
d = filtered_baselines(strcmp(filtered_baselines.dtw_mode, 'joint_states'), :);
[dims, comp] = calcDimComposite(d, 6, cfg.w);
h3 = plot(dims, comp, 'o--', 'Color', c_motion*0.6, 'MarkerFaceColor', c_motion*0.6, 'LineWidth', 2, 'MarkerSize', 8);

% Segment - Shape
d = filtered_baselines(strcmpi(filtered_baselines.level, 'segment') & strcmp(filtered_baselines.dtw_mode, 'position'), :);
[dims, comp] = calcDimComposite(d, 3, cfg.w);
h4 = plot(dims, comp, 's--', 'Color', c_shape*0.6, 'MarkerFaceColor', c_shape*0.6, 'LineWidth', 2, 'MarkerSize', 8);

%ylim([0.299,0.801])
ylim([0.8599, 0.9601])

% Styling
set(gca, 'XScale', 'log', 'XTick', [6,15,30,60,150,300,600,1200,3600], ...
    'FontName', 'Courier New', 'FontWeight', 'bold', 'FontSize', 18);
set(gca, 'YTick', [0.86,0.88,0.90,0.92,0.94,0.96],'FontName', 'Courier New', 'FontWeight', 'bold', 'FontSize', 18);
%set(gca, 'YTick', [0.3,0.4,0.5,0.6,0.7,0.8],'FontName', 'Courier New', 'FontWeight', 'bold', 'FontSize', 18);
xlabel('Total embedding dimensions', 'FontWeight', 'bold', 'FontSize', 20);
ylabel('NDCG@50_{GT}', 'FontWeight', 'bold', 'FontSize', 20);
%ylabel('\rho_{500}', 'FontWeight', 'bold', 'FontSize', 20);
legend([h1 h2 h3 h4], {'Motion (T)', 'Shape (T)', 'Motion (S)', 'Shape (S)'}, 'Location', [0.741333335558573,0.166500002622604,0.217333328882853,0.186999994754791],'FontSize', 14);
%legend([h1 h2 h3 h4], {'Motion (T)', 'Shape (T)', 'Motion (S)', 'Shape (S)'}, 'Location', 'best','FontSize', 14);
ax = gca; ax.Position = [0.15 0.14 0.84 0.84]; ax.YGrid = 'on';
hold off;

%exportgraphics(gca, fullfile(cfg.figure_folder, 'dimensionality_ndcg_baseline.pdf'));

%% FIGURE 1: DIMENSIONALITY (ALL)
dtw_modes = {'joint_states', 'position'};
baseline_weights = {'Joint + All', 'Pos + All'};

filtered_baselines = [];
for m = 1:length(dtw_modes)
    mask = strcmp(data_table.dtw_mode, dtw_modes{m}) & ...
           strcmp(data_table.weight_mode, baseline_weights{m}) & ...
           data_table.top_k == cfg.top_k;
    filtered_baselines = [filtered_baselines; data_table(mask, :)];
end
fprintf('✓ Baselines: %d experiments\n\n', height(filtered_baselines));

exclude_configs = {'Multi-15', 'Multi-25'};
% Falls deine Spalte anders heißt (z.B. 'emb_name'), hier anpassen:
keep_mask = ~ismember(filtered_baselines.emb_config, exclude_configs);
filtered_baselines = filtered_baselines(keep_mask, :);

fprintf('Creating Figure 1...\n');

fig1 = figure('Position', [200, 100, 750, 500], 'Color', 'w');
hold on;

% Farben
c_motion = [0.86, 0.13, 0.15];
c_shape  = [0.15, 0.39, 0.91];

% Trajectory - Motion
d = filtered_baselines(strcmpi(filtered_baselines.level, 'bahn') & strcmp(filtered_baselines.dtw_mode, 'joint_states'), :);
[dims, comp] = calcDimComposite(d, 6, cfg.w);
h1 = plot(dims, comp, 'o-', 'Color', c_motion, 'MarkerFaceColor', c_motion, 'LineWidth', 3, 'MarkerSize', 11);

% Trajectory - Shape
d = filtered_baselines(strcmpi(filtered_baselines.level, 'bahn') & strcmp(filtered_baselines.dtw_mode, 'position'), :);
[dims, comp] = calcDimComposite(d, 3, cfg.w);
h2 = plot(dims, comp, 's-', 'Color', c_shape, 'MarkerFaceColor', c_shape, 'LineWidth', 3, 'MarkerSize', 11);

% Segment - Motion
d = filtered_baselines(strcmp(filtered_baselines.dtw_mode, 'joint_states'), :);
[dims, comp] = calcDimComposite(d, 6, cfg.w);
h3 = plot(dims, comp, 'o--', 'Color', c_motion*0.6, 'MarkerFaceColor', c_motion*0.6, 'LineWidth', 2, 'MarkerSize', 8);

% Segment - Shape
d = filtered_baselines(strcmpi(filtered_baselines.level, 'segment') & strcmp(filtered_baselines.dtw_mode, 'position'), :);
[dims, comp] = calcDimComposite(d, 3, cfg.w);
h4 = plot(dims, comp, 's--', 'Color', c_shape*0.6, 'MarkerFaceColor', c_shape*0.6, 'LineWidth', 2, 'MarkerSize', 8);


ylim([0.299,0.801])
%ylim([0.8599, 0.9601])

% Styling
set(gca, 'XScale', 'log', 'XTick', [6,15,30,60,150,300,600,1200,3600], ...
    'FontName', 'Courier New', 'FontWeight', 'bold', 'FontSize', 18);
%set(gca, 'YTick', [0.86,0.88,0.90,0.92,0.94,0.96], 'FontName', 'Courier New', 'FontWeight', 'bold', 'FontSize', 18);
set(gca, 'YTick', [0.3,0.4,0.5,0.6,0.7,0.8],'FontName', 'Courier New', 'FontWeight', 'bold', 'FontSize', 18);
xlabel('Total embedding dimensions', 'FontWeight', 'bold', 'FontSize', 20);
%ylabel('NDCG@50_{GT}', 'FontWeight', 'bold', 'FontSize', 20);
ylabel('\rho_{500}', 'FontWeight', 'bold', 'FontSize', 20);
%legend([h1 h2 h3 h4], {'Motion (T)', 'Shape (T)', 'Motion (S)', 'Shape (S)'}, 'Location', [0.755333335558573,0.256500002622604,0.217333328882853,0.186999994754791], 'FontSize', 14);
legend([h1 h2 h3 h4], {'Motion (T)', 'Shape (T)', 'Motion (S)', 'Shape (S)'}, 'Location', 'best', 'FontSize', 14);
ax = gca; ax.Position = [0.15 0.14 0.84 0.84]; ax.YGrid = 'on';
hold off;

%exportgraphics(gca, fullfile(cfg.figure_folder, 'dimensionality_spearman_all.pdf'));


%% === TABLE 1 ===
table1_data = {};
levels = {'bahn', 'segment'};
level_labels = {'Trajectory', 'Segment'};
dtw_modes = {'joint_states', 'position'};
mode_labels = {'MOTION (Joint)', 'SHAPE (Cartesian)'};
baseline_weights = {'Joint only', 'Position only'};
dim_multiplier = [6, 3];

for lvl = 1:length(levels)
    level_data = filtered_baselines(strcmpi(filtered_baselines.level, levels{lvl}), :);
    
    for m = 1:length(dtw_modes)
        mode_data = level_data(strcmp(level_data.dtw_mode, dtw_modes{m}), :);
        mode_data.total_dims = mode_data.n_coarse + mode_data.n_fine;
        dims_unique = sort(unique(mode_data.total_dims));
        
        for i = 1:length(dims_unique)
            d = mode_data(mode_data.total_dims == dims_unique(i), :);
            metrics = calcMetrics(d, cfg.w);
            table1_data(end+1, :) = {level_labels{lvl}, mode_labels{m}, baseline_weights{m}, ...
                dims_unique(i), dims_unique(i) * dim_multiplier(m), height(d), ...
                metrics.spearman, metrics.ndcg50_dtw, metrics.ndcg50_gt,  ...
                metrics.ndcg10_dtw, metrics.ndcg10_gt, ...
                metrics.r50_dtw , metrics.r50_gt, ...
                metrics.r10_dtw, metrics.r10_gt, ...
                metrics.composite};
        end
    end
end

T1 = cell2table(table1_data, 'VariableNames', ...
    {'Level', 'Mode', 'Config', 'Dims_Per_Param', 'Actual_Dims', 'N', ...
    'Spearman', 'NDCG50_DTW', 'NDCG50_GT', ...
    'NDCG10_DTW', 'NDCG10_GT', ...
    'R50_DTW', 'R50_GT', ...
    'R10_DTW', 'R10_GT', ...
    'Composite'});
T1 = sortrows(T1, {'Level', 'Mode', 'Composite'}, {'ascend', 'ascend', 'descend'});
writetable(T1, fullfile('results/', 'table1_dimensionality.csv'));
fprintf('✓ Saved: table1_dimensionality.csv\n\n');


%% FIGURE 2: WEIGHT MODE CONTRIBUTION
fprintf('Creating Figure 2...\n');

% Filter
filtered_wm = data_table(data_table.database_size == 5000 & ...
                       data_table.top_k == 50 & ...
                       strcmpi(data_table.emb_config, 'Single-10') & ...
                       data_table.dtw_normalize == 0, :);

% Farben
c_motion = [0.86, 0.13, 0.15];
c_shape  = [0.15, 0.39, 0.91];

fig2 = figure('Position', [100, 100, 750, 500], 'Color', 'w');
hold on;

% Feature order für X-Achse
feature_order_joint = {'Joint only', 'Joint + Meta', 'Joint + Velocity', 'Joint + Orient', 'Joint + Position', 'Joint + All'};
feature_order_pos = {'Position only', 'Pos + Meta', 'Pos + Velocity', 'Pos + Orient', 'Pos + Joint', 'Pos + All'};

% Trajectory - Motion
d = filtered_wm(strcmpi(filtered_wm.level, 'bahn') & strcmp(filtered_wm.dtw_mode, 'joint_states'), :);
comp = calcWeightModeComposite(d, feature_order_joint, cfg.w);
h1 = plot(1:length(comp), comp, 'o-', 'Color', c_motion, 'MarkerFaceColor', c_motion, 'LineWidth', 3, 'MarkerSize', 11);

% Trajectory - Shape
d = filtered_wm(strcmpi(filtered_wm.level, 'bahn') & strcmp(filtered_wm.dtw_mode, 'position'), :);
comp = calcWeightModeComposite(d, feature_order_pos, cfg.w);
h2 = plot(1:length(comp), comp, 's-', 'Color', c_shape, 'MarkerFaceColor', c_shape, 'LineWidth', 3, 'MarkerSize', 11);

% Segment - Motion
d = filtered_wm(strcmpi(filtered_wm.level, 'segment') & strcmp(filtered_wm.dtw_mode, 'joint_states'), :);
comp = calcWeightModeComposite(d, feature_order_joint, cfg.w);
h3 = plot(1:length(comp), comp, 'o--', 'Color', c_motion*0.6, 'MarkerFaceColor', c_motion*0.6, 'LineWidth', 2, 'MarkerSize', 8);

% Segment - Shape
d = filtered_wm(strcmpi(filtered_wm.level, 'segment') & strcmp(filtered_wm.dtw_mode, 'position'), :);
comp = calcWeightModeComposite(d, feature_order_pos, cfg.w);
h4 = plot(1:length(comp), comp, 's--', 'Color', c_shape*0.6, 'MarkerFaceColor', c_shape*0.6, 'LineWidth', 2, 'MarkerSize', 8);

% Styling
set(gca, 'XTick', 1:6, 'XTickLabel', {'Baseline', '+Meta', '+Vel', '+Orient', '+Cross', '+All'});
xlabel('Incremental feature addition', 'FontWeight', 'bold', 'FontSize', 20);
ylabel('NDCG@50 (GT)', 'FontWeight', 'bold', 'FontSize', 20);
set(gca, 'FontName', 'Courier New', 'FontWeight', 'bold', 'FontSize', 18);
legend([h1 h2 h3 h4], {'Motion (Traj.)', 'Shape (Traj.)', 'Motion (Seg.)', 'Shape (Seg.)'}, 'Location', 'best', 'FontSize', 14);
ax = gca; ax.Position = [0.14 0.14 0.82 0.84]; ax.YGrid = 'on';
hold off;

if cfg.save_figures
    exportgraphics(gca, fullfile(figure_folder, 'incremental_feature.pdf'));
end

%% TABLE 2: WEIGHT MODES
table2_data = {};
levels = {'bahn', 'segment'};
level_labels = {'Trajectory', 'Segment'};
dtw_modes = {'joint_states', 'position'};
mode_labels = {'Motion', 'Shape'};

for lvl = 1:length(levels)
    level_data = filtered_wm(strcmpi(filtered_wm.level, levels{lvl}), :);
    
    for m = 1:length(dtw_modes)
        mode_data = level_data(strcmp(level_data.dtw_mode, dtw_modes{m}), :);
        weight_modes = unique(mode_data.weight_mode);
        
        for w = 1:length(weight_modes)
            wm = weight_modes{w};
            d = mode_data(strcmp(mode_data.weight_mode, wm), :);
            if height(d) < 3, continue; end
            
            metrics = calcMetrics(d, cfg.w);
            table2_data(end+1, :) = {level_labels{lvl}, mode_labels{m}, wm, height(d), ...
                metrics.spearman, metrics.ndcg50_dtw, metrics.ndcg50_gt};
        end
    end
end

T2 = cell2table(table2_data, 'VariableNames', ...
    {'Level', 'Mode', 'Weight_Mode', 'N', 'Spearman', 'NDCG50_DTW', 'NDCG50_GT'});
T2 = sortrows(T2, 'NDCG50_GT', 'descend');
writetable(T2, fullfile('results/', 'table2_weight_modes.csv'));
fprintf('✓ Saved: table2_weight_modes.csv\n\n');

% ========================================================================
%% COMPOSITE SCORE CONFIGURATION
% ========================================================================

function [dims, composite] = calcDimComposite(data, multiplier, w)
    data.total_dims = data.n_coarse + data.n_fine;
    dims_unique = sort(unique(data.total_dims));
    dims = dims_unique * multiplier;
    composite = zeros(length(dims_unique), 1);
    
    for i = 1:length(dims_unique)
        d = data(data.total_dims == dims_unique(i), :);
        m = calcMetrics(d, w);
        composite(i) = m.composite;
    end
end

function metrics = calcMetrics(data, w)
    metrics.spearman = mean(data.spearman_dtw_eb);
    metrics.ndcg50_dtw = mean(data.ndcg_50_dtw_eb, 'omitnan');
    metrics.ndcg50_gt = mean(data.ndcg_50_gt_eb, 'omitnan');
    metrics.r50_dtw = mean(data.r50_dtw_eb);
    metrics.r50_gt = mean(data.r50_gt_eb);
    metrics.r10_dtw = mean(data.r10_dtw_eb);
    metrics.r10_gt = mean(data.r10_gt_eb);
    metrics.ndcg10_dtw = mean(data.ndcg_10_dtw_eb);
    metrics.ndcg10_gt = mean(data.ndcg_10_gt_eb);
    metrics.composite = w.spearman * metrics.spearman + ...
                        w.ndcg50_dtw * metrics.ndcg50_dtw + ...
                        w.ndcg50_gt * metrics.ndcg50_gt + ...
                        w.r50_dtw * metrics.r50_dtw + ...
                        w.r50_gt * metrics.r50_gt;
end

function composite = calcWeightModeComposite(data, feature_order, w)
    composite = zeros(length(feature_order), 1);
    for i = 1:length(feature_order)
        d = data(strcmp(data.weight_mode, feature_order{i}), :);
        if height(d) > 0
            m = calcMetrics(d, w);
            composite(i) = m.composite;
        else
            composite(i) = NaN;
        end
    end
end
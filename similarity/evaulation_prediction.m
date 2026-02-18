%% EXCEL EXPORT
export_filename = 'results/experiment_data.xlsx';

% Embedding Validation
appendCSVtoExcel('results/embedding_validation_*.csv', ...
                    export_filename, 'embedding_validation', 'timestamp');

% Similarity Search
appendCSVtoExcel('results/similarity_search_*.csv', ...
                    export_filename, 'similarity_search', 'Timestamp');

%% SCATTER PLOT SECTION 5.1 OVERALL PERFORMANCE (BASELINE)
% === 1. LOAD & FILTER ===
data2_table = readtable('experiment_data.xlsx', 'Sheet', 'similarity_search', 'VariableNamingRule', 'preserve');
fprintf('Loaded total: %d rows\n', height(data2_table));

% --- FILTER: Level 'bahn', entsprechende Weights, K=50 ---
filter_idx_motion = strcmpi(data2_table.level, 'bahn') & ...
             strcmp(data2_table.weights, '0,1,0,0,0') & ...
             data2_table.K == 50;
filter_idx_shape = strcmpi(data2_table.level, 'bahn') & ...
             strcmp(data2_table.weights, '1,0,0,0,0') & ...
             data2_table.K == 50; 

data2_table_motion = data2_table(filter_idx_motion, :);
data2_table_shape = data2_table(filter_idx_shape, :);

% --- TIMESTAMPS ---
ts_motion = "20260201_230706"; 
ts_shape  = "20260201_231314"; 

% DATEN EXTRAHIEREN
idx_m = string(data2_table_motion.Timestamp) == ts_motion;
gt_motion = data2_table_motion.ground_truth(idx_m);
pred1_motion = data2_table_motion.fromsegs_s1_weighted(idx_m);
pred2_motion = data2_table_motion.fromsegs_s2_weighted(idx_m);

idx_s = string(data2_table_shape.Timestamp) == ts_shape;
gt_shape = data2_table_shape.ground_truth(idx_s);
pred1_shape = data2_table_shape.fromsegs_s1_weighted(idx_s);
pred2_shape = data2_table_shape.fromsegs_s2_weighted(idx_s);

fprintf('Final Selection -> Motion: %d samples, Shape: %d samples\n', length(gt_motion), length(gt_shape));

% Style-Definitionen
c_motion = [0.86, 0.13, 0.15]; % Rot
c_shape  = [0.15, 0.39, 0.91]; % Blau
c_ideal  = 'black';    % Grau

%%%%%%%%%%%%%%%%%%%%%%%%%% === FIGURE A: STAGE 1 (Embedding) ===
fprintf('Creating Plot Stage 2...\n');

% 1. Statistik berechnen (Stage 2)
% Motion
R_m1 = corrcoef(gt_motion, pred1_motion); r_m1 = R_m1(1,2);
rmse_m1 = sqrt(mean((pred1_motion - gt_motion).^2));
str_leg_m1 = sprintf('Motion');

% Shape
R_s1 = corrcoef(gt_shape, pred1_shape); r_s1 = R_s1(1,2);
rmse_s1 = sqrt(mean((pred1_shape - gt_shape).^2));
str_leg_s1 = sprintf('Shape');

% 2. Plotten
fig1 = figure('Position', [100, 100, 750, 500], 'Color', 'w'); % Fenster verschoben damit man beide sieht
hold on;

% Limits basierend auf Stage 2 (oder global max, wie du willst)
max_val1 = max([gt_motion; gt_shape; pred1_motion; pred1_shape]);
limit_max1 = ceil(max_val1); 

h_ideal = plot([0 limit_max1], [0 limit_max1], '--', 'Color', c_ideal, 'LineWidth', 2);
h_motion = scatter(gt_motion, pred1_motion, 90, c_motion, 'filled', 'o', 'MarkerFaceAlpha', 0.8);
h_shape  = scatter(gt_shape, pred1_shape, 90, c_shape, 'filled', 's', 'MarkerFaceAlpha', 0.8);

grid on; box on; axis square;
xlim([0 limit_max1]); ylim([0 limit_max1]);

xlabel('Measured TCP error [mm]', 'FontWeight', 'bold', 'FontSize', 16);
ylabel('Predicted TCP error [mm]', 'FontWeight', 'bold', 'FontSize', 16);
yticks([0,0.2,0.4,0.6,0.8, 1.0])
yticklabels([0,0.2,0.4,0.6,0.8, 1.0])
set(gca, 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', 14);

% Legende mit Stats
[~, objh] = legend([h_ideal, h_motion, h_shape], ...
    {'Ideal', str_leg_m1, str_leg_s1}, ...
    'Location', 'northwest', 'FontSize', 14);

objhl = findobj(objh, 'type', 'patch'); %// objects of legend of type line
set(objhl, 'Markersize', 10); %// set marker size as desired

str_stats = { ...
    sprintf('Motion: r = %.2f, RMSE = %.3f', r_m1, rmse_m1), ...
    sprintf('Shape:  r = %.2f, RMSE = %.3f', r_s1, rmse_s1)};

text(0.96, 0.04, str_stats, ...
    'Units', 'normalized', ...           % Position relativ zur Achse (0 bis 1)
    'HorizontalAlignment', 'right', ...  % Rechtsbündig an 0.96
    'VerticalAlignment', 'bottom', ...   % Unten andocken an 0.04
    'BackgroundColor', 'w', ...          % Weißer Hintergrund
    'EdgeColor', 'k', ...                % Schwarzer Rahmen
    'LineWidth', 0.5, ...                  % Rahmen-Dicke
    'Margin', 4, ...
    'FontWeight', 'bold', ...% Innenabstand
    'FontSize', 15, ...
    'FontName', 'Times New Roman');


hold off;
exportgraphics(gca, fullfile('figs', 'prediction_scatter_s1_baseline.pdf'));


% === FIGURE B: STAGE 2 (DTW Refinement) ===
fprintf('Creating Plot Stage 2...\n');

% 1. Statistik berechnen (Stage 2)
% Motion
R_m2 = corrcoef(gt_motion, pred2_motion); r_m2 = R_m2(1,2);
rmse_m2 = sqrt(mean((pred2_motion - gt_motion).^2));
str_leg_m2 = sprintf('Motion');

% Shape
R_s2 = corrcoef(gt_shape, pred2_shape); r_s2 = R_s2(1,2);
rmse_s2 = sqrt(mean((pred2_shape - gt_shape).^2));
str_leg_s2 = sprintf('Shape');

% 2. Plotten
fig2 = figure('Position', [100, 100, 750, 500], 'Color', 'w'); % Fenster verschoben damit man beide sieht
hold on;

% Limits basierend auf Stage 2 (oder global max, wie du willst)
max_val2 = max([gt_motion; gt_shape; pred2_motion; pred2_shape]);
limit_max2 = ceil(max_val2); 

h_ideal = plot([0 limit_max2], [0 limit_max2], '--', 'Color', c_ideal, 'LineWidth', 2);
h_motion = scatter(gt_motion, pred2_motion, 90, c_motion, 'filled', 'o', 'MarkerFaceAlpha', 0.8);
h_shape  = scatter(gt_shape, pred2_shape, 90, c_shape, 'filled', 's', 'MarkerFaceAlpha', 0.8);

grid on; box on; axis square;
xlim([0 limit_max2]); ylim([0 limit_max2]);

xlabel('Measured error [mm]', 'FontWeight', 'bold', 'FontSize', 16);
ylabel('Predicted error [mm]', 'FontWeight', 'bold', 'FontSize', 16);
yticks([0,0.2,0.4,0.6,0.8, 1.0]);
yticklabels([0,0.2,0.4,0.6,0.8, 1.0]);
set(gca, 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', 14);

% Legende mit Stats
[~, objh] = legend([h_ideal, h_motion, h_shape], ...
    {'Ideal', str_leg_m1, str_leg_s1}, ...
    'Location', 'northwest', 'FontSize', 14);

objhl = findobj(objh, 'type', 'patch'); %// objects of legend of type line
set(objhl, 'Markersize', 10); %// set marker size as desired

str_stats = { ...
    sprintf('Motion: r = %.2f, RMSE = %.3f', r_m2, rmse_m2), ...
    sprintf('Shape:  r = %.2f, RMSE = %.3f', r_s2, rmse_s2)};

text(0.96, 0.04, str_stats, ...
    'Units', 'normalized', ...           % Position relativ zur Achse (0 bis 1)
    'HorizontalAlignment', 'right', ...  % Rechtsbündig an 0.96
    'VerticalAlignment', 'bottom', ...   % Unten andocken an 0.04
    'BackgroundColor', 'w', ...          % Weißer Hintergrund
    'EdgeColor', 'k', ...                % Schwarzer Rahmen
    'LineWidth', 0.5, ...                  % Rahmen-Dicke
    'Margin', 4, ...
    'FontWeight', 'bold', ...% Innenabstand
    'FontSize', 15, ...
    'FontName', 'Times New Roman');

hold off;

exportgraphics(gca, fullfile('figs', 'prediction_scatter_s2_baseline.pdf'));

%% SCATTER PLOT SECTION 5.1 OVERALL PERFORMANCE (ALL)
% === 1. LOAD & FILTER ===
data2_table = readtable('experiment_data.xlsx', 'Sheet', 'similarity_search', 'VariableNamingRule', 'preserve');
fprintf('Loaded total: %d rows\n', height(data2_table));

% --- FILTER: Level 'bahn', Weights '1,1,1,1,1', K=50 ---
filter_idx = strcmpi(data2_table.level, 'bahn') & ...
             strcmp(data2_table.weights, '1,1,1,1,1') & ...
             data2_table.K == 50; 
data2_table = data2_table(filter_idx, :);
fprintf('Filtered (Bahn & Weights [1,1,1,1,1], K=50): %d rows\n', height(data2_table));
%
% --- TIMESTAMPS ---
ts_motion = "20260201_230706"; 
ts_shape  = "20260201_231314"; 

% DATEN EXTRAHIEREN
idx_m = string(data2_table.Timestamp) == ts_motion;
gt_motion = data2_table.ground_truth(idx_m);
pred1_motion = data2_table.fromsegs_s1_weighted(idx_m);
pred2_motion = data2_table.fromsegs_s2_weighted(idx_m);

idx_s = string(data2_table.Timestamp) == ts_shape;
gt_shape = data2_table.ground_truth(idx_s);
pred1_shape = data2_table.fromsegs_s1_weighted(idx_s);
pred2_shape = data2_table.fromsegs_s2_weighted(idx_s);

fprintf('Final Selection -> Motion: %d samples, Shape: %d samples\n', length(gt_motion), length(gt_shape));

% Style-Definitionen
c_motion = [0.86, 0.13, 0.15]; % Rot
c_shape  = [0.15, 0.39, 0.91]; % Blau
c_ideal  = 'black';    % Grau

%%%%%%%%%%%%%%%%%%%%%%%%%% === FIGURE A: STAGE 1 (Embedding) ===
fprintf('Creating Plot Stage 2...\n');

% 1. Statistik berechnen (Stage 2)
% Motion
R_m1 = corrcoef(gt_motion, pred1_motion); r_m1 = R_m1(1,2);
rmse_m1 = sqrt(mean((pred1_motion - gt_motion).^2));
str_leg_m1 = sprintf('Motion');

% Shape
R_s1 = corrcoef(gt_shape, pred1_shape); r_s1 = R_s1(1,2);
rmse_s1 = sqrt(mean((pred1_shape - gt_shape).^2));
str_leg_s1 = sprintf('Shape');

% 2. Plotten
fig1 = figure('Position', [100, 100, 750, 500], 'Color', 'w'); % Fenster verschoben damit man beide sieht
hold on;

% Limits basierend auf Stage 2 (oder global max, wie du willst)
max_val1 = max([gt_motion; gt_shape; pred1_motion; pred1_shape]);
limit_max1 = ceil(max_val1); 

h_ideal = plot([0 limit_max1], [0 limit_max1], '--', 'Color', c_ideal, 'LineWidth', 2);
h_motion = scatter(gt_motion, pred1_motion, 90, c_motion, 'filled', 'o', 'MarkerFaceAlpha', 0.8);
h_shape  = scatter(gt_shape, pred1_shape, 90, c_shape, 'filled', 's', 'MarkerFaceAlpha', 0.8);

grid on; box on; axis square;
xlim([0 limit_max1]); ylim([0 limit_max1]);

xlabel('Measured error [mm]', 'FontWeight', 'bold', 'FontSize', 16);
ylabel('Predicted error [mm]', 'FontWeight', 'bold', 'FontSize', 16);
yticks([0,0.2,0.4,0.6,0.8, 1.0]);
yticklabels([0,0.2,0.4,0.6,0.8, 1.0]);
set(gca, 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', 14);

% Legende mit Stats
[~, objh] = legend([h_ideal, h_motion, h_shape], ...
    {'Ideal', str_leg_m1, str_leg_s1}, ...
    'Location', 'northwest', 'FontSize', 14);

objhl = findobj(objh, 'type', 'patch'); %// objects of legend of type line
set(objhl, 'Markersize', 10); %// set marker size as desired

str_stats = { ...
    sprintf('Motion: r = %.2f, RMSE = %.3f', r_m1, rmse_m1), ...
    sprintf('Shape:  r = %.2f, RMSE = %.3f', r_s1, rmse_s1)};

text(0.96, 0.04, str_stats, ...
    'Units', 'normalized', ...           % Position relativ zur Achse (0 bis 1)
    'HorizontalAlignment', 'right', ...  % Rechtsbündig an 0.96
    'VerticalAlignment', 'bottom', ...   % Unten andocken an 0.04
    'BackgroundColor', 'w', ...          % Weißer Hintergrund
    'EdgeColor', 'k', ...                % Schwarzer Rahmen
    'LineWidth', 0.5, ...                  % Rahmen-Dicke
    'Margin', 4, ...
    'FontWeight', 'bold', ...% Innenabstand
    'FontSize', 15, ...
    'FontName', 'Times New Roman');


hold off;
exportgraphics(gca, fullfile('figs', 'prediction_scatter_s1_all.pdf'));


% === FIGURE B: STAGE 2 (DTW Refinement) ===
fprintf('Creating Plot Stage 2...\n');

% 1. Statistik berechnen (Stage 2)
% Motion
R_m2 = corrcoef(gt_motion, pred2_motion); r_m2 = R_m2(1,2);
rmse_m2 = sqrt(mean((pred2_motion - gt_motion).^2));
str_leg_m2 = sprintf('Motion');

% Shape
R_s2 = corrcoef(gt_shape, pred2_shape); r_s2 = R_s2(1,2);
rmse_s2 = sqrt(mean((pred2_shape - gt_shape).^2));
str_leg_s2 = sprintf('Shape');

% 2. Plotten
fig2 = figure('Position', [100, 100, 750, 500], 'Color', 'w'); % Fenster verschoben damit man beide sieht
hold on;

% Limits basierend auf Stage 2 (oder global max, wie du willst)
max_val2 = max([gt_motion; gt_shape; pred2_motion; pred2_shape]);
limit_max2 = ceil(max_val2); 

h_ideal = plot([0 limit_max2], [0 limit_max2], '--', 'Color', c_ideal, 'LineWidth', 2);
h_motion = scatter(gt_motion, pred2_motion, 90, c_motion, 'filled', 'o', 'MarkerFaceAlpha', 0.8);
h_shape  = scatter(gt_shape, pred2_shape, 90, c_shape, 'filled', 's', 'MarkerFaceAlpha', 0.8);

grid on; box on; axis square;
xlim([0 limit_max2]); ylim([0 limit_max2]);

xlabel('Measured error [mm]', 'FontWeight', 'bold', 'FontSize', 16);
ylabel('Predicted error [mm]', 'FontWeight', 'bold', 'FontSize', 16);
yticks([0,0.2,0.4,0.6,0.8, 1.0]);
yticklabels([0,0.2,0.4,0.6,0.8, 1.0]);
set(gca, 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', 14);

% Legende mit Stats
[~, objh] = legend([h_ideal, h_motion, h_shape], ...
    {'Ideal', str_leg_m1, str_leg_s1}, ...
    'Location', 'northwest', 'FontSize', 14);

objhl = findobj(objh, 'type', 'patch'); %// objects of legend of type line
set(objhl, 'Markersize', 10); %// set marker size as desired

str_stats = { ...
    sprintf('Motion: r = %.2f, RMSE = %.3f', r_m2, rmse_m2), ...
    sprintf('Shape:  r = %.2f, RMSE = %.3f', r_s2, rmse_s2)};

text(0.96, 0.04, str_stats, ...
    'Units', 'normalized', ...           % Position relativ zur Achse (0 bis 1)
    'HorizontalAlignment', 'right', ...  % Rechtsbündig an 0.96
    'VerticalAlignment', 'bottom', ...   % Unten andocken an 0.04
    'BackgroundColor', 'w', ...          % Weißer Hintergrund
    'EdgeColor', 'k', ...                % Schwarzer Rahmen
    'LineWidth', 0.5, ...                  % Rahmen-Dicke
    'Margin', 4, ...
    'FontWeight', 'bold', ...% Innenabstand
    'FontSize', 15, ...
    'FontName', 'Times New Roman');

hold off;

exportgraphics(gca, fullfile('figs', 'prediction_scatter_s2_all.pdf'));

%% FIGURE C: RESIDUAL HISTOGRAM // OVERALL PERFORMANCE
fprintf('Creating Residual Histogram...\n');

resid_motion = pred2_motion - gt_motion;
resid_shape  = pred2_shape - gt_shape;

c_motion = [0.86, 0.13, 0.15]; % Rot
c_shape  = [0.15, 0.39, 0.91]; % Blau
c_zero   = [0.4, 0.4, 0.4];    % Grau für die Nulllinie

fig4 = figure('Position', [150, 150, 700, 500], 'Color', 'w');
hold on;

% 1. Nulllinie zeichnen (Referenz)
xline(0, '--', 'Color', c_zero, 'LineWidth', 2, 'DisplayName', 'Ideal Prediction');

% 2. Histogramme plotten
% Wir nutzen 'FaceAlpha' für Transparenz, damit man Überlappungen sieht
h_hist_m = histogram(resid_motion, 'BinWidth', 0.02, 'FaceColor', c_motion, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.6);

h_hist_s = histogram(resid_shape, 'BinWidth', 0.02, 'FaceColor', c_shape, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.6);

% 3. Styling
grid on;
box on;
xlabel('Prediction Residual ($MAE_{pred} - MAE_{true}$) [mm]', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 16);
ylabel('Count (Number of Trajectories)', 'FontWeight', 'bold', 'FontSize', 16);
title('Error Distribution & Safety Bias', 'FontWeight', 'bold', 'FontSize', 16);

% Achsen hübsch machen
set(gca, 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', 14);
xlim([-0.2 0.2]); % Fokus auf den relevanten Bereich (anpassbar!)

% 4. Legende
legend([h_hist_m, h_hist_s], ...
    {'Motion Resid.', 'Shape Resid.'}, ...
    'Location', 'northeast', 'FontSize', 12);

% 5. Text-Annotationen für "Safe" vs "Risky" (Optional, aber cool)
yl = ylim;
text(-0.15, yl(2)*0.9, 'Underestimation', 'Color', [0.8 0 0], 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
text(-0.15, yl(2)*0.85, '(Risky)', 'Color', [0.8 0 0], 'FontSize', 10, 'HorizontalAlignment', 'center');

text(0.15, yl(2)*0.9, 'Overestimation', 'Color', [0 0.6 0], 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
text(0.15, yl(2)*0.85, '(Safe)', 'Color', [0 0.6 0], 'FontSize', 10, 'HorizontalAlignment', 'center');

hold off;

% Speichern
% exportgraphics(gca, fullfile(figure_folder, 'prediction_histogram.pdf'));


%% SCATTER PLOT SECTION 5.2 INFLUENCE GRANULARITY
% === 1. LOAD & FILTER ===
data2_table = readtable('experiment_data.xlsx', 'Sheet', 'similarity_search', 'VariableNamingRule', 'preserve');
fprintf('Loaded total: %d rows\n', height(data2_table));

% Filter: Bahn-Level, Weights 1,1,1,1,1, K=50
filter_idx = strcmpi(data2_table.level, 'bahn') & ...
             strcmp(data2_table.weights, '1,1,1,1,1') & ...
             data2_table.K == 50 & ...
             strcmp(data2_table.dtw_mode, 'position');
data2_table = data2_table(filter_idx, :);

% Timestamps
ts_motion = "20260201_230706"; 
ts_shape  = "20260201_231314"; 

% DATEN EXTRAHIEREN (Wir kombinieren Motion & Shape für den Granularity-Vergleich)
idx_m = string(data2_table.Timestamp) == ts_motion;
idx_s = string(data2_table.Timestamp) == ts_shape;

% Ground Truth (Alle zusammen)
gt_all = [data2_table.ground_truth(idx_m); data2_table.ground_truth(idx_s)];

% STAGE 1 Predictions
pred1_direct = [data2_table.direct_s1_weighted(idx_m); data2_table.direct_s1_weighted(idx_s)];
pred1_seg    = [data2_table.fromsegs_s1_weighted(idx_m); data2_table.fromsegs_s1_weighted(idx_s)];

% STAGE 2 Predictions
pred2_direct = [data2_table.direct_s2_weighted(idx_m); data2_table.direct_s2_weighted(idx_s)];
pred2_seg    = [data2_table.fromsegs_s2_weighted(idx_m); data2_table.fromsegs_s2_weighted(idx_s)];

fprintf('Data Points: %d (Combined Motion+Shape)\n', length(gt_all));

% Style-Definitionen für Granularität
c_direct = [0.85, 0.33, 0.10]; % Orange/Rot (Trajectory)
c_seg    = [0.15, 0.39, 0.91]; % Royal Blau (Segment)
c_ideal  = [0.4, 0.4, 0.4];    % Grau

% === FIGURE A: STAGE 1 (Granularity Comparison) ===
fprintf('Creating Plot Stage 1 (Direct vs Segment)...\n');

% Statistik berechnen
[R_d1, ~] = corrcoef(gt_all, pred1_direct); r_d1 = R_d1(1,2);
rmse_d1 = sqrt(mean((pred1_direct - gt_all).^2));

[R_s1, ~] = corrcoef(gt_all, pred1_seg); r_s1 = R_s1(1,2);
rmse_s1 = sqrt(mean((pred1_seg - gt_all).^2));

% Plotten
fig1 = figure('Position', [100, 100, 600, 600], 'Color', 'w');
hold on;

max_val = max([gt_all; pred1_direct; pred1_seg]);
limit_max = ceil(max_val * 10) / 10 + 0.05; 

h_ideal = plot([0 limit_max], [0 limit_max], '--', 'Color', c_ideal, 'LineWidth', 2);

% Scatter Plots (Direct zuerst, damit Segment darüber liegt)
h_direct = scatter(gt_all, pred1_direct, 50, c_direct, 'filled', 'v', 'MarkerFaceAlpha', 0.5);
h_seg    = scatter(gt_all, pred1_seg, 50, c_seg, 'filled', 'o', 'MarkerFaceAlpha', 0.6);

grid on; box on; axis square;
xlim([0 limit_max]); ylim([0 limit_max]);

xlabel('Measured TCP error [mm]', 'FontWeight', 'bold', 'FontSize', 16);
ylabel('Predicted TCP error [mm]', 'FontWeight', 'bold', 'FontSize', 16);
set(gca, 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', 14);

% Legende
legend([h_ideal, h_direct, h_seg], ...
    {'Ideal', 'Trajectory (Direct)', 'Segment (Decomposed)'}, ...
    'Location', 'northwest', 'FontSize', 14);

% Statistik Box
str_stats1 = { ...
    sprintf('Direct:  r = %.2f, RMSE = %.3f', r_d1, rmse_d1), ...
    sprintf('Decomposed: r = %.2f, RMSE = %.3f', r_s1, rmse_s1)};

text(0.96, 0.04, str_stats1, 'Units', 'normalized', 'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'bottom', 'BackgroundColor', 'w', 'EdgeColor', 'k', ...
    'LineWidth', 1, 'Margin', 6, 'FontSize', 15, 'FontName', 'Times New Roman', 'FontWeight', 'bold');

hold off;
exportgraphics(gca, fullfile('figs', 'granularity_stage1.pdf'));


% === FIGURE B: STAGE 2 (Granularity Comparison) ===
fprintf('Creating Plot Stage 2 (Direct vs Segment)...\n');

% Statistik berechnen
[R_d2, ~] = corrcoef(gt_all, pred2_direct); r_d2 = R_d2(1,2);
rmse_d2 = sqrt(mean((pred2_direct - gt_all).^2));

[R_s2, ~] = corrcoef(gt_all, pred2_seg); r_s2 = R_s2(1,2);
rmse_s2 = sqrt(mean((pred2_seg - gt_all).^2));

% Plotten
fig2 = figure('Position', [750, 100, 600, 600], 'Color', 'w');
hold on;

max_val2 = max([gt_all; pred2_direct; pred2_seg]);
limit_max2 = ceil(max_val2 * 10) / 10 + 0.05; 

h_ideal = plot([0 limit_max2], [0 limit_max2], '--', 'Color', c_ideal, 'LineWidth', 2);

h_direct = scatter(gt_all, pred2_direct, 50, c_direct, 'filled', 'v', 'MarkerFaceAlpha', 0.5);
h_seg    = scatter(gt_all, pred2_seg, 50, c_seg, 'filled', 'o', 'MarkerFaceAlpha', 0.6);

grid on; box on; axis square;
xlim([0 limit_max2]); ylim([0 limit_max2]);

xlabel('Measured TCP error [mm]', 'FontWeight', 'bold', 'FontSize', 16);
ylabel('Predicted TCP error [mm]', 'FontWeight', 'bold', 'FontSize', 16);
set(gca, 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', 14);

% Legende
legend([h_ideal, h_direct, h_seg], ...
    {'Ideal', 'Trajectory Direct', 'Decomposed'}, ...
    'Location', 'northwest', 'FontSize', 14);

% Statistik Box
str_stats2 = { ...
    sprintf('Direct:  r = %.2f, RMSE = %.3f', r_d2, rmse_d2), ...
    sprintf('Decomposed: r = %.2f, RMSE = %.3f', r_s2, rmse_s2)};

text(0.96, 0.04, str_stats2, 'Units', 'normalized', 'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'bottom', 'BackgroundColor', 'w', 'EdgeColor', 'k', ...
    'LineWidth', 1, 'Margin', 6, 'FontSize', 15, 'FontName', 'Times New Roman', 'FontWeight', 'bold');

hold off;
exportgraphics(gca, fullfile('figs', 'granularity_stage2.pdf'));


%% === RMSE CALCULATION: 5.3 INFLUENCE OF WEIGHTS ===
% Load Data
data2_table = readtable('experiment_data.xlsx', 'Sheet', 'similarity_search', 'VariableNamingRule', 'preserve');
fprintf('Loaded total: %d rows\n', height(data2_table));

% 1. BASIS-FILTER (Alles außer Weights)
% Wir filtern auf: Level='bahn', K=50, Mode='joint_states'
% ACHTUNG: Wir filtern HIER NOCH NICHT auf bestimmte Weights!
base_filter = strcmpi(data2_table.level, 'bahn') & ...
              data2_table.K == 50 & ...
              strcmp(data2_table.dtw_mode, 'joint_states');

data_filtered = data2_table(base_filter, :);

% 2. TIMESTAMPS (Motion & Shape)
ts_motion = "20260201_230706"; 
ts_shape  = "20260201_231314"; 

% Nur Daten von diesen beiden Experimenten behalten
time_filter = (string(data_filtered.Timestamp) == ts_motion) | ...
              (string(data_filtered.Timestamp) == ts_shape);
data_final = data_filtered(time_filter, :);

% 3. IDENTIFIZIERE ALLE VORHANDENEN GEWICHTE
unique_weights = unique(data_final.weights);
fprintf('\nFound %d unique weight configurations.\n', length(unique_weights));
fprintf('-----------------------------------------------------------------\n');
fprintf('%-20s | %-12s | %-12s | %-12s\n', 'Weights', 'N (Samples)', 'RMSE (S1)', 'RMSE (S2)');
fprintf('-----------------------------------------------------------------\n');

% 4. LOOP ÜBER ALLE GEWICHTE & BERECHNUNG
for i = 1:length(unique_weights)
    w_str = unique_weights{i};
    
    % Daten für dieses spezifische Gewicht extrahieren
    idx_w = strcmp(data_final.weights, w_str);
    sub_data = data_final(idx_w, :);
    
    % Ground Truth & Predictions (Wir nutzen NUR 'fromsegs' = Decomposed)
    gt   = sub_data.ground_truth;
    pred_s1 = sub_data.fromsegs_s1_weighted;
    pred_s2_dec = sub_data.fromsegs_s2_weighted;
    
    % RMSE Berechnung
    rmse_s1 = sqrt(mean((pred_s1 - gt).^2));
    rmse_s2 = sqrt(mean((pred_s2_dec - gt).^2));
    
    % Ausgabe in Tabelle
    fprintf('%-20s | %5d        | %.5f      | %.5f\n', w_str, length(gt), rmse_s1, rmse_s2);
end
fprintf('-----------------------------------------------------------------\n');

%% === RMSE ANALYSIS: INFLUENCE OF NEIGHBORHOOD SIZE (K) ===
% Load Data
data2_table = readtable('experiment_data.xlsx', 'Sheet', 'similarity_search', 'VariableNamingRule', 'preserve');
fprintf('Loaded total: %d rows\n', height(data2_table));

target_timestamps = [ ...
    "20260201_231314", ...
    "20260202_232739", ...
    "20260203_143040", ...
    "20260201_230706", ...
    "20260202_181345", ...
    "20260203_143127", ...
    "20260203_172314", ...
    "20260203_172348" ...
];

% 1. FILTERUNG (Alles in einem Schritt)
% Level = 'bahn'
% Mode  = 'position' (SHAPE)
% Weights = '1,1,1,1,1'
% CONSTRAINT 1: K muss gleich dtw_calls sein (Volle Pipeline)
% CONSTRAINT 2: Timestamp muss in der Whitelist sein
%%
filter_idx = strcmpi(data2_table.level, 'bahn') & ...
             strcmp(data2_table.dtw_mode, 'position') & ... 
             strcmp(data2_table.weights, '1,1,1,1,1') & ...
             (data2_table.K == data2_table.dtw_calls) & ...
             ismember(string(data2_table.Timestamp), target_timestamps);

data_filtered = data2_table(filter_idx, :);
fprintf('Filtered rows (Valid Timestamps, Position, K == dtw_calls): %d\n', height(data_filtered));

% 2. GRUPPIEREN NACH K
unique_K = unique(data_filtered.K);

fprintf('\n=== RESULTS: Influence of Neighborhood Size (K) - SHAPE ===\n');
fprintf('-------------------------------------------------------------\n');
fprintf('%-10s | %-12s | %-12s | %-12s\n', 'Size (K)', 'Samples', 'MAE (S1)', 'MAE (S2)','RMSE (S1)', 'RMSE (S2)');
fprintf('-------------------------------------------------------------\n');

% 3. LOOP & BERECHNUNG
for i = 1:length(unique_K)
    k_val = unique_K(i);
    
    % Daten für dieses K
    sub_data = data_filtered(data_filtered.K == k_val, :);
    
    % Ground Truth & Predictions (Decomposed)
    gt   = sub_data.ground_truth;
    pred_s1 = sub_data.fromsegs_s1_weighted;
    pred_s2_dec = sub_data.fromsegs_s2_weighted;
    
    % RMSE Berechnung
    mae_s1 = mean(abs(pred_s1 - gt));
    mae_2 = mean(abs(pred_s2_dec - gt));
    rmse_s1 = sqrt(mean((pred_s1 - gt).^2));
    rmse_s2 = sqrt(mean((pred_s2_dec - gt).^2));
    
    % Ausgabe
    fprintf('K = %-6d | %5d        | %.5f      | %.5f      | %.5f      | %.5f\n', k_val, length(gt), mae_s1, mae_2, rmse_s1, rmse_s2);
end
fprintf('-------------------------------------------------------------\n');

%% === EFFICIENCY + ACCURACY ANALYSIS: RUNTIME (ms), MAE, RMSE ===
% Load Data
data2_table = readtable('experiment_data.xlsx', 'Sheet', 'similarity_search', 'VariableNamingRule', 'preserve');
fprintf('Loaded total: %d rows\n', height(data2_table));

% 1. TIMESTAMPS
target_timestamps = [ ...
"20260201_231314", ...
"20260202_232739", ...
"20260203_143040", ...
"20260201_230706", ...
"20260202_181345", ...
"20260203_143127", ...
"20260203_172314", ...
"20260203_172348" ...
];

% 2. FILTERUNG
% Level = 'bahn'
% Mode  = 'position' (Shape) - Shape und Motion sind zeitlich fast identisch
% Weights = '1,1,1,1,1' (Wir messen das finale Fused-Modell!)
% K == dtw_calls (Volle Pipeline)
filter_idx = strcmpi(data2_table.level, 'bahn') & ...
             strcmp(data2_table.weights, '1,1,1,1,1') & ...
             (data2_table.K == data2_table.dtw_calls) & ...
             ismember(string(data2_table.Timestamp), target_timestamps);

data_filtered = data2_table(filter_idx, :);
fprintf('Filtered rows for Efficiency Analysis: %d\n', height(data_filtered));

% 3. ANALYSE NACH K
unique_K = unique(data_filtered.K);
fprintf('\n=== RESULTS: Runtime (ms), MAE & RMSE ===\n');
fprintf('%-10s | %-12s | %-12s | %-12s | %-12s | %-12s | %-12s\n | %-12s | %-12s | %-12s | %-12s | %-12s | %-12s\n', ...
    'Size (K)', 'Stage 1 (ms)', 'Stage 2 (ms)', 'Loading (ms)', 'Total (ms)', 'MAE-S1-DIR', 'RMSE-S1-DIR', ...
    'MAE-S1-DEC', 'RMSE-S1-DEC', 'MAE-S2-DIR', 'RMSE-S2-DIR', 'MAE-S2-DEC', 'RMSE-S2-DEC');


for i = 1:length(unique_K)
    k_val = unique_K(i);
    sub_data = data_filtered(data_filtered.K == k_val, :);
    
    % Zeiten berechnen (Mittelwert * 1000 für ms)
    t_s1 = mean(sub_data.stage1_time_sec) * 1000;
    t_s2 = mean(sub_data.stage2_time_sec) * 1000;
    t_load = mean(sub_data.loading_time_sec) * 1000;
    t_total = mean(sub_data.total_time_sec) * 1000;
    
    % Ground Truth & Predictions
    gt = sub_data.ground_truth;
    pred_s1_dir = sub_data.direct_s1_weighted;
    pred_s1_dec = sub_data.fromsegs_s1_weighted;
    pred_s2_dir = sub_data.direct_s2_weighted;
    pred_s2_dec = sub_data.fromsegs_s2_weighted;
    
    % MAE & RMSE Berechnung
    mae_val_s1_dec = mean(abs(pred_s1_dec - gt));
    rmse_val_s1_dec = sqrt(mean((pred_s1_dec - gt).^2));
    mae_val_s1_dir = mean(abs(pred_s1_dir - gt));
    rmse_val_s1_dir = sqrt(mean((pred_s1_dir - gt).^2));
    mae_val_s2_dir = mean(abs(pred_s2_dir - gt));
    rmse_val_s2_dir = sqrt(mean((pred_s2_dir - gt).^2));
    mae_val_s2_dec = mean(abs(pred_s2_dec - gt));
    rmse_val_s2_dec = sqrt(mean((pred_s2_dec - gt).^2));
    
    fprintf('K = %-6d | %8.2f ms   | %8.2f ms   | %8.2f ms   | %8.2f ms   | %8.5f mm  | %8.5f mm | %8.5f mm | %8.5f mm | %8.5f mm | %8.5f mm | %8.5f mm | %8.5f mm\n', ...
        k_val, t_s1, t_s2, t_load, t_total, mae_val_s1_dir, rmse_val_s1_dir ,mae_val_s1_dec, rmse_val_s1_dec, mae_val_s2_dir, rmse_val_s2_dir, mae_val_s2_dec, rmse_val_s2_dec);
end
fprintf('---------------------------------------------------------------------------------------------\n');


%% === FINAL PLOT: ACCURACY VS EFFICIENCY (CLEAN LOG SCALE) ===
clear; clc; close all;
% 1. DATA LOADING
filename = 'experiment_data.xlsx';
data_all = readtable(filename, 'Sheet', 'similarity_search', 'VariableNamingRule', 'preserve');
target_timestamps = [ ...
    "20260201_231314", "20260202_232739", "20260203_143040", ...
    "20260201_230706", "20260202_181345", "20260203_143127", "20260203_172314", ...
    "20260203_172348"];
% 2. BERECHNUNG
K_vals = [5, 10, 25, 50];
p_s1_dir = zeros(3, 2); p_s1_dec = zeros(3, 2);
p_s2_dir = zeros(3, 2); p_s2_dec = zeros(3, 2);

fprintf('Calculating metrics...\n');
for i = 1:4
    k = K_vals(i);
    idx_base = strcmp(data_all.weights, '1,1,1,1,1') & ...
               (data_all.K == k) & (data_all.K == data_all.dtw_calls) & ...
               ismember(string(data_all.Timestamp), target_timestamps);
    rows_bahn = data_all(idx_base & strcmpi(data_all.level, 'bahn'), :);
    if isempty(rows_bahn), continue; end
    
    % Metrics
    p_s1_dir(i, 1) = mean(rows_bahn.stage1_time_sec) * 1000;
    p_s1_dir(i, 2) = mean(rows_bahn.err_direct_s1_weighted); 
    p_s1_dec(i, 1) = mean(rows_bahn.stage1_time_sec) * 1000;
    p_s1_dec(i, 2) = mean(rows_bahn.err_fromsegs_s1_weighted); 
    p_s2_dir(i, 1) = mean(rows_bahn.total_time_sec) * 1000;
    p_s2_dir(i, 2) = mean(rows_bahn.err_direct_s2_weighted);
    p_s2_dec(i, 2) = mean(rows_bahn.err_fromsegs_s2_weighted);
    
    % Segment Time Summation
    rows_seg_all = data_all(idx_base & contains(lower(data_all.level), 'seg'), :);
    total_times_sum = zeros(height(rows_bahn), 1);
    for r = 1:height(rows_bahn)
        ids_str = string(rows_bahn.segment_id{r});
        if ismissing(ids_str) || strlength(ids_str)==0, continue; end
        mask_segs = ismember(rows_seg_all.query_id, ids_str);
        if any(mask_segs), total_times_sum(r) = sum(rows_seg_all.total_time_sec(mask_segs)); end
    end
    p_s2_dec(i, 1) = mean(total_times_sum) * 1000;
end

%% === PLOTTING ===
fig = figure('Position', [100, 100, 750, 500], 'Color', 'w');
ax = axes('XScale', 'log', 'YScale', 'linear', 'XMinorGrid', 'off');
hold on; grid on; box off;

c_dir = [0.86, 0.13, 0.15]; c_dec = [0.15, 0.39, 0.91];
sz = 120;

% Verbindungslinien
for i = 1:4
    plot([p_s1_dir(i,1), p_s2_dir(i,1)], [p_s1_dir(i,2), p_s2_dir(i,2)], '-', 'Color', [c_dir, 0.5], 'LineWidth', 2.5); 
    plot([p_s1_dec(i,1), p_s2_dec(i,1)], [p_s1_dec(i,2), p_s2_dec(i,2)], '-', 'Color', [c_dec, 0.5], 'LineWidth', 2.5);
end

% Marker
h_s1_dir = scatter(p_s1_dir(:,1), p_s1_dir(:,2), sz, c_dir, 'filled', '^', 'MarkerEdgeColor', 'k');
h_s2_dir = scatter(p_s2_dir(:,1), p_s2_dir(:,2), sz, c_dir, 'filled', 'o', 'MarkerEdgeColor', 'k');
h_s1_dec = scatter(p_s1_dec(:,1), p_s1_dec(:,2), sz, c_dec, 'filled', '^', 'MarkerEdgeColor', 'k');
h_s2_dec = scatter(p_s2_dec(:,1), p_s2_dec(:,2), sz, c_dec, 'filled', 'o', 'MarkerEdgeColor', 'k');

% Labels
font_args = {'FontName', 'Times New Roman', 'FontSize', 18, 'FontWeight', 'bold'};
for i = 1:4
    text(p_s2_dir(i,1)*1.15, p_s2_dir(i,2), sprintf('K=%d  ', K_vals(i)), 'Color', c_dir, 'HorizontalAlignment', 'left', font_args{:});
    text(p_s2_dec(i,1)*1.15, p_s2_dec(i,2), sprintf('K=%d  ', K_vals(i)), 'Color', c_dec, 'HorizontalAlignment', 'left', font_args{:});
end

xlabel('Inference time [ms]', 'FontWeight', 'bold', 'FontSize', 16, 'FontName', 'Times New Roman');
ylabel('MAE [mm]', 'FontWeight', 'bold', 'FontSize', 16, 'FontName', 'Times New Roman');
yticks([0.0, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 1.0])
ylim([0.02, 0.09])
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, 'LineWidth', 1, 'FontWeight', 'bold');
xticks([50 100 200 500 1000 2000 5000 10000]);
xticklabels({'50', '100', '200', '500', '1000', '2000', '5000', '10000'});
xlim([40, 15000]);

[~, objh] = legend([h_s1_dir,  h_s2_dir, h_s1_dec, h_s2_dec], ...
    {'Stage 1 (Direct)', 'Stage 2 (Direct)', 'Stage 1 (Decomposed)',  'Stage 2 (Decomposed)'}, ...
    'Location', 'northeast', 'FontSize', 15, 'FontWeight', 'bold');

objhl = findobj(objh, 'type', 'patch'); %// objects of legend of type line
set(objhl, 'Markersize', 13); %// set marker size as desired

exportgraphics(gca, fullfile('figs', 'efficiency_tradeoff_log_mae.pdf'));
fprintf('Plot saved: efficiency_tradeoff_log_clean.pdf\n');

%% === TRADE-OFF PLOT: S1 vs DIRECT vs DECOMPOSED (RMSE CALCULATED) ===
clear; clc; close all;

% 1. DATA LOADING
filename = 'experiment_data.xlsx';
data_all = readtable(filename, 'Sheet', 'similarity_search', 'VariableNamingRule', 'preserve');
target_timestamps = [ ...
    "20260201_231314", "20260202_232739", "20260203_143040", ...
    "20260201_230706", "20260202_181345", "20260203_143127", "20260203_172314", ...
    "20260203_172348"];

% 2. BERECHNUNG
K_vals = [5, 10, 25, 50];
p_s1_dir = zeros(3, 2); p_s1_dec = zeros(3, 2);
p_s2_dir = zeros(3, 2); p_s2_dec = zeros(3, 2);

fprintf('Calculating RMSE from deviations...\n');

for i = 1:4
    k = K_vals(i);
    
    idx_base = strcmp(data_all.weights, '1,1,1,1,1') & ...
               (data_all.K == k) & (data_all.K == data_all.dtw_calls) & ...
               ismember(string(data_all.Timestamp), target_timestamps);
           
    rows_bahn = data_all(idx_base & strcmpi(data_all.level, 'bahn'), :);
    if isempty(rows_bahn), continue; end
    
    % --- RMSE BERECHNUNG (sqrt(mean(x^2))) ---
    
    % 1. STAGE 1 (Baseline)
    p_s1_dir(i, 1) = mean(rows_bahn.stage1_time_sec) * 1000;
    p_s1_dir(i, 2) = sqrt(mean(rows_bahn.err_direct_s1_weighted .^ 2)); % RMSE berechnen
    
    p_s1_dec(i, 1) = mean(rows_bahn.stage1_time_sec) * 1000;
    p_s1_dec(i, 2) = sqrt(mean(rows_bahn.err_fromsegs_s1_weighted .^ 2)); % RMSE berechnen
    
    % 2. STAGE 2 DIRECT
    p_s2_dir(i, 1) = mean(rows_bahn.total_time_sec) * 1000;
    p_s2_dir(i, 2) = sqrt(mean(rows_bahn.err_direct_s2_weighted .^ 2)); % RMSE berechnen
    
    % 3. STAGE 2 DECOMPOSED
    p_s2_dec(i, 2) = sqrt(mean(rows_bahn.err_fromsegs_s2_weighted .^ 2)); % RMSE berechnen
    
    % Zeit: ECHTE Berechnung über Segment-IDs (Wie gehabt)
    rows_seg_all = data_all(idx_base & contains(lower(data_all.level), 'seg'), :);
    total_times_sum = zeros(height(rows_bahn), 1);
    for r = 1:height(rows_bahn)
        ids_str = string(rows_bahn.segment_id{r});
        if ismissing(ids_str) || strlength(ids_str)==0, continue; end
        mask_segs = ismember(rows_seg_all.query_id, ids_str);
        if any(mask_segs), total_times_sum(r) = sum(rows_seg_all.total_time_sec(mask_segs)); end
    end
    p_s2_dec(i, 1) = mean(total_times_sum) * 1000;
end

%% === PLOTTING (PERFECT STYLE - RMSE) ===
fig = figure('Position', [100, 100, 750, 500], 'Color', 'w');
ax = axes('XScale', 'log', 'YScale', 'linear', 'XMinorGrid', 'off');
hold on; grid on; box off;

c_dir = [0.86, 0.13, 0.15]; 
c_dec = [0.15, 0.39, 0.91];
sz = 120;

% Verbindungslinien
for i = 1:4
    plot([p_s1_dir(i,1), p_s2_dir(i,1)], [p_s1_dir(i,2), p_s2_dir(i,2)], '-', 'Color', [c_dir, 0.5], 'LineWidth', 2.5); 
    plot([p_s1_dec(i,1), p_s2_dec(i,1)], [p_s1_dec(i,2), p_s2_dec(i,2)], '-', 'Color', [c_dec, 0.5], 'LineWidth', 2.5);
end

% Marker
h_s1_dir = scatter(p_s1_dir(:,1), p_s1_dir(:,2), sz, c_dir, 'filled', '^', 'MarkerEdgeColor', 'k');
h_s2_dir = scatter(p_s2_dir(:,1), p_s2_dir(:,2), sz, c_dir, 'filled', 'o', 'MarkerEdgeColor', 'k');
h_s1_dec = scatter(p_s1_dec(:,1), p_s1_dec(:,2), sz, c_dec, 'filled', '^', 'MarkerEdgeColor', 'k');
h_s2_dec = scatter(p_s2_dec(:,1), p_s2_dec(:,2), sz, c_dec, 'filled', 'o', 'MarkerEdgeColor', 'k');

% Labels
font_args = {'FontName', 'Times New Roman', 'FontSize', 18, 'FontWeight', 'bold'};
for i = 1:4
    text(p_s2_dir(i,1)*1.15, p_s2_dir(i,2), sprintf('K=%d  ', K_vals(i)), 'Color', c_dir, 'HorizontalAlignment', 'left', font_args{:});
    text(p_s2_dec(i,1)*1.15, p_s2_dec(i,2), sprintf('K=%d  ', K_vals(i)), 'Color', c_dec, 'HorizontalAlignment', 'left', font_args{:});
end

xlabel('Inference time [ms]', 'FontWeight', 'bold', 'FontSize', 16, 'FontName', 'Times New Roman');
ylabel('RMSE [mm]', 'FontWeight', 'bold', 'FontSize', 16, 'FontName', 'Times New Roman');

% ACHTUNG: Limits angepasst für RMSE (Werte sind höher als MAE)
yticks([0.0, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 1.0])
ylim([0.02, 0.09])

set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, 'LineWidth', 1, 'FontWeight', 'bold');
xticks([50 100 200 500 1000 2000 5000 10000]);
xticklabels({'50', '100', '200', '500', '1000', '2000', '5000', '10000'});
xlim([40, 15000]);

[~, objh] = legend([h_s1_dir,  h_s2_dir, h_s1_dec, h_s2_dec], ...
    {'Stage 1 (Direct)', 'Stage 2 (Direct)', 'Stage 1 (Decomposed)',  'Stage 2 (Decomposed)'}, ...
    'Location', 'southeast', 'FontSize', 14, 'FontWeight', 'bold');
objhl = findobj(objh, 'type', 'patch'); 
set(objhl, 'Markersize', 13); 

%exportgraphics(gca, fullfile('figs', 'efficiency_tradeoff_log_rmse.pdf'));
fprintf('Plot saved: efficiency_tradeoff_log_rmse.pdf\n');
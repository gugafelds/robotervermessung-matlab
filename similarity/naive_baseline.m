%  NAIVE MEAN BASELINE
%  ========================================================================
%  Uses the global mean SIDTW path error as a constant prediction
%  for all 500 query trajectories. Simplest possible baseline.
%  ========================================================================

clear; clc;

addpath(genpath(pwd));
addpath(genpath('../main'));

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║  NAIVE MEAN BASELINE                                           ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

% ========================================================================
%% SECTION 1: CONFIGURATION
% ========================================================================

queries_quantity = 20000;
random_seed = 21;

exclude_ids = {
    '1765989370'; '1765989294'; '1765988821'; '1765988920'; '1765989411';
    '1765990630'; '1765990747'; '1765990822'; '1765991047'; '1765991234';
    '1765991190'; '1765991445'; '1765991515'; '1765991949'; '1765991743';
    '1769770498'; '1769770684'; '1769770935'; '1769771107'; '1769771447';
    '1769773928'; '1769772060'; '1769772213'; '1769773985'; '1769774278';
    '1769772609'; '1769773593'; '1769772776'; '1769772900'; '1769773333';
    '1769774581';
};

fprintf('Queries:        %d\n', queries_quantity);
fprintf('Random seed:    %d\n\n', random_seed);

% ========================================================================
%% SECTION 2: LOAD DATA
% ========================================================================

fprintf('═══ LOADING DATA ═══\n\n');

conn = connectingToPostgres();
schema = 'bewegungsdaten';

% --- All bahn-level ground truth (for global mean) ---
fprintf('Loading all bahn-level SIDTW values...\n');
all_gt_sql = sprintf([...
    'SELECT s.segment_id, s.sidtw_average_distance ' ...
    'FROM auswertung.info_sidtw s ' ...
    'INNER JOIN %s.bahn_metadata m ON s.segment_id = m.bahn_id ' ...
    'WHERE m.segment_id = m.bahn_id'], schema);
all_gt = fetch(conn, all_gt_sql);
fprintf('  ✓ Loaded ground truth for %d trajectories\n', height(all_gt));

% --- Select query IDs (same as main evaluation) ---
ids_str_list = sprintf('''%s'',', exclude_ids{:});
ids_str_list = ids_str_list(1:end-1);

query_sql = sprintf([...
    'SELECT bahn_id FROM (' ...
    'SELECT DISTINCT b.bahn_id FROM %s.bahn_metadata b ' ...
    'INNER JOIN auswertung.info_sidtw s ON b.bahn_id = s.segment_id ' ...
    'WHERE b.segment_id = b.bahn_id ' ...
    'AND b.bahn_id NOT IN (' ...
        'SELECT t_all.bahn_id ' ...
        'FROM %s.bahn_info t_exclude ' ...
        'JOIN %s.bahn_info t_all ON t_exclude.record_filename = t_all.record_filename ' ...
        'WHERE t_exclude.bahn_id IN (%s) ' ...
    ') ' ...
    ') AS distinct_bahnen ' ...
    'ORDER BY md5(bahn_id || ''%d'') LIMIT %d'], ...
    schema, schema, schema, ids_str_list, random_seed, queries_quantity);
query_results = fetch(conn, query_sql);
query_ids = query_results.bahn_id;

close(conn);

fprintf('  ✓ Selected %d query trajectories\n\n', length(query_ids));

% ========================================================================
%% SECTION 3: COMPUTE GLOBAL MEAN
% ========================================================================

global_mean = mean(all_gt.sidtw_average_distance);
global_median = median(all_gt.sidtw_average_distance);
global_std = std(all_gt.sidtw_average_distance);

fprintf('═══ DATABASE STATISTICS ═══\n\n');
fprintf('  Total trajectories:  %d\n', height(all_gt));
fprintf('  Global mean:         %.4f mm\n', global_mean);
fprintf('  Global median:       %.4f mm\n', global_median);
fprintf('  Global std:          %.4f mm\n\n', global_std);

% ========================================================================
%% SECTION 4: EVALUATE ON QUERY SET
% ========================================================================

fprintf('═══ NAIVE BASELINE PREDICTION ═══\n\n');

% Get ground truth for query trajectories
query_gt = nan(length(query_ids), 1);
for i = 1:length(query_ids)
    idx = find(strcmp(all_gt.segment_id, query_ids{i}), 1);
    if ~isempty(idx)
        query_gt(i) = all_gt.sidtw_average_distance(idx);
    end
end

valid = ~isnan(query_gt);
query_gt = query_gt(valid);
n_valid = sum(valid);

fprintf('  Valid queries: %d / %d\n\n', n_valid, length(query_ids));

% --- Mean predictor ---
errors_mean = abs(query_gt - global_mean);
mae_mean = mean(errors_mean);
rmse_mean = sqrt(mean(errors_mean.^2));
max_err_mean = max(errors_mean);

% --- Median predictor ---
errors_median = abs(query_gt - global_median);
mae_median = mean(errors_median);
rmse_median = sqrt(mean(errors_median.^2));
max_err_median = max(errors_median);

fprintf('  ┌─────────────────┬──────────────┬──────────────┐\n');
fprintf('  │ Metric          │ Mean Pred.   │ Median Pred. │\n');
fprintf('  ├─────────────────┼──────────────┼──────────────┤\n');
fprintf('  │ Prediction      │ %.4f mm   │ %.4f mm   │\n', global_mean, global_median);
fprintf('  │ MAE             │ %.4f mm   │ %.4f mm   │\n', mae_mean, mae_median);
fprintf('  │ RMSE            │ %.4f mm   │ %.4f mm   │\n', rmse_mean, rmse_median);
fprintf('  │ Max Error       │ %.4f mm   │ %.4f mm   │\n', max_err_mean, max_err_median);
fprintf('  └─────────────────┴──────────────┴──────────────┘\n\n');

% ========================================================================
%% SECTION 5: DISTRIBUTION ANALYSIS
% ========================================================================

fprintf('═══ QUERY SET STATISTICS ═══\n\n');
fprintf('  Query mean:          %.4f mm\n', mean(query_gt));
fprintf('  Query median:        %.4f mm\n', median(query_gt));
fprintf('  Query std:           %.4f mm\n', std(query_gt));
fprintf('  Query min:           %.4f mm\n', min(query_gt));
fprintf('  Query max:           %.4f mm\n\n', max(query_gt));

% Error percentiles
fprintf('  Error percentiles (Mean Predictor):\n');
fprintf('    25th: %.4f mm\n', prctile(errors_mean, 25));
fprintf('    50th: %.4f mm\n', prctile(errors_mean, 50));
fprintf('    75th: %.4f mm\n', prctile(errors_mean, 75));
fprintf('    90th: %.4f mm\n', prctile(errors_mean, 90));
fprintf('    95th: %.4f mm\n\n', prctile(errors_mean, 95));

% ========================================================================
%% SECTION 6: SAVE
% ========================================================================

timestamp = datestr(now, 'yyyymmdd_HHMMSS');
if ~exist('results', 'dir'), mkdir('results'); end

results = table(query_ids(valid), query_gt, ...
    repmat(global_mean, n_valid, 1), errors_mean, ...
    repmat(global_median, n_valid, 1), errors_median, ...
    'VariableNames', {'query_id', 'ground_truth', ...
    'pred_mean', 'err_mean', 'pred_median', 'err_median'});

filename = sprintf('results/naive_baseline_%s.csv', timestamp);
writetable(results, filename);

fprintf('✓ Saved: %s\n', filename);
fprintf('✓ Done!\n');
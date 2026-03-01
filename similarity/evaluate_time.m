%% === EFFICIENCY + ACCURACY ANALYSIS: RUNTIME (ms), MAE, RMSE === // DIRECT + DECOMPOSED
% Load Data
data2_table = readtable('experiment_data.xlsx', 'Sheet', 'similarity_search', 'VariableNamingRule', 'preserve');
%fprintf('Loaded total: %d rows\n', height(data2_table));

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
%%
% 2. FILTERUNG - Bahn level (Direct)
filter_idx = strcmpi(data2_table.level, 'bahn') & ...
             strcmp(data2_table.weights, '1,1,1,1,1') & ...
             (data2_table.K == data2_table.dtw_calls) & ...
             ismember(string(data2_table.Timestamp), target_timestamps);
data_filtered = data2_table(filter_idx, :);
fprintf('Filtered rows (bahn): %d\n', height(data_filtered));

% Pre-filter segment rows (for Decomposed timing)
filter_seg_idx = contains(lower(data2_table.level), 'seg') & ...
                 strcmp(data2_table.weights, '1,1,1,1,1') & ...
                 (data2_table.K == data2_table.dtw_calls) & ...
                 ismember(string(data2_table.Timestamp), target_timestamps);
data_seg_all = data2_table(filter_seg_idx, :);
fprintf('Filtered rows (seg):  %d\n', height(data_seg_all));

% 3. ANALYSE NACH K
unique_K = unique(data_filtered.K);

fprintf('\n========================================================================\n');
fprintf('=== TIMING STATISTICS (ms) - Direct vs Decomposed ===\n');
fprintf('========================================================================\n');

for i = 1:length(unique_K)
    k_val = unique_K(i);
    filter_join = strcmpi(data_filtered.dtw_mode, 'position') & ...
                 (data_filtered.K == k_val);
    filter_joint = strcmp(data_filtered.dtw_mode, 'joint_states') & ...
                 (data_filtered.K == k_val);
    sub_data_pos = data_filtered(filter_join, :);
    sub_data_joint = data_filtered(filter_joint, :);


    % --- Direct Zeiten Pos (in ms) ---
    t_s1_all_pos     = sub_data_pos.stage1_time_sec * 1000;
    t_s2_dir_all_pos = sub_data_pos.stage2_time_sec * 1000;
    t_total_dir_all_pos = sub_data_pos.total_time_sec * 1000;

    % --- Direct Zeiten Joint (in ms) ---
    t_s1_all_joint     = sub_data_joint.stage1_time_sec * 1000;
    t_s2_dir_all_joint = sub_data_joint.stage2_time_sec * 1000;
    t_total_dir_all_joint = sub_data_joint.total_time_sec * 1000;

    t_s1_all = (t_s1_all_pos+t_s1_all_joint)/2;
    t_s2_dir_all = (t_s2_dir_all_pos+t_s2_dir_all_joint)/2;
    t_total_dir_all = (t_total_dir_all_pos+t_total_dir_all_joint)/2;

    % --- Decomposed Zeiten (Summe der total_time_sec der Segmente pro Bahn) ---
    filter_join = strcmpi(data_seg_all.dtw_mode, 'position') & ...
                 (data_seg_all.K == k_val);
    seg_rows_k_pos = data_seg_all(filter_join, :);

    t_total_dec_all_pos = zeros(height(sub_data_pos), 1);
    t_s2_dec_all_pos    = zeros(height(sub_data_pos), 1);

    for r = 1:height(sub_data_pos)
        ids_str = string(sub_data_pos.segment_id{r});
        
        segs_filter = startsWith(seg_rows_k_pos.segment_id, ids_str);
        segs = seg_rows_k_pos(segs_filter, :);

        t_total_dec_all_pos(r) = sum(seg_rows_k_pos.total_time_sec(segs_filter)) * 1000;
        t_s2_dec_all_pos(r)    = sum(seg_rows_k_pos.stage2_time_sec(segs_filter)) * 1000;
    end

    % --- Decomposed Zeiten (Summe der total_time_sec der Segmente pro Bahn) ---
    filter_join = strcmpi(data_seg_all.dtw_mode, 'joint_states') & ...
                 (data_seg_all.K == k_val);
    seg_rows_k_joint = data_seg_all(filter_join, :);

    t_total_dec_all_joint = zeros(height(sub_data_joint), 1);
    t_s2_dec_all_joint    = zeros(height(sub_data_joint), 1);

    for r = 1:height(sub_data_joint)
        ids_str = string(sub_data_joint.segment_id{r});
        
        segs_filter = startsWith(seg_rows_k_joint.segment_id, ids_str);
        segs = seg_rows_k_joint(segs_filter, :);

        t_total_dec_all_joint(r) = sum(seg_rows_k_joint.total_time_sec(segs_filter)) * 1000;
        t_s2_dec_all_joint(r)    = sum(seg_rows_k_joint.stage2_time_sec(segs_filter)) * 1000;
    end

    t_s2_dec_all = (t_s2_dec_all_pos+t_s2_dec_all_joint)/2;
    t_total_dec_all = (t_total_dec_all_pos+t_total_dec_all_joint)/2;
    

    fprintf('\n--- K = %d ---\n', k_val);
    fprintf('%-20s | %8s | %8s | %8s | %8s | %8s\n', ...
        '', 'Mean', 'Median', 'Std', 'P90', 'Max');
    fprintf('%s\n', repmat('-', 1, 75));

    % Stage 1
    fprintf('%-20s | %8.1f | %8.1f | %8.1f | %8.1f | %8.1f\n', ...
        'Stage 1', ...
        mean(t_s1_all, 'omitnan'), median(t_s1_all, 'omitnan'), ...
        std(t_s1_all, 'omitnan'), prctile(t_s1_all, 90), max(t_s1_all));

    % Stage 2 Direct
    fprintf('%-20s | %8.1f | %8.1f | %8.1f | %8.1f | %8.1f\n', ...
        'S2 Direct', ...
        mean(t_s2_dir_all, 'omitnan'), median(t_s2_dir_all, 'omitnan'), ...
        std(t_s2_dir_all, 'omitnan'), prctile(t_s2_dir_all, 90), max(t_s2_dir_all));

    % Stage 2 Decomposed
    fprintf('%-20s | %8.1f | %8.1f | %8.1f | %8.1f | %8.1f\n', ...
        'S2 Decomposed', ...
        mean(t_s2_dec_all, 'omitnan'), median(t_s2_dec_all, 'omitnan'), ...
        std(t_s2_dec_all, 'omitnan'), prctile(t_s2_dec_all, 90), max(t_s2_dec_all));

    fprintf('%s\n', repmat('-', 1, 75));

    % Total Direct
    fprintf('%-20s | %8.1f | %8.1f | %8.1f | %8.1f | %8.1f\n', ...
        'Total Direct', ...
        mean(t_total_dir_all, 'omitnan'), median(t_total_dir_all, 'omitnan'), ...
        std(t_total_dir_all, 'omitnan'), prctile(t_total_dir_all, 90), max(t_total_dir_all));

    % Total Decomposed
    fprintf('%-20s | %8.1f | %8.1f | %8.1f | %8.1f | %8.1f\n', ...
        'Total Decomposed', ...
        mean(t_total_dec_all, 'omitnan'), median(t_total_dec_all, 'omitnan'), ...
        std(t_total_dec_all, 'omitnan'), prctile(t_total_dec_all, 90), max(t_total_dec_all));

    % --- Ground Truth & Predictions ---
    gt = data_filtered.ground_truth;
    pred_s1_dir = data_filtered.direct_s1_weighted;
    pred_s1_dec = data_filtered.fromsegs_s1_weighted;
    pred_s2_dir = data_filtered.direct_s2_weighted;
    pred_s2_dec = data_filtered.fromsegs_s2_weighted;

    % --- MAE & RMSE ---
    mae_s1_dir  = mean(abs(pred_s1_dir - gt));
    mae_s1_dec  = mean(abs(pred_s1_dec - gt));
    mae_s2_dir  = mean(abs(pred_s2_dir - gt));
    mae_s2_dec  = mean(abs(pred_s2_dec - gt));
    rmse_s1_dir = sqrt(mean((pred_s1_dir - gt).^2));
    rmse_s1_dec = sqrt(mean((pred_s1_dec - gt).^2));
    rmse_s2_dir = sqrt(mean((pred_s2_dir - gt).^2));
    rmse_s2_dec = sqrt(mean((pred_s2_dec - gt).^2));

    fprintf('%s\n', repmat('-', 1, 75));
    fprintf('%-20s | %8s | %8s\n', 'Prediction', 'MAE', 'RMSE');
    fprintf('%s\n', repmat('-', 1, 45));
    fprintf('%-20s | %8.4f | %8.4f\n', 'S1 Direct',      mae_s1_dir, rmse_s1_dir);
    fprintf('%-20s | %8.4f | %8.4f\n', 'S1 Decomposed',  mae_s1_dec, rmse_s1_dec);
    fprintf('%-20s | %8.4f | %8.4f\n', 'S2 Direct',      mae_s2_dir, rmse_s2_dir);
    fprintf('%-20s | %8.4f | %8.4f\n', 'S2 Decomposed',  mae_s2_dec, rmse_s2_dec);
end

fprintf('\n========================================================================\n');
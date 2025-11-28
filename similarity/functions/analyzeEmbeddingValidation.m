function analyzeEmbeddingValidation(csv_file)
    % ANALYZEEMBEDDINGVALIDATION - Comprehensive analysis of embedding validation results
    %
    % Usage:
    %   analyzeEmbeddingValidation('similarity/results/embedding_validation_2025-11-28.csv')
    %
    % Input:
    %   csv_file - Path to CSV file with experiment results
    %
    % Generates:
    %   - Console output with 5 comprehensive analyses
    %   - Recommendations for best configurations
    
    fprintf('\n========================================\n');
    fprintf('EMBEDDING VALIDATION ANALYSIS\n');
    fprintf('========================================\n\n');
    
    % Load results
    fprintf('Reading: %s\n\n', csv_file);
    results_table = readtable(csv_file);
    
    % Extract arrays for analysis
    embedding_names = cellstr(results_table.Embedding_Config);
    query_bahn_ids = string(results_table.Query_Bahn_ID);  % ← DIESER
    weight_modes = cellstr(results_table.Weight_Mode);
    dtw_modes = cellstr(results_table.DTW_Mode);
    levels = cellstr(results_table.Level);
    
    % Dimensions
    total_dims_arr = results_table.Total_Dims;
    
    % Separate trajectory and segment results
    traj_mask = strcmp(levels, 'Trajectory');
    seg_mask = strcmp(levels, 'Segment');
    
    % Metrics columns
    metric_cols = {'Spearman', 'P_50', 'P_10', 'P_5', 'P_3', 'P_1'};
    
    % Build results matrices
    traj_results_matrix = zeros(sum(traj_mask), length(metric_cols));
    seg_results_matrix = zeros(sum(seg_mask), length(metric_cols));
    
    for m = 1:length(metric_cols)
        if ismember(metric_cols{m}, results_table.Properties.VariableNames)
            traj_results_matrix(:, m) = results_table.(metric_cols{m})(traj_mask);
            seg_results_matrix(:, m) = results_table.(metric_cols{m})(seg_mask);
        end
    end
    
    % Extract names for trajectory-level only (segment has same structure)
    embedding_names_traj = embedding_names(traj_mask);
    query_bahn_ids_traj = query_bahn_ids(traj_mask);
    weight_modes_traj = weight_modes(traj_mask);
    dtw_modes_traj = dtw_modes(traj_mask);
    total_dims_arr_traj = total_dims_arr(traj_mask);
    
    embedding_names_seg = embedding_names(seg_mask);
    query_bahn_ids_seg = query_bahn_ids(seg_mask);
    weight_modes_seg = weight_modes(seg_mask);
    
    total_experiments_traj = sum(traj_mask);
    total_experiments_seg = sum(seg_mask);
    
    fprintf('Total experiments: %d trajectory + %d segment = %d\n\n', ...
        total_experiments_traj, total_experiments_seg, height(results_table));
    
    %% --- ANALYSIS 1: Performance by Embedding Architecture ---
    fprintf('========================================\n');
    fprintf('ANALYSIS 1: Performance by Embedding Architecture\n');
    fprintf('========================================\n\n');
    
    unique_embeddings = unique(embedding_names_traj);
    
    for e = 1:length(unique_embeddings)
        emb_name = unique_embeddings{e};
        emb_idx_traj = strcmp(embedding_names_traj, emb_name);
        emb_idx_seg = strcmp(embedding_names_seg, emb_name);
        
        % Trajectory level
        avg_spearman_traj = mean(traj_results_matrix(emb_idx_traj, 1));
        avg_p10_traj = mean(traj_results_matrix(emb_idx_traj, 3));
        avg_p1_traj = mean(traj_results_matrix(emb_idx_traj, 6));
        
        % Segment level
        avg_spearman_seg = mean(seg_results_matrix(emb_idx_seg, 1), 'omitnan');
        avg_p10_seg = mean(seg_results_matrix(emb_idx_seg, 3), 'omitnan');
        avg_p1_seg = mean(seg_results_matrix(emb_idx_seg, 6), 'omitnan');
        
        % Get dims
        dims = unique(total_dims_arr_traj(emb_idx_traj));
        
        fprintf('%s (Total Dims: %d):\n', emb_name, dims);
        fprintf('  Trajectory: ρ=%.4f, P@10=%.3f, P@1=%.3f\n', ...
            avg_spearman_traj, avg_p10_traj, avg_p1_traj);
        fprintf('  Segment:    ρ=%.4f, P@10=%.3f, P@1=%.3f\n\n', ...
            avg_spearman_seg, avg_p10_seg, avg_p1_seg);
    end
    
    %% --- ANALYSIS 2: Performance by Query Trajectory ---
    fprintf('========================================\n');
    fprintf('ANALYSIS 2: Performance by Query Trajectory\n');
    fprintf('========================================\n\n');
    
    unique_queries = unique(query_bahn_ids_traj);
    
    for q = 1:length(unique_queries)
        query_id = unique_queries{q};
        query_idx_traj = strcmp(query_bahn_ids_traj, query_id);
        query_idx_seg = strcmp(query_bahn_ids_seg, query_id);
        
        % Trajectory level
        avg_spearman_traj = mean(traj_results_matrix(query_idx_traj, 1));
        avg_p10_traj = mean(traj_results_matrix(query_idx_traj, 3));
        avg_p1_traj = mean(traj_results_matrix(query_idx_traj, 6));
        
        % Segment level
        avg_spearman_seg = mean(seg_results_matrix(query_idx_seg, 1), 'omitnan');
        avg_p10_seg = mean(seg_results_matrix(query_idx_seg, 3), 'omitnan');
        avg_p1_seg = mean(seg_results_matrix(query_idx_seg, 6), 'omitnan');
        
        fprintf('Query %s:\n', query_id);
        fprintf('  Trajectory: ρ=%.4f, P@10=%.3f, P@1=%.3f\n', ...
            avg_spearman_traj, avg_p10_traj, avg_p1_traj);
        fprintf('  Segment:    ρ=%.4f, P@10=%.3f, P@1=%.3f\n\n', ...
            avg_spearman_seg, avg_p10_seg, avg_p1_seg);
    end
    
    %% --- ANALYSIS 3: Performance by Weight-Mode Configuration ---
    fprintf('========================================\n');
    fprintf('ANALYSIS 3: Performance by Weight-Mode Configuration\n');
    fprintf('========================================\n\n');
    
    unique_weight_modes = unique(weight_modes_traj);
    
    % Separate into Joint-based and Position-based
    joint_modes = unique_weight_modes(contains(unique_weight_modes, 'Joint'));
    pos_modes = unique_weight_modes(contains(unique_weight_modes, {'Position', 'Pos'}));
    
    fprintf('JOINT-BASED CONFIGURATIONS (Trajectory Level):\n');
    for w = 1:length(joint_modes)
        wm_name = joint_modes{w};
        wm_idx = strcmp(weight_modes_traj, wm_name);
        
        avg_spearman_traj = mean(traj_results_matrix(wm_idx, 1));
        avg_p10_traj = mean(traj_results_matrix(wm_idx, 3));
        
        fprintf('  %s: ρ=%.4f, P@10=%.3f\n', wm_name, avg_spearman_traj, avg_p10_traj);
    end
    
    fprintf('\nPOSITION-BASED CONFIGURATIONS (Trajectory Level):\n');
    for w = 1:length(pos_modes)
        wm_name = pos_modes{w};
        wm_idx = strcmp(weight_modes_traj, wm_name);
        
        avg_spearman_traj = mean(traj_results_matrix(wm_idx, 1));
        avg_p10_traj = mean(traj_results_matrix(wm_idx, 3));
        
        fprintf('  %s: ρ=%.4f, P@10=%.3f\n', wm_name, avg_spearman_traj, avg_p10_traj);
    end
    
    %% --- ANALYSIS 4: Top 5 Overall Configurations ---
    fprintf('\n========================================\n');
    fprintf('ANALYSIS 4: Top 5 Overall Configurations\n');
    fprintf('========================================\n\n');
    
    % Trajectory level
    [sorted_spearman_traj, sort_idx_traj] = sort(traj_results_matrix(:, 1), 'descend');
    
    fprintf('TRAJECTORY LEVEL (by Spearman ρ):\n');
    for i = 1:min(5, total_experiments_traj)
        idx = sort_idx_traj(i);
        fprintf('%d. ρ=%.4f | %s | Query_%s | %s\n', ...
            i, sorted_spearman_traj(i), ...
            embedding_names_traj{idx}, query_bahn_ids_traj{idx}, weight_modes_traj{idx});
    end
    
    % Segment level
    valid_seg_idx = ~isnan(seg_results_matrix(:, 1));
    valid_seg_spearman = seg_results_matrix(valid_seg_idx, 1);
    valid_seg_indices = find(valid_seg_idx);
    [sorted_spearman_seg, sort_idx_seg_rel] = sort(valid_seg_spearman, 'descend');
    sort_idx_seg = valid_seg_indices(sort_idx_seg_rel);
    
    fprintf('\nSEGMENT LEVEL (by Spearman ρ):\n');
    for i = 1:min(5, length(sort_idx_seg))
        idx = sort_idx_seg(i);
        fprintf('%d. ρ=%.4f | %s | Query_%s | %s\n', ...
            i, sorted_spearman_seg(i), ...
            embedding_names_seg{idx}, query_bahn_ids_seg{idx}, weight_modes_seg{idx});
    end
    
    %% --- ANALYSIS 5: Performance vs Dimensionality Trade-off ---
    fprintf('\n========================================\n');
    fprintf('ANALYSIS 5: Performance vs. Dimensionality Trade-off\n');
    fprintf('========================================\n\n');
    
    unique_dims = unique(total_dims_arr_traj);
    unique_dims = sort(unique_dims);
    
    fprintf('Average Performance by Total Dimensions (Trajectory Level):\n');
    for d = 1:length(unique_dims)
        dims = unique_dims(d);
        dims_idx = (total_dims_arr_traj == dims);
        
        avg_spearman = mean(traj_results_matrix(dims_idx, 1));
        avg_p10 = mean(traj_results_matrix(dims_idx, 3));
        
        % Find embedding config with these dims
        emb_config = unique(embedding_names_traj(dims_idx));
        
        fprintf('  %d dims (%s): ρ=%.4f, P@10=%.3f\n', ...
            dims, emb_config{1}, avg_spearman, avg_p10);
    end
    
    %% --- FINAL RECOMMENDATIONS ---
    fprintf('\n========================================\n');
    fprintf('KEY RECOMMENDATIONS\n');
    fprintf('========================================\n\n');
    
    % Best embedding overall
    emb_avg_spearman = zeros(length(unique_embeddings), 1);
    for e = 1:length(unique_embeddings)
        emb_idx = strcmp(embedding_names_traj, unique_embeddings{e});
        emb_avg_spearman(e) = mean(traj_results_matrix(emb_idx, 1));
    end
    [best_emb_spearman, best_emb_idx] = max(emb_avg_spearman);
    best_embedding = unique_embeddings{best_emb_idx};
    
    fprintf('1. BEST EMBEDDING ARCHITECTURE (averaged across all scenarios):\n');
    fprintf('   %s - Average ρ=%.4f\n\n', best_embedding, best_emb_spearman);
    
    % Best single configuration
    [best_single_spearman, best_single_idx] = max(traj_results_matrix(:, 1));
    fprintf('2. BEST SINGLE CONFIGURATION:\n');
    fprintf('   %s | Query_%s | %s\n', ...
        embedding_names_traj{best_single_idx}, query_bahn_ids_traj(best_single_idx), ...
        weight_modes_traj{best_single_idx});
    fprintf('   Trajectory: ρ=%.4f, P@12=%.3f, P@10=%.3f, P@5=%.3f, P@3=%.3f, P@1=%.3f\n\n', ...
        best_single_spearman, ...
        traj_results_matrix(best_single_idx, 2), ...
        traj_results_matrix(best_single_idx, 3), ...
        traj_results_matrix(best_single_idx, 4), ...
        traj_results_matrix(best_single_idx, 5), ...
        traj_results_matrix(best_single_idx, 6));
    
    % DTW Mode comparison
    fprintf('4. DTW MODE COMPARISON:\n');
    joint_mask = strcmp(dtw_modes_traj, 'joint_states');
    pos_mask = strcmp(dtw_modes_traj, 'position');
    
    avg_joint = mean(traj_results_matrix(joint_mask, 1));
    avg_pos = mean(traj_results_matrix(pos_mask, 1));
    
    fprintf('   Joint States: Average ρ=%.4f\n', avg_joint);
    fprintf('   Position:     Average ρ=%.4f\n', avg_pos);
    
    if avg_joint > avg_pos
        fprintf('   → Joint States mode performs better on average\n\n');
    else
        fprintf('   → Position mode performs better on average\n\n');
    end
    
    fprintf('========================================\n');
    fprintf('ANALYSIS COMPLETE\n');
    fprintf('========================================\n\n');
    
end
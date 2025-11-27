function results = runExperiment(config)
    % RUNEXPERIMENT - Universal wrapper for dtw_baseline.m
    fprintf('\n>>> Running Experiment: %s <<<\n', config.exp_name);

    % ====================================================================
    % STEP 1: Set workspace variables (will be used by dtw_baseline.m)
    % ====================================================================
    
    % Query
    query_bahn_id = config.query_bahn_id;
    
    % Database
    use_database_sampling = true;
    database_sample_size = config.database_sample_size;
    random_seed = config.random_seed;
    
    % DTW Mode
    if isfield(config, 'dtw_mode')
        dtw_mode = config.dtw_mode;
    else
        dtw_mode = 'joint_states';  % Default
    end
    
    % DTW Config
    normalize_dtw = true;
    use_rotation_alignment = false;
    cdtw_window = 0.10;
    
    % Lower Bounds
    lb_kim_keep_ratio = 0.40;
    lb_keogh_candidates = 100;
    
    % Embedding Architecture
    use_multi_scale = true;
    n_coarse = config.n_coarse;
    n_medium = 0;
    n_fine = config.n_fine;
    norm_strategy = 'max_extent';
    
    % RRF Config
    rrf_k = 60;
    
    % ⭐ Weights - Keep as array for dtw_baseline.m
    weights = config.weights;  % [pos, joint, orient, vel, meta]
    
    % Output
    if isfield(config, 'top_k_trajectories')
        top_k_trajectories = config.top_k_trajectories;
    else
        top_k_trajectories = 50;
    end
    
    % Database
    schema = 'bewegungsdaten';
    db_name = 'robotervermessung';
    
    % ====================================================================
    % STEP 2: Run main pipeline
    % ====================================================================
    
    fprintf('  Config: Mode=%s, N=%d, Weights=[%.1f,%.1f,%.1f,%.1f,%.1f], Dims=%d+%d\n', ...
        dtw_mode, database_sample_size, ...
        weights(1), weights(2), weights(3), weights(4), weights(5), ...
        n_coarse, n_fine);
    
    % Run dtw_baseline.m
    run('dtw_baseline.m');
    
    % ====================================================================
    % STEP 3: Collect results
    % ====================================================================
    
    results = struct();
    
    % Config
    results.exp_name = config.exp_name;
    results.dtw_mode = dtw_mode;
    results.query_id = query_bahn_id;
    results.database_size = num_candidates;
    results.sample_size = database_sample_size;
    results.random_seed = random_seed;
    
    % Weights
    results.weight_pos = weights(1);
    results.weight_joint = weights(2);
    results.weight_orient = weights(3);
    results.weight_vel = weights(4);
    results.weight_meta = weights(5);
    
    % Dimensions
    results.n_coarse = n_coarse;
    results.n_fine = n_fine;
    results.total_dims = (n_coarse + n_fine) * 3;
    
    % Metrics - Trajectory Level
    results.spearman = rho_spearman;
    results.p_at_k = prec_k;
    results.p_at_10 = prec_10;
    results.p_at_5 = prec_5;
    results.p_at_3 = prec_3;
    results.p_at_1 = prec_1;
    
    % ⭐ Metrics - Segment Level (ALL metrics if available)
    if exist('seg_rho_all', 'var') && ~isempty(seg_rho_all)
        results.seg_spearman = mean(seg_rho_all);
        
        % Check which P@k metrics exist for segments
        if exist('seg_prec_k_all', 'var') && ~isempty(seg_prec_k_all)
            results.seg_p_at_k = mean(seg_prec_k_all);
        else
            results.seg_p_at_k = NaN;
        end
        
        if exist('seg_prec_10_all', 'var') && ~isempty(seg_prec_10_all)
            results.seg_p_at_10 = mean(seg_prec_10_all);
        else
            results.seg_p_at_10 = NaN;
        end
        
        if exist('seg_prec_5_all', 'var') && ~isempty(seg_prec_5_all)
            results.seg_p_at_5 = mean(seg_prec_5_all);
        else
            results.seg_p_at_5 = NaN;
        end
        
        if exist('seg_prec_3_all', 'var') && ~isempty(seg_prec_3_all)
            results.seg_p_at_3 = mean(seg_prec_3_all);
        else
            results.seg_p_at_3 = NaN;
        end
        
        if exist('seg_prec_1_all', 'var') && ~isempty(seg_prec_1_all)
            results.seg_p_at_1 = mean(seg_prec_1_all);
        else
            results.seg_p_at_1 = NaN;
        end
    else
        % No segment data available
        results.seg_spearman = NaN;
        results.seg_p_at_k = NaN;
        results.seg_p_at_10 = NaN;
        results.seg_p_at_5 = NaN;
        results.seg_p_at_3 = NaN;
        results.seg_p_at_1 = NaN;
    end
    
    % Timing
    results.dtw_time = dtw_trajectory_time;
    
    fprintf('  >>> Results: Mode=%s, ρ=%.4f, P@10=%.3f, P@1=%.0f, DTW=%.2fs <<<\n\n', ...
        dtw_mode, results.spearman, results.p_at_10, results.p_at_1, results.dtw_time);
end
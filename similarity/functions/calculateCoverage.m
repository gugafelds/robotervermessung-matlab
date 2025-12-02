function coverage_results = calculateCoverage(dtw_ranking, embedding_ranking, k_values, id_field)
% CALCULATECOVERAGE - Calculate coverage metrics for embedding vs DTW rankings
%
%   Analyzes how many embedding top-X results are needed to cover all
%   DTW top-K results. This quantifies the reliability and practical
%   deployment requirements of embedding-based retrieval.
%
%   INPUTS:
%       dtw_ranking        - Table with DTW rankings (sorted by dtw_distance)
%                            Must contain ID field (e.g., 'bahn_id' or 'segment_id')
%       embedding_ranking  - Table with embedding rankings (sorted by rrf_score)
%                            Must contain same ID field
%       k_values           - Array of K values to analyze (e.g., [10, 20, 50])
%       id_field           - Name of ID field (e.g., 'bahn_id' or 'segment_id')
%
%   OUTPUT:
%       coverage_results - Struct with fields:
%           .k_values              - Input K values
%           .coverage_points       - Coverage point for each K (min X for 100% recall)
%           .expansion_ratios      - Expansion ratio for each K (coverage_point / K)
%           .recall_at_k           - Recall when using embedding top-K
%           .worst_ranks           - Worst embedding rank for each DTW top-K
%           .mean_rank_shift       - Average rank difference
%           .median_rank_shift     - Median rank difference
%
%   EXAMPLE:
%       coverage = calculateCoverage(trajectory_table, fused_ranking, [10, 50], 'bahn_id');
%       fprintf('Need top-%d embeddings for all DTW top-50\n', coverage.coverage_points(2));
%       fprintf('Expansion ratio: %.2fx\n', coverage.expansion_ratios(2));
%
%   INTERPRETATION:
%       Expansion Ratio:
%           1.0     - Perfect! Same ranking
%           < 1.5   - Excellent reliability
%           < 2.0   - Good reliability
%           < 3.0   - Moderate reliability
%           >= 3.0  - Poor reliability
%
%   Author: Gustavo Barros
%   Date: 02.12.2025

%% ========================================================================
%  VALIDATION
%  ========================================================================

% Check inputs
if nargin < 4
    error('calculateCoverage requires 4 inputs: dtw_ranking, embedding_ranking, k_values, id_field');
end

% Validate tables
if ~istable(dtw_ranking) || ~istable(embedding_ranking)
    error('dtw_ranking and embedding_ranking must be tables');
end

% Check if ID field exists
if ~ismember(id_field, dtw_ranking.Properties.VariableNames)
    error('ID field "%s" not found in dtw_ranking', id_field);
end

if ~ismember(id_field, embedding_ranking.Properties.VariableNames)
    error('ID field "%s" not found in embedding_ranking', id_field);
end

% Validate k_values
if ~isnumeric(k_values) || any(k_values <= 0)
    error('k_values must be positive numbers');
end

k_values = sort(k_values);  % Ensure ascending order
num_k = length(k_values);

% Get total number of candidates
num_dtw_candidates = height(dtw_ranking);
num_emb_candidates = height(embedding_ranking);

%% ========================================================================
%  CREATE LOOKUP MAP FOR EMBEDDING RANKS
%  ========================================================================

% Create a map: ID → Embedding Rank
embedding_rank_map = containers.Map('KeyType', 'char', 'ValueType', 'double');

for i = 1:num_emb_candidates
    id = embedding_ranking.(id_field){i};
    embedding_rank_map(id) = i;  % Store embedding rank
end

%% ========================================================================
%  ANALYZE COVERAGE FOR EACH K
%  ========================================================================

coverage_points = zeros(num_k, 1);
expansion_ratios = zeros(num_k, 1);
recall_at_k = zeros(num_k, 1);
worst_ranks = zeros(num_k, 1);
mean_rank_shifts = zeros(num_k, 1);
median_rank_shifts = zeros(num_k, 1);

for k_idx = 1:num_k
    K = k_values(k_idx);
    
    % Check if K is valid
    if K > num_dtw_candidates
        warning('K=%d exceeds DTW ranking size (%d). Using max size.', K, num_dtw_candidates);
        K = num_dtw_candidates;
    end
    
    % Get DTW top-K IDs
    dtw_top_k_ids = dtw_ranking.(id_field)(1:K);
    
    % Find embedding ranks for these DTW top-K
    embedding_ranks_of_dtw_top_k = zeros(K, 1);
    
    for i = 1:K
        dtw_id = dtw_top_k_ids{i};
        
        if isKey(embedding_rank_map, dtw_id)
            embedding_ranks_of_dtw_top_k(i) = embedding_rank_map(dtw_id);
        else
            % ID not found in embeddings (should not happen, but handle gracefully)
            warning('DTW top-%d ID "%s" not found in embedding ranking', i, dtw_id);
            embedding_ranks_of_dtw_top_k(i) = inf;  % Mark as not found
        end
    end
    
    % Remove any inf values (not found)
    valid_ranks = embedding_ranks_of_dtw_top_k(~isinf(embedding_ranks_of_dtw_top_k));
    
    if isempty(valid_ranks)
        warning('No valid embedding ranks found for DTW top-%d', K);
        coverage_points(k_idx) = inf;
        expansion_ratios(k_idx) = inf;
        recall_at_k(k_idx) = 0;
        worst_ranks(k_idx) = inf;
        mean_rank_shifts(k_idx) = inf;
        median_rank_shifts(k_idx) = inf;
        continue;
    end
    
    % ====================================================================
    % METRIC 1: Coverage Point (worst case rank)
    % ====================================================================
    % "How many embedding top-X do I need to get all DTW top-K?"
    coverage_point = max(valid_ranks);
    coverage_points(k_idx) = coverage_point;
    
    % ====================================================================
    % METRIC 2: Expansion Ratio
    % ====================================================================
    % "By what factor do I need to expand?"
    expansion_ratio = coverage_point / K;
    expansion_ratios(k_idx) = expansion_ratio;
    
    % ====================================================================
    % METRIC 3: Recall@K
    % ====================================================================
    % "How many DTW top-K are in embedding top-K?"
    num_in_top_k = sum(valid_ranks <= K);
    recall = num_in_top_k / K;
    recall_at_k(k_idx) = recall;
    
    % ====================================================================
    % METRIC 4: Worst Rank
    % ====================================================================
    worst_ranks(k_idx) = coverage_point;  % Same as coverage point
    
    % ====================================================================
    % METRIC 5: Rank Shifts
    % ====================================================================
    % "How much do ranks shift on average?"
    % DTW rank i should ideally be at embedding rank i
    dtw_ranks = (1:length(valid_ranks))';
    rank_shifts = abs(valid_ranks - dtw_ranks);
    
    mean_rank_shifts(k_idx) = mean(rank_shifts);
    median_rank_shifts(k_idx) = median(rank_shifts);
end

%% ========================================================================
%  BUILD OUTPUT STRUCT
%  ========================================================================

coverage_results = struct();

% Input parameters
coverage_results.k_values = k_values;

% Core metrics
coverage_results.coverage_points = coverage_points;
coverage_results.expansion_ratios = expansion_ratios;

% Additional metrics
coverage_results.recall_at_k = recall_at_k;
coverage_results.worst_ranks = worst_ranks;
coverage_results.mean_rank_shift = mean_rank_shifts;
coverage_results.median_rank_shift = median_rank_shifts;

% Summary statistics
coverage_results.mean_expansion_ratio = mean(expansion_ratios(~isinf(expansion_ratios)));
coverage_results.max_expansion_ratio = max(expansion_ratios(~isinf(expansion_ratios)));
coverage_results.min_expansion_ratio = min(expansion_ratios(~isinf(expansion_ratios)));

%% ========================================================================
%  DISPLAY SUMMARY (OPTIONAL)
%  ========================================================================

if nargout == 0
    % If no output requested, display results
    fprintf('\n=== Coverage Analysis Summary ===\n');
    fprintf('ID Field: %s\n', id_field);
    fprintf('DTW Candidates: %d\n', num_dtw_candidates);
    fprintf('Embedding Candidates: %d\n\n', num_emb_candidates);
    
    fprintf('%-10s | %-15s | %-15s | %-12s | %-12s\n', ...
        'K', 'Coverage Point', 'Expansion', 'Recall@K', 'Worst Rank');
    fprintf('%s\n', repmat('-', 75, 1));
    
    for k_idx = 1:num_k
        fprintf('%-10d | %-15d | %-15.3f | %-12.2f | %-12d\n', ...
            k_values(k_idx), ...
            coverage_points(k_idx), ...
            expansion_ratios(k_idx), ...
            recall_at_k(k_idx), ...
            worst_ranks(k_idx));
    end
    
    fprintf('\n--- Interpretation ---\n');
    fprintf('Coverage Point: Minimum embedding top-X needed for all DTW top-K\n');
    fprintf('Expansion Ratio: Coverage Point / K\n');
    fprintf('  < 1.5:  Excellent reliability ✅\n');
    fprintf('  < 2.0:  Good reliability ✅\n');
    fprintf('  < 3.0:  Moderate reliability ⚠️\n');
    fprintf('  >= 3.0: Poor reliability ❌\n');
    fprintf('Recall@K: Fraction of DTW top-K found in embedding top-K\n');
    fprintf('Worst Rank: Worst embedding rank among DTW top-K\n\n');
end

end
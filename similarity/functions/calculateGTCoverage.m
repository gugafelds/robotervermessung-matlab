function gt_coverage = calculateGTCoverage(embedding_ranking, ground_truth_ids, k_values, id_field)
% CALCULATEGTCOVERAGE - Calculate GT-specific coverage metrics for embedding ranking
%
%   Analyzes how well the embedding ranking captures ground truth (GT)
%   trajectories. Unlike calculateCoverage() which compares against DTW
%   top-K, this function measures coverage of KNOWN similar trajectories.
%
%   INPUTS:
%       embedding_ranking  - Table with embedding rankings (sorted by rrf_score)
%                            Must contain ID field (e.g., 'bahn_id' or 'segment_id')
%       ground_truth_ids   - Cell array of GT trajectory IDs that are known
%                            to be similar to the query (e.g., from same recording)
%       k_values           - Array of K values to analyze (e.g., [10, 50])
%       id_field           - Name of ID field (e.g., 'bahn_id' or 'segment_id')
%
%   OUTPUT:
%       gt_coverage - Struct with fields:
%           .k_values              - Input K values
%           .recall_at_k           - Recall@K: fraction of GT found in top-K
%           .expansion_ratios      - Expansion ratio: coverage_point / num_gt
%           .coverage_points       - Min embedding top-X needed for all GT
%           .num_gt                - Number of ground truth trajectories
%           .gt_ranks              - Embedding rank of each GT
%           .mean_gt_rank          - Average embedding rank of GT
%           .median_gt_rank        - Median embedding rank of GT
%
%   EXAMPLE:
%       % Query has 9 GT trajectories from same recording
%       query_gt_ids = {'1763567278', '1763567279', ..., '1763567286'};
%       gt_cov = calculateGTCoverage(fused_ranking, query_gt_ids, [10, 50], 'bahn_id');
%       
%       fprintf('Recall@50(GT) = %.2f\n', gt_cov.recall_at_k(2));
%       fprintf('Need top-%d embeddings for all GT\n', gt_cov.coverage_points(2));
%
%   INTERPRETATION:
%       Recall@K:
%           1.0     - Perfect! All GT found in top-K
%           >= 0.9  - Excellent (max 1 GT lost)
%           >= 0.8  - Good (1-2 GT lost)
%           < 0.8   - Poor (3+ GT lost)
%
%       Expansion Ratio (Coverage Point / num_gt):
%           < 2x    - Excellent efficiency
%           < 5x    - Good efficiency
%           < 10x   - Acceptable efficiency
%           >= 10x  - Poor efficiency
%
%   DIFFERENCE FROM calculateCoverage():
%       - calculateCoverage: "How many Emb top-X for DTW top-K?"
%         → Measures consistency between two ranking methods
%       
%       - calculateGTCoverage: "How many Emb top-X for ALL known GT?"
%         → Measures recall of known-similar trajectories
%
%   Author: Gustavo Barros
%   Date: 03.12.2025

%% ========================================================================
%  VALIDATION
%  ========================================================================

% Check inputs
if nargin < 4
    error('calculateGTCoverage requires 4 inputs: embedding_ranking, ground_truth_ids, k_values, id_field');
end

% Validate tables
if ~istable(embedding_ranking)
    error('embedding_ranking must be a table');
end

% Validate ground_truth_ids
if ~iscell(ground_truth_ids) || isempty(ground_truth_ids)
    error('ground_truth_ids must be a non-empty cell array');
end

% Check if ID field exists
if ~ismember(id_field, embedding_ranking.Properties.VariableNames)
    error('ID field "%s" not found in embedding_ranking', id_field);
end

% Validate k_values
if ~isnumeric(k_values) || any(k_values <= 0)
    error('k_values must be positive numbers');
end

k_values = sort(k_values);  % Ensure ascending order
num_k = length(k_values);
num_gt = length(ground_truth_ids);

% Get total number of candidates in embedding ranking
num_emb_candidates = height(embedding_ranking);

%% ========================================================================
%  FIND EMBEDDING RANKS FOR ALL GT
%  ========================================================================

% Create a map: ID → Embedding Rank
embedding_rank_map = containers.Map('KeyType', 'char', 'ValueType', 'double');

for i = 1:num_emb_candidates
    id = embedding_ranking.(id_field){i};
    embedding_rank_map(id) = i;  % Store embedding rank
end

% Find embedding ranks for all GT
gt_ranks = zeros(num_gt, 1);
gt_found = false(num_gt, 1);

for gt_idx = 1:num_gt
    gt_id = ground_truth_ids{gt_idx};
    
    if isKey(embedding_rank_map, gt_id)
        gt_ranks(gt_idx) = embedding_rank_map(gt_id);
        gt_found(gt_idx) = true;
    else
        % GT not found in embedding ranking (should rarely happen)
        warning('GT "%s" not found in embedding ranking', gt_id);
        gt_ranks(gt_idx) = inf;  % Mark as not found
        gt_found(gt_idx) = false;
    end
end

% Remove any inf values (not found)
valid_gt_ranks = gt_ranks(~isinf(gt_ranks));
num_gt_found = sum(gt_found);

if num_gt_found == 0
    warning('No GT trajectories found in embedding ranking!');
    % Return empty result
    gt_coverage = struct();
    gt_coverage.k_values = k_values;
    gt_coverage.recall_at_k = zeros(num_k, 1);
    gt_coverage.expansion_ratios = inf(num_k, 1);
    gt_coverage.coverage_points = inf(num_k, 1);
    gt_coverage.num_gt = num_gt;
    gt_coverage.num_gt_found = 0;
    gt_coverage.gt_ranks = gt_ranks;
    gt_coverage.mean_gt_rank = inf;
    gt_coverage.median_gt_rank = inf;
    return;
end

%% ========================================================================
%  CALCULATE METRICS FOR EACH K
%  ========================================================================

recall_at_k = zeros(num_k, 1);
coverage_points = zeros(num_k, 1);
expansion_ratios = zeros(num_k, 1);

for k_idx = 1:num_k
    K = k_values(k_idx);
    
    % Check if K is valid
    if K > num_emb_candidates
        warning('K=%d exceeds embedding ranking size (%d). Using max size.', K, num_emb_candidates);
        K = num_emb_candidates;
    end
    
    % ====================================================================
    % METRIC 1: Recall@K
    % ====================================================================
    % "How many GT are found in embedding top-K?"
    num_gt_in_top_k = sum(valid_gt_ranks <= K);
    recall_at_k(k_idx) = num_gt_in_top_k / num_gt;
    
    % ====================================================================
    % METRIC 2: Coverage Point
    % ====================================================================
    % "What's the minimum embedding top-X needed for all GT?"
    coverage_points(k_idx) = max(valid_gt_ranks);
    
    % ====================================================================
    % METRIC 3: Expansion Ratio
    % ====================================================================
    % "By what factor do we need to expand from ideal?"
    % Ideal: All num_gt would be in ranks 1:num_gt (expansion = 1x)
    expansion_ratios(k_idx) = coverage_points(k_idx) / num_gt;
end

%% ========================================================================
%  BUILD OUTPUT STRUCT
%  ========================================================================

gt_coverage = struct();

% Input parameters
gt_coverage.k_values = k_values;
gt_coverage.num_gt = num_gt;
gt_coverage.num_gt_found = num_gt_found;

% Core metrics
gt_coverage.recall_at_k = recall_at_k;
gt_coverage.coverage_points = coverage_points;
gt_coverage.expansion_ratios = expansion_ratios;

% GT rank statistics
gt_coverage.gt_ranks = gt_ranks;
gt_coverage.mean_gt_rank = mean(valid_gt_ranks);
gt_coverage.median_gt_rank = median(valid_gt_ranks);
gt_coverage.min_gt_rank = min(valid_gt_ranks);
gt_coverage.max_gt_rank = max(valid_gt_ranks);

% Overall coverage point (for all GT)
gt_coverage.overall_coverage_point = max(valid_gt_ranks);
gt_coverage.overall_expansion_ratio = gt_coverage.overall_coverage_point / num_gt;

%% ========================================================================
%  DISPLAY SUMMARY (OPTIONAL)
%  ========================================================================

if nargout == 0
    % If no output requested, display results
    fprintf('\n=== Ground Truth Coverage Analysis ===\n');
    fprintf('ID Field: %s\n', id_field);
    fprintf('Total GT: %d\n', num_gt);
    fprintf('GT Found in Ranking: %d\n', num_gt_found);
    fprintf('Embedding Candidates: %d\n\n', num_emb_candidates);
    
    fprintf('--- GT Rank Statistics ---\n');
    fprintf('Mean GT Rank: %.1f\n', gt_coverage.mean_gt_rank);
    fprintf('Median GT Rank: %.1f\n', gt_coverage.median_gt_rank);
    fprintf('Min GT Rank: %d\n', gt_coverage.min_gt_rank);
    fprintf('Max GT Rank: %d\n', gt_coverage.max_gt_rank);
    fprintf('Overall Coverage Point: %d\n', gt_coverage.overall_coverage_point);
    fprintf('Overall Expansion Ratio: %.1fx\n\n', gt_coverage.overall_expansion_ratio);
    
    fprintf('%-10s | %-12s | %-15s | %-15s\n', ...
        'K', 'Recall@K', 'Coverage Point', 'Expansion');
    fprintf('%s\n', repmat('-', 60, 1));
    
    for k_idx = 1:num_k
        fprintf('%-10d | %-12.2f | %-15d | %-15.2fx\n', ...
            k_values(k_idx), ...
            recall_at_k(k_idx), ...
            coverage_points(k_idx), ...
            expansion_ratios(k_idx));
    end
    
    fprintf('\n--- Interpretation ---\n');
    fprintf('Recall@K: Fraction of GT found in embedding top-K\n');
    fprintf('  1.0:    Perfect! All GT found ✅\n');
    fprintf('  >= 0.9: Excellent (max 1 GT lost) ✅\n');
    fprintf('  >= 0.8: Good (1-2 GT lost) ⚠️\n');
    fprintf('  < 0.8:  Poor (3+ GT lost) ❌\n\n');
    
    fprintf('Expansion Ratio: Coverage Point / num_gt\n');
    fprintf('  < 2x:   Excellent efficiency ✅\n');
    fprintf('  < 5x:   Good efficiency ✅\n');
    fprintf('  < 10x:  Acceptable efficiency ⚠️\n');
    fprintf('  >= 10x: Poor efficiency ❌\n\n');
    
    fprintf('Coverage Point: Min embedding top-X needed for all GT\n');
    fprintf('  Practical use: Set pgvector_top_k >= Coverage Point\n\n');
end

end
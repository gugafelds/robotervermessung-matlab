function avg_coverage = averageGTCoverage(coverage_cells)
% AVERAGEGTCOVERAGE - Average GT coverage metrics over multiple segments
%
%   Computes the average GT coverage metrics across multiple segments.
%   This is useful for segment-level analysis where each query segment
%   has its own set of GT segments (constructed from GT trajectories).
%
%   INPUTS:
%       coverage_cells - Cell array of GT coverage structs from 
%                        calculateGTCoverage(). Each cell corresponds
%                        to one query segment.
%                        Empty cells are automatically filtered out.
%
%   OUTPUT:
%       avg_coverage - Averaged GT coverage struct with fields:
%           .k_values              - K values analyzed
%           .recall_at_k           - Average recall@K across segments
%           .expansion_ratios      - Average expansion ratio across segments
%           .coverage_points       - Average coverage point across segments
%           .mean_gt_rank          - Average mean GT rank across segments
%           .median_gt_rank        - Average median GT rank across segments
%           .num_segments_analyzed - Number of non-empty segments
%
%   EXAMPLE:
%       % Calculate GT coverage for each segment
%       gt_coverage_segments = cell(3, 1);
%       for seg_idx = 1:3
%           gt_coverage_segments{seg_idx} = calculateGTCoverage(...);
%       end
%       
%       % Average across segments
%       avg_cov = averageGTCoverage(gt_coverage_segments);
%       fprintf('Average Recall@50(GT) = %.2f\n', avg_cov.recall_at_k(2));
%
%   HANDLING MISSING DATA:
%       - Empty cells are ignored
%       - NaN values are excluded from averaging
%       - Inf values are excluded from averaging
%       - If no valid data for a K, result is NaN
%
%   Author: Gustavo Barros
%   Date: 03.12.2025

%% ========================================================================
%  VALIDATION
%  ========================================================================

% Check input
if nargin < 1
    error('averageGTCoverage requires 1 input: coverage_cells');
end

if ~iscell(coverage_cells)
    error('coverage_cells must be a cell array');
end

% Filter out empty cells
valid_coverages = coverage_cells(~cellfun(@isempty, coverage_cells));

if isempty(valid_coverages)
    warning('No valid coverage data to average');
    avg_coverage = [];
    return;
end

num_valid_segments = length(valid_coverages);

%% ========================================================================
%  GET K VALUES FROM FIRST VALID COVERAGE
%  ========================================================================

k_values = valid_coverages{1}.k_values;
num_k = length(k_values);

% Validate that all segments have same k_values
for i = 2:num_valid_segments
    if ~isequal(valid_coverages{i}.k_values, k_values)
        warning('Segment %d has different k_values than segment 1. Using segment 1 k_values.', i);
    end
end

%% ========================================================================
%  INITIALIZE ACCUMULATORS
%  ========================================================================

% For each K value
recall_sum = zeros(num_k, 1);
expansion_sum = zeros(num_k, 1);
coverage_sum = zeros(num_k, 1);
counts = zeros(num_k, 1);  % How many valid values for each K

% For rank statistics
mean_gt_rank_sum = 0;
median_gt_rank_sum = 0;
min_gt_rank_sum = 0;
max_gt_rank_sum = 0;
overall_coverage_sum = 0;
overall_expansion_sum = 0;
rank_stats_count = 0;

%% ========================================================================
%  ACCUMULATE OVER ALL SEGMENTS
%  ========================================================================

for seg_idx = 1:num_valid_segments
    cov = valid_coverages{seg_idx};
    
    % ====================================================================
    % Accumulate K-specific metrics
    % ====================================================================
    for k_idx = 1:num_k
        % Only accumulate if value is valid (not NaN, not Inf)
        if isfield(cov, 'recall_at_k') && ...
           ~isnan(cov.recall_at_k(k_idx)) && ...
           isfinite(cov.recall_at_k(k_idx))
            
            recall_sum(k_idx) = recall_sum(k_idx) + cov.recall_at_k(k_idx);
            counts(k_idx) = counts(k_idx) + 1;
        end
        
        if isfield(cov, 'expansion_ratios') && ...
           ~isnan(cov.expansion_ratios(k_idx)) && ...
           isfinite(cov.expansion_ratios(k_idx))
            
            expansion_sum(k_idx) = expansion_sum(k_idx) + cov.expansion_ratios(k_idx);
        end
        
        if isfield(cov, 'coverage_points') && ...
           ~isnan(cov.coverage_points(k_idx)) && ...
           isfinite(cov.coverage_points(k_idx))
            
            coverage_sum(k_idx) = coverage_sum(k_idx) + cov.coverage_points(k_idx);
        end
    end
    
    % ====================================================================
    % Accumulate rank statistics
    % ====================================================================
    if isfield(cov, 'mean_gt_rank') && ...
       ~isnan(cov.mean_gt_rank) && ...
       isfinite(cov.mean_gt_rank)
        
        mean_gt_rank_sum = mean_gt_rank_sum + cov.mean_gt_rank;
        rank_stats_count = rank_stats_count + 1;
    end
    
    if isfield(cov, 'median_gt_rank') && ...
       ~isnan(cov.median_gt_rank) && ...
       isfinite(cov.median_gt_rank)
        
        median_gt_rank_sum = median_gt_rank_sum + cov.median_gt_rank;
    end
    
    if isfield(cov, 'min_gt_rank') && ...
       ~isnan(cov.min_gt_rank) && ...
       isfinite(cov.min_gt_rank)
        
        min_gt_rank_sum = min_gt_rank_sum + cov.min_gt_rank;
    end
    
    if isfield(cov, 'max_gt_rank') && ...
       ~isnan(cov.max_gt_rank) && ...
       isfinite(cov.max_gt_rank)
        
        max_gt_rank_sum = max_gt_rank_sum + cov.max_gt_rank;
    end
    
    if isfield(cov, 'overall_coverage_point') && ...
       ~isnan(cov.overall_coverage_point) && ...
       isfinite(cov.overall_coverage_point)
        
        overall_coverage_sum = overall_coverage_sum + cov.overall_coverage_point;
    end
    
    if isfield(cov, 'overall_expansion_ratio') && ...
       ~isnan(cov.overall_expansion_ratio) && ...
       isfinite(cov.overall_expansion_ratio)
        
        overall_expansion_sum = overall_expansion_sum + cov.overall_expansion_ratio;
    end
end

%% ========================================================================
%  COMPUTE AVERAGES
%  ========================================================================

avg_coverage = struct();

% Input parameters
avg_coverage.k_values = k_values;
avg_coverage.num_segments_analyzed = num_valid_segments;

% K-specific metrics
avg_coverage.recall_at_k = recall_sum ./ max(counts, 1);
avg_coverage.expansion_ratios = expansion_sum ./ max(counts, 1);
avg_coverage.coverage_points = coverage_sum ./ max(counts, 1);

% Set to NaN where no valid data
avg_coverage.recall_at_k(counts == 0) = NaN;
avg_coverage.expansion_ratios(counts == 0) = NaN;
avg_coverage.coverage_points(counts == 0) = NaN;

% Rank statistics
if rank_stats_count > 0
    avg_coverage.mean_gt_rank = mean_gt_rank_sum / rank_stats_count;
    avg_coverage.median_gt_rank = median_gt_rank_sum / rank_stats_count;
    avg_coverage.min_gt_rank = min_gt_rank_sum / rank_stats_count;
    avg_coverage.max_gt_rank = max_gt_rank_sum / rank_stats_count;
    avg_coverage.overall_coverage_point = overall_coverage_sum / rank_stats_count;
    avg_coverage.overall_expansion_ratio = overall_expansion_sum / rank_stats_count;
else
    avg_coverage.mean_gt_rank = NaN;
    avg_coverage.median_gt_rank = NaN;
    avg_coverage.min_gt_rank = NaN;
    avg_coverage.max_gt_rank = NaN;
    avg_coverage.overall_coverage_point = NaN;
    avg_coverage.overall_expansion_ratio = NaN;
end

%% ========================================================================
%  DISPLAY SUMMARY (OPTIONAL)
%  ========================================================================

if nargout == 0
    % If no output requested, display results
    fprintf('\n=== Average GT Coverage Across Segments ===\n');
    fprintf('Number of segments analyzed: %d\n\n', num_valid_segments);
    
    fprintf('--- Average GT Rank Statistics ---\n');
    fprintf('Mean GT Rank: %.1f\n', avg_coverage.mean_gt_rank);
    fprintf('Median GT Rank: %.1f\n', avg_coverage.median_gt_rank);
    fprintf('Min GT Rank: %.1f\n', avg_coverage.min_gt_rank);
    fprintf('Max GT Rank: %.1f\n', avg_coverage.max_gt_rank);
    fprintf('Overall Coverage Point: %.1f\n', avg_coverage.overall_coverage_point);
    fprintf('Overall Expansion Ratio: %.1fx\n\n', avg_coverage.overall_expansion_ratio);
    
    fprintf('%-10s | %-12s | %-15s | %-15s\n', ...
        'K', 'Recall@K', 'Coverage Point', 'Expansion');
    fprintf('%s\n', repmat('-', 60, 1));
    
    for k_idx = 1:num_k
        if ~isnan(avg_coverage.recall_at_k(k_idx))
            fprintf('%-10d | %-12.2f | %-15.1f | %-15.2fx\n', ...
                k_values(k_idx), ...
                avg_coverage.recall_at_k(k_idx), ...
                avg_coverage.coverage_points(k_idx), ...
                avg_coverage.expansion_ratios(k_idx));
        else
            fprintf('%-10d | %-12s | %-15s | %-15s\n', ...
                k_values(k_idx), 'N/A', 'N/A', 'N/A');
        end
    end
    
    fprintf('\n--- Notes ---\n');
    fprintf('These are AVERAGES across %d segments\n', num_valid_segments);
    fprintf('Segments with missing/invalid data were excluded\n\n');
end

end
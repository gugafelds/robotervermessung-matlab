function fused_ranking = fuseRankingsRRF(rankings, weights, k, id_column)
    % fuseRankingsRRF - Fuse multiple rankings using Reciprocal Rank Fusion (RRF)
    %
    % Syntax: fused_ranking = fuseRankingsRRF(rankings, weights, k, id_column)
    %
    % Inputs:
    %   rankings - Struct with fields for each modality (e.g., rankings.position, rankings.joint)
    %   weights - Array [pos, joint, orient, vel, meta] OR Struct with weight fields
    %   k - RRF constant (default: 60)
    %   id_column - Name of ID column (default: 'bahn_id')
    %
    % Output:
    %   fused_ranking - Table sorted by RRF score
    
    if nargin < 4
    id_column = 'bahn_id';
end
if nargin < 3
    k = 60;
end

% ========================================================================
% STEP 1: Convert weights Array → Struct if needed
% ========================================================================
if ~isstruct(weights)
    w = weights;
    weights = struct();
    weights.position = w(1);
    weights.joint = w(2);
    weights.orientation = w(3);
    weights.velocity = w(4);
    weights.metadata = w(5);
end

% ========================================================================
% STEP 2: Get all ranking fields (skip empty rankings)
% ========================================================================
ranking_fields = fieldnames(rankings);

% Filter out empty rankings
valid_fields = {};
for i = 1:length(ranking_fields)
    field = ranking_fields{i};
    if ~isempty(rankings.(field))
        valid_fields{end+1} = field;
    end
end

if isempty(valid_fields)
    error('No valid rankings provided');
end

% ========================================================================
% STEP 3: Collect all unique IDs
% ========================================================================
all_ids = {};
for i = 1:length(valid_fields)
    field = valid_fields{i};
    ranking_table = rankings.(field);
    all_ids = [all_ids; ranking_table.(id_column)];
end
unique_ids = unique(all_ids);

% ========================================================================
% STEP 4: Compute RRF scores
% ========================================================================
rrf_scores = zeros(length(unique_ids), 1);

for i = 1:length(valid_fields)
    field = valid_fields{i};
    ranking_table = rankings.(field);
    
    % Get weight for this modality
    if isfield(weights, field)
        weight = weights.(field);
    else
        weight = 1.0;  % Default weight
    end
    
    % Skip if weight is 0
    if weight == 0
        continue;
    end
    
    % For each ID, find its rank and compute RRF contribution
    for j = 1:length(unique_ids)
        id = unique_ids{j};
        
        % Find rank of this ID in current ranking
        idx = find(strcmp(ranking_table.(id_column), id));
        
        if ~isempty(idx)
            rank = idx(1);  % Take first occurrence
            rrf_contribution = weight / (k + rank);
            rrf_scores(j) = rrf_scores(j) + rrf_contribution;
        end
        % If ID not in this ranking, contribution = 0 (skip)
    end
end

% ========================================================================
% STEP 5: Create fused ranking table
% ========================================================================
fused_ranking = table(unique_ids, rrf_scores, ...
    'VariableNames', {id_column, 'rrf_score'});

% Sort by RRF score (descending)
fused_ranking = sortrows(fused_ranking, 'rrf_score', 'descend');

% ========================================================================
% STEP 6: Add metadata from first valid ranking
% ========================================================================
first_ranking = rankings.(valid_fields{1});

% Get all column names except id_column and distance
meta_columns = first_ranking.Properties.VariableNames;
meta_columns = meta_columns(~ismember(meta_columns, {id_column, 'distance', 'cosine_distance'}));

% Add metadata columns
for col_idx = 1:length(meta_columns)
    col_name = meta_columns{col_idx};
    
    % Initialize column with appropriate type
    sample_value = first_ranking.(col_name)(1);
    
    if isnumeric(sample_value)
        % Numeric column
        new_col = zeros(height(fused_ranking), 1);
        
    elseif isstring(sample_value)
        % String column → convert to cell array of strings
        new_col = cell(height(fused_ranking), 1);
        for i = 1:height(fused_ranking)
            new_col{i} = "";
        end
        
    elseif iscell(sample_value)
        % Cell column
        new_col = cell(height(fused_ranking), 1);
        
    elseif ischar(sample_value)
        % Char column
        new_col = cell(height(fused_ranking), 1);
        
    else
        % Unknown type → skip this column
        warning('Skipping column %s (unsupported type: %s)', col_name, class(sample_value));
        continue;
    end
    
    % Fill values with proper type handling
    for i = 1:height(fused_ranking)
        id = fused_ranking.(id_column){i};
        idx = find(strcmp(first_ranking.(id_column), id));
        
        if ~isempty(idx)
            value = first_ranking.(col_name)(idx(1));
            
            if isnumeric(sample_value)
                new_col(i) = value;
                
            elseif isstring(sample_value)
                % String → store as string in cell
                new_col{i} = value;
                
            elseif iscell(sample_value)
                new_col{i} = value;
                
            elseif ischar(sample_value)
                new_col{i} = value;
            end
        end
    end
    
    % Convert cell array back to appropriate type for table
    if isstring(sample_value)
        % Convert cell of strings back to string array
        fused_ranking.(col_name) = string(new_col);
    else
        fused_ranking.(col_name) = new_col;
    end
end
%% ========================================================================
%  DTW BASELINE - TRAJECTORY SIMILARITY SEARCH
%  ========================================================================
%  Hierarchical similarity search using Dynamic Time Warping (DTW)
%  
%  Author: [Your Name]
%  Date: [Date]
%  
%  Description:
%  Finds similar trajectories to a query trajectory using DTW distance.
%  Works hierarchically: first on trajectory level, then on segment level.
%  
%  Modes:
%  - 'position': DTW on 3D positions (x, y, z)
%  - 'joint_states': DTW on 6D joint angles (joint_1 to joint_6)
%  ========================================================================

%% SECTION 1: CONFIGURATION & PARAMETERS
%  ========================================================================
%  Set all parameters here

clear; clc;

% === Query Parameters ===
query_bahn_id = '1762443169';  % ID der gesuchten Bahn
mode = 'joint_states';                % 'position' oder 'joint_states'

% === DTW Method ===
use_fast_dtw = false;             % false = standard DTW, true = FastDTW
fast_dtw_radius = 1;              % Radius für FastDTW (wenn aktiviert)

% === Pre-filtering Parameters ===
length_tolerance = 0.01;          % 10% Toleranz für Bahnlänge
use_duration_filter = false;      % Optional: Dauer-Filter
duration_tolerance = 0.10;        % 10% Toleranz für Dauer (wenn aktiviert)

% === Output Parameters ===
top_k_trajectories = 5;          % Anzahl Top ähnlichste Bahnen

% === Database Configuration ===
schema = 'bewegungsdaten';
db_name = 'robotervermessung';


%% SECTION 2: DATABASE CONNECTION
%  ========================================================================

fprintf('\n=== Connecting to Database ===\n');

conn = connectingToPostgres;

if isopen(conn)
    fprintf('✓ Database connection successful\n');
else
    error('✗ Database connection failed');
end


%% SECTION 3: LOAD QUERY TRAJECTORY DATA
%  ========================================================================

fprintf('\n=== Loading Query Trajectory: %s ===\n', query_bahn_id);

% Load trajectory data based on mode
if strcmp(mode, 'position')
    query_data = loadTrajectoryPosition(conn, schema, query_bahn_id);
    query_start = query_data(1, :);  % Erster Punkt
    query_data = query_data - query_start;  % Alle Punkte verschieben
    fprintf('✓ Loaded position data: %d points (x, y, z) - normalized to origin\n', size(query_data, 1));
    
elseif strcmp(mode, 'joint_states')
    query_data = loadTrajectoryJointStates(conn, schema, query_bahn_id);
    query_start = query_data(1, :);
    query_data = query_data - query_start;
    fprintf('✓ Loaded joint states data: %d points (6 joints) - normalized to start position\n', size(query_data, 1));
    
else
    error('Invalid mode. Use "position" or "joint_states"');
end

% Load query trajectory metadata
query_metadata = loadTrajectoryMetadata(conn, schema, query_bahn_id);
fprintf('  Query length: %.2f mm\n', query_metadata.length);
fprintf('  Query duration: %.2f s\n', query_metadata.duration);
fprintf('  Query segments: %d\n', query_metadata.num_segments);


%% SECTION 4: PRE-FILTERING CANDIDATES
%  ========================================================================

fprintf('\n=== Pre-filtering Candidate Trajectories ===\n');

% Calculate filter bounds
min_length = query_metadata.length * (1 - length_tolerance);
max_length = query_metadata.length * (1 + length_tolerance);

fprintf('  Length filter: %.2f - %.2f mm (±%.0f%%)\n', ...
    min_length, max_length, length_tolerance * 100);

% Build SQL query for filtering
filter_query = sprintf(...
    ['SELECT bahn_id, length, duration ' ...
     'FROM robotervermessung.%s.bahn_metadata ' ...
     'WHERE bahn_id = segment_id ' ...                     
     'AND length BETWEEN %.2f AND %.2f ' ...
     'AND bahn_id != ''%s'''], ...
    schema, min_length, max_length, query_bahn_id);

% Add duration filter if enabled
if use_duration_filter
    min_duration = query_metadata.duration * (1 - duration_tolerance);
    max_duration = query_metadata.duration * (1 + duration_tolerance);
    filter_query = sprintf('%s AND duration BETWEEN %.2f AND %.2f', ...
    filter_query, min_duration, max_duration);
    fprintf('  Duration filter: %.2f - %.2f s (±%.0f%%)\n', ...
        min_duration, max_duration, duration_tolerance * 100);
end

% Execute filtering
fprintf('DEBUG - SQL Query:\n%s\n', filter_query);
candidate_metadata = fetch(conn, filter_query);
num_candidates = height(candidate_metadata);

fprintf('✓ Found %d candidate trajectories after filtering\n', num_candidates);

if num_candidates == 0
    warning('No candidates found. Try increasing tolerance.');
    close(conn);
    return;
end


%% SECTION 5: DTW COMPUTATION - TRAJECTORY LEVEL
%  ========================================================================

fprintf('\n=== Computing DTW Distances (Trajectory Level) ===\n');
fprintf('  Method: %s\n', iif(use_fast_dtw, 'FastDTW', 'Standard DTW'));

% Initialize results structure
trajectory_results = struct();
trajectory_results.bahn_id = cell(num_candidates, 1);
trajectory_results.dtw_distance = zeros(num_candidates, 1);
trajectory_results.length = zeros(num_candidates, 1);
trajectory_results.duration = zeros(num_candidates, 1);

% Compute DTW for each candidate
% Compute DTW for each candidate
for i = 1:num_candidates
    candidate_id = candidate_metadata.bahn_id{i};
    
    % Load candidate data
    if strcmp(mode, 'position')
        candidate_data = loadTrajectoryPosition(conn, schema, candidate_id);
    else
        candidate_data = loadTrajectoryJointStates(conn, schema, candidate_id);
    end
    
    % === NORMALISIERUNG: Translation auf Startpunkt ===
    if ~isempty(candidate_data)
        candidate_start = candidate_data(1, :);
        candidate_data = candidate_data - candidate_start;
    end
    % === ENDE NORMALISIERUNG ===
    
    % Compute DTW distance
    if use_fast_dtw
        dist = fastDTW(query_data, candidate_data, fast_dtw_radius);
    else
        dist = standardDTW(query_data, candidate_data);
    end
    
    % Store results
    trajectory_results.bahn_id{i} = candidate_id;
    trajectory_results.dtw_distance(i) = dist;
    trajectory_results.length(i) = candidate_metadata.length(i);
    trajectory_results.duration(i) = candidate_metadata.duration(i);
    
    % Progress indicator
    if mod(i, 10) == 0 || i == num_candidates
        fprintf('  Progress: %d/%d trajectories\n', i, num_candidates);
    end
end

% Convert to table and sort by distance
trajectory_table = struct2table(trajectory_results);
trajectory_table = sortrows(trajectory_table, 'dtw_distance', 'ascend');

fprintf('✓ DTW computation completed\n');


%% SECTION 6: DTW COMPUTATION - SEGMENT LEVEL
%  ========================================================================

fprintf('\n=== Computing DTW Distances (Segment Level) ===\n');

% Load query segments
query_segments = loadSegments(conn, schema, query_bahn_id, mode);
num_query_segments = length(query_segments);

fprintf('  Query segments: %d\n', num_query_segments);

% ===== FÜR JEDES QUERY-SEGMENT: UNABHÄNGIGE SUCHE =====
segment_results = cell(num_query_segments, 1);

for seg_idx = 1:num_query_segments
    query_segment_id = sprintf('%s_%d', query_bahn_id, seg_idx);
    query_segment_data = query_segments{seg_idx};
    
    % === NORMALISIERUNG des Query-Segments ===
    if ~isempty(query_segment_data)
        query_seg_start = query_segment_data(1, :);
        query_segment_data = query_segment_data - query_seg_start;
    end
    
    fprintf('\n  Processing Query Segment %d/%d (ID: %s)...\n', seg_idx, num_query_segments, query_segment_id);
    
    % === VORFILTERUNG: Hole Metadata für dieses Query-Segment ===
    query_seg_meta = loadSegmentMetadata(conn, schema, query_segment_id);
    
    seg_min_length = query_seg_meta.length * (1 - length_tolerance);
    seg_max_length = query_seg_meta.length * (1 + length_tolerance);
    
    fprintf('    Query segment length: %.2f mm\n', query_seg_meta.length);
    fprintf('    Length filter: %.2f - %.2f mm (±%.0f%%)\n', ...
        seg_min_length, seg_max_length, length_tolerance * 100);
    
    % === HOLE ALLE KANDIDATEN-SEGMENTE (aus ALLEN Bahnen, nicht nur Top-K!) ===
    query_all_segments = sprintf(...
        ['SELECT segment_id, bahn_id, length, duration ' ...
         'FROM robotervermessung.%s.bahn_metadata ' ...
         'WHERE segment_id != bahn_id ' ...                    % Nur Segmente, keine Bahnen
         'AND length BETWEEN %.2f AND %.2f ' ...               % Längenfilter
         'AND segment_id != ''%s'''], ...                      % Nicht Query-Segment selbst
        schema, seg_min_length, seg_max_length, query_segment_id);
    
    % Optional: Dauer-Filter auch für Segmente
    if use_duration_filter
        seg_min_duration = query_seg_meta.duration * (1 - duration_tolerance);
        seg_max_duration = query_seg_meta.duration * (1 + duration_tolerance);
        query_all_segments = sprintf('%s AND duration BETWEEN %.2f AND %.2f', ...
            query_all_segments, seg_min_duration, seg_max_duration);
    end
    
    candidate_segments = fetch(conn, query_all_segments);
    num_candidate_segments = height(candidate_segments);
    
    fprintf('    Found %d candidate segments after filtering\n', num_candidate_segments);
    
    if num_candidate_segments == 0
        fprintf('    WARNING: No candidate segments found!\n');
        segment_results{seg_idx} = table();
        continue;
    end
    
    % === DTW BERECHNUNG für alle Kandidaten-Segmente ===
    segment_dtw_distances = zeros(num_candidate_segments, 1);
    
    for cand_seg_idx = 1:num_candidate_segments
        candidate_segment_id = candidate_segments.segment_id{cand_seg_idx};
        
        % Load candidate segment data
        if strcmp(mode, 'position')
            candidate_segment_data = loadSegmentData(conn, schema, candidate_segment_id, 'position');
        else
            candidate_segment_data = loadSegmentData(conn, schema, candidate_segment_id, 'joint_states');
        end
        
        % === NORMALISIERUNG ===
        if ~isempty(candidate_segment_data)
            candidate_seg_start = candidate_segment_data(1, :);
            candidate_segment_data = candidate_segment_data - candidate_seg_start;
        end
        
        % Compute DTW
        if use_fast_dtw
            dist = fastDTW(query_segment_data, candidate_segment_data, fast_dtw_radius);
        else
            dist = standardDTW(query_segment_data, candidate_segment_data);
        end
        
        segment_dtw_distances(cand_seg_idx) = dist;
        
        % Progress indicator
        if mod(cand_seg_idx, 50) == 0 || cand_seg_idx == num_candidate_segments
            fprintf('      Progress: %d/%d segments\n', cand_seg_idx, num_candidate_segments);
        end
    end
    
    % === ERSTELLE ERGEBNIS-TABELLE (gleich wie Trajectory-Level) ===
    segment_results{seg_idx} = table(...
        candidate_segments.segment_id, ...
        candidate_segments.bahn_id, ...
        segment_dtw_distances, ...
        candidate_segments.length, ...
        candidate_segments.duration, ...
        'VariableNames', {'segment_id', 'bahn_id', 'dtw_distance', 'length', 'duration'});
    
    % Sortieren nach DTW-Distanz
    segment_results{seg_idx} = sortrows(segment_results{seg_idx}, 'dtw_distance', 'ascend');
    
    fprintf('    ✓ Segment %d: Compared against %d candidates\n', seg_idx, num_candidate_segments);
end

fprintf('\n✓ Segment-level DTW computation completed\n');

%% SECTION 7: RESULTS & RANKING
%  ========================================================================

fprintf('\n========================================\n');
fprintf('RESULTS: TRAJECTORY SIMILARITY RANKING\n');
fprintf('========================================\n\n');

fprintf('Query Trajectory: %s\n', query_bahn_id);
fprintf('Mode: %s\n', mode);
fprintf('DTW Method: %s\n', iif(use_fast_dtw, 'FastDTW', 'Standard DTW'));
fprintf('\n');

% === Trajectory-Level Ranking ===
fprintf('--- TOP %d SIMILAR TRAJECTORIES (Complete Paths) ---\n\n', top_k_trajectories);
fprintf('%-5s %-20s %-15s %-12s %-12s\n', ...
    'Rank', 'Bahn ID', 'DTW Distance', 'Length (mm)', 'Duration (s)');
fprintf('%s\n', repmat('-', 1, 70));

for i = 1:min(top_k_trajectories, height(trajectory_table))
    fprintf('%-5d %-20s %-15.2f %-12.2f %-12.2f\n', ...
        i, ...
        trajectory_table.bahn_id{i}, ...
        trajectory_table.dtw_distance(i), ...
        trajectory_table.length(i), ...
        trajectory_table.duration(i));
end

fprintf('\n\n');

% === Segment-Level Ranking ===
fprintf('--- SEGMENT-LEVEL MATCHES (Individual Segments) ---\n\n');

for seg_idx = 1:num_query_segments
    segment_table = segment_results{seg_idx};
    
    if isempty(segment_table)
        fprintf('Query Segment: %s_%d\n', query_bahn_id, seg_idx);
        fprintf('  No matching segments found!\n\n');
        continue;
    end
    
    fprintf('Query Segment: %s_%d\n', query_bahn_id, seg_idx);
    fprintf('%-5s %-25s %-20s %-15s %-12s %-12s\n', ...
        'Rank', 'Segment ID', 'Parent Bahn ID', 'DTW Distance', 'Length (mm)', 'Duration (s)');
    fprintf('%s\n', repmat('-', 1, 95));
    
    % Zeige Top-K Segmente (gleicher Parameter wie bei Bahnen)
    for i = 1:min(top_k_trajectories, height(segment_table))
        fprintf('%-5d %-25s %-20s %-15.2f %-12.2f %-12.2f\n', ...
            i, ...
            segment_table.segment_id{i}, ...
            segment_table.bahn_id{i}, ...
            segment_table.dtw_distance(i), ...
            segment_table.length(i), ...
            segment_table.duration(i));
    end
    
    fprintf('\n');
end


%% SECTION 8: CLEANUP
%  ========================================================================

fprintf('\n=== Cleanup ===\n');
close(conn);
fprintf('✓ Database connection closed\n');
fprintf('\nAnalysis complete!\n');


%% ========================================================================
%  FUNCTION DEFINITIONS
%  ========================================================================
%  All functions must be at the end of the file in MATLAB

function data = loadTrajectoryPosition(conn, schema, bahn_id)
    % Load position data (x, y, z) for a trajectory
    query = sprintf(...
        ['SELECT x_soll, y_soll, z_soll ' ...
         'FROM robotervermessung.%s.bahn_position_soll ' ...
         'WHERE bahn_id = ''%s'' ' ...
         'ORDER BY timestamp'], ...
        schema, bahn_id);
    
    result = fetch(conn, query);
    data = table2array(result);
end

function data = loadTrajectoryJointStates(conn, schema, bahn_id)
    % Load joint states data (6 joints) for a trajectory
    query = sprintf(...
        ['SELECT joint_1, joint_2, joint_3, joint_4, joint_5, joint_6 ' ...
         'FROM robotervermessung.%s.bahn_joint_states ' ...
         'WHERE bahn_id = ''%s'' ' ...
         'ORDER BY timestamp'], ...
        schema, bahn_id);
    
    result = fetch(conn, query);
    data = table2array(result);
end

function metadata = loadTrajectoryMetadata(conn, schema, bahn_id)
    % Load metadata for trajectory (only the row where bahn_id = segment_id)
    query = sprintf(...
        ['SELECT length, duration ' ...
         'FROM robotervermessung.%s.bahn_metadata ' ...
         'WHERE bahn_id = ''%s'' AND segment_id = ''%s'''], ...    % <-- bahn_id = segment_id
        schema, bahn_id, bahn_id);
    
    result = fetch(conn, query);
    
    metadata = struct();
    metadata.length = result.length;
    metadata.duration = result.duration;
    
    % Count segments: number of rows with same bahn_id
    query_segments = sprintf(...
        ['SELECT COUNT(*) as num_segments ' ...
         'FROM robotervermessung.%s.bahn_metadata ' ...
         'WHERE bahn_id = ''%s'' AND segment_id != bahn_id'], ...  % <-- Nur Segment-Zeilen zählen
        schema, bahn_id);
    
    result_segments = fetch(conn, query_segments);
    metadata.num_segments = result_segments.num_segments;
end

function segments = loadSegments(conn, schema, bahn_id, mode)
    % Load all segments for a trajectory
    query_count = sprintf(...
        ['SELECT COUNT(*) as num_segments ' ...
         'FROM robotervermessung.%s.bahn_metadata ' ...
         'WHERE bahn_id = ''%s'' AND segment_id != bahn_id'], ...  % <-- GEÄNDERT
        schema, bahn_id);
    
    result = fetch(conn, query_count);
    num_segments = result.num_segments;
    
    segments = cell(num_segments, 1);
    
    for i = 1:num_segments
        segment_id = sprintf('%s_%d', bahn_id, i);
        
        if strcmp(mode, 'position')
            query = sprintf(...
                ['SELECT x_soll, y_soll, z_soll ' ...
                 'FROM robotervermessung.%s.bahn_position_soll ' ...
                 'WHERE segment_id = ''%s'' ' ...
                 'ORDER BY timestamp'], ...
                schema, segment_id);
        else  % joint_states
            query = sprintf(...
                ['SELECT joint_1, joint_2, joint_3, joint_4, joint_5, joint_6 ' ...
                 'FROM robotervermessung.%s.bahn_joint_states ' ...
                 'WHERE segment_id = ''%s'' ' ...
                 'ORDER BY timestamp'], ...
                schema, segment_id);
        end
        
        result = fetch(conn, query);
        segments{i} = table2array(result);
    end
end

function dist = standardDTW(seq1, seq2)
    % Standard Dynamic Time Warping algorithm
    % 
    % Inputs:
    %   seq1: N x D matrix (N points, D dimensions)
    %   seq2: M x D matrix (M points, D dimensions)
    % Output:
    %   dist: DTW distance

    if isempty(seq1) || isempty(seq2)
        warning('standardDTW: Empty sequence detected!');
        dist = Inf;
        return;
    end

    if any(isnan(seq1(:))) || any(isnan(seq2(:)))
        warning('standardDTW: NaN values detected in sequences!');
        dist = Inf;
        return;
    end
    
    n = size(seq1, 1);
    m = size(seq2, 1);
    
    % Initialize cost matrix
    dtw_matrix = inf(n+1, m+1);
    dtw_matrix(1, 1) = 0;
    
    % Fill cost matrix
    for i = 1:n
        for j = 1:m
            % Euclidean distance between points
            cost = norm(seq1(i, :) - seq2(j, :));
            
            % DTW recurrence relation
            dtw_matrix(i+1, j+1) = cost + min([...
                dtw_matrix(i, j+1), ...    % insertion
                dtw_matrix(i+1, j), ...    % deletion
                dtw_matrix(i, j)]);        % match
        end
    end
    
    % Return final DTW distance
    dist = dtw_matrix(n+1, m+1);
end

function dist = fastDTW(seq1, seq2, radius)
    % FastDTW - Approximate DTW with linear time complexity
    % 
    % Inputs:
    %   seq1: N x D matrix
    %   seq2: M x D matrix
    %   radius: constraint radius (typically 1)
    % Output:
    %   dist: Approximate DTW distance
    
    min_time_size = radius + 2;
    
    % Base case: use standard DTW for small sequences
    if size(seq1, 1) < min_time_size || size(seq2, 1) < min_time_size
        dist = standardDTW(seq1, seq2);
        return;
    end
    
    % Recursive case: shrink and compute
    % Downsample sequences
    seq1_shrunk = downsampleSequence(seq1);
    seq2_shrunk = downsampleSequence(seq2);
    
    % Recursively compute DTW on smaller sequences
    dist_shrunk = fastDTW(seq1_shrunk, seq2_shrunk, radius);
    
    % Get path from shrunk sequences (simplified for this baseline)
    % Use constrained DTW with radius on original sequences
    dist = constrainedDTW(seq1, seq2, radius);
end

function seq_down = downsampleSequence(seq)
    % Downsample sequence by factor of 2 (simple averaging)
    n = size(seq, 1);
    d = size(seq, 2);
    
    n_down = floor(n / 2);
    seq_down = zeros(n_down, d);
    
    for i = 1:n_down
        idx1 = 2*i - 1;
        idx2 = min(2*i, n);
        seq_down(i, :) = mean(seq(idx1:idx2, :), 1);
    end
end

function dist = constrainedDTW(seq1, seq2, radius)
    % DTW with Sakoe-Chiba band constraint
    n = size(seq1, 1);
    m = size(seq2, 1);
    
    % Initialize cost matrix
    dtw_matrix = inf(n+1, m+1);
    dtw_matrix(1, 1) = 0;
    
    % Fill cost matrix with constraint
    for i = 1:n
        for j = max(1, i-radius):min(m, i+radius)
            % Euclidean distance
            cost = norm(seq1(i, :) - seq2(j, :));
            
            % DTW recurrence
            dtw_matrix(i+1, j+1) = cost + min([...
                dtw_matrix(i, j+1), ...
                dtw_matrix(i+1, j), ...
                dtw_matrix(i, j)]);
        end
    end
    
    dist = dtw_matrix(n+1, m+1);
end

function result = iif(condition, true_val, false_val)
    % Inline if-function (ternary operator)
    if condition
        result = true_val;
    else
        result = false_val;
    end
end

function metadata = loadSegmentMetadata(conn, schema, segment_id)
    % Load metadata for a specific segment
    query = sprintf(...
        ['SELECT length, duration ' ...
         'FROM robotervermessung.%s.bahn_metadata ' ...
         'WHERE segment_id = ''%s'''], ...
        schema, segment_id);
    
    result = fetch(conn, query);
    
    metadata = struct();
    metadata.length = result.length;
    metadata.duration = result.duration;
end

function data = loadSegmentData(conn, schema, segment_id, mode)
    % Load data for a specific segment
    if strcmp(mode, 'position')
        query = sprintf(...
            ['SELECT x_soll, y_soll, z_soll ' ...
             'FROM robotervermessung.%s.bahn_position_soll ' ...
             'WHERE segment_id = ''%s'' ' ...
             'ORDER BY timestamp'], ...
            schema, segment_id);
    else  % joint_states
        query = sprintf(...
            ['SELECT joint_1, joint_2, joint_3, joint_4, joint_5, joint_6 ' ...
             'FROM robotervermessung.%s.bahn_joint_states ' ...
             'WHERE segment_id = ''%s'' ' ...
             'ORDER BY timestamp'], ...
            schema, segment_id);
    end
    
    result = fetch(conn, query);
    data = table2array(result);
end
function embedding = computeMetadataEmbedding(metadata)
    % Movement Type
    movement_str = lower(metadata.movement_type);
    linear_count = sum(movement_str == 'l') + contains(movement_str, 'linear');
    circular_count = sum(movement_str == 'c') + contains(movement_str, 'circular');
    total = linear_count + circular_count;
    
    if total > 0
        linear_ratio = linear_count / total;
        circular_ratio = circular_count / total;
    else
        linear_ratio = 0;
        circular_ratio = 0;
    end
    
    % Length
    length_norm = min(metadata.length / 9000, 1);

    duration_norm = min(metadata.duration / 25, 1);
    
    % Twist Stats (5 features)
    twist_max = 3100;
    min_twist_norm = (metadata.min_twist + twist_max) / (2 * twist_max);
    max_twist_norm = (metadata.max_twist + twist_max) / (2 * twist_max);
    mean_twist_norm = (metadata.mean_twist + twist_max) / (2 * twist_max);
    median_twist_norm = (metadata.median_twist + twist_max) / (2 * twist_max);
    std_twist_norm = metadata.std_twist / twist_max;
    
    % Acceleration Stats (5 features)
    accel_max = 10200;
    min_accel_norm = (metadata.min_acceleration + accel_max) / (2 * accel_max);
    max_accel_norm = (metadata.max_acceleration + accel_max) / (2 * accel_max);
    mean_accel_norm = (metadata.mean_acceleration + accel_max) / (2 * accel_max);
    median_accel_norm = (metadata.median_acceleration + accel_max) / (2 * accel_max);
    std_accel_norm = metadata.std_acceleration / accel_max;
    
    % Zusammensetzen (13D)
    embedding = [linear_ratio, circular_ratio, length_norm, duration_norm, ...
                 min_twist_norm, max_twist_norm, mean_twist_norm, ...
                 median_twist_norm, std_twist_norm, ...
                 min_accel_norm, max_accel_norm, mean_accel_norm, ...
                 median_accel_norm, std_accel_norm];
    
    % L2 normalize
    embedding = embedding / norm(embedding);
end
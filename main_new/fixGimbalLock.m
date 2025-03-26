function euler_fixed = fixGimbalLock(euler_angles)
    euler_fixed = euler_angles;
    
    for i = 1:3  % Check each angle component
        angle_data = euler_angles(:,i);
        
        % Check if we have values close to ±180
        near_180 = abs(abs(angle_data) - 180) < 5;
        
        if any(near_180)
            % If we have values near 180, fix sign flips
            mask_neg = angle_data < 0;
            angle_data(mask_neg) = angle_data(mask_neg) + 360;
            euler_fixed(:,i) = angle_data;
        end
    end
end

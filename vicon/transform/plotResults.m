function plotResults(data_ist, data_ist_trafo, data_orientation_soll, position_soll, euler_soll, euler_trans)
    % Colors
    c1 = [0 0.4470 0.7410];    % Blue - SOLL
    c2 = [0.8500 0.3250 0.0980]; % Orange - IST
    c3 = [0.9290 0.6940 0.1250]; % Yellow - Euler Trans

    % Timestamps in seconds
    time_ist = str2double(data_ist.timestamp);
    time_soll = str2double(data_orientation_soll.timestamp);
    timestamps_ist = (time_ist(:,1)- time_soll(1,1))/1e9;
    timestamps_soll = (time_soll(:,1)- time_soll(1,1))/1e9;

    % Get IST Euler angles
    q_ist = table2array(data_ist(:,8:11));
    q_ist = [q_ist(:,4), q_ist(:,3), q_ist(:,2), q_ist(:,1)];
    euler_ist = rad2deg(quat2eul(q_ist));

    % Plot angles
    figure('Color','white','Name','Euler Angles Comparison', 'Position', [100 100 1200 800])
    
    % Roll
    subplot(3,1,1)
    hold on
    plot(timestamps_soll, euler_soll(:,1), '--', 'Color', c1, 'LineWidth', 5, 'DisplayName', 'SOLL')
    plot(timestamps_ist, euler_ist(:,1), ':', 'Color', c2, 'LineWidth', 1, 'DisplayName', 'IST')
    plot(timestamps_ist, euler_trans(:,1), '-', 'Color', c3, 'LineWidth', 1.5, 'DisplayName', 'Euler Trans')
    title('Yaw Angle')
    ylabel('Angle [°]')
    legend('Location', 'best')
    grid on
    hold off

    % Pitch
    subplot(3,1,2)
    hold on
    plot(timestamps_soll, euler_soll(:,2), '--', 'Color', c1, 'LineWidth', 5, 'DisplayName', 'SOLL')
    plot(timestamps_ist, euler_ist(:,2), ':', 'Color', c2, 'LineWidth', 1, 'DisplayName', 'IST')
    plot(timestamps_ist, euler_trans(:,2), '-', 'Color', c3, 'LineWidth', 1.5, 'DisplayName', 'Euler Trans')
    title('Pitch Angle')
    ylabel('Angle [°]')
    legend('Location', 'best')
    grid on
    hold off

    % Yaw
    subplot(3,1,3)
    hold on
    plot(timestamps_soll, euler_soll(:,3), '--', 'Color', c1, 'LineWidth', 5, 'DisplayName', 'SOLL')
    plot(timestamps_ist, euler_ist(:,3), ':', 'Color', c2, 'LineWidth', 1, 'DisplayName', 'IST')
    plot(timestamps_ist, euler_trans(:,3), '-', 'Color', c3, 'LineWidth', 1.5, 'DisplayName', 'Euler Trans')
    title('Roll Angle')
    xlabel('Time [s]')
    ylabel('Angle [°]')
    legend('Location', 'best')
    grid on
    hold off

    % Plot Position
    figure('Color', 'white', 'Name', 'Position Comparison');
    hold on
    plot3(position_soll(:,1), position_soll(:,2), position_soll(:,3), '--', 'Color', c1, 'LineWidth', 1.5, 'DisplayName', 'SOLL')
    %plot3(data_ist(:,5), data_ist(:,6), data_ist(:,7), ':', 'Color', c2, 'LineWidth', 1.5, 'DisplayName', 'IST')
    plot3(data_ist_trafo(:,1), data_ist_trafo(:,2), data_ist_trafo(:,3), '-', 'Color', c3, 'LineWidth', 2, 'DisplayName', 'Transformed')
    grid on
    xlabel('X [mm]')
    ylabel('Y [mm]')
    zlabel('Z [mm]')
    legend('Location', 'best')
    view(3)
    hold off
end
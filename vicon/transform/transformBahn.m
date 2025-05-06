function [status, message] = transformBahn(bahn_id, plots, upload_single, upload_all, transform_only, schema)
% TRANSFORM_BAHN Transforms trajectory data for a given bahn_id or all trajectories
%   Inputs:
%   bahn_id - ID of the path to transform (or empty if upload_all is true)
%   plots - boolean for plotting (default false)
%   upload_single - boolean for single upload (default true)
%   upload_all - boolean for uploading all (default false)
%   transform_only - boolean for transform only (default false)
%   schema - database schema name (default 'bewegungsdaten')

try
    % Database connection
    datasource = "RobotervermessungMATLAB";
    username = "felixthomas";
    password = "manager";
    conn = postgresql(datasource,username,password);

    % Check connection
    if isopen(conn)
        disp('Connection successful');
    else
        error('Database connection failed');
    end

    % Initialize results
    messages = {};
    all_results = struct('bahn_id', {}, 'status', {}, 'message', {});

    % Get bahn_ids to process
    if upload_all
        % Get all available bahn_ids
        query = ['SELECT bahn_id FROM robotervermessung.' schema '.bahn_info'];
        all_bahn_ids = fetch(conn, query);
        bahn_ids_to_process = all_bahn_ids.bahn_id;
        
        % Get existing transformations
        query = ['SELECT DISTINCT bahn_id FROM robotervermessung.' schema '.bahn_pose_trans'];
        existing_bahn_ids = fetch(conn, query);
        existing_bahn_ids = existing_bahn_ids.bahn_id;
    else
        bahn_ids_to_process = {bahn_id};
        
        % Check if transformation exists for single bahn_id
        query = ['SELECT DISTINCT bahn_id FROM robotervermessung.' schema '.bahn_pose_trans WHERE bahn_id = ''' bahn_id ''''];
        existing = fetch(conn, query);
        existing_bahn_ids = existing.bahn_id;
    end

    % Process each bahn_id
    for i = 1:length(bahn_ids_to_process)
        current_bahn_id = bahn_ids_to_process{i};
        
        % Skip if already exists
        if ismember(current_bahn_id, existing_bahn_ids)
            result.bahn_id = current_bahn_id;
            result.status = false;
            result.message = 'Already transformed - skipping';
            all_results(end+1) = result;
            continue;
        end

        try
            % Search for calibration run
            [calibration_id, is_calibration_run] = findCalibrationRun(conn, current_bahn_id, schema);
            
            % Extract calibration data for position
            tablename_cal = ['robotervermessung.' schema '.bahn_pose_ist'];
            opts_cal = databaseImportOptions(conn,tablename_cal);
            opts_cal.RowFilter = opts_cal.RowFilter.bahn_id == calibration_id;
            data_cal_ist = sqlread(conn,tablename_cal,opts_cal);
            data_cal_ist = sortrows(data_cal_ist,'timestamp');
            
            tablename_cal = ['robotervermessung.' schema '.bahn_events'];
            opts_cal = databaseImportOptions(conn,tablename_cal);
            opts_cal.RowFilter = opts_cal.RowFilter.bahn_id == calibration_id;
            data_cal_soll = sqlread(conn,tablename_cal,opts_cal);
            data_cal_soll = sortrows(data_cal_soll,'timestamp');
            
            % Position data for coordinate transformation
            [trafo_rot, trafo_trans, error_metrics] = calibration(data_cal_ist, data_cal_soll, plots);
            
            % Extract calibration data for orientation
            tablename_cal = ['robotervermessung.' schema '.bahn_orientation_soll'];
            opts_cal = databaseImportOptions(conn,tablename_cal);
            opts_cal.RowFilter = opts_cal.RowFilter.bahn_id == calibration_id;
            data_cal_soll = sqlread(conn,tablename_cal,opts_cal);
            data_cal_soll = sortrows(data_cal_soll,'timestamp');

            % Calculate relative rotation matrix for orientation
            q_transform = calibrateQuaternion(data_cal_ist, data_cal_soll);
            
            clear data_cal opts_cal tablename_cal
            
            % Read complete IST data of the path to be transformed
            query = ['SELECT * FROM robotervermessung.' schema '.bahn_pose_ist ' ...
                    'WHERE robotervermessung.' schema '.bahn_pose_ist.bahn_id = ''' current_bahn_id ''''];
            data_ist = fetch(conn, query);
            data_ist = sortrows(data_ist,'timestamp');
            
            % Read complete SOLL data
            query = ['SELECT * FROM robotervermessung.' schema '.bahn_orientation_soll ' ...
                    'WHERE robotervermessung.' schema '.bahn_orientation_soll.bahn_id = ''' current_bahn_id ''''];
            data_orientation_soll = fetch(conn, query);
            data_orientation_soll = sortrows(data_orientation_soll,'timestamp');
            
            query = ['SELECT * FROM robotervermessung.' schema '.bahn_position_soll ' ...
                    'WHERE robotervermessung.' schema '.bahn_position_soll.bahn_id = ''' current_bahn_id ''''];
            data_position_soll = fetch(conn, query);
            data_position_soll = sortrows(data_position_soll,'timestamp');
            
            % Transform data
            position_ist = table2array(data_ist(:,5:7));
            position_soll = table2array(data_position_soll(:,5:7));
            data_ist_trafo = coordTransformation(position_ist, trafo_rot, trafo_trans);
            q_transformed = transformQuaternion(data_ist, data_orientation_soll, q_transform, trafo_rot);

            % Visualization if requested
            if plots
                euler_trans = quat2eul(q_transformed);
                q_soll = [data_orientation_soll.qw_soll, data_orientation_soll.qx_soll, data_orientation_soll.qy_soll, data_orientation_soll.qz_soll];
                euler_soll = quat2eul(q_soll);
                euler_trans_deg = rad2deg(euler_trans);
                euler_soll_deg = rad2deg(euler_soll);
                euler_trans_deg_fixed = fixGimbalLock(euler_trans_deg);
                euler_soll_deg_fixed = fixGimbalLock(euler_soll_deg);
                plotResults(data_ist, data_ist_trafo, data_orientation_soll, position_soll, euler_soll_deg_fixed, euler_trans_deg_fixed);
            end

            % Save to database if not transform_only
            if ~transform_only
                uploadTransformedData(conn, current_bahn_id, calibration_id, data_ist, data_ist_trafo, q_transformed, schema);
            end
            
            result.bahn_id = current_bahn_id;
            result.status = true;
            result.message = sprintf('Successfully processed bahn_id %s (Calibration: %s)', current_bahn_id, calibration_id);
            
        catch ME
            result.bahn_id = current_bahn_id;
            result.status = false;
            result.message = ['Error: ' ME.message];
        end
        
        all_results(end+1) = result;
    end

    % Prepare final status and message
    successful = sum([all_results.status]);
    total = length(all_results);
    
    if upload_all
        status = successful > 0;
        message = sprintf('Processed %d trajectories. Success: %d, Failed: %d', ...
            total, successful, total - successful);
    else
        status = all_results(1).status;
        message = all_results(1).message;
    end

catch ME
    status = false;
    message = ['Database error: ' ME.message];
end

% Clean up
if exist('conn', 'var')
    close(conn);
end
clear datasource username password
end
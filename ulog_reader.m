function plotdata = ulog_reader(filename)
    %% Initial setup
    % clear; clc;
    misc.filename = filename;
    ulog = ulogreader(misc.filename);
    misc.msgs = readTopicMsgs(ulog);
    misc.t0 = ulog.StartTime;
    misc.tf = ulog.EndTime;
    %% vehicle_local_position
    topics_raw.vehicle_local_position_data = readTopicMsgs(ulog,'TopicNames',{'vehicle_local_position'}, ... 
    'InstanceID',{0},'Time',[misc.t0 misc.tf]);
    plotdata.est_drone_states.p.timestamp = topics_raw.vehicle_local_position_data.TopicMessages{1,1}.timestamp;
    % convert to seconds
    plotdata.est_drone_states.p.vector_time = datevec(plotdata.est_drone_states.p.timestamp);
    plotdata.est_drone_states.p.seconds = (plotdata.est_drone_states.p.vector_time(:,5).*60) + plotdata.est_drone_states.p.vector_time(:,6);
    plotdata.est_drone_states.p.Time = plotdata.est_drone_states.p.seconds - plotdata.est_drone_states.p.seconds(1);
    % read data
    plotdata.est_drone_states.p.Data(:,1) = topics_raw.vehicle_local_position_data.TopicMessages{1,1}{:,5};
    plotdata.est_drone_states.p.Data(:,2) = topics_raw.vehicle_local_position_data.TopicMessages{1,1}{:,6};
    plotdata.est_drone_states.p.Data(:,3) = topics_raw.vehicle_local_position_data.TopicMessages{1,1}{:,7};
    plotdata.Time = plotdata.est_drone_states.p.Time;
    %% vehicle_local_velocity
    plotdata.est_drone_states.v.Time = plotdata.est_drone_states.p.Time;
    % read data
    plotdata.est_drone_states.v.Data(:,1) = topics_raw.vehicle_local_position_data.TopicMessages{1,1}{:,10};
    plotdata.est_drone_states.v.Data(:,2) = topics_raw.vehicle_local_position_data.TopicMessages{1,1}{:,11};
    plotdata.est_drone_states.v.Data(:,3) = topics_raw.vehicle_local_position_data.TopicMessages{1,1}{:,12};
    %% vehicle_local_acceleration
    plotdata.est_drone_states.a.Time = plotdata.est_drone_states.p.Time;
    % read data
    plotdata.est_drone_states.a.Data(:,1) = topics_raw.vehicle_local_position_data.TopicMessages{1,1}{:,16};
    plotdata.est_drone_states.a.Data(:,2) = topics_raw.vehicle_local_position_data.TopicMessages{1,1}{:,17};
    plotdata.est_drone_states.a.Data(:,3) = topics_raw.vehicle_local_position_data.TopicMessages{1,1}{:,18};
    %% vehicle_attitude
    topics_raw.vehicle_attitude_data = readTopicMsgs(ulog,'TopicNames',{'vehicle_attitude'}, ... 
    'InstanceID',{0},'Time',[misc.t0 misc.tf]);
    plotdata.est_drone_states.eta.timestamp = topics_raw.vehicle_attitude_data.TopicMessages{1,1}.timestamp;
    % convert to seconds
    plotdata.est_drone_states.eta.vector_time = datevec(plotdata.est_drone_states.eta.timestamp);
    plotdata.est_drone_states.eta.seconds = (plotdata.est_drone_states.eta.vector_time(:,5).*60) + plotdata.est_drone_states.eta.vector_time(:,6);
    plotdata.est_drone_states.eta.Time = plotdata.est_drone_states.eta.seconds - plotdata.est_drone_states.eta.seconds(1);
    % read data
    plotdata.est_drone_states.eta.Data_q = topics_raw.vehicle_attitude_data.TopicMessages{1,1}{:,2};
    plotdata.est_drone_states.eta.Data = rad2deg(quat2eul( plotdata.est_drone_states.eta.Data_q,'XYZ'));
    %% vehicle_angular_velocity
    topics_raw.vehicle_angular_velocity_data = readTopicMsgs(ulog,'TopicNames',{'vehicle_angular_velocity'}, ... 
    'InstanceID',{0},'Time',[misc.t0 misc.tf]);
    plotdata.est_drone_states.omega.timestamp = topics_raw.vehicle_angular_velocity_data.TopicMessages{1,1}.timestamp;
    % convert to seconds
    plotdata.est_drone_states.omega.vector_time = datevec(plotdata.est_drone_states.omega.timestamp);
    plotdata.est_drone_states.omega.seconds = (plotdata.est_drone_states.omega.vector_time(:,5).*60) + plotdata.est_drone_states.omega.vector_time(:,6);
    plotdata.est_drone_states.omega.Time = plotdata.est_drone_states.omega.seconds - plotdata.est_drone_states.omega.seconds(1);
    % read data
    plotdata.est_drone_states.omega.Data = topics_raw.vehicle_angular_velocity_data.TopicMessages{1,1}{:,2};
    %% vehicle_local_position_setpoint
    topics_raw.vehicle_local_position_setpoint_data = readTopicMsgs(ulog,'TopicNames',{'vehicle_local_position_setpoint'}, ... 
    'InstanceID',{0},'Time',[misc.t0 misc.tf]);
    plotdata.references.p.timestamp = topics_raw.vehicle_local_position_setpoint_data.TopicMessages{1,1}.timestamp;
    % convert to seconds
    plotdata.references.p.vector_time = datevec(plotdata.references.p.timestamp);
    plotdata.references.p.seconds = (plotdata.references.p.vector_time(:,5).*60) + plotdata.references.p.vector_time(:,6);
    plotdata.references.p.Time = plotdata.references.p.seconds - plotdata.references.p.seconds(1);
    % read data
    plotdata.references.p.Data(:,1) = topics_raw.vehicle_local_position_setpoint_data.TopicMessages{1,1}{:,1};
    plotdata.references.p.Data(:,2) = topics_raw.vehicle_local_position_setpoint_data.TopicMessages{1,1}{:,2};
    plotdata.references.p.Data(:,3) = topics_raw.vehicle_local_position_setpoint_data.TopicMessages{1,1}{:,3};
    %% vehicle_local_velocity_setpoint
    plotdata.references.v.Time = plotdata.references.p.Time;
    % read data
    plotdata.references.v.Data(:,1) = topics_raw.vehicle_local_position_setpoint_data.TopicMessages{1,1}{:,6};
    plotdata.references.v.Data(:,2) = topics_raw.vehicle_local_position_setpoint_data.TopicMessages{1,1}{:,7};
    plotdata.references.v.Data(:,3) = topics_raw.vehicle_local_position_setpoint_data.TopicMessages{1,1}{:,8};
    %% vehicle_local_acceleration_setpoint
    plotdata.references.a.Time = plotdata.references.p.Time;
    % read data
    plotdata.references.a.Data = topics_raw.vehicle_local_position_setpoint_data.TopicMessages{1,1}{:,9};
    %% vehicle_attitude sp
    topics_raw.vehicle_attitude_setpoint_data = readTopicMsgs(ulog,'TopicNames',{'vehicle_attitude_setpoint'}, ... 
    'InstanceID',{0},'Time',[misc.t0 misc.tf]);
    plotdata.references.eta.timestamp = topics_raw.vehicle_attitude_setpoint_data.TopicMessages{1,1}.timestamp;
    % convert to seconds
    plotdata.references.eta.vector_time = datevec(plotdata.references.eta.timestamp);
    plotdata.references.eta.seconds = (plotdata.references.eta.vector_time(:,5).*60) + plotdata.references.eta.vector_time(:,6);
    plotdata.references.eta.Time = plotdata.references.eta.seconds - plotdata.references.eta.seconds(1);
    % read data
    plotdata.references.eta.Data_q = topics_raw.vehicle_attitude_setpoint_data.TopicMessages{1,1}{:,5};
    plotdata.references.eta.Data = rad2deg(quat2eul( plotdata.references.eta.Data_q,'XYZ'));
    %% vehicle_rates_setpoint
    topics_raw.vehicle_rates_setpoint_data = readTopicMsgs(ulog,'TopicNames',{'vehicle_rates_setpoint'}, ... 
    'InstanceID',{0},'Time',[misc.t0 misc.tf]);
    plotdata.references.omega.timestamp = topics_raw.vehicle_rates_setpoint_data.TopicMessages{1,1}.timestamp;
    % convert to seconds
    plotdata.references.omega.vector_time = datevec(plotdata.references.omega.timestamp);
    plotdata.references.omega.seconds = (plotdata.references.omega.vector_time(:,5).*60) + plotdata.references.omega.vector_time(:,6);
    plotdata.references.omega.Time = plotdata.references.omega.seconds - plotdata.references.omega.seconds(1);
    % read data
    plotdata.references.omega.Data(:,1) = topics_raw.vehicle_rates_setpoint_data.TopicMessages{1,1}{:,1};
    plotdata.references.omega.Data(:,2) = topics_raw.vehicle_rates_setpoint_data.TopicMessages{1,1}{:,2};
    plotdata.references.omega.Data(:,3) = topics_raw.vehicle_rates_setpoint_data.TopicMessages{1,1}{:,3};
    %% actuators outputs
    topics_raw.actuator_outputs_data = readTopicMsgs(ulog,'TopicNames',{'actuator_outputs'}, ... 
    'InstanceID',{2},'Time',[misc.t0 misc.tf]);
    misc.n_of_actuators = topics_raw.actuator_outputs_data.TopicMessages{1,1}{:,1}(1);
    plotdata.est_drone_states.actuator_outputs.timestamp = topics_raw.actuator_outputs_data.TopicMessages{1,1}.timestamp;
    % convert to seconds
    plotdata.est_drone_states.actuator_outputs.vector_time = datevec(plotdata.est_drone_states.actuator_outputs.timestamp);
    plotdata.est_drone_states.actuator_outputs.seconds = (plotdata.est_drone_states.actuator_outputs.vector_time(:,5).*60) + plotdata.est_drone_states.actuator_outputs.vector_time(:,6);
    plotdata.est_drone_states.actuator_outputs.Time= plotdata.est_drone_states.actuator_outputs.seconds - plotdata.est_drone_states.actuator_outputs.seconds(1);
    % read data
    plotdata.est_drone_states.actuator_outputs.Data = topics_raw.actuator_outputs_data.TopicMessages{1,1}{:,2}(:,1:misc.n_of_actuators);
    %% battery_status
    topics_raw.battery_status_data = readTopicMsgs(ulog,'TopicNames',{'battery_status'}, ... 
    'InstanceID',{0},'Time',[misc.t0 misc.tf]);
    plotdata.battery.timestamp = topics_raw.battery_status_data.TopicMessages{1,1}.timestamp;
    % convert to seconds
    plotdata.battery.vector_time = datevec(plotdata.battery.timestamp);
    plotdata.battery.seconds = ...
        (plotdata.battery.vector_time(:,5).*60)...
        + plotdata.battery.vector_time(:,6);
    plotdata.battery.Time= ...
        plotdata.battery.seconds - plotdata.battery.seconds(1);
    % read data
    plotdata.battery.voltage = topics_raw.battery_status_data.TopicMessages{1,1}{:,1};
    plotdata.battery.voltage_filtered = topics_raw.battery_status_data.TopicMessages{1,1}{:,2};
    plotdata.battery.current = topics_raw.battery_status_data.TopicMessages{1,1}{:,3};
    plotdata.battery.current_filtered = topics_raw.battery_status_data.TopicMessages{1,1}{:,4};
    plotdata.battery.current_avg = topics_raw.battery_status_data.TopicMessages{1,1}{:,5};
    %% sensors
    topics_raw.sensor_combined_data = readTopicMsgs(ulog,'TopicNames',{'sensor_combined'}, ... 
    'InstanceID',{0},'Time',[misc.t0 misc.tf]);
    plotdata.sensors_raw.timestamp = topics_raw.sensor_combined_data.TopicMessages{1,1}.timestamp;
    % convert to seconds
    plotdata.sensors_raw.vector_time = datevec(plotdata.sensors_raw.timestamp);
    plotdata.sensors_raw.seconds = ...
        (plotdata.sensors_raw.vector_time(:,5).*60)...
        + plotdata.sensors_raw.vector_time(:,6);
    plotdata.sensors_raw.Time= ...
        plotdata.sensors_raw.seconds - plotdata.sensors_raw.seconds(1);
    % Gyro
    plotdata.sensors_raw.gyro = topics_raw.sensor_combined_data.TopicMessages{1,1}{:,1};
    % acc
    plotdata.sensors_raw.acc = topics_raw.sensor_combined_data.TopicMessages{1,1}{:,4};
    %% Extra
    % Read all system information.
    extra.systeminfo = readSystemInformation(ulog);
    % Read all initial parameter values.
    extra.params = readParameters(ulog);
    % Read all logged output messages.
    extra.loggedoutput = readLoggedOutput(ulog);
end
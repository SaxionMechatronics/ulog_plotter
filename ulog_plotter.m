close all;

%% Save already-existenting variables
workspaceVars = who; % get worspace variables' names



%% Load plotdata
if ~exist('plotdata','var')
    filename = strcat(paths.output_dir, filesep, 'plotdata');
    
    fprintf('\nLooking for data to plot into: %s', filename);
    
    load(filename);
    
    disp('done');
    disp('Plotting...');
end

%% Plotting options
display.opt.position = [200 200 800 300];
display.opt.fontsize = 20.0;
display.opt.linewidth = 1.5;
display.opt.offset = 15;
display.opt.end = 90;

%% Choose what to display
display.plot.position = 1; % --- p
display.plot.velocity = 1; % --- v
display.plot.attitude = 1; % --- eta
display.plot.errors = 1;
display.plot.omega = 1; % --- omega
display.plot.propellers_thrusts = 0; % --- gamma (x_a)
display.plot.pwms = 1;
display.plot.inputs = 0; % propellers_thrusts_der --- gamma_dot (u)
display.plot.propellers_speeds = 0; % --- w
display.plot.propellers_accels = 0; % --- w_dot
display.plot.errors_boxplot = 1;

%% Handlers
display.pfig = []; % it will contain the handlers to all figures plotted

%% Position tracking
if display.plot.position
    title_str = 'Position tracking';
    
    fh = figure('position', display.opt.position, 'name', title_str);
    set(gca, 'defaulttextinterpreter', 'latex');
    % set background to white
    set(gcf,'color','w');
    hold on;
    plot(plotdata.references.p.Time - display.opt.offset, plotdata.references.p.Data(:, 1),'--','Color','#0072BD','LineWidth',display.opt.linewidth);
    plot(plotdata.references.p.Time - display.opt.offset, plotdata.references.p.Data(:, 2), '--','Color','#D95319', 'LineWidth', display.opt.linewidth);
    plot(plotdata.references.p.Time - display.opt.offset, plotdata.references.p.Data(:, 3), '--','Color','#EDB120', 'LineWidth', display.opt.linewidth);
    plot(plotdata.est_drone_states.p.Time - display.opt.offset  , plotdata.est_drone_states.p.Data(:, 1),'Color','#0072BD','LineWidth', display.opt.linewidth);
    plot(plotdata.est_drone_states.p.Time - display.opt.offset  , plotdata.est_drone_states.p.Data(:, 2),'Color','#D95319','LineWidth', display.opt.linewidth);
    plot(plotdata.est_drone_states.p.Time - display.opt.offset  , plotdata.est_drone_states.p.Data(:, 3),'Color','#EDB120','LineWidth', display.opt.linewidth);
    hold off;
    set(gca, 'FontSize', (display.opt.fontsize - 6))
    title(title_str, 'interpreter', 'latex', 'fontsize', display.opt.fontsize);
    xlim([0 display.opt.end]);
    lh = legend('$p_{r,x}$', '$p_{r,y}$', '$p_{r,z}$', ...
        '$p_{x}$', '$p_{y}$', '$p_{z}$', ...
        'Location', 'southeast', ...
        'Orientation', 'horizontal', 'interpreter', 'latex');
    xlabel('time [s]', 'fontsize', display.opt.fontsize, 'interpreter', 'latex');
    ylabel('[m]', 'fontsize', display.opt.fontsize, 'interpreter', 'latex');
    set(lh, 'Interpreter', 'latex', 'fontsize', display.opt.fontsize);
    grid on; box on
    display.pfig = [display.pfig fh];
end

%% Linear velocity tracking
if display.plot.velocity
    title_str = 'Linear velocity tracking';
    
    fh = figure('position',display.opt.position, 'name', title_str);
    set(gca, 'defaulttextinterpreter', 'latex');   
    % set background to white     
    set(gcf,'color','w');
    hold on;
    plot(plotdata.references.v.Time - display.opt.offset, plotdata.references.v.Data(:, :, 1), '--', 'LineWidth', display.opt.linewidth);
    plot(plotdata.est_drone_states.v.Time - display.opt.offset, plotdata.est_drone_states.v.Data(:, :, 1), 'LineWidth', display.opt.linewidth);
    hold off;
    set(gca, 'FontSize', (display.opt.fontsize - 6));
    title(title_str, 'interpreter', 'latex', 'fontsize', display.opt.fontsize);
    xlim([0 display.opt.end]);
    ylim([-3.7 3]);
    ylabel('[m/s]', 'fontsize', display.opt.fontsize, 'interpreter', 'latex');
    lh = legend('$\dot{p}_{r,x}$', '$\dot{p}_{r,y}$', '$\dot{p}_{r,z}$', ...
        '$\dot{p}_{x}$', '$\dot{p}_{y}$', '$\dot{p}_{z}$', ...
        'Location', 'southeast', ...
        'Orientation', 'horizontal', 'interpreter', 'latex');
    xlabel('time [s]', 'fontsize', display.opt.fontsize, 'interpreter', 'latex');
    set(lh, 'Interpreter', 'latex', 'fontsize', display.opt.fontsize);
    grid on; box on
    display.pfig = [display.pfig fh];
end

%% Attitude (RPY) tracking
if display.plot.attitude
    title_str = 'Orientation tracking';
    
    fh = figure('position', display.opt.position, 'name', 'Orientation tracking');
    set(gca, 'defaulttextinterpreter', 'latex');     
    % set background to white     
    set(gcf,'color','w');
    hold on;
    plot(plotdata.references.eta.Time - display.opt.offset, plotdata.references.eta.Data(:, 1),'--','Color','#0072BD','LineWidth',display.opt.linewidth);
    plot(plotdata.references.eta.Time - display.opt.offset, plotdata.references.eta.Data(:, 2), '--','Color','#D95319', 'LineWidth', display.opt.linewidth);
    plot(plotdata.references.eta.Time - display.opt.offset, plotdata.references.eta.Data(:, 3), '--','Color','#EDB120', 'LineWidth', display.opt.linewidth);
    plot(plotdata.est_drone_states.eta.Time - display.opt.offset  , plotdata.est_drone_states.eta.Data(:, 1),'Color','#0072BD','LineWidth', display.opt.linewidth);
    plot(plotdata.est_drone_states.eta.Time - display.opt.offset  , plotdata.est_drone_states.eta.Data(:, 2),'Color','#D95319','LineWidth', display.opt.linewidth);
    plot(plotdata.est_drone_states.eta.Time - display.opt.offset  , plotdata.est_drone_states.eta.Data(:, 3),'Color','#EDB120','LineWidth', display.opt.linewidth);
    hold off;
    set(gca,'FontSize', (display.opt.fontsize - 6))
    title(title_str, 'interpreter', 'latex', 'fontsize', display.opt.fontsize);
    xlim([0 display.opt.end]);
    ylabel('[deg]', 'fontsize', display.opt.fontsize, 'interpreter', 'latex');
    lh = legend('$\phi_r$', '$\theta_r$', '$\psi_r$', ...
        '$\phi$', '$\theta$', '$\psi$', ...
        'Location', 'southeast', ...
        'Orientation', 'horizontal', 'interpreter', 'latex');
    xlabel('time [s]', 'fontsize', display.opt.fontsize, 'interpreter', 'latex');
    set(lh,'Interpreter', 'latex', 'fontsize', display.opt.fontsize);
    grid on; box on
    display.pfig = [display.pfig fh];
end

%% Attitude tracking plot
if (display.plot.attitude)
    title_str = 'Individual attitude tracking';
    fh = figure('position', [1 1 1 3].*display.opt.position, 'name', title_str);
    set(gca, 'defaulttextinterpreter', 'latex');
    % set background to white
    set(gcf,'color','w');
    sub1 = subplot(3,1,1);
    hold on;
    plot(plotdata.references.eta.Time - display.opt.offset, plotdata.references.eta.Data(:, 1),'--','LineWidth',display.opt.linewidth);
    plot(plotdata.est_drone_states.eta.Time - display.opt.offset  , plotdata.est_drone_states.eta.Data(:, 1),'LineWidth', display.opt.linewidth);
    hold off;
    box on
    set(gca, 'FontSize', (display.opt.fontsize - 6))
    xlim([0 display.opt.end]);
    lh = legend('$\phi_r$', ...
        '$\phi$', ...
        'Location', 'southeast', ...
        'Orientation', 'horizontal', 'interpreter', 'latex');
    ylabel('[deg]', 'fontsize', display.opt.fontsize, 'interpreter', 'latex');
    set(lh, 'Interpreter', 'latex', 'fontsize', display.opt.fontsize);
    grid on;
    
    sub2 = subplot(3,1,2);
    hold on;
    plot(plotdata.references.eta.Time - display.opt.offset, plotdata.references.eta.Data(:, 2), '--', 'LineWidth', display.opt.linewidth);
    plot(plotdata.est_drone_states.eta.Time - display.opt.offset  , plotdata.est_drone_states.eta.Data(:, 2),'LineWidth', display.opt.linewidth);
    hold off;
    box on
    set(gca, 'FontSize', (display.opt.fontsize - 6))
    xlim([0 display.opt.end]);
    lh = legend('$\theta_r$', ...
        '$\theta$', ...
        'Location', 'southeast', ...
        'Orientation', 'horizontal', 'interpreter', 'latex');
    ylabel('[deg]', 'fontsize', display.opt.fontsize, 'interpreter', 'latex');
    set(lh, 'Interpreter', 'latex', 'fontsize', display.opt.fontsize);
    grid on;
    
    sub3 = subplot(3,1,3);
    hold on;
    plot(plotdata.references.eta.Time - display.opt.offset, plotdata.references.eta.Data(:, 3), '--', 'LineWidth', display.opt.linewidth);
    plot(plotdata.est_drone_states.eta.Time - display.opt.offset  , plotdata.est_drone_states.eta.Data(:, 3),'LineWidth', display.opt.linewidth);
    hold off;
    box on
    set(gca, 'FontSize', (display.opt.fontsize - 6))
    xlim([0 display.opt.end]);
    lh = legend('$\psi_r$', ...
        '$\psi$', ...
        'Location', 'southeast', ...
        'Orientation', 'horizontal', 'interpreter', 'latex');
    xlabel('Time [s]', 'fontsize', display.opt.fontsize, 'interpreter', 'latex');
    ylabel('[deg]', 'fontsize', display.opt.fontsize, 'interpreter', 'latex');
    set(lh, 'Interpreter', 'latex', 'fontsize', display.opt.fontsize);
    grid on;
    sgtitle(title_str, 'interpreter', 'latex', 'fontsize', display.opt.fontsize)
    display.pfig = [display.pfig fh];
end

%% Angular velocity tracking
if display.plot.omega
    title_str = 'Angular velocity tracking';
    
    fh = figure('position', display.opt.position, 'name', 'Angular velocity tracking');
    set(gca, 'defaulttextinterpreter', 'latex');     
    % set background to white     
    set(gcf,'color','w');
    hold on;
    plot(plotdata.references.omega.Time - display.opt.offset, plotdata.references.omega.Data(:, :, 1), '--', 'LineWidth', display.opt.linewidth);
    plot(plotdata.est_drone_states.omega.Time - display.opt.offset, plotdata.est_drone_states.omega.Data(:, :, 1), 'LineWidth', display.opt.linewidth);
    hold off;
    set(gca, 'FontSize', (display.opt.fontsize - 6));
    title(title_str, 'interpreter', 'latex','fontsize', display.opt.fontsize);
    xlim([0 display.opt.end]);
    ylim([-0.7 0.5]);
    ylabel('[rad/s]', 'fontsize', display.opt.fontsize, 'interpreter', 'latex');
    lh = legend('$\omega_{r,x}$', '$\omega_{r,y}$', '$\omega_{r,z}$', ...
        '$\omega_{x}$', '$\omega_{y}$', '$\omega_{z}$', ...
        'Location', 'southeast', ...
        'Orientation','horizontal', 'interpreter', 'latex');
    xlabel('time [s]', 'fontsize', display.opt.fontsize, 'interpreter', 'latex');
    set(lh, 'Interpreter', 'latex', 'fontsize', display.opt.fontsize);
    grid on; box on
    display.pfig = [display.pfig fh];
end

%% Angular velocity tracking plot
if (display.plot.omega)
    title_str = 'Individual angular velocity  tracking';
    fh = figure('position', [1 1 1 3].*display.opt.position, 'name', title_str);
    set(gca, 'defaulttextinterpreter', 'latex');
    % set background to white
    set(gcf,'color','w');
    sub1 = subplot(3,1,1);
    hold on;
    plot(plotdata.references.omega.Time - display.opt.offset, plotdata.references.omega.Data(:, 1),'--','LineWidth',display.opt.linewidth);
    plot(plotdata.est_drone_states.omega.Time - display.opt.offset  , plotdata.est_drone_states.omega.Data(:, 1),'LineWidth', display.opt.linewidth);
    hold off;
    box on
    set(gca, 'FontSize', (display.opt.fontsize - 6))
    xlim([0 display.opt.end]);
    lh = legend('$\omega_{r,x}$', ...
        '$\omega_{x}$', ...
        'Location', 'southeast', ...
        'Orientation', 'horizontal', 'interpreter', 'latex');
    ylabel('[deg]', 'fontsize', display.opt.fontsize, 'interpreter', 'latex');
    set(lh, 'Interpreter', 'latex', 'fontsize', display.opt.fontsize);
    grid on;
    
    sub2 = subplot(3,1,2);
    hold on;
    plot(plotdata.references.omega.Time - display.opt.offset, plotdata.references.omega.Data(:, 2), '--', 'LineWidth', display.opt.linewidth);
    plot(plotdata.est_drone_states.omega.Time - display.opt.offset  , plotdata.est_drone_states.omega.Data(:, 2),'LineWidth', display.opt.linewidth);
    hold off;
    box on
    set(gca, 'FontSize', (display.opt.fontsize - 6))
    xlim([0 display.opt.end]);
    lh = legend('$\omega_{r,y}$', ...
        '$\omega_{y}$', ...
        'Location', 'southeast', ...
        'Orientation', 'horizontal', 'interpreter', 'latex');
    ylabel('[deg]', 'fontsize', display.opt.fontsize, 'interpreter', 'latex');
    set(lh, 'Interpreter', 'latex', 'fontsize', display.opt.fontsize);
    grid on;
    
    sub3 = subplot(3,1,3);
    hold on;
    plot(plotdata.references.omega.Time - display.opt.offset, plotdata.references.omega.Data(:, 3), '--', 'LineWidth', display.opt.linewidth);
    plot(plotdata.est_drone_states.omega.Time - display.opt.offset  , plotdata.est_drone_states.omega.Data(:, 3),'LineWidth', display.opt.linewidth);
    hold off;
    box on
    set(gca, 'FontSize', (display.opt.fontsize - 6))
    xlim([0 display.opt.end]);
    lh = legend('$\omega_{r,z}$', ...
        '$\omega_{z}$', ...
        'Location', 'southeast', ...
        'Orientation', 'horizontal', 'interpreter', 'latex');
    xlabel('Time [s]', 'fontsize', display.opt.fontsize, 'interpreter', 'latex');
    ylabel('[deg]', 'fontsize', display.opt.fontsize, 'interpreter', 'latex');
    set(lh, 'Interpreter', 'latex', 'fontsize', display.opt.fontsize);
    grid on;
    sgtitle(title_str, 'interpreter', 'latex', 'fontsize', display.opt.fontsize)
    display.pfig = [display.pfig fh];
end

%% Position and Attitude tracking error
if display.plot.errors
    % Computing the position error
    plotdata.errors.p.Data = plotdata.references.p.Data - plotdata.est_drone_states.p.Data;
    plotdata.errors.p.Time = plotdata.references.p.Time; % Real and desired position should have the same timeline!
    
    title_str = 'Position error';
    
    fh = figure('position', display.opt.position, 'name', 'position error');
    set(gca, 'defaulttextinterpreter', 'latex');     
    % set background to white     
    set(gcf,'color','w');
    plot(plotdata.errors.p.Time-display.opt.offset, plotdata.errors.p.Data, 'LineWidth', display.opt.linewidth);
    set(gca, 'FontSize', (display.opt.fontsize - 6));
    title(title_str, 'interpreter', 'latex', 'fontsize', display.opt.fontsize);
    xlim([0 display.opt.end]);
    ylabel('[m]', 'fontsize', display.opt.fontsize, 'interpreter', 'latex');
    xlabel('time [s]', 'fontsize', display.opt.fontsize, 'interpreter', 'latex');
    lh = legend('$e_{p_x}$', '$e_{p_y}$', '$e_{p_z}$', ...
        'Location', 'southeast', ...
        'Orientation', 'horizontal', 'interpreter', 'latex');
    set(lh, 'Interpreter', 'latex', 'fontsize', display.opt.fontsize);
    grid on; box on
    display.pfig = [display.pfig fh];
    
    % Computing the attitude error
    plotdata.errors.attitude.Data = plotdata.references.eta.Data - plotdata.est_drone_states.eta.Data;
    plotdata.errors.attitude.Time = plotdata.references.eta.Time; % Real and desired attitude should have the same timeline!
    
    title_str = 'Attitude error';
    
    fh = figure('position', display.opt.position, 'name', 'attitude error');
    set(gca, 'defaulttextinterpreter', 'latex');     
    % set background to white     
    set(gcf,'color','w');
    plot(plotdata.errors.attitude.Time -display.opt.offset, plotdata.errors.attitude.Data, ...
        'LineWidth', display.opt.linewidth);
    set(gca, 'FontSize', (display.opt.fontsize - 6));
    title(title_str, 'interpreter', 'latex', 'fontsize', display.opt.fontsize);
    xlim([0 display.opt.end]);
    ylabel('[deg]', 'fontsize', display.opt.fontsize, 'interpreter', 'latex');
    xlabel('time [s]', 'fontsize', display.opt.fontsize, 'interpreter', 'latex');
    lh = legend('$e_{\phi}$', '$e_{\theta}$', '$e_{\psi}$', ...
        'Location', 'southeast', ...
        'Orientation', 'horizontal', 'interpreter', 'latex');
    set(lh, 'Interpreter', 'latex', 'fontsize', display.opt.fontsize);
    grid on; box on
    display.pfig = [display.pfig fh];
end
%% boxplots for attitude erros
if (display.plot.errors_boxplot)
    title_str = 'Attitude errors';
    fh = figure('position', display.opt.position, 'name', title_str);
    set(gca, 'defaulttextinterpreter', 'latex');
    % set background to white
    set(gcf,'color','w');
    boxplot(plotdata.errors.attitude.Data(:, :),'Labels',{'Roll','Pitch','Yaw'});grid;
    set(gca, 'FontSize', display.opt.fontsize-6)
    title(title_str, 'interpreter', 'latex', 'fontsize', display.opt.fontsize);
    ylabel('[deg]', 'fontsize', display.opt.fontsize, 'interpreter', 'latex');


    title_str = 'Position errors';
    fh = figure('position', display.opt.position, 'name', title_str);
    set(gca, 'defaulttextinterpreter', 'latex');
    % set background to white
    set(gcf,'color','w');
    boxplot(plotdata.errors.p.Data(:, :),'Labels',{'x','y','z'});grid;
    set(gca, 'FontSize', display.opt.fontsize-6)
    title(title_str, 'interpreter', 'latex', 'fontsize', display.opt.fontsize);
    ylabel('[m]', 'fontsize', display.opt.fontsize, 'interpreter', 'latex');

end
%% Propellers - Thrusts (Gamma)
if display.plot.propellers_thrusts
    title_str = 'Propellers - Thrusts';
    fh = figure('position', display.opt.position, 'name', title_str);
    set(gca, 'defaulttextinterpreter', 'latex');     
    % set background to white     
    set(gcf,'color','w');
    
    thrust_index_1 = find(nbx_idx==13);
    hold on;
    plot(plotdata.est_drone_states.u.Time, plotdata.est_drone_states.gamma.Data, 'LineWidth', display.opt.linewidth);
    % I can use the same bound of the propellers thrusts since it is equal
    % for all the motors
    % Bogdan's modification in order to have it for the general plotter; in
    % the case that the state constraints change, we will need the index of
    % the position of the first propeller in the state constraints array
    % (i.e. number 13), which is thrust_index from above
    plot(plotdata.Time-display.opt.offset, plotdata.constraints.lbx.Data(:,thrust_index_1:thrust_index_1+nu_lambda-1,1), 'LineWidth', display.opt.linewidth, 'color', 'r');
    plot(plotdata.Time-display.opt.offset, plotdata.constraints.ubx.Data(:,thrust_index_1:thrust_index_1+nu_lambda-1,1), 'LineWidth', display.opt.linewidth, 'color', 'r');
    hold off;
    set(gca, 'FontSize', (display.opt.fontsize - 6));
    title(title_str, 'interpreter', 'latex', 'fontsize', display.opt.fontsize);
    xlim([0 display.opt.end]);
    ylabel('[N]', 'fontsize', display.opt.fontsize, 'interpreter', 'latex');
    
    label_legend = strings(settings.nu_lambda,1);
    for i=1:indata.model.nu_lambda
        label_legend(i,1) = "f$_" + num2str(i) + "$";
    end
    
    lh = legend(label_legend', ...
        'Location', 'southeast', ...
        'Orientation', 'horizontal', 'interpreter', 'latex');
    xlabel('time [s]', 'fontsize', display.opt.fontsize, 'interpreter', 'latex');
    set(lh, 'Interpreter', 'latex', 'fontsize', display.opt.fontsize);
    ymax = max(max(plotdata.constraints.ubx.Data(:,thrust_index_1:thrust_index_1+nu_lambda-1,1)));
    ymin = min(min(plotdata.constraints.lbx.Data(:,thrust_index_1:thrust_index_1+nu_lambda-1,1)));
    ymax = ymax + 0.25 * ymax;
    ymin = ymin - 0.25 * ymin;
    ylim([ymin ymax]);
    grid on; box on
    display.pfig = [display.pfig fh];
end

%% Propellers - Thrusts derivatives (Gamma_dot)
% - dot_omega_prop_MPC is the propellers' angular acceleration dictated by the
% controller.
% Is computed using the controller's output (i.e. gamma_dot vector) and its
% integrated counterpart (i.e. gamma vector).
%
% - dot_omega_prop_(max/min)_MPC are the bounds (max/min) for the propellers'
% angular acceleration that the controller can dictates.
if display.plot.inputs
    
    gamma_dot = plotdata.controls.Data(:, 1:settings.nu_lambda,1);
    gamma_dot_max = plotdata.constraints.ubu.Data(:, 1:settings.nu_lambda,1);
    gamma_dot_min = plotdata.constraints.lbu.Data(:, 1:settings.nu_lambda,1);
    
    title_str = 'Propellers - Thrusts derivatives';
    fh = figure('position', display.opt.position, 'name', title_str);
    % changed the previous code with the one below where are plotted the
    % force derivative for each propeller
    for i=1:1:settings.nu_lambda
        subplot(round(settings.nu_lambda/2),round(settings.nu_lambda/2),i);
        set(gca, 'defaulttextinterpreter', 'latex');
        hold on;
        plot(plotdata.Time-display.opt.offset, gamma_dot(:, i), 'LineWidth', display.opt.linewidth);
        plot(plotdata.Time-display.opt.offset, gamma_dot_max(:, i), 'LineWidth', display.opt.linewidth, 'color','r');
        plot(plotdata.Time-display.opt.offset, gamma_dot_min(:, i), 'LineWidth', display.opt.linewidth, 'color','r');
        hold off;
        xlabel('time [s]', 'fontsize', display.opt.fontsize, 'interpreter', 'latex');
        ylabel('[N/s]', 'fontsize', display.opt.fontsize, 'interpreter', 'latex');
        title("${\dot{f}}_"+num2str(i)+"$", 'fontsize', display.opt.fontsize, 'interpreter', 'latex')
        xlim([0 display.opt.end]);
        ymax = max(gamma_dot_max(:, 1));
        ymin = min(gamma_dot_min(:, 1));
        ymax = ymax + 0.25 * ymax;
        ymin = ymin + 0.25 * ymin;
        ylim([ymin ymax]);
        grid on; box on
    end
    
    display.pfig = [display.pfig fh];
end

%% Propellers - Angular speeds (w)
if display.plot.propellers_speeds
    % omega_prop_MPC is the propellers' speed dictated by the controller.
    % It is computed using the integrated counterpart of the controller's output
    %(i.e. gamma vector).
    
    c_f = indata.model.c_f;
    % Deal with the case in which all the propellers have the same thrust
    % coefficients
    if length(c_f) == 1
        c_f = c_f*ones(settings.nu_lambda,1);
    end
    % Compute omega_prop_MPC
    plotdata.omega_prop_MPC.Data = zeros(length(plotdata.est_drone_states.gamma.Data), settings.nu_lambda);
    for i = 1 : settings.nu_lambda
        plotdata.omega_prop_MPC.Data(:, i) = sqrt(plotdata.est_drone_states.gamma.Data(:, i) / c_f(i));
    end
    plotdata.omega_prop_MPC.Time = plotdata.est_drone_states.u.Time;
    
    title_str = 'Propellers - Angular speeds';
    
    fh = figure('position',display.opt.position, 'name', 'Propellers - Angular speeds');
    set(gca,'defaulttextinterpreter','latex');
    hold on;
    plot(plotdata.Time-display.opt.offset, plotdata.omega_prop_MPC.Data,'LineWidth', display.opt.linewidth);
    % I can use the same bound of the propellers speeds since it is equal
    % for all the motors
    % Bogdan -> changed "6" with thrust_index in the code below for the generality of
    % the code
    plot(plotdata.Time-display.opt.offset, sqrt(plotdata.constraints.lbx.Data(:,thrust_index_1:thrust_index_1+nu_lambda-1,1)./(repmat(c_f',[length(plotdata.omega_prop_MPC.Time),1]))),'LineWidth', display.opt.linewidth, 'color', 'r');
    plot(plotdata.Time-display.opt.offset, sqrt(plotdata.constraints.ubx.Data(:,thrust_index_1:thrust_index_1+nu_lambda-1,1)./(repmat(c_f',[length(plotdata.omega_prop_MPC.Time),1]))),'LineWidth', display.opt.linewidth, 'color', 'r');
    hold off;
    set(gca, 'FontSize', (display.opt.fontsize - 6));
    title(title_str, 'interpreter', 'latex', 'fontsize', display.opt.fontsize);
    xlim([0 display.opt.end]);
    ylim([5 120])
    ylabel('[Hz]', 'fontsize', display.opt.fontsize, 'interpreter', 'latex');
    
    label_legend = strings(settings.nu_lambda,1);
    for i=1:indata.model.nu_lambda
        label_legend(i,1) = "w$_" + num2str(i) + "$";
    end
    lh = legend(label_legend', ...
        'Location', 'southeast', ...
        'Orientation', 'horizontal',  'interpreter', 'latex');
    
    xlabel('time [s]', 'fontsize', display.opt.fontsize, 'interpreter', 'latex');
    set(lh, 'Interpreter', 'latex', 'fontsize', display.opt.fontsize);
    grid on; box on
    display.pfig = [display.pfig fh];
end

%% Propellers - Angular accelerations (w_dot)
% - dot_omega_prop_MPC is the propellers' angular acceleration dictated by the
% controller.
% Is computed using the controller's output (i.e. gamma_dot vector) and its
% integrated counterpart (i.e. gamma vector).
%
% - dot_omega_prop_(max/min)_MPC are the bounds (max/min) for the propellers'
% angular accelerations.

if display.plot.propellers_accels
    
    c_f = indata.model.c_f;
    % Deal with the case in which all the propellers have the same thrust
    % coefficients
    if length(c_f) == 1
        c_f = c_f*ones(settings.nu_lambda,1);
    end
    
    gamma_dot = plotdata.controls;
    
    % Compute omega_prop_MPC
    plotdata.omega_prop_MPC.Data = zeros(length(plotdata.est_drone_states.gamma.Data), settings.nu_lambda);
    for i = 1 : settings.nu_lambda
        plotdata.omega_prop_MPC.Data(:, i) = sqrt(plotdata.est_drone_states.gamma.Data(:, i) / c_f(i));
    end
    plotdata.omega_prop_MPC.Time = plotdata.est_drone_states.u.Time;
    
    % Compute dot_omega_prop_MPC
    omega_dot = zeros(length(gamma_dot.Time), settings.nu_lambda);
    for i = 1 : settings.nu_lambda
        omega_dot(:, i) = (0.5 / c_f(i)) * (gamma_dot.Data(:, i) ./ plotdata.omega_prop_MPC.Data(:, i));
    end
    
    % Compute dot_omega_prop_(max/min)_MPC
    omega_dot_max = plotdata.constraints.ubu.Data(:, 1:settings.nu_lambda, 1) ./ (2* sqrt(plotdata.est_drone_states.gamma.Data * c_f));
    omega_dot_min = plotdata.constraints.lbu.Data(:, 1:settings.nu_lambda, 1) ./ (2* sqrt(plotdata.est_drone_states.gamma.Data * c_f));
    fh = figure('position', display.opt.position, 'name', 'Propellers - Angular accs');
    
    for i=1:indata.model.nu_lambda
        subplot(round(settings.nu_lambda/2),round(settings.nu_lambda/2),i);
        set(gca, 'defaulttextinterpreter', 'latex');
        plot(plotdata.Time-display.opt.offset, omega_dot(:, i), 'LineWidth', display.opt.linewidth);
        hold;
        plot(plotdata.Time-display.opt.offset, omega_dot_max(:, i), 'LineWidth', display.opt.linewidth, 'color','r');
        plot(plotdata.Time-display.opt.offset, omega_dot_min(:, i), 'LineWidth', display.opt.linewidth, 'color','r');
        xlabel('time [s]', 'fontsize', display.opt.fontsize, 'interpreter', 'latex');
        ylabel('[Hz/s]', 'fontsize', display.opt.fontsize, 'interpreter', 'latex');
        
        title(["$\dot{w}_" + num2str(i) + "$"], 'fontsize', display.opt.fontsize, 'interpreter', 'latex');
        
        xlim([0 display.opt.end]);
        ymax = max(omega_dot_max(:, 1));
        ymin = min(omega_dot_min(:, 1));
        ymax = ymax + 0.25 * ymax;
        ymin = ymin + 0.25 * ymin;
        ylim([ymin ymax]);
        grid on; box on
    end
    
    display.pfig = [display.pfig fh];
end

%% PWMs
if display.plot.pwms
    title_str = 'Rotors PWMs';
    
    fh = figure('position', display.opt.position, 'name', title_str);
    set(gca, 'defaulttextinterpreter', 'latex');
    % set background to white
    set(gcf,'color','w');
    hold on;
    for i=1:1:size(plotdata.est_drone_states.actuator_outputs.Data,2)
        plot(plotdata.est_drone_states.actuator_outputs.Time(1:end) - display.opt.offset, plotdata.est_drone_states.actuator_outputs.Data(1:end,i),'LineWidth',display.opt.linewidth);
    end
    hold off;
    set(gca, 'FontSize', (display.opt.fontsize - 6))
    title(title_str, 'interpreter', 'latex', 'fontsize', display.opt.fontsize);
    label_legend = strings(size(plotdata.est_drone_states.actuator_outputs.Data,2),1);
    for i=1:size(plotdata.est_drone_states.actuator_outputs.Data,2)
        label_legend(i,1) = "w$_" + num2str(i) + "$";
    end
    lh = legend(label_legend', ...
        'Location', 'southeast', ...
        'Orientation', 'horizontal',  'interpreter', 'latex');
    xlabel('time [s]', 'fontsize', display.opt.fontsize, 'interpreter', 'latex');
    ylabel('[us]', 'fontsize', display.opt.fontsize, 'interpreter', 'latex');
    ylim([900 2100]);
    xlim([0 display.opt.end]);
    set(lh, 'Interpreter', 'latex', 'fontsize', display.opt.fontsize);
    grid on; box on
    display.pfig = [display.pfig fh];

    fh = figure('position', display.opt.position, 'name', title_str);
    set(gca, 'defaulttextinterpreter', 'latex');
    % set background to white
    set(gcf,'color','w');
    histogram(plotdata.est_drone_states.actuator_outputs.Data)
    set(gca, 'FontSize', (display.opt.fontsize - 6))
    title(title_str, 'interpreter', 'latex', 'fontsize', display.opt.fontsize);
   
    xlabel('[us]', 'fontsize', display.opt.fontsize, 'interpreter', 'latex');
    ylabel('Histogram', 'fontsize', display.opt.fontsize, 'interpreter', 'latex');
    set(lh, 'Interpreter', 'latex', 'fontsize', display.opt.fontsize);
    grid on; box on
    display.pfig = [display.pfig fh];
end
%% Specify variables of this script to keep
vars2keep = {'display'};

%% Workspace cleanup
for i= 1:length(vars2keep)
    if ~any(strcmp(workspaceVars, vars2keep{i}))
        workspaceVars = union(workspaceVars, vars2keep{i});
    end
end
variables2remove = setdiff(who, workspaceVars);
clear(variables2remove{:}); % clean all unecessary variables created here
clear variables2remove
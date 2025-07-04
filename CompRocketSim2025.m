%% 3-DOF Simulation Environment for 2025 IREC

clc;
clear;
close all;

addpath("RASAero\");

%% Variables

M_dry = 23.7;                         % [kg] Dry mass of rocket

motor_fname = "AeroTech_M2500T.rse";    %'thrust_curves/AeroTech_M2500T.rse'
motor_wet_mass = 8.108;                 % [kg] Mass with no fuel
motor_prop_mass = 4.766;                % [kg] Mass of prop
motor_dry_mass = motor_wet_mass - motor_prop_mass;

g = 9.81;                               % [m/s^2] Gravity
h_fins = 0.105;                         % [m] Fin Height
t_fins = 0.0047625;                     % [m] Fin Thickness
N_fins = 4;                             % Number of Fins

pad_altitude = 2871;                    % [m] Midland Air & Space Port Pad Altitude ASL

frontal_area = 0.013439; % [m^2] from RASAero, units converted

rail_length = 5.6; % [m]
rail_angle = 5; % [degrees] angle off of vertical
rail_angle_horizontal = deg2rad(90 - rail_angle); % [rad] angle off of horizontal

%% Drag
% drag_curve = readmatrix("drag_curve.csv"); % col 1 - mach #, col 2 - Cd

rasaero_data = readmatrix("RASAero\RASAero_Non_Extended.xlsx");

drag_curve = approx_drag_curve(rasaero_data, frontal_area);

%% Simulation Initial Conditions + Parameters
dT = 0.005;       % [s]

% position vector [z, s] represents altitude, downrange
% initial conditions
position = [0, 0]; % [m] AGL, downrange
velocity = [0, 0]; % [m/s]
acceleration = [0, 0]; % [m/s]

sim_end_time = 60;
t = 0;

%% Motor Function
motor = motor_generator(dT, motor_fname);
time_lookup = motor.time;
thrust_lookup = motor.thrust_lookup;
prop_mass_lookup = motor.prop_mass_lookup;

motor_burn_time = max(time_lookup); % [s] Motor burn time

%% Recorder Setup
preallocate_iter = sim_end_time / dT;

r_time = zeros(1,preallocate_iter);
r_position = zeros(preallocate_iter, 2);
r_velocity = zeros(preallocate_iter, 2);
r_acceleration = zeros(preallocate_iter, 2);
r_T = zeros(1,preallocate_iter);
r_a = zeros(1,preallocate_iter);
r_P = zeros(1,preallocate_iter);
r_M = zeros(1,preallocate_iter);
r_motor_mass = zeros(1,preallocate_iter);
r_Th = zeros(1,preallocate_iter);
r_Cd = zeros(1,preallocate_iter);
r_Fd = zeros(1,preallocate_iter);
r_Mach = zeros(1,preallocate_iter);

%% Simulation driver

iter = 0;
bool_cont = true;

% Event booleans
event_launch_rail_departure = false;
event_burnout = false;
event_apogee = false;
event_main_deploy = false;

disp("Sim started");

while bool_cont
   % Calcuate Properties at Simulation Time
   t = t + dT;
   iter = iter + 1;

   % if(mod(iter, 200) == 0)
   %     %disp("Iter: " + iter);
   %     disp("t: " + t);
   % end


  % Uses International Standard Atmosphere based on launch pad height
  % Returns Temperature (T), speed of sound (a), pressure (P), density (rho)
  % May want to eventually replace; keep for now

  % [T, a, P, rho] = atmosisa(z);
  % disp([T, a, P, rho]);
  [rho, a, T, P] = stdatmo(position(1) + pad_altitude);
  % disp([T, a, P, rho]);

  speed = norm(velocity);

  if event_launch_rail_departure == true
    direction = velocity ./ speed;
  else
    direction = [sin(rail_angle_horizontal),cos(rail_angle_horizontal)];
  end

  

%% Calculate forces and z_dot_dot
    % Atmospheric Drag
    Fd = drag_force(drag_curve, speed, frontal_area, a, rho);% gvign nan
    %disp(Fd)

    % drag is opposite the velocity vector
    Fd = -1 * direction * Fd;

    % Motor Thrust and Mass
    if t >= 0 && t <= motor_burn_time
        Th = thrust_lookup(iter);
        motor_mass = prop_mass_lookup(iter) + motor_dry_mass;
    else
        Th = 0;
        motor_mass = motor_dry_mass; %[kg]
    end

    Th = Th*direction;

    % Weight
    M = M_dry + motor_mass;
    W = (M*g)*[1,0];

    % Solve governing eqn for z_dot_dot
    acceleration = (Th + Fd - W)/M;

    % Calculate any other additional parameters
    mach = speed / a;

  %% Log Current Values to the Recorders
    
    r_time(iter) = t;
    r_position(iter,1:2) = position;
    r_velocity(iter,1:2) = velocity;
    r_acceleration(iter,1:2) = acceleration;
    r_T(iter) = T;
    r_a(iter) = a;
    r_P(iter) = P;
    r_M(iter) = M;
    r_motor_mass(iter) = motor_mass;
    r_Th(iter) = norm(Th);
    r_Fd(iter) = norm(Fd);
    r_Mach(iter) = mach;
   

  %% Calculate z and z_dot for the next timestep
    %z_dot = z_dot + z_dot_dot * dT;
    %z = z + z_dot * dT;
    
    velocity = velocity + acceleration * dT;
    position = position + velocity * dT;


  %% Check for simulation events
    if (t < sim_end_time && position(1) >= 0) || (iter < 10) % simulation continues
        % when sim is good to continue
        bool_cont = true;
    
        % launch rail departure
        if (norm(position) >= rail_length) && ~event_launch_rail_departure
            event_launch_rail_departure = true;
            disp("Launch rail departure at t = " + t);
        end

        % burnout
        if (t >= motor_burn_time) && ~event_burnout
            event_burnout = true;
            disp("Motor burnout at t = " + t);
        end

        % apogee
        if (velocity(1) <= 0) && (~event_apogee) && (event_launch_rail_departure)
            event_apogee = true;
            disp("Apogee at t = " + t);
        end
        
        % main deployment
        if (event_apogee && (position(1) - pad_altitude) <= 304.8) && ~event_main_deploy
            event_main_deploy = true;
            disp("Main deployment at t=" + t);
        end

    else
        % end sim
        disp("Sim ended at iter = " + iter + " t = " + t + " z = " + position(1) + " vs. pad alt: " + pad_altitude);
        bool_cont = false;
    end
end

% Copy z recorder, modify for ASL altitude
r_z_asl = r_position(:,1) + pad_altitude;

% trim preallocation
r_time = r_time(1:iter);
r_position = r_position(1:iter, :);
r_velocity = r_velocity(1:iter, :);
r_acceleration = r_acceleration(1:iter, :);
r_T = r_T(1:iter);
r_a = r_a(1:iter);
r_P = r_P(1:iter);
r_M = r_M(1:iter);
r_motor_mass = r_motor_mass(1:iter);
r_Th = r_Th(1:iter);
r_Fd = r_Fd(1:iter);
r_Mach = r_Mach(1:iter);
r_z_asl = r_z_asl(1:iter);

%% Plots

fig_1_transform = figure;
set(fig_1_transform, 'Units', 'Normalized', 'OuterPosition', [0, 0.25, 0.5, 0.75]);

subplot(3,1,1);
plot(r_time, r_position(:,1)); % AGL altitude
% plot(r_time, r_z_asl); % ASL altitude
title("Altitude vs Time");
xlabel("Time [s]");
ylabel("Altitude [m]");
grid on;

subplot(3,1,2); 
plot(r_time, r_velocity(:,1));
title("Vertical Velocity vs Time");
xlabel("Time [s]");
ylabel("Velocity [m/s]");
yline(0,'k');
grid on;

subplot(3,1,3); 
plot(r_time, r_acceleration(:,1));
title("Vertical Acceleration vs Time");
xlabel("Time [s]");
ylabel("Acceleration [m/s^2]");
yline(0,'k');
grid on;

fig_2_transform = figure;
set(fig_2_transform, 'Units', 'Normalized', 'OuterPosition', [0.5, 0.25, 0.5, 0.75]);

subplot(3,1,1);
plot(r_time, r_position(:,2)); % downrange
title("Downrange Distance vs Time");
xlabel("Time [s]");
ylabel("Altitude [m]");
grid on;

subplot(3,1,2); 
plot(r_time, r_velocity(:,2));
title("Horiziontal Velocity vs Time");
xlabel("Time [s]");
ylabel("Velocity [m/s]");
yline(0,'k');
grid on;

subplot(3,1,3); 
plot(r_time, r_acceleration(:,2));
title("Horiziontal Acceleration vs Time");
xlabel("Time [s]");
ylabel("Acceleration [m/s^2]");
yline(0,'k');
grid on;

%% Flight summary

disp("Flight summary")

% launch rail departure (t, z_dot)
past_rail = r_z_asl(:) >= rail_length;
past_rail_idx = find(past_rail==1,1,"first");
past_rail_time = r_time(past_rail_idx);
past_rail_z_dot = r_velocity(past_rail_idx,1);
disp("Clears Rail at " + past_rail_time + " [s] at speed " + past_rail_z_dot + " [m/s]")

% burnout (t, z, z_agl, z_dot)
burnout_idx = find(r_time >= motor_burn_time, 1, "first");
burnout_t = r_time(burnout_idx);
burnout_z = r_position(burnout_idx,1);
burnout_z_agl = r_z_asl(burnout_idx);
burnout_z_dot = r_velocity(burnout_idx,1);
disp("Burnout occurs at " + burnout_t + " [s] when Altitude = " + burnout_z + " [m] ASL (" + burnout_z_agl + " [m] AGL) and Velocity = " + burnout_z_dot + " [m/s]");

% apogee (t, z, z_agl)
z_max = max(r_position(:,1));
z_max_idx = find(r_position(:,1)==z_max);
z_max_t = r_time(z_max_idx);
z_asl_max = max(r_z_asl);
disp("Apogee of " + z_asl_max + " [m] ASL (" + z_max + " [m] AGL) at t = " + z_max_t + " [s] ");

% maximums (z_dot, mach)
z_dot_max = max(r_velocity(:,1));
z_dot_max_idx = find(r_velocity(:,1)==z_dot_max);
z_dot_max_t = r_time(z_dot_max_idx);
disp("Max Velocity of " + z_dot_max + " [m/s] at t = " + z_dot_max_t + " [s] ");

%% Load Reference Data and Plot to Compare

% These are all copy pastes from last year's sim- just for sintax reference

or_alt_csv = readmatrix('OR_alt.csv');
or_t = or_alt_csv(:,1);
or_z = or_alt_csv(:,2);

%or_data = readtable(fullfile('or_sim_data', 'all_data_4.csv'));
%or_time = table2array(or_data(:,"x_Time_s_"))';
%or_z = table2array(or_data(:, "Altitude_m_"))';
%or_z_dot = table2array(or_data(:, "TotalVelocity_m_s_"))';
%or_z_dot_dot = table2array(or_data(:, "TotalAcceleration_m_s__"))';
%or_mass = table2array(or_data(:, "Mass_g_"))';
%or_thrust = table2array(or_data(:, "Thrust_N_"))';
%or_drag_coefficient = table2array(or_data(:, "DragCoefficient___"))';
%or_drag_force = table2array(or_data(:, "DragForce_N_"))';
%or_motor_mass = table2array(or_data(:, "MotorMass_g_"))';
%or_speed_of_sound = table2array(or_data(:,"SpeedOfSound_m_s_"))';
%or_mach = table2array(or_data(:,"MachNumber___"))';
%or_pressure = table2array(or_data(:,"AirPressure_mbar_"))';

%% Plot Our values and O.R. Values over each other
    % position
    figure;
    hold on;
    plot(r_time, r_position(:,1));
    plot(or_t, or_z);
    title('Position AGL (m)');
    ylabel('Position (m)');
    xlabel('Time (s)');
    legend("3 DoF", "OpenRocket");

    % velocity
    %figure(2)
    %plot(time, r_z_dot, or_time, or_z_dot)
    %title('Velocity (m/s)')
    %ylabel('Velocity (m/s)')
    %xlabel('Time (s)')
    %legend("1 DoF", "OpenRocket")

    % acceleration
    %figure(3)
    %plot(time, abs(r_z_dot_dot), or_time, or_z_dot_dot)
    %title('Acceleration (m/s^2)')
    %legend("1 DoF", "OpenRocket")
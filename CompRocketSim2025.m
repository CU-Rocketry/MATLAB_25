%% 3-DOF Simulation Environment for 2025 IREC

clc;
clear;
close all;

addpath("RASAero\");

%% Variables

M_dry = 23;                         % [kg] Dry mass of rocket

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

%% Drag
% drag_curve = readmatrix("drag_curve.csv"); % col 1 - mach #, col 2 - Cd

rasaero_data = readmatrix("RASAero\RASAero_Non_Extended.xlsx");

drag_curve = approx_drag_curve(rasaero_data, frontal_area);

%% Simulation Initial Conditions + Parameters
dT = 0.005;       % [s]
z = pad_altitude; % [m]
z_dot = 0;        % [m/s]
z_dot_dot = 0;    % [m/s^2]

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
r_z = zeros(1,preallocate_iter);
r_z_dot = zeros(1,preallocate_iter);
r_z_dot_dot = zeros(1,preallocate_iter);
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
apogee_reached = false;

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
  [rho, a, T, P] = stdatmo(z);
  % disp([T, a, P, rho]);

%% Calculate forces and z_dot_dot
    % Atmospheric Drag
    Fd = drag_force(drag_curve, z_dot, frontal_area, a, rho);% gvign nan
    %disp(Fd)

    if z_dot > 0
        % drag is opposite the velocity vector
        Fd = -1 * Fd;
    end

    % Motor Thrust and Mass
    if t >= 0 && t <= motor_burn_time
        Th = thrust_lookup(iter);
        motor_mass = prop_mass_lookup(iter) + motor_dry_mass;
    else
        Th = 0;
        motor_mass = motor_dry_mass; %[kg]
    end

    % Weight
    M = M_dry + motor_mass;
    W = M*g;

    % Solve governing eqn for z_dot_dot
    z_dot_dot = (Th + Fd - W)/M;

    % Calculate any other additional parameters
    mach = z_dot / a;

  %% Log Current Values to the Recorders
    
    r_time(iter) = t;
    r_z(iter) = z;
    r_z_dot(iter) = z_dot;
    r_z_dot_dot(iter) = z_dot_dot;
    r_T(iter) = T;
    r_a(iter) = a;
    r_P(iter) = P;
    r_M(iter) = M;
    r_motor_mass(iter) = motor_mass;
    r_Th(iter) = Th;
    r_Fd(iter) = Fd;
    r_Mach(iter) = mach;
   

  %% Calculate z and z_dot for the next timestep
    z_dot = z_dot + z_dot_dot * dT;
    z = z + z_dot * dT;
 
  %% Check for Simulation Events
        % if statement for when apogee_reached = true;

  %% Evaluate if sim continues
    if (t < sim_end_time && z >= pad_altitude) || (iter < 10)
        % when sim is good to continue
        bool_cont = true;
    else
        % end sim
        disp("Sim ended at iter=" + iter + " t=" + t + " z=" + z + " vs. pad alt: " + pad_altitude);
        bool_cont = false;
    end
end

% Create recorder for AGL altitude
r_z_agl = r_z - pad_altitude;

% trim preallocation
r_time = r_time(1:iter);
r_z = r_z(1:iter);
r_z_dot = r_z_dot(1:iter);
r_z_dot_dot = r_z_dot_dot(1:iter);
r_T = r_T(1:iter);
r_a = r_a(1:iter);
r_P = r_P(1:iter);
r_M = r_M(1:iter);
r_motor_mass = r_motor_mass(1:iter);
r_Th = r_Th(1:iter);
r_Fd = r_Fd(1:iter);
r_Mach = r_Mach(1:iter);
r_z_agl = r_z_agl(1:iter);

%% Plots

fig_transform = figure;
set(fig_transform, 'Units', 'Normalized', 'OuterPosition', [0, 0.25, 0.5, 0.75]);

subplot(3,1,1);
% plot(r_time, r_z); % ASL altitude
plot(r_time, r_z_agl);
title("Altitude vs Time");
xlabel("Time [s]");
ylabel("Altitude [m]");
grid on;

subplot(3,1,2); 
plot(r_time, r_z_dot);
title("Vertical Velocity vs Time");
xlabel("Time [s]");
ylabel("Velocity [m/s]");
yline(0,'k');
grid on;

subplot(3,1,3); 
plot(r_time, r_z_dot_dot);
title("Vertical Acceleration vs Time");
xlabel("Time [s]");
ylabel("Acceleration [m/s^2]");
yline(0,'k');
grid on;

%% Flight summary

% max altitude
z_agl_max = max(r_z_agl);
disp("Apogee (AGL) = " + z_agl_max + " [m]");

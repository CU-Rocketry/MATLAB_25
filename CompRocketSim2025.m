%% 3-DOF Simulation Environment for 2025 IREC

%% Variables

motor_fname = "AeroTech_M2500.eng";     %'thrust_curves/AeroTech_M2500T.rse'
motor_wet_mass = 8.108;                 % [kg] Mass with no fuel
motor_prop_mass = 4.766;                % [kg] Mass of prop
motor_dry_mass = motor_wet_mass - motor_prop_mass;

g = 9.81;                               % [m/s^2] Gravity
h_fins = 0.105;                         % [m] Fin Height
t_fins = 0.0047625;                     % [m] Fin Thickness
N_fins = 4;                             % Number of Fins

pad_altitude = 2871;                    % [m] Midland Air & Space Port Pad Altitude ASL

% Drag
% drag_curve = readmatrix("drag_curve.csv"); % col 1 - mach #, col 2 - Cd

frontal_area = 0.013439; % [m^2] from RASAero, units converted
rasaero_data = readmatrix("RASAero\RASAero_Non_Extended.xlsx");

max_mach = 3;

cd_curve = [rasaero_data(1:max_mach*100,1),rasaero_data(1:max_mach*100,3)];

figure;
plot(cd_curve(:,1), cd_curve(:,2));
xlabel("Mach number");
ylabel("Drag Coefficient");

coeffs = polyfit(cd_curve(:,1), cd_curve(:,2),10);
y = polyval(coeffs,cd_curve(:,1));
figure;
plot(cd_curve(:,1), y);
xlabel("Mach number");
ylabel("Drag Coefficient");

%% Simulation Initial Conditions + Parameters
dT = 0.005;       % [s]
z = pad_altitude; % [m]
z_dot = 0;        % [m/s]
z_dot_dot = 0;    % [m/s^2]

sim_end_time = 60;
t = 0;

%% Functions

%% Drag Function
drag_curve = readmatrix("drag_curve.csv");
drag_curve_mach = drag_curve(:,1);
drag_curve_drag = drag_curve(:,2);

frontal_area = 0.013439; % [m^2] from RASAero, units converted

figure;
plot(drag_curve_mach, drag_curve_drag);
xlabel("Mach number");
ylabel("Drag Coefficient");

% D = Cd*0.5*ro*v^2*A

%% Motor Function
motor = motor_generator(dT, motor_fname);

%% Recorder Setup

time = [];
r_z = [];
r_z_dot = [];
r_z_dot_dot = [];
r_T = [];
r_a = [];
r_P = [];
r_M = [];
r_motor_mass = [];
r_Th = [];
r_Cd = [];
r_Fd = [];
r_Mach = [];

%% Simulation driver

iter = 0;
bool_cont = true;
apogee_reached = false;

while bool_cont
   % Calcuate Properties at Simulation Time
   t = t + dT;
   iter = iter + 1;

  % Uses International Standard Atmosphere based on launch pad height
  % Returns Temperature (T), speed of sound (a), pressure (P), density (rho)
    [T, a, P, rho] = atmosisa(z); 

% Calculate forces and z_dot_dot
    % Atmospheric drag
    % Motor Thrust and Mass
    % Weight
    % Solve governing eqn for z_dot_dot
    % Calculate any other additional parameters

  % Log Current Values to the Recorders
   %{
   r_z(iter) = ;
   r_z_dot(iter) = ;
   r_z_dot_dot(iter) = ;
   time(iter) = ;
   r_T(iter) = ;
   r_a(iter) = ;
   r_P(iter) = ;
   r_M(iter) = ;
   r_motor_mass(iter) = ;
   r_Th(iter) = ;
   r_Cd(iter) = ;
   r_Fd(iter) = ;
   r_Mach(iter) = ;
   %}
  
  % Check for Simulation Events
        % if statement for when apogee_reached = true;

  % Evaluate if sim continues
    if (t < sim_end_time && z > 0) || (iter < 5)
        % when sim is good to continue
        bool_cont = true;
    else
        % end sim
        bool_cont = false;
    end
end

%% Plots

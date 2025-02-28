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

%% Functions

%% Motor Function


t = 0;
%% Simulation driver

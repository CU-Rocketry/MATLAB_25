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

%% fitting
% index 1 to 4 quad
% index 5 to 91 quad
% index 92 to 105 linear
% index 106 to end linear

fit1 = polyfit(cd_curve(1:4,1), cd_curve(1:4,2),2);
fit2 = polyfit(cd_curve(5:91,1), cd_curve(5:91,2),2);
fit3 = polyfit(cd_curve(92:105,1), cd_curve(92:105,2),1);
fit4 = polyfit(cd_curve(106:150,1), cd_curve(106:150,2),1);

mach = linspace(0,3,3*1000);
% 
% y1 = polyval(fit1,mach);
% y2 = polyval(fit2,mach);
% y3 = polyval(fit3,mach);
% y4 = polyval(fit4,mach);

p1 = poly2sym(fit1);
p2 = poly2sym(fit2);
p3 = poly2sym(fit3);
p4 = poly2sym(fit4);

int1 = min(solve(p1 == p2)); % solve where p1=p2, take lower of two
int2 = min(solve(p2 == p3));
int3 = solve(p3 == p4);

y = piecewise((x>=0) & (x<int1),p1, (x>=int1) & (x<int2),p2, (x>=int2) & (x<int3),p3, (x>=int3) & (x<2),p4);

figure;
hold on
fplot(y)
xlim([0,2])
ylim([0.4,0.65])
plot(cd_curve(:,1), cd_curve(:,2), "o")

% coeffs = polyfit(cd_curve(:,1), cd_curve(:,2),10);
% y = polyval(coeffs,cd_curve(:,1));
figure;
hold on;
%plot(cd_curve(:,1), y);
plot(mach,y1);
plot(mach,y2);
plot(mach,y3);
plot(mach,y4);

xlabel("Mach number");
ylabel("Drag Coefficient");

xlim([0,1.5])
ylim([0.4,0.65])

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

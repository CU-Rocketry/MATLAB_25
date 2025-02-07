%% 3-DOF Simulation Environment for 2025 IREC

%% Variables

motor_fname = "AeroTech_M2500.eng"; %'thrust_curves/AeroTech_M2500T.rse'
motor_wet_mass = 8.108;  % [kg] Mass with no fuel
motor_prop_mass = 4.766; % [kg] Mass of prop
motor_dry_mass = motor_wet_mass - motor_prop_mass;

g = 9.81;                % [m/s^2] Gravity
%% Functions
% Drag Function
drag_curve = readmatrix("drag_curve.csv");
drag_curve_mach = drag_curve(:,1);
drag_curve_drag = drag_curve(:,2);

figure;
plot(drag_curve_mach, drag_curve_drag);
xlabel("Mach number");
ylabel("Maybe the Cd??")

% Motor Function


% Simulation driver
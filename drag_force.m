function f_drag = drag_force(drag_curve, rho, velocity, temperature, frontal_area)


% Calculating Cd
    a = sqrt(1.4*286*temperature); % speed of sound [m/s], assuming constant gamma = 1.4, R = 286

    mach = velocity / a; % convert [m/s] to mach number
    
    Cd = drag_curve(mach);

% Calculating drag
    f_drag = Cd*0.5*rho*velocity^2*frontal_area;

end
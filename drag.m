function f_drag = drag(drag_curve, density_curve, altitude, velocity)
    % ignore temp_curve for now plz :'(
    mach = velocity / 343.3; % convert [m/s] to mach number (constant temperature)
    
    % D=Cd*0.5*ro*v^2*A
    f_drag = drag_curve

    f_drag = 0;
end
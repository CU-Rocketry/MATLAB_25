% Args:
%   - drag_curve: a symbolic function with respect to x of Cd versus mach
%   - velocity: rocket's velocity, assuming into airstream
%   - frontal area: wetted area of the rocket based on RASAero
%   - altitude: meters
function f_drag = drag_force(drag_curve, velocity, frontal_area, a, rho)
    mach = velocity / a; % convert [m/s] to mach number based on atmosphere
    %disp("Mach: " + mach);
    
    % syms x;
    % Cd = double(subs(drag_curve, x, mach));
    Cd = drag_curve(mach);
    % disp("Cd: " + Cd);

    % Calculate and return drag force (D=Cd*(1/2)*rho*v^2*A)
    f_drag = Cd*0.5*rho*velocity^2*frontal_area;
end
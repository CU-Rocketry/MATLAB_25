function drag_curve = approx_drag_curve(rasaero_data, frontal_area)
    max_mach = 3;

    cd_curve = [rasaero_data(1:max_mach*100,1),rasaero_data(1:max_mach*100,3)];
    
    figure;
    plot(cd_curve(:,1), cd_curve(:,2));
    xlabel("Mach number");
    ylabel("Drag Coefficient");
    
    % Fitting
    % index 1 to 4 quad
    % index 5 to 91 quad
    % index 92 to 105 linear
    % index 106 to end linear
    
    fit1 = polyfit(cd_curve(1:4,1), cd_curve(1:4,2),2);
    fit2 = polyfit(cd_curve(5:91,1), cd_curve(5:91,2),2);
    fit3 = polyfit(cd_curve(92:105,1), cd_curve(92:105,2),1);
    fit4 = polyfit(cd_curve(106:150,1), cd_curve(106:150,2),1);
    
    % 
    % y1 = polyval(fit1,mach);
    % y2 = polyval(fit2,mach);
    % y3 = polyval(fit3,mach);
    % y4 = polyval(fit4,mach);
        
    p1 = poly2sym(fit1);
    p2 = poly2sym(fit2);
    p3 = poly2sym(fit3);
    p4 = poly2sym(fit4);
    
    int1 = double(min(solve(p1 == p2))); % solve where p1=p2, take lower of two
    int2 = double(min(solve(p2 == p3)));
    int3 = double(solve(p3 == p4));
    
    %y = piecewise((x<0),0, (x>=0) & (x<int1),p1, (x>=int1) & (x<int2),p2, (x>=int2) & (x<int3),p3, (x>=int3) & (x<2),p4);

    drag_curve = @(x) (x < 0) .* 0 + ...
    ((x >= 0) & (x < int1)) .* polyval(fit1, x) + ...
    ((x >= int1) & (x < int2)) .* polyval(fit2, x) + ...
    ((x >= int2) & (x < int3)) .* polyval(fit3, x) + ...
    ((x >= int3) & (x < 2)) .* polyval(fit4, x);
    
    % example usage of drag_force:
    % f_drag = drag_force(y,200,frontal_area, a, rho);
    % disp(f_drag);

    
    
    % uncomment to plot approximated drag curve
    hold on;
    %plot(cd_curve(:,1), cd_curve(:,2), "r.", "MarkerSize",10); % plot original
    plot(cd_curve(:,1), drag_curve(cd_curve(:,1)), "r");
    legend("Original", "Regression");
    xlabel("Mach number");
    ylabel("Drag Coefficient");
    grid on;
    xlim([0,2]);
    ylim([0.4,0.65]);

    % for comparison to OR
    or_drag_csv = readmatrix('OR_drag.csv');
    or_drag_mach = or_drag_csv(:,1);
    or_drag_cd = or_drag_csv(:,2);

    plot(or_drag_mach, or_drag_cd);
    legend("Original", "Regression", "OpenRocket");
    hold off;
end

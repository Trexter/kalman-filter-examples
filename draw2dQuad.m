function draw2dQuad(pos_function, theta_function, synthetic_z_func, mu, t_eval)
    
    syms t

    theta = double(subs(theta_function, t, t_eval));
    pos = double(subs(pos_function, t, t_eval));
    
    R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    
    z = double(subs(synthetic_z_func, t, t_eval));
    
    
    hold on

    axis([-8 8 0 8 0 1])    
    view(0,90)
    pbaspect([2 1 1])
    
    
    
    quiver(pos(1), pos(2), -z(1)*sin(theta), -z(1)*cos(theta));
    
    quiver(pos(1), pos(2), cos(theta), -sin(theta));
    quiver(pos(1), pos(2), sin(theta), cos(theta));
    
    text(pos(1) -z(1)/2 * sin(theta), pos(2) -z(1)/2 * cos(theta), num2str(z(1), 'Range Sensor: %.2f'))
    
    text(-7.5, 0.5, num2str(t_eval, 'Time: %.2f'))
    
    
    quiver(mu(1), mu(2), 0.5*cos(mu(3)), -0.5*sin(mu(3)));
    quiver(mu(1), mu(2), 0.5*sin(mu(3)), 0.5*cos(mu(3)));
    
    
    hold off
    
end


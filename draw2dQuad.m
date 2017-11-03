function draw2dQuad(pos_function, theta_function, synthetic_z_func, t_eval)
    
    syms t

    theta = double(subs(theta_function, t, t_eval));
    pos = double(subs(pos_function, t, t_eval));
    
    R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    
    z = double(subs(synthetic_z_func, t, t_eval))
    
    clf;
    hold on
    
    axis([-8 8 0 8 0 1])    
    view(0,90)
    
    
    
    quiver(pos(1), pos(2), cos(theta), -sin(theta));
    quiver(pos(1), pos(2), sin(theta), cos(theta));
    
    quiver(pos(1), pos(2), -z(1)*sin(theta), -z(1)*cos(theta));
    
    hold off
    drawnow
end


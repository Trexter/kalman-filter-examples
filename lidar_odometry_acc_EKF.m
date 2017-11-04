% a simulated EKF implementation where we estimate the state of a 2d
% quadrotor flying over a planar surface measuring its angle, bodyframe
% velocities, body frame accelerations, and range measurements from a rigidly mounted 1-d lidar

% this is a very realistic/common EKF implementation
syms x z theta b_dx b_dz b_dtheta b_ax b_az bias_ax bias_az t

G = 9.8;

dt = 0.02; % this is the time increments between measurements

gt_accel_bias_x = 0.1;
gt_accel_bias_z = 0.3;

symbolic_state = [x; z; theta; b_dx; b_dz; b_dtheta; b_ax; b_az; bias_ax; bias_az];

% form a parametric equation of motion to compute the synthetic sensor
% measurments from
loop_time = 3;

ground_truth_pos = [4*cos((t/loop_time)*2*pi);
    sin((t/loop_time)*2*pi) + 2];

% the angle of this fake quad can be theoretically computed from the
% direction of its acceleration
ground_truth_vel = diff(ground_truth_pos, t);
% add gravity to the y accel
ground_truth_accel = [diff(ground_truth_vel(1), t);
                       diff(ground_truth_vel(2), t) + G];
                   
% finally compute this vector's angle
u = [ground_truth_accel(1);0;ground_truth_accel(2)];
u = u/norm(u);
v = [0; 0;1];
ground_truth_theta = sign(u(1)-v(1))*acos(dot(u, v));


% PRECOMPUTE THE KALMAN FILTER FUNCTIONS/JACOBIANS

% precompute the kalman filter update and process functions
f_func = [x + cos(theta)*(dt*b_dx + 0.5*dt^2*b_ax) + sin(theta)*(dt*b_dz + 0.5*dt^2*b_az);
        z - sin(theta)*(dt*b_dx + 0.5*dt^2*b_ax) + cos(theta)*(dt*b_dz + 0.5*dt^2*b_az);
        theta + dt*b_dtheta;
        b_dx + dt*b_ax;
        b_dz + dt*b_az;
        b_dtheta;
        b_ax;
        b_az;
        bias_ax;
        bias_az]; % predict the state forward in time
    
F_jaco = jacobian(f_func, symbolic_state); % describes how uncertainty propagates

Q = eye(10); % prediction/innovation uncertainty
Q(1, 1) = 2*dt;
Q(2, 2) = 2*dt;
Q(3, 3) = 2*dt;
Q(4, 4) = 2*dt;
Q(5, 5) = 2*dt;
Q(6, 6) = 2*dt;
Q(7, 7) = 4*dt;
Q(8, 8) = 4*dt;
Q(9, 9) = 0.001*dt;
Q(10, 10) = 0.001*dt;

% synthetic measurement functions
%[range, b_dx, b_dz, b_ax, b_az]
% range is measured from the origin in the -z direction

% ground truth range function
lidar_dir_gt = ([cos(ground_truth_theta), -sin(ground_truth_theta); sin(ground_truth_theta), cos(ground_truth_theta)]*[0;1]);
range_gt = norm((ground_truth_pos(2)/lidar_dir_gt(2)) * lidar_dir_gt);

accelerometer_gt = [cos(ground_truth_theta), -sin(ground_truth_theta); sin(ground_truth_theta), cos(ground_truth_theta)] * ground_truth_accel + [gt_accel_bias_x;gt_accel_bias_z];

odometer_gt = [cos(ground_truth_theta), -sin(ground_truth_theta); sin(ground_truth_theta), cos(ground_truth_theta)] * ground_truth_vel;
odometer_gt(2, 1) = odometer_gt(2, 1)
odometer_gt(3, 1) = diff(ground_truth_theta, t)

synthetic_z_func = [range_gt;
                    accelerometer_gt(1);
                    accelerometer_gt(2);
                    odometer_gt];
                
% observation function - this can be played with
% given your predicted state what measurement would you expect from your
% sensors
h_func = []

% RUN THE FILTER

%initialize the state
mu = [double(subs(ground_truth_pos(1), t, 0));
    double(subs(ground_truth_pos(2), t, 0));
    double(subs(ground_truth_theta, t, 0));
    0;
    10;
    0;
    0;
    0;
    gt_accel_bias_x;
    gt_accel_bias_z];

Sigma = eye(10);
Sigma(1, 1) = 0.25;
Sigma(2, 2) = 0.25;
Sigma(3, 3) = 0.25;
Sigma(4, 4) = 100;
Sigma(5, 5) = 100;
Sigma(6, 6) = 100;
Sigma(7, 7) = 100;
Sigma(8, 8) = 100;
Sigma(9, 9) = 0.1;
Sigma(10, 10) = 0.1;

for time = (0:dt:20)
    %clear the figure
    clf;
    %draw ground truth of quad
    subplot(2, 4, [1, 2])
    draw2dQuad(ground_truth_pos, ground_truth_theta, synthetic_z_func, time)
    
    %compute the F jacobian
    F = double(subs(F_jaco, symbolic_state, mu));
    Sigma = F*Sigma*F' + Q; % propagate the uncertainty
    
    mu = double(subs(f_func, symbolic_state, mu)) % predict forward
    
    %draw a gaussian of the prediction
    subplot(2, 4, [5, 6])
    drawGaussian2D(Sigma(1:2, 1:2), mu(1:2))
    
    %draw our sigma with an image
    subplot(2, 4, 3)
    image(Sigma)
    
    drawnow;
end


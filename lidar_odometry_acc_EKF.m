% a simulated EKF implementation where we estimate the state of a 2d
% quadrotor flying over a planar surface measuring its angle, bodyframe
% velocities, body frame accelerations, and range measurements from a rigidly mounted 1-d lidar

% this is a very realistic/common EKF implementation
syms x z theta b_dx b_dz b_dtheta b_ax b_az bias_ax bias_az t

G = 9.8;

dt = 0.02; % this is the time increments between measurements

gt_accel_bias_x = 0.1;
gt_accel_bias_z = 0.3;

mu = [x; z; theta; b_dx; b_dz; b_dtheta; b_ax; b_az; bias_ax; bias_az];

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
f_func = [x + cos(theta)*(dt*b_dx + 0.5*dt^2*b_ax) - sin(theta)*(dt*b_dz + 0.5*dt^2*b_az);
        z + sin(theta)*(dt*b_dx + 0.5*dt^2*b_ax) + cos(theta)*(dt*b_dz + 0.5*dt^2*b_az);
        theta + dt*b_dtheta;
        b_dx + dt*b_ax;
        b_dz + dt*b_az;
        b_dtheta;
        b_ax;
        b_az;
        bias_ax;
        bias_az];
    
F_jaco = jacobian(f_func, mu);



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

% RUN THE FILTER

for time = (0:dt:20)
    draw2dQuad(ground_truth_pos, ground_truth_theta, synthetic_z_func, time)
end


% a simulated EKF implementation where we estimate the state of a 2d
% quadrotor flying over a planar surface measuring its angle, bodyframe
% velocities, body frame accelerations, and range measurements from a rigidly mounted 1-d lidar

% this is a very realistic/common EKF implementation
syms x z theta b_dx b_dz b_dtheta b_ax b_az bias_ax bias_az t

G = 9.8;

dt = 0.1; % this is the time increments between measurements

gt_accel_bias_x = 0.1;
gt_accel_bias_z = 0.3;

%landmark pos at orgin

symbolic_state = [x; z; theta; b_dx; b_dz; b_dtheta; b_ax; b_az; bias_ax; bias_az];

% form a parametric equation of motion to compute the synthetic sensor
% measurments from
loop_time = 7;

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
        bias_az] % predict the state forward in time
    
F_jaco = jacobian(f_func, symbolic_state) % describes how uncertainty propagates

Q = eye(10); % prediction/innovation uncertainty
Q(1, 1) = 0.003*dt;
Q(2, 2) = 0.003*dt;
Q(3, 3) = 0.003*dt;
Q(4, 4) = 2*dt;
Q(5, 5) = 2*dt;
Q(6, 6) = 2*dt;
Q(7, 7) = 2*dt;
Q(8, 8) = 2*dt;
Q(9, 9) = 0.005*dt;
Q(10, 10) = 0.005*dt;

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

% create the gt of the image sensor detecting the landmark
%TODO

synthetic_z_func = [range_gt;
                    accelerometer_gt(1);
                    accelerometer_gt(2);
                    odometer_gt];
                
% observation function - this can be played with
% given your predicted state what measurement would you expect from your
% sensors
lidar_dir = ([cos(theta), -sin(theta); sin(theta), cos(theta)]*[0;1]);
range_obs_function = norm((z/lidar_dir(2)) * lidar_dir);
accel_obs_function = [b_ax + bias_ax;
                      b_az + bias_az;] + [cos(theta), -sin(theta); sin(theta), cos(theta)] * [0;G];

h_func = [range_obs_function;
         accel_obs_function;
         b_dx;
         b_dz;
         b_dtheta] % a function which maps the state into an expected measurement
     
H_jaco = jacobian(h_func, symbolic_state) % describes how a measurement can update the current state

% this matrix describes the uncertainties in each measurement
R = eye(6);
R(1, 1) = 0.05;
R(2, 2) = 0.02;
R(3, 3) = 0.02;
R(4, 4) = 0.003;
R(5, 5) = 0.003;
R(6, 6) = 0.003;

% RUN THE FILTER

%initialize the state
mu = [double(subs(ground_truth_pos(1), t, 0));
    1;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0];

Sigma = eye(10);
Sigma(1, 1) = 0.25;
Sigma(2, 2) = 10000;
Sigma(3, 3) = 10000;
Sigma(4, 4) = 10000;
Sigma(5, 5) = 10000;
Sigma(6, 6) = 10000;
Sigma(7, 7) = 10000;
Sigma(8, 8) = 10000;
Sigma(9, 9) = 100000;
Sigma(10, 10) = 100000;

for time = (0:dt:200)
    %clear the figure
    clf;
    
    %compute the F jacobian
    F = double(subs(F_jaco, symbolic_state, mu));
    Sigma = F*Sigma*F' + Q; % propagate the uncertainty
    
    mu = double(subs(f_func, symbolic_state, mu)); % predict forward
    
    % start the update procedure
    
    %compute the H jacobian
    H = double(subs(H_jaco, symbolic_state, mu));
    
    measurement = double(subs(synthetic_z_func, t, time));
    
    % add noise to measurement
    noise = diag(R) .* randn(size(measurement));
    measurement = measurement + noise;
    
    expected_measurement = double(subs(h_func, symbolic_state, mu));
    y = measurement - expected_measurement;
    S = R + H*Sigma*H';
    K = Sigma*H'/(S);
    
    mu = mu + K*y;
    
    % use Josephs form: P = ( I ? KH) P (I ? KH)' + KRK'
    I_KH = (eye(size(Sigma)) - K*H);
    Sigma = I_KH*Sigma*I_KH' + K*R*K'
    
    %draw a gaussian of the position
    subplot(2, 4, [5, 6])
    title('Position Estimate')
    hold on;
    plot_gaussian_ellipsoid(Sigma(1:2, 1:2), mu(1:2), 0.5)
    plot_gaussian_ellipsoid(Sigma(1:2, 1:2), mu(1:2), 1)
    plot_gaussian_ellipsoid(Sigma(1:2, 1:2), mu(1:2), 2)
    plot_gaussian_ellipsoid(Sigma(1:2, 1:2), mu(1:2), 10)
    xlabel('X')
    ylabel('Z')
    
    %draw a gaussian of the body frame acceleration
    subplot(2, 4, 7)
    title('Body Frame Acceleration Estimate')
    hold on;
    plot_gaussian_ellipsoid(Sigma(7:8, 7:8), mu(7:8), 0.5)
    plot_gaussian_ellipsoid(Sigma(7:8, 7:8), mu(7:8), 1)
    plot_gaussian_ellipsoid(Sigma(7:8, 7:8), mu(7:8), 2)
    plot_gaussian_ellipsoid(Sigma(7:8, 7:8), mu(7:8), 10)
    xlabel('ax_b_o_d_y (m/s^2)')
    ylabel('az_b_o_d_y (m/s^2)')
    
    axis([-10 10 -10 10 0 1])
    view(0,90)
    pbaspect([1 1 1])
    
    %draw a gaussian of the bias estimate
    subplot(2, 4, 8)
    title('Accelerometer Bias Estimate')
    hold on;
    plot_gaussian_ellipsoid(Sigma(9:10, 9:10), mu(9:10), 0.5)
    plot_gaussian_ellipsoid(Sigma(9:10, 9:10), mu(9:10), 1)
    plot_gaussian_ellipsoid(Sigma(9:10, 9:10), mu(9:10), 2)
    plot_gaussian_ellipsoid(Sigma(9:10, 9:10), mu(9:10), 10)
    xlabel('Bias_x (m/s^2)')
    ylabel('Bias_z (m/s^2)')
    
    axis([-3 3 -3 3 0 1])
    view(0,90)
    pbaspect([1 1 1])
    
    %draw our sigma with an image
    subplot(2, 4, 3)
    imagesc(Sigma)
    title('Full Covariance Matrix')
    colorbar;
    
    %draw the Kalman gain
    subplot(2, 4, 4)
    imagesc(K)
    title('Kalman Gain')
    colorbar;
    
     %draw ground truth of quad
    subplot(2, 4, [1, 2])
    title('Ground Truth vs EKF Estimate')
    draw2dQuad(ground_truth_pos, ground_truth_theta, synthetic_z_func, mu, time)
    
    % plot the landmark
    plot(0,0,'r*')
    
    %plot the image plane with detection
    %TODO
    
    xlabel('X')
    ylabel('Z')
    
    drawnow;
end


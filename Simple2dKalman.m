
% [x, y, dx, dy]
x = [0; 0; 0; 0;]

Sigma = [0.1, 0, 0, 0; 0, 0.1, 0, 0; 0, 0, 2, 0; 0, 0, 0, 2]

dt = 0.5;

% State Transistion
F = [1, 0, dt, 0;
    0, 1, 0, dt;
    0, 0, 1, 0;
    0, 0, 0, 1]

% Measurement function
H = [1, 0, 0, 0;
    0, 0, 0, 1];

% prediction uncertainty 
Q = [(2*dt)^2, 0, 0, 0;
    0, (2*dt)^2, 0, 0;
    0, 0, dt, 0;
    0, 0, 0, dt]

% measurement uncertainty
R = [0.1, 0;
    0, 0.01]

% [x, dy]
z = [1, 2, 3, 4, 5, 6, 7, 8;
    0.5, 0.5, 0.5, 1, 1, 2, 2, 2] % measurement 

for index = (1:8)
% run through the prediction and update
x = F*x;
Sigma = F*Sigma*F' + Q;

y = z(1:2, index) - H*x;
S = H*Sigma*H' + R;
K = Sigma*H'*inv(S);
x = x + K*y
Sigma = (eye(4, 4) - K*H)*Sigma

drawGaussian2D(Sigma(1:2, 1:2), x(1:2))

pause(1)

end



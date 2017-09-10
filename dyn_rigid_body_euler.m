function [dX] = dyn_rigid_body_euler(t, X)
%DYN_RIGI_BODY3 Calculate the values of w_dot with respect to to time for
%Euler's 3-1-3 rotation sequence.

global JP;

% Get input parameters
W = X(1:3,1);

% Get Angles.
psi = X(4,1);
theta = X(5,1);
phi = X(6,1);

Wx = [0 -W(3) W(2);...
    W(3) 0 -W(1);...
    -W(2) W(1) 0];

dW = -JP^-1*(Wx*JP*W);

psi_dot = (W(1) * sin(phi) + W(2) * cos(phi)) / sin(theta);
theta_dot = W(1) * cos(phi) - W(2) * sin(phi);
%phi_dot = W(3) - (sin(phi) * cos(theta) * W(1) + cos(phi) * cos(theta) * W(2)) / sin(theta);
phi_dot = W(3) - (W(1) * sin(phi) + W(2) * cos(phi)) / tan(theta);

dX = [dW; psi_dot; theta_dot; phi_dot];

end


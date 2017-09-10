function [dX] = dyn_rigid_bryant(t,X)
%DYN_RIGI_BODY3 Calculate the values of w_dot with respect to to time for
%Bryant's 3-2-1 rotation sequence.

global JP;

W = X(1:3,1);

% Get Bryant's angles.
psi_b = X(4,1);
theta_b = X(5,1);
phi_b = X(6,1);

Wx = [0 -W(3) W(2);...
    W(3) 0 -W(1);...
    -W(2) W(1) 0];

dW = -JP^-1*(Wx*JP*W);

% For Bryant: WHYWHWYWHWYHWYWHWYWHYWHWYWYHWY
%phi_dot_b   = W(1) + (W(2) * sin(phi_b) + W(2) * cos(phi_b)) * tan(theta_b);
%theta_dot_b = W(2) * cos(phi_b) - W(3) * sin(phi_b);
%psi_dot_b   = (W(2) * sin(phi_b) + W(3) * cos(phi_b)) / cos(theta_b);

psi_dot_b   = (W(2) * sin(phi_b) + W(3) * cos(phi_b)) / cos(theta_b);
theta_dot_b = W(2) * cos(phi_b) - W(3) * sin(phi_b);
phi_dot_b   = W(1) + ((W(2) * sin(phi_b) + W(2) * cos(phi_b)) * tan (theta_b));
%phi_dot_b = W(1)+W(2)*sin(phi_b)*sin(theta_b)+W(3)*cos(phi_b)*sin(theta_b);

dX = [dW; psi_dot_b; theta_dot_b; phi_dot_b];

end
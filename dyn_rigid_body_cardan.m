function [dX] = dyn_rigid_body_cardan(t,X)
%DYN_RIGI_BODY3 Calculate the values of w_dot with respect to to time for
%Cardan's 1-2-4 rotation sequence.

global JP;

W = X(1:3,1);

% Get Cardan's angles.
psi_c = X(4,1);
theta_c = X(5,1);
phi_c = X(6,1);

Wx = [0 -W(3) W(2);...
    W(3) 0 -W(1);...
    -W(2) W(1) 0];

dW = -JP^-1*(Wx*JP*W);

% For Cardan:
psi_dot_c = W(3) - (W(1) * cos(psi_c) - W(2) * sin(psi_c)) * tan(phi_c); % Spin
theta_dot_c = (W(1) * cos(psi_c) - W(2) * sin(psi_c)) / cos(phi_c); % Nutation
phi_dot_c = W(1) * sin(psi_c) - W(2) * cos(psi_c); % Precession


dX= [dW; psi_dot_c; theta_dot_c; phi_dot_c];

end
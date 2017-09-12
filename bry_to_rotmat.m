function [ mb ] = bry_to_rotmat(psi_b, theta_b, phi_b)
%BRY_TO_ROTMAT converts Bryant Angles to a rotation matrix.

% First row of the rotation matrix.
mb1 = cos(theta_b) * cos(psi_b);
mb2 = cos(theta_b) * sin (psi_b);
mb3 = -sin(theta_b);

% Second row of the rotation matrix.
mb4 = -cos(phi_b) * sin (psi_b) + sin(phi_b) * sin (theta_b) * cos(psi_b);
mb5 = cos(psi_b) * cos (phi_b) + sin(psi_b) * sin(theta_b) * sin(phi_b);
mb6 = sin(phi_b) * cos(theta_b);

% Third row of the rotation matrix.
mb7 = sin(psi_b) * sin (phi_b) + cos(psi_b) * sin(theta_b) * cos(phi_b);
mb8 = sin(phi_b) * cos(psi_b) + cos(theta_b) * sin(theta_b) * sin(psi_b);
mb9 = cos(phi_b) * cos(theta_b);

mb = [mb1, mb2, mb3;...
    mb4, mb5, mb6;...
    mb7, mb8, mb9];

% mb = [cos(theta_b) * cos(psi_b), cos(theta_b) * sin (psi_b), -sin(theta_b);...
%    -cos(phi_b) * sin (psi_b) + sin(phi_b) * sin (theta_b) * cos(psi_b), cos(psi_b) * cos (phi_b) + sin(psi_b) * sin(theta_b) * sin(phi_b), sin(phi_b) * cos(theta_b);...
%    sin(psi_b) * sin (phi_b) + cos(psi_b) * sin(theta_b) * cos(phi_b), sin(phi_b) * cos(psi_b) + cos(theta_b) * sin(theta_b) * sin(psi_b), cos(phi_b) * cos(theta_b) ]
% 

end



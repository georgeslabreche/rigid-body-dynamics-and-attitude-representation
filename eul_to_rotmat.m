function [ m ] = eul_to_rotmat(psi, theta, phi)
%EUL_TO_ROTMAT converts Euler Angles to a rotation matrix.

% First row of the rotation matrix.
m1 = cos(psi) * cos(phi) - cos(theta) * sin(phi) * sin(psi);
m2 = cos(psi) * sin(phi) + cos(theta) * cos(phi) * sin(psi);
m3 = sin(theta) * sin(psi);

% Second row of the rotation matrix.
m4 = -sin(psi) * cos(phi) - cos(theta) * sin(phi) * cos(psi);
m5 = -sin(psi) * sin(phi) + cos(theta) * cos(phi) * cos(psi);
m6 = sin(theta) * cos(psi);

% Third row of the rotation matrix.
m7 = sin(phi) * sin(theta);
m8 = -cos(phi) * sin(theta);
m9 = cos(theta);

m = [m1, m2, m3;...
    m4, m5, m6;...
    m7, m8, m9];

end


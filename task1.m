clear all 
close all
clc

global JJ

% Test inertia matrix
JJ = [100, -20, -100;...
      -20, 300, -50;...
      -100, -50, 500];
  
  
% Values for the inertia matrix.
J11 = 4.011760901487551;
J12 = 0.003212192682296;
J13 = 0.107760217621921;
J22 = 4.000877329160449;
J23 = 0.029431977034595;
J33 = 4.987361769351999;

% Inertia Matrix.
JJ = [J11, -J12, J13;...
      -J12, J22, -J23;...
      J13, -J23, J33];
 
JJ = [100, -20, -100;...
      -20, 300, -50;...
      -100, -50, 500];

% The principal moment of inertia of the body is represented by the eigan 
% values. First we determine what those eigen valyes are by solving lambda:
sol_lambda = get_principal_moments_of_inertia(JJ);

% For each lambda, calculate the unit vector and extract its 
% coefficients/components:
%  - A unit vector defines a principal director of the intertia tensor.
%  - The component of a unit vector, will define the rows
%    of the orthogonal transformation q from the xyz system into the x'y'z' 
%    system aligned along the three principal directions.
%
% Curtis, Howard. "Moments of Interia" in Orbital Mechanics: For
% Engineering Students, 516.
unit_vect1 = get_unit_vect(JJ, sol_lambda(1)); % For lambda1 532.052.
coeffs1 = coeffs(unit_vect1);

unit_vect2 = get_unit_vect(JJ, sol_lambda(2)); % For lambda2 72.1083.
coeffs2 = coeffs(unit_vect2);

unit_vect3 = get_unit_vect(JJ, sol_lambda(3)); % For lambda3 295.840.
coeffs3 = coeffs(unit_vect3);

% TODOD: check our work?

% Build the orthogonal transformation q from the xyz system into the x'y'z' 
% system aligned along the three principal directions.
q = [coeffs1(1,3), coeffs1(1,2), coeffs1(1,1);... % Coefficients of lambda1 532.052
    coeffs2(1,3), coeffs2(1,2), coeffs2(1,1);... % Coefficients of lambda2 72.1083
    coeffs3(1,3), coeffs3(1,2), coeffs3(1,1)]; % Coefficients of lambda3 295.840
single(q)

% Apply the transformation/
result = q*JJ*transpose(q);
single(result)

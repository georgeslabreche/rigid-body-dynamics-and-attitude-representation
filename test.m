clear all 
close all
clc

global JJ

JJ = [100, -20, -100;...
      -20, 300, -50;...
      -100, -50, 500];
  
  
% Values for the inertia matrix
J11 = 4.011760901487551;
J12 = 0.003212192682296;
J13 = 0.107760217621921;
J22 = 4.000877329160449;
J23 = 0.029431977034595;
J33 = 4.987361769351999;

% Inertia Matrix
JJ = [J11, -J12, J13;...
      -J12, J22, -J23;...
      J13, -J23, J33];

JJ = [100, -20, -100;...
  -20, 300, -50;...
  -100, -50, 500];

% The principal moment of inertia of the body is represented by the eigan 
% values. First we determine what those eigen valyes are by solving lambda:
syms lambda
jj_lambda_eye = JJ-lambda*eye(3);
det_fn = det(jj_lambda_eye) == 0;
sol_lambda = solve(det_fn, lambda);

% Print the resulting eigan values
lambda1 = sol_lambda(1); % 532.05
lambda2 = sol_lambda(3); % 295.84
lambda3 = sol_lambda(2); % 72.11

% Subsitute these back in JJ-lambda*eye(3)
lambda1_subbed = subs(jj_lambda_eye, lambda, lambda1);
unit_vect1 = get_unit_vect(lambda1_subbed);
coeffs1 = coeffs(unit_vect1);

testo = [100 - lambda1, -20, -100;...
        -20, 300 - lambda1, -50;...
        -100, -50, 500 - lambda1],
 

% Let's check our work:
%check1 = (JJ-lambda1*eye(3))*[coeffs1(1,3); coeffs1(1,2); coeffs1(1,1)];
%double(check1)

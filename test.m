clear all 
close all
clc

global JJ

%% Input
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
      J13, -J23, J33]

% The principal moment of inertia of the body
% is represented by the eigan values and
% Z is the matrix that diagonalizes JG in JP:
[Z JP] = eig(JJ)

% Verify the result:
check = JJ==Z*JP*transpose(Z)

syms lambda
detFn = det(JJ-lambda*eye(3)) == 0;
solLambda = solve(detFn, lambda);
solLambda

fzero(@(lambda) det(JJ-lambda*eye(3)), 1)



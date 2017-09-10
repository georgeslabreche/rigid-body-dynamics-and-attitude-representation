%%-------------------------------------------------------------------------
%     R0008R  Introduction to Space Mechanics and Electronics
%                   Lulea University of Technology
%                       Teacher: Leonard Felicetti  
%--------------------------------------------------------------------------
%  Assignment 1: Rigid Body Dynamics and Attitude Representation
%%-------------------------------------------------------------------------
%                       Students: XXX and YYY
%%-------------------------------------------------------------------------


clear all 
close all
clc

global JJ JP

%% Input

% Inertia Matrix
 J11 = 4.011760901487551;
  J22 = 4.000877329160449;
  J33 = 4.987361769351999;
  J12 = 0.003212192682296;
  J13 = -0.107760217621921;
  J23 = 0.029431977034595;
  

JJ = [J11 -J12 -J13;...
      -J12 J22 -J23;...
      -J13 -J23 J33]; %kg.

  [Z, JP] = eig(JJ);
  
% Initial Angular Velocity
w_1 = [0.1; 0.1; 1.0]; %rpm
w_2 = [0.1; 1.0; 0.1]; %rpm
w_3 = [1.0; 0.1; 0.1]; %rpm

W_0 = [w_1]*(2*pi/60); %rad/s

% Inertial reference frame with e3=h0
H_0 = JP*W_0;

% Initial condtions for Euler's angles (3-1-3 rotation sequence).
psi_e_0 = 0.0; 
theta_e_0 = acos(H_0(3) / norm(H_0)); 
phi_e_0 = atan2(H_0(1) / (norm(H_0) * sin(theta_e_0)), H_0(2) / (norm(H_0) * sin(theta_e_0))); 
euler_0 = [psi_e_0; theta_e_0; phi_e_0]; % rad/s

% Initial conditions for Bryant's angles (3-2-1 rotation sequence).
psi_b_0 = 0.0;
theta_b_0 = asin(-H_0(1) / norm(H_0)); 
phi_b_0 = atan2(H_0(2) / (norm(H_0) * cos(theta_b_0)), H_0(3) / (norm(H_0) * cos(theta_b_0)));
bryant_0 = [psi_b_0; theta_b_0; phi_b_0]; % rad/s

% Initial conditions for Cardan's angles (1-2-3 rotation sequence)
theta_c_0 = asin(H_0(1) / norm(H_0));
psi_c_0 = atan2(-H_0(2) / (norm(H_0) * cos(theta_c_0)), H_0(3) / (norm(H_0) * cos(theta_c_0))); 
phi_c_0 = 0.0;
cardan_0 = [psi_c_0; theta_c_0; phi_c_0]; % rad/s

% Initial conditions vector
X0_e = [W_0; euler_0];
X0_b = [W_0; bryant_0];
X0_c = [W_0; cardan_0];

% Time
t_0 = 0;
t_f = 5*60; % minutes 
dt = 0.1; % ? 

tspan = [t_0 : dt : t_f];

% Set options for ode45
options=odeset('AbsTol', 1.e-10,'RelTol', 1.e-10);

%% EULER: 3-1-3 Rotation Sequence
[time,X]=ode45('dyn_rigid_body_euler', tspan, X0_e, options);
angle_plot('Euler Angles: 3-1-3 Rotation Sequence', time, X);
W = X(:,1:3);

%% BRYANT
[time,X]=ode45('dyn_rigid_body_bryant', tspan, X0_b, options);
angle_plot('Bryant Angles: 3-2-1 Rotation Sequence', time, X);

%% CARDAN
[time,X]=ode45('dyn_rigid_body_cardan', tspan, X0_c, options);
angle_plot('Cardan Angles: 1-2-3 Rotation Sequence', time, X);


%% ANGULAR VELOCITIES
figure('Name', 'Angular Velocities');
plot(time, W * 60/(2*pi));
xlabel('time [s]');
ylabel('angular velocity [rpm]');
legend('w_1', 'w_2', 'w_3');


 
 
 
 

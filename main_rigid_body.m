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

  %[Z, JP] = eig(JJ);
  JP = [3 0 0;...
      0 4 0;...
      0 0 5]
  
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
dt = 1; % ? 

tspan = [t_0 : dt : t_f];

% Set options for ode45
options=odeset('AbsTol', 1.e-10,'RelTol', 1.e-10);

%% EULER: 3-1-3 Rotation Sequence
[time,X]=ode45('dyn_rigid_body_euler', tspan, X0_e, options);
angle_plot('Euler Angles: 3-1-3 Rotation Sequence', time, X);

% save this for later when we will plot the angular velocities.
W = X(:,1:3);

% Save these for later when we will plot the cones.
psi = X(:,4);
theta = X(:,5);
phi = X(:,6);


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


%% task5

% acc = w * (w * r)
ra = [0.5 0 0]';
rb = [0 0.5 0]';
rc = [0 0 0.5]';
va = zeros(length(time),3);
vb = zeros(length(time),3);
vc = zeros(length(time),3);
aa = zeros(length(time),3);
ab = zeros(length(time),3);
ac = zeros(length(time),3);

for t_time = 1 : length(time)
    va(t_time,:) = cross(W(t_time,:),ra);
    vb(t_time,:) = cross(W(t_time,:),rb);
    vc(t_time,:) = cross(W(t_time,:),rc);
end

for t_time = 1 : length(time)
    aa(t_time,:) = cross(W(t_time,:),va(t_time,:));
    ab(t_time,:) = cross(W(t_time,:),vb(t_time,:));
    ac(t_time,:) = cross(W(t_time,:),vc(t_time,:));
end
    
% accelerations
figure('Name','Acceleration')
subplot(1,3,1)
plot(time,aa)
xlabel('time [s]')
ylabel('accerlation')

subplot(1,3,2)
plot(time,ab)
xlabel('time [s]')
ylabel('accerlation')

subplot(1,3,3)
plot(time,ac)
xlabel('time [s]')
ylabel('accerlation')

    
%% BODY AND SPACE CONES
% other quantities
RR = zeros(3,3,length(time));
T = zeros(length(time),1);

b1_b = [1 0 0]';
b2_b = [0 1 0]';
b3_b = [0 0 1]';
e1_e = [1 0 0]';
e2_e = [0 1 0]';
e3_e = [0 0 1]';
    
b1_e = zeros(3,length(time));
b2_e = zeros(3,length(time));
b3_e = zeros(3,length(time));

for t_time = 1 : length(time)
   RR(:,:,t_time) = eul_to_rotmat(psi(t_time), theta(t_time), phi(t_time)); 
   T(t_time,1) = 1/2 * W(t_time,:) * JP * W(t_time,:)';
   
   % The body reference frame wrt inertial frame in time:
   b1_e(:,t_time) = RR(:,:,t_time) * b1_b;
   b2_e(:,t_time) = RR(:,:,t_time) * b2_b;
   b3_e(:,t_time) = RR(:,:,t_time) * b3_b;
   
end


figure('Name','Kinetic Energy');
plot(time,T);
xlabel('time [s]');
ylabel('Kinetic Energy [T]');


%  figure('Name','trajectory of   b3   wrt inertial reference frame ')
%
%  figure('Name','trajectory of   W  wrt body reference frame')
%
%  figure('Name','trajectory of   W  wrt inertial reference frame')

%  animation
figure('Name','Movie')


for t_time =  1:length(time)
    hold on
    
    %  body axes wrt inertial reference frame
    quiver3(0,0,0,b1_e(1,t_time),b1_e(2,t_time),b1_e(3,t_time),'r');
    quiver3(0,0,0,b2_e(1,t_time),b2_e(2,t_time),b2_e(3,t_time),'g');
    quiver3(0,0,0,b3_e(1,t_time),b3_e(2,t_time),b3_e(3,t_time),'b');
    
    %  angular momentum wrt inertial reference frame
    
    xlim([-1.,1.5]);
    ylim([-1.,1.5]);
    zlim([-0.5,1.5]);
    
    if t_time ==1
        grid
        quiver3(0,0,0,1.5*e3_e(1,1),1.5*e3_e(2,1),1.5*e3_e(3,1),'k');
        quiver3(0,0,0,1.5*e1_e(1,1),1.5*e1_e(2,1),1.5*e1_e(3,1),'k');
        quiver3(0,0,0,1.5*e2_e(1,1),1.5*e2_e(2,1),1.5*e2_e(3,1),'k');
        keyboard
    end
    
    pause(0.1);
     % hold off
     % plot3(0,0,0,'b');
      
end
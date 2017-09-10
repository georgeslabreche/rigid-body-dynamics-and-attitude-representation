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
JJ = [J11, -J12, -J13;...
      -J12, J22, -J23;...
      J13, -J23, J33];

% Initial Angular Velocity

W_0 = 

% Inertial reference frame with e3=h0
H_0 = JJ*W_0;

psi_0 = 0.0; 
theta_0 = acos(H_0(3)/norm(H_0)); 
phi_0 = atan2(H_0(1)/(norm(H_0)*sin(theta_0)),H_0(2)/(norm(H_0)*sin(theta_0))); 
Euler_0 = [psi_0; theta_0; phi_0]; % rad

% Initial conditions vector

X0 = 


% Time
t_0 = 
t_f =  
dt =  

tspan = [t_0 : dt : t_f];

%% Processing

options=odeset('AbsTol',1.e-10,'RelTol',1.e-10);
[time,X]=ode45('dyn_rigid_body', tspan, X0, options);

%% Post-Processing

W = 
Euler = 

psi = 
theta = 
phi = 


%% Figures and Results

% angular velocity
figure('Name','Angular Velocities')
plot(time,W * 60/(2*pi))
xlabel('time [s]')
ylabel('angular velocity [rpm]')
legend('w_1','w_2','w_3')


% euler angles 313
figure('Name','Euler Angles 313')
subplot(1,3,1)
plot(time,psi * 180/pi)
xlabel('time [s]')
ylabel('precession angle [deg]')

subplot(1,3,2)
plot(time,theta * 180/pi)
xlabel('time [s]')
ylabel('nutation angle [deg]')

subplot(1,3,3)
plot(time,phi * 180/pi)
xlabel('time [s]')
ylabel('spin angle [deg]')


% other quantities
RR = zeros(3,3,length(time));
T = zeros(length(time),1);



b1_b = [1 0 0]';
b2_b = [0 1 0]';
b3_b = [0 0 1]';

b1_e = zeros(3,length(time));
b2_e = zeros(3,length(time));
b3_e = zeros(3,length(time));


for t_time = 1 : length(time)
   RR(:,:,t_time) = Eul_2_RR(psi(t_time), theta(t_time), phi(t_time)); 
   T(t_time,1) = ;
   
   
   b1_e(:,t_time) = 
   b2_e(:,t_time) = 
   b3_e(:,t_time) = 
   
end


figure('Name','Kinetic Energy')


% figure('Name','trajectory of b3 wrt inertial reference frame ')
% 
% figure('Name','trajectory of W wrt body reference frame')
% 
% figure('Name','trajectory of W wrt inertial reference frame')

% animation
figure('Name','Movie')


for t_time = 1:length(time)
    hold on
    
    % body axes wrt inertial reference frame
    quiver3(0,0,0,b1_e(1,t_time),b1_e(2,t_time),b1_e(3,t_time),'r');
    quiver3(0,0,0,b2_e(1,t_time),b2_e(2,t_time),b2_e(3,t_time),'g');
    quiver3(0,0,0,b3_e(1,t_time),b3_e(2,t_time),b3_e(3,t_time),'b');
    
    % angular momentum wrt inertial reference frame
    
    
    xlim([-1.,1.5]);
    ylim([-1.,1.5]);
    zlim([-0.5,1.5]);
    
    if t_time ==1
        
        grid
        quiver3(0,0,0,1.5*e3_e(1,1),1.5*e3_e(2,1),1.5*e3_e(3,1),'k')
        quiver3(0,0,0,1.5*e1_e(1,1),1.5*e1_e(2,1),1.5*e1_e(3,1),'k')
        quiver3(0,0,0,1.5*e2_e(1,1),1.5*e2_e(2,1),1.5*e2_e(3,1),'k')
        keyboard
        
    end
    
    pause(0.1)
%     hold off
%     plot3(0,0,0,'b')
%     
    
    
end


function [] = angle_plot(title, time, X)
%ANGLE_PLOT Summary of this function goes here
%   Detailed explanation goes here
psi = wrapTo2Pi(X(:,4));
theta = wrapTo2Pi(X(:,5));
phi = wrapTo2Pi(X(:,6));

figure('Name', title)
subplot(1,3,1)
plot(time, psi * 180/pi)
xlabel('time [s]')
ylabel('Precession angle [deg]')

subplot(1,3,2)
plot(time, theta * 180/pi)
xlabel('time [s]')
ylabel('Nutation angle [deg]')

subplot(1,3,3)
plot(time, phi * 180/pi)
xlabel('time [s]')
ylabel('Spin angle [deg]')

end


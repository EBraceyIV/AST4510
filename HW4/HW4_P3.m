%Homework 4 Problem 3
clc; clear;
% %Parts A-D
%Define initial position & velocity vectors
rPCI0 = [68524.298 -17345.863 -51486.409]'; %km
vPCI0 = [-0.578936 0.957665 0.357759]'; %km*s^-1
%Define initial and terminal times (convert to seconds)
t0_m = 1329.16; t0 = t0_m*60; %s
tf_m = 3885.73; tf = tf_m*60; %s
%Define gravitational paramater
mu = 398600; %km^3*s^-2
%Define inputs to ode113
tspan = [t0 tf];
Z = [rPCI0; vPCI0]; options = odeset('RelTol',1e-6);
%Calculate the position and velocity along the orbit with ode113
[t_ode,Y] = ode113(@rv_orbit,tspan,Z,options);
rPCIf_ode = Y(length(Y),1:3)'; vPCIf_ode = Y(length(Y),4:6)';
plot3(Y(:,1),Y(:,2),Y(:,3),'-')
hold on

% %Parts E-I
%Separate time span into enough intervals for a smooth curve
time = linspace(t0,tf,400)';
%Preallocate space for the position and velocity vectors
rPCI = [rPCI0'; zeros(length(time)-1,3)];
vPCI = [vPCI0'; zeros(length(time)-1,3)];
%For each time value, calculate the position & velocity using Kepler solver
for ii = 2:length(time)
    [rPCIf,vPCIf,E0,nu0,E,nu] = propagateKepler_BraceyIV_Edgar(t0,time(ii),rPCI0,vPCI0,mu);
    %Store new position & velocity
    rPCI(ii,:) = rPCIf';
    vPCI(ii,:) = vPCIf';
end
plot3(rPCI(:,1),rPCI(:,2),rPCI(:,3),'.')
title('Part (i)'); xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]'); legend('ode113 - r(t)','Kepler - r(t)');

% %Printouts
fprintf('Part (a)\n')
fprintf('Final eccentric anomaly = %.3f degrees\n',rad2deg(E))
fprintf('\nPart (b)\n')
fprintf('Final true anomaly = %.3f degrees\n',rad2deg(nu))
fprintf('\nPart (c)\n')
fprintf('Final position (km) = [%.3f %.3f %.3f]''\n',rPCIf)
fprintf('Final velocity (km*s^-1) = [%.3f %.3f %.3f]''\n',vPCIf)
fprintf('\nPart (d)\n')
fprintf('Kepler vs ODE position [km]:\n')
disp([rPCIf rPCIf_ode])
fprintf('Kepler vs ODE velcoity [km*s^-1]:\n')
disp([vPCIf vPCIf_ode])
fprintf('\nPart (h)\n')
fprintf('Kepler position [km] & velocity [km*s^-1] array:\n')
format short e
disp([rPCI vPCI])
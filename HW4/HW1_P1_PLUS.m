%Homework 4 Problem 1
clc; clear;
options = odeset('RelTol',1e-6);
%Set t_o & t_f span
tspan_min = [1329.16 3885.73]; tspan = tspan_min/60;
%Define inital conditions, r & v, in a column vector
Z = [68524.298 -17345.863 -51486.409 -0.578936 0.957665 0.357759]';
t_all = zeros(20,1);
Y_all = zeros(20,6);
%Solve the system of ODEs using the ICs
[t,Y] = ode113(@rv_orbit,tspan,Z,options);
t_all = t; Y = Y_all;
%Plot the points of the orbit
plot3(Y(:,1),Y(:,2),Y(:,3));
title('Position on orbit'); xlabel('x'); ylabel('y'); zlabel('z'); legend('r(t)','v(t)');

figure
tspan_min = [3.93 1771.58]; tspan = tspan_min/60;
Z = [2721.965 3522.863 5267.244 9.572396 -0.474701 -2.725664]';
[t,Y] = ode113(@rv_orbit,tspan,Z,options);
%Plot the points of the orbit
plot3(Y(:,1),Y(:,2),Y(:,3));
t_all(1:length(t),2) = t; 
Y_all(1:length(Y),2) = Y;

% 
% figure
% tspan_min = [242.82 612.69]; tspan = tspan_min/60;
% Z = [6997.56 -34108.00 20765.49 0.15599 0.25517 1.80763]';
% [t,Y] = ode113(@rv_orbit,tspan,Z,options);
% %Plot the points of the orbit
% plot3(Y(:,1),Y(:,2),Y(:,3));
% 
% figure
% tspan_min = [616.79 1880.41]; tspan = tspan_min/60;
% Z = [1882.725 9864.690 4086.088 -5.565367 5.451548 2.258105]';
% [t,Y] = ode113(@rv_orbit,tspan,Z,options);
% %Plot the points of the orbit
% plot3(Y(:,1),Y(:,2),Y(:,3));
% 
% figure
% tspan_min = [21.02 1913.38]; tspan = tspan_min/60;
% Z = [1882.725 9864.690 4086.088 -5.565367 5.451548 2.258105]';
% [t,Y] = ode113(@rv_orbit,tspan,Z,options);
% %Plot the points of the orbit
% plot3(Y(:,1),Y(:,2),Y(:,3));
% 
% figure
% tspan_min = [27 57]; tspan = tspan_min/60;
% Z = [-10515.45 -5235.37 49.1700 -2.10305 -4.18146 5.56329]';
% [t,Y] = ode113(@rv_orbit,tspan,Z,options);
% %Plot the points of the orbit
% plot3(Y(:,1),Y(:,2),Y(:,3));
% 

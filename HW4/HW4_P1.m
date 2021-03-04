%Homework 4 Problem 1
clc; clear;
%Case 1
options = odeset('RelTol',1e-6);
%Set time interval as [t0 tf]  %Convert from minutes to seconds
tspan_min = [1329.16 3885.73]; tspan = tspan_min*60;
%Define inital position and velocity in a column vector
rPCI = [68524.298 -17345.863 -51486.409]';
vPCI = [ -0.578936 0.957665 0.357759]';
Z = [rPCI; vPCI];
%Solve the system of ODEs using the inital conditions
[t,Y] = ode113(@rv_orbit,tspan,Z,options);
%Plot the points of the orbit
plot3(Y(:,1),Y(:,2),Y(:,3)); hold on
title('Orbital path (Case 1)'); xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]'); legend('r(t) [km]');
% %Above comments also apply to following cases:
%Case 2
figure
tspan_min = [3.93 1771.58]; tspan = tspan_min*60;
Z = [2721.965 3522.863 5267.244 9.572396 -0.474701 -2.725664]';
[t,Y] = ode113(@rv_orbit,tspan,Z,options);
plot3(Y(:,1),Y(:,2),Y(:,3)); hold on
title('Orbital path (Case 2)'); xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]'); legend('r(t) [km]');
%Case 3
figure
tspan_min = [242.82 612.69]; tspan = tspan_min*60;
Z = [6997.56 -34108.00 20765.49 0.15599 0.25517 1.80763]';
[t,Y] = ode113(@rv_orbit,tspan,Z,options);
plot3(Y(:,1),Y(:,2),Y(:,3)); hold on
title('Orbital path (Case 3)'); xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]'); legend('r(t) [km]');
%Case 4
figure
tspan_min = [616.79 1880.41]; tspan = tspan_min*60;
Z = [1882.725 9864.690 4086.088 -5.565367 5.451548 2.258105]';
[t,Y] = ode113(@rv_orbit,tspan,Z,options);
plot3(Y(:,1),Y(:,2),Y(:,3)); hold on
title('Orbital path (Case 4)'); xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]'); legend('r(t) [km]');
%Case 5
figure
tspan_min = [21.02 1913.38]; tspan = tspan_min*60;
Z = [-664.699 8112.75 4479.81 -0.87036 -0.068046 -8.290459]';
[t,Y] = ode113(@rv_orbit,tspan,Z,options);
plot3(Y(:,1),Y(:,2),Y(:,3)); hold on
title('Orbital path (Case 5)'); xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]'); legend('r(t) [km]');
%Case 6
figure
tspan_min = [27 57]; tspan = tspan_min*60;
Z = [-10515.45 -5235.37 49.1700 -2.10305 -4.18146 5.56329]';
[t,Y] = ode113(@rv_orbit,tspan,Z,options);
plot3(Y(:,1),Y(:,2),Y(:,3));
title('Orbital path (Case 6)'); xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]'); legend('r(t) [km]');

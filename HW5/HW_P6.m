%Homework 5 Problem 6
clc; clear;
%Define initial values
Earth_r = 6378.145; %km
geoAlt = 35786; %km
mu = 398600; %km^3*s^-2
alt = 300; %km
r1 = alt + Earth_r; %km
tau2_hr = 24; %hours
tau2 = tau2_hr*3600; %seconds
r2 = (tau2^2*mu/(4*pi^2))^(1/3);
e1 = 0;
e2 = 0;
Omega = 0;
inc1 = deg2rad(28.5); %radians
inc2 = 0; %radians
omega = 0;
tau2 = 24*3600; %s
nu1 = linspace(0,2*pi,200);
%Define range for f
f = 0:0.01:1;
%Calculate a&e of transfer orbit
a = (r1+r2)/2; %semi-major axis of transfer orbit
%Initialize position & velocity vectors
rv1 = zeros(length(nu1),3);
vv1 = zeros(length(nu1),3);
rv2 = zeros(length(nu1),3);
vv2 = zeros(length(nu1),3);
%Calculate position & velocity for orbits 1 & 2
for jj = 1:length(nu1)
    nu = nu1(jj);
    oe1 = [r1 e1 Omega inc1 omega nu]';
    [rv1(jj,:), vv1(jj,:)] = oe2rv_BraceyIV_Edgar(oe1,mu);
    oe2 = [r2 e2 Omega inc2 omega nu]';
    [rv2(jj,:), vv2(jj,:)] = oe2rv_BraceyIV_Edgar(oe2,mu);
end
%Calculate cranking angle
hv1 = cross(rv1(1,:),vv1(1,:));
hv2 = cross(rv2(1,:),vv2(1,:));
theta = acos(dot(hv1,hv2)/(norm(hv1,2)*norm(hv2,2)));
% %Delta V calculations
v1minus = sqrt(mu/r1);
v1plus = sqrt((2*mu)/r1-(mu/a));
v2minus = sqrt((2*mu)/r2-(mu/a));
v2plus = sqrt(mu/r2);
%Impulse calculations for parts A & B
%Inclination change @ transfer apoapsis
dv1 = sqrt(v1minus^2+v1plus^2-2*v1minus*v1plus*cos(f(1)*theta));
dv2 = sqrt(v2minus^2+v2plus^2-2*v2minus*v2plus*cos(f(end)*theta));
fprintf('Part (a)\n')
fprintf('Total impulse when inclination changed at apoapsis\n')
fprintf('of transfer orbit: %.3f km*s^-1\n',(dv1+dv2))
%Inclination change @ initial orbit
dv1 = sqrt(v1minus^2+v1plus^2-2*v1minus*v1plus*cos(f(end)*theta));
dv2 = sqrt(v2minus^2+v2plus^2-2*v2minus*v2plus*cos(f(1)*theta));
fprintf('\nPart (b)\n')
fprintf('Total impulse when inclination changed at initial\n')
fprintf('orbit: %.3f km*s^-1\n',(dv1+dv2))
%Normalize impulses
vc1 = sqrt(mu/r1);
v1minus = v1minus/vc1;
v1plus = v1plus/vc1;
v2minus = v2minus/vc1;
v2plus = v2plus/vc1;
%Initialize impulse storage
dv_hold = zeros(size(f));
%Calculate impulses with ratio of inclination change
for jj = 1:length(f)
    dv1 = sqrt(v1minus^2+v1plus^2-2*v1minus*v1plus*cos(f(jj)*theta));
    dv2 = sqrt(v2minus^2+v2plus^2-2*v2minus*v2plus*cos((1-f(jj))*theta));
    dv_hold(jj) = dv1+dv2;
end
plot(f,dv_hold,'LineWidth',2)
title('Total \Deltav depending on inclination change ratio')
xlabel('f') %Amount of inclination change applied at initial orbit
ylabel('Total impulse')
%Find minimum impulse & value of f used to obtain
[dv_min, dv_min_ind] = min(dv_hold);
f_min = f(dv_min_ind);
fprintf('\nPart (d)\n')
fprintf('Lowest impulse = %.3f when f = %.3f.\n',dv_min,f_min)


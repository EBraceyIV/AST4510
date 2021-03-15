%Homework 5 Problem 4
clc; clear;
%Define inital values
mu = 398600; %km^3*s^-2
Earth_r = 6378.145; %km
alt = 300; %km
r1 = Earth_r + alt; %km
tau2_hr = 23.934; %hours
tau2 = 23.934*3600; %seconds
a2 = (tau2^2*mu/(4*pi^2))^(1/3);
r2 = a2;
R = r2/r1;
e1 = 0;
e2 = 0;
Omega = deg2rad(60); %radians
inc1 = deg2rad(57); %radians
inc2 = 0; %radians
omega = 0;
nu1 = linspace(0,2*pi,200);
nu2 = linspace(0,pi,200);
%Calculate a & e of transfer orbit
a = (r1+r2)/2; %semi-major axis of transfer orbit
e = (r2-r1)/(r2+r1); %eccentricity of transfer orbit
%Initialze position & velocity matrices
rv1 = zeros(length(nu1),3);
vv1 = zeros(length(nu1),3);
rv2 = zeros(length(nu1),3);
vv2 = zeros(length(nu1),3);
rvt = zeros(length(nu2),3);
vvt = zeros(length(nu2),3);
%Calculate position & velocity for orbits 1 & 2
for ii = 1:length(nu1)
    nu = nu1(ii);
    oe1 = [r1 e1 Omega inc1 omega nu]';
    [rv1(ii,:), vv1(ii,:)] = oe2rv_BraceyIV_Edgar(oe1,mu);
    oe2 = [a2 e2 Omega inc2 omega nu]';
    [rv2(ii,:), vv2(ii,:)] = oe2rv_BraceyIV_Edgar(oe2,mu);
end
%Calculate position & velocity for trasnfer orbit
for jj = 1:length(nu2)
    nu = nu2(jj);
    oe_t = [a e Omega inc1 omega nu]';
    [rvt(jj,:), vvt(jj,:)] = oe2rv_BraceyIV_Edgar(oe_t,mu);
end
%Plot orbits & transfer
plot3(rv1(:,1),rv1(:,2),rv1(:,3),rv2(:,1),rv2(:,2),rv2(:,3),'LineWidth',2)
hold on
plot3(rvt(:,1),rvt(:,2),rvt(:,3),'LineWidth',2)
%Plot Line of Intersections
rvtf = rvt(end,:).';
LoI = [-rvtf rvtf].';
plot3(LoI(:,1),LoI(:,2),LoI(:,3),'-m','LineWidth',2);
%Impulse calculations
Dv1 = vvt(1,:) - vv1(1,:);
Dv2 = vvt(end,:) - vv2(end,:);
%Plot the impulses (magnified for clarity)
quiver3(rvt(1,1),rvt(1,2),rvt(1,3),Dv1(1,1),Dv1(1,2),Dv1(1,3),3000,'LineWidth',2)
quiver3(rvt(end,1),rvt(end,2),rvt(end,3),Dv2(1,1),Dv2(1,2),Dv2(1,3),3000,'LineWidth',2)
legend('Initial Orbit','Final Orbit','Transfer Orbit','Line of Intersection','Impulse 1 (x3000)','Impulse 2 (x3000)')
xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]');
title('Orbit Transfer')
%Calculate impulse magnitudes
Dv1_mag = norm(Dv1,2);
Dv2_mag = norm(Dv2,2);
Dv_mag = Dv1_mag + Dv2_mag;
fprintf('Part (a)\n')
fprintf('Impulse 1 = %g km*s^-1\n',Dv1_mag)
fprintf('Impulse 2 = %g km*s^-1\n',Dv2_mag)
fprintf('\nPart (b)\n')
fprintf('Total Delta V required = %g km*s^-1\n',Dv_mag)
%Calculate time to complete transfer
fprintf('\nPart (c)\n')
taut = 2*pi*sqrt(oe_t(1)^3/mu); %s
taut_hr = taut/3600; %hr
transfer_time = taut_hr/2; %hr
fprintf('Time to complete transfer = %.3f hours.\n',transfer_time)
%Calculate mass ratios
fprintf('\nPart (d)\n')
g0 = 0.00980665; %km*s^-2
mRatio1 = exp(Dv1_mag/(g0*320));
mRatio2 = exp(Dv2_mag/(g0*320));
fprintf('Mass ratio 1 = %g\n',mRatio1)
fprintf('Mass ratio 2 = %g\n',mRatio2)
fprintf('\nPart (f)\n')
fprintf('Changing the longitude of the ascending node does not\n')
fprintf('change where the transfer orbit starts on the inital orbit.\n')


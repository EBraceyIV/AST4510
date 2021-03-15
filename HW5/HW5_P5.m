%Homework 5 Problem 5
clc; clear;
%Define initial values
alt = 350; %km
Earth_r = 6378.145; %km
mu = 398600; %km^3*s^-2
r1 = Earth_r + alt; %km
r2 = 26558; %km
R = r2/r1;
e1 = 0;
e2 = 0;
Omega = 0;
inc1 = deg2rad(28); %radians
inc2 = deg2rad(55); %radians
omega = 0;
nu1 = linspace(0,2*pi,200);
nu2 = linspace(0,pi,200);
%Calculate a & e of transfew orbit
a = (r1+r2)/2; %semi-major axis of transfer orbit
e = (r2-r1)/(r2+r1); %eccentricity of transfer orbit
%Initialze position & velocity vectors
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
    oe2 = [r2 e2 Omega inc2 omega nu]';
    [rv2(ii,:), vv2(ii,:)] = oe2rv_BraceyIV_Edgar(oe2,mu);
end
%Calculate position & velocity for transfer orbit
for jj = 1:length(nu2)
    nu = nu2(jj);
    oe_t = [a e Omega inc1 omega nu]';
    [rvt(jj,:), vvt(jj,:)] = oe2rv_BraceyIV_Edgar(oe_t,mu);
end
%Calculate how long transfer takes
taut = 2*pi*sqrt(oe_t(1)^3/mu); %s
taut_hr = taut/3600; %hr
transfer_time = taut_hr/2; %hr
%Calculate cranking angle
hv1 = cross(rv1(1,:),vv1(1,:));
hv2 = cross(rv2(1,:),vv2(1,:));
theta = acos(dot(hv1,hv2)/(norm(hv1,2)*norm(hv2,2)));
%Line of Intersection calculations
lv = cross(hv1,hv2)/(norm(cross(hv1,hv2),2));
rvft = rvt(end,:).';
LoI = [-rvft rvft].';
plot3(LoI(:,1),LoI(:,2),LoI(:,3),'-m','LineWidth',2);
hold on
%Plot orbits & transfer
plot3(rv1(:,1),rv1(:,2),rv1(:,3),rv2(:,1),rv2(:,2),rv2(:,3),'LineWidth',2)
hold on
plot3(rvt(:,1),rvt(:,2),rvt(:,3),'LineWidth',2)
%Impulse calculations
Dv1 = vvt(1,:) - vv1(1,:);
Dv2 = vvt(end,:) - vv2(end,:);
Dv_mag = norm(Dv1,2)+norm(Dv2,2);
%Plot the impulses (magnified for clarity)
quiver3(rvt(1,1),rvt(1,2),rvt(1,3),Dv1(1,1),Dv1(1,2),Dv1(1,3),3000,'LineWidth',2)
quiver3(rvt(end,1),rvt(end,2),rvt(end,3),Dv2(1,1),Dv2(1,2),Dv2(1,3),3000,'LineWidth',2)
legend('Line of Intersection','Initial Orbit','Final Orbit','Transfer Orbit','Impulse 1 (x3000)','Impulse 2 (x3000)')
xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]');
title('Orbit Transfer')
fprintf('Part (a)\n')
fprintf('Line of Intersection = [%.3f %.3f %.3f] (km)\n',lv)
fprintf('\nPart (b)\n')
fprintf('Impulse 1 applied at [%.3f %.3f %.3f] (km)\n',rvt(1,:))
fprintf('Impulse 2 applied at [%.3f %.3f %.3f] (km)\n',rvt(end,:))
fprintf('\nPart (c)\n')
fprintf('Total Delta v required = %.3f (km*s^-1)\n',Dv_mag)
fprintf('\nPart (d)\n')
fprintf('Total time to complete transfer = %.3f hours.\n',transfer_time)
fprintf('\nPart (e)\n')
fprintf('Eccentricity of transfer orbit = %.3f\n',oe_t(2))
fprintf('\nPart (f)\n')
fprintf('Theta = %.3f degrees\n',rad2deg(theta))




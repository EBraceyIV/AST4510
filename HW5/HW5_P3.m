%Homework 5 Problem 3
clc; clear;
%Define inital values
mu = 1;
vc1 = 1;
vc2 = 0.5;
r1 = mu/vc1^2;
r2 = mu/vc2^2;
e1 = 0;
e2 = 0;
Omega = 0;
inc = 0;
omega = 0;
nu1 = linspace(0,2*pi,200)';
nu2 = linspace(0,pi,200)';
%Calculate a & e of transfer orbit
a = (r1+r2)/2;
e = (r2-r1)/(r2+r1);
%Initialize matrices for r&v values
rv1 = zeros(length(nu1),3);
vv1 = zeros(length(nu1),3);
rv2 = zeros(length(nu1),3);
vv2 = zeros(length(nu1),3);
rvt = zeros(length(nu2),3);
vvt = zeros(length(nu2),3);
%Calculate inital and final orbits
for ii = 1:length(nu1)
    nu = nu1(ii);
    oe1 = [r1 e1 Omega inc omega nu]';
    [rv1(ii,:), vv1(ii,:)] = oe2rv_BraceyIV_Edgar(oe1,mu);
    oe2 = [r2 e2 Omega inc omega nu]';
    [rv2(ii,:), vv2(ii,:)] = oe2rv_BraceyIV_Edgar(oe2,mu);
end
%Calculate transfer orbit
for jj = 1:length(nu2)
    nu = nu2(jj);
    oe_t = [a e Omega inc omega nu]';
    [rvt(jj,:), vvt(jj,:)] = oe2rv_BraceyIV_Edgar(oe_t,mu);
end
%Calculate impulses
Dv1 = vvt(1,:) - vv1(1,:);
Dv2 = vvt(end,:) - vv2(end,:);
%Plot orbits & impulses
plot(rv1(:,1),rv1(:,2),rv2(:,1),rv2(:,2),rvt(:,1),rvt(:,2),'LineWidth',2)
hold on
quiver(rvt(1,1),rvt(1,2),Dv1(1,1),Dv1(1,2),2,'LineWidth',2)
quiver(rvt(end,1),rvt(end,2),Dv2(1,1),Dv2(1,2),2,'LineWidth',2)
axis equal
xlabel('x'); ylabel('y');
title('Orbit Transfer')
legend('Inital Oribt','Final Orbit','Transfer Orbit','Impulse 1 (x2)','Impulse 2 (x2)')

%Homework 3 Problem 8
clc; clear;
mu = 398600;
a = 19133.333;
e = 0.5;
Omega = 30*pi/180;
inc = 45*pi/180;
omega = 45*pi/180;
nu = 0*pi/180;
oe = [a e Omega inc omega nu]';
[rPCI, vPCI] = oe2rv_BraceyIV_Edgar(oe,mu);
fprintf('position vector [km] = \n')
for ii = 1:length(rPCI)
    fprintf('\t\t%.8f\n',rPCI(ii))
end
fprintf('\n')
fprintf('velocity vector [km*s^-1] = \n')
for ii = 1:length(vPCI)
    fprintf('\t\t%.8f\n',vPCI(ii))
end

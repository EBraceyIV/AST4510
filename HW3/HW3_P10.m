%Homework 3 Problem 10
clc; clear;
mu = 1;
a = 1.6;
e = 0.4;
Omega = 287*pi/180;
inc = 46*pi/180;
omega = 28*pi/180;
nu = 139*pi/180;
oe = [a e Omega inc omega nu]';
[rPCI, vPCI] = oe2rv_BraceyIV_Edgar(oe,mu);
fprintf('position vector = \n')
for ii = 1:length(rPCI)
    fprintf('\t\t%.8f\n',rPCI(ii))
end
fprintf('\n')
fprintf('velocity vector = \n')
for ii = 1:length(vPCI)
    fprintf('\t\t%.8f\n',vPCI(ii))
end


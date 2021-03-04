%Homework 3 Problem 7
clc; clear;
mu = 398600;
a = 15307.548;
e = 0.7;
Omega = 194*pi/180;
inc = 39*pi/180;
omega = 85*pi/180;
nu = 48*pi/180;
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

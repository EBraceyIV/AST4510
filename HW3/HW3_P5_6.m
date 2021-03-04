%Homework 3 Problem 5 & 6
clc; clear;

%%Problem 5
fprintf('Problem 5\n')
%Define given values
mu = 1;
rPCI = [0.7 0.6 0.3]'; 
vPCI = [-0.8 0.8 0]';
%Part A
fprintf('Part A\n')
%Calculate & define orbital elements
oe = rv2oe_BraceyIV_Edgar(rPCI,vPCI,mu);
a = oe(1); fprintf('semi-major axis = %.8f\n',a)
e = oe(2); fprintf('eccentricity =  %.8f\n', e)
Omega = oe(3); fprintf('longitude of the ascending node = %.8f\n',Omega)
inc = oe(4); fprintf('orbital inclination = %.8f\n', inc)
omega = oe(5); fprintf('argument of the periapsis = %.8f\n',omega)
nu = oe(6); fprintf('true anomaly = %.8f\n',nu)


%Part B
fprintf('\nPart B\n')
tau = 2*pi*sqrt(a^3/mu); %Calculate orbital period
fprintf('orbital period = %.8f\n',tau)

%Part C
fprintf('\nPart C\n')
p = a*(1-e^2); %Calculate semi-latus rectum
fprintf('semi latus rectum = %.8f\n',p)

%Part D
fprintf('\nPart D\n')
hv = cross(rPCI,vPCI); %Calculate specific angular momentum vector
h = norm(hv,2); %Calculate magnitude of specific angular momentum
fprintf('specific angular momentum =\n')
for ii = 1:length(hv)
    fprintf('\t\t\t%.8f\n',hv(ii))
end
fprintf('magnitude of specific angular momentum = %.8f\n',h)

%Part E
fprintf('\nPart E\n')
epsilon = -mu/(2*a); %Calculate specific mechanical energy
fprintf('specific mechanical energy= %.8f\n',epsilon)

fprintf('\n\n')

%%Problem 6
fprintf('Problem 6\n')
%Calculate & define the position and velocity vectors
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

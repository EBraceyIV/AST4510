function [rPCI, vPCI] = oe2rv_BraceyIV_Edgar(oe,mu)
%%%INPUTS:
%oe = orbital elements
%mu = gravitational constant

%%%OUTPUTS:
%rv = position column vector
%vv = velocity column vector

%%%CALCULATIONS:
%Define each orbital element from matrix
a = oe(1);
e = oe(2);
Omega = oe(3);
inc = oe(4);
omega = oe(5);
nu = oe(6);
%Calculate semi-latus rectum
p = a*(1-e^2);
%Calculate radius @ position
r = p/(1+e*cos(nu));
%Calculate the position & velocity vectors in P
rv = [r*cos(nu) r*sin(nu) 0]';
vv = sqrt(mu/p).*[-sin(nu) e+cos(nu) 0]';
%Calculate the transformation matrices
T_NI = [cos(Omega) -sin(Omega) 0; sin(Omega) cos(Omega) 0; 0 0 1];
T_QN = [1 0 0; 0 cos(inc) -sin(inc); 0 sin(inc) cos(inc)];
T_PQ = [cos(omega) -sin(omega) 0; sin(omega) cos(omega) 0; 0 0 1];
T_PI = T_NI*T_QN*T_PQ;
%Transform the position and velocity vectors from P to I
rPCI = T_PI*rv;
vPCI = T_PI*vv;
end

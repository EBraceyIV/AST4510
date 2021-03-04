function oe = rv2oe_BraceyIV_Edgar(rPCI,vPCI,mu)
%%%INPUTS:
%rPCI = position column vector
%vPCI = velocity column vector
%mu = gravitational constant

%%%OUPUT:
%oe = [a e Omega i omega nu]'
    %a = semi-latus rectum
    %e = eccentricity
    %Omega = longitude of the ascending node
    %i = inclination
    %omega = argument of the periapsis
    %nu = true anomaly
    
%%%CALCULATIONS:
%Calculate specific angular momentum vector & magnitude
hv = cross(rPCI,vPCI); h = norm(hv,2);
%Calculate eccentricity vector & magnitude
ev = cross(vPCI,hv)/mu-rPCI/norm(rPCI,2); e = norm(ev,2);
%Calculate semi-latus rectum
p = h^2/mu;
%Calculate semi-major axis
a = p/(1-e^2);
%Initialize orbital elements matrix
oe = zeros(6,1);
oe(1) = a; %Add semi-major axis to orbital elements list
oe(2) = e; %Add eccentricity to orbital elements list
%Define inertial unit vectors
Ix = [1 0 0]'; Iy = [0 1 0]'; Iz = [0 0 1];
%Calculate line of nodes & magnitude
nv = cross(Iz,hv);
n = norm(nv,2);
%Calculate longitude of the ascending node
Omega = atan2(-dot(nv,Iy),-dot(nv,Ix)) + pi;
%Add longitude of the ascending node to orbital elements list
oe(3) = Omega; 
%Calculate orbital inclination
i = atan2(dot(hv,cross(nv,Iz)),n*dot(hv,Iz));
%Add orbital inclination to orbital elements list
oe(4) = i; 
%Calculate arugment of the periapsis
omega = atan2(-dot(ev,cross(hv,nv)),-h*dot(ev,nv)) + pi;
%Add argument of the periapsis to orbital elements list
oe(5) = omega; 
%Calculate eccentricity
nu = atan2(-dot(rPCI,cross(hv,ev)),-h*dot(rPCI,ev)) + pi;
%Add true anomaly to orbital elements list
oe(6) = nu; 
%Output orbital elements as a column vector
oe = oe';
end
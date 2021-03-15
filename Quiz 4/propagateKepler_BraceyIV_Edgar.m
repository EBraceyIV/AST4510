function [rPCIf,vPCIf,E0,nu0,E,nu] = propagateKepler_BraceyIV_Edgar(t0,t,rPCI0,vPCI0,mu)
% % ------------------------------------------------------------------ %
% % ----------------------- propagateKepler.m ------------------------ %
% % ------- Propagate Spacecraft Orbit Using Kepler's Equation ------- %
% % ------------------------------------------------------------------ %
% % Given a position and inertial velocity at a time t0 expressed in   %
% % planet-centered inertial (PCI) coordinates, determine the position %
% % and inertial velocity at a later time t on an elliptic orbit by    %
% % solving Kepler's equation.                                         %
% % ------------------------------------------------------------------ %
% % ------------- Inputs (Supplied Data) for Test Cases -------------- %
% % ------------------------------------------------------------------ %
% %    t0 = initial time                                               %
% %     t = terminal time                                              %
% % rPCI0 = Initial PCI Position                                       %
% % vPCI0 = Initial PCI Inertial Velocity                              %
% %    mu = planet gravitational parameter                             %
% % ------------------------------------------------------------------ %
% % --------- Output (Computed Quantities) from Test Cases ----------- %
% % ------------------------------------------------------------------ %
% % rPCIf = Terminal PCI Position                                      %
% % vPCIf = Terminal PCI Inertial Velocity                             %
% %    E0 = Eccentric Anomaly at Time t0                               %
% %   nu0 = True Anomaly at Time t0                                    %
% %     E = Eccentric Anomaly at Time t                                %
% %    nu = True Anomaly at Time t                                     %
% % ------------------------------------------------------------------ %
% % -------- Note: all quantities must be in consistent units -------- %
% % ------------------------------------------------------------------ %

%Calculate orbital elements given initial conditions
oe = rv2oe_BraceyIV_Edgar(rPCI0,vPCI0,mu);
%Define semi-major axis, eccentricity, and initial true anomaly
a = oe(1); e = oe(2); nu0 = oe(6);
%Calculate the initial eccentric anomaly using inital true anomaly
E0 = 2*atan2(sqrt(1-e)*sin(nu0/2),sqrt(1+e)*cos(nu0/2));
%Use Kepler Solver to find the terminal eccentric and true anomaly
[E, nu] = KeplerSolver_BraceyIV_Edgar(a,e,t0,t,nu0,mu);
%Redfine the true anomaly in the orbital elements as the terminal value
oe(6) = nu;
%Calculate terminal position and velocity vectors using orbital elements
[rPCIf, vPCIf] = oe2rv_BraceyIV_Edgar(oe,mu);
end


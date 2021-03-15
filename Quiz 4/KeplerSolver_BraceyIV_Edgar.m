function [E,nu] = KeplerSolver_BraceyIV_Edgar(a,e,t0,t,nu0,mu)    
% % ---------------------------------------------------------------------- %
% % --------- Kepler Solver for Propagation Along Elliptic Orbits -------- %
% % Solve Kepler's equation for the eccentric anomaly and the true anomaly %
% % given an initial time, t0, the initial true anomaly, nu0, the terminal %
% % time, t, the semi-major axis, a, eccentricity, e, and the planet       %
% % gravitational parameter, mu.  This function takes into account the     %
% % number of periapsis crossings en route from t0 to t and assumes that   %
% % the orbit is elliptic (that is, is is assumed that the eccentricity is %
% % strictly less than unity.                                              %
% % ---------------------------------------------------------------------- %
% % Inputs: a, e, mu, t, t0, and nu0                                       %
% %   a     = semi-major axis                                              %
% %   e     = eccentricity                                                 %
% %   t0    = initial time                                                 %
% %   t     = terminal time                                                %
% %   nu0   = true anomaly at t0                                           %
% %   mu    = gravitational parameter of planet                            %
% % ---------------------------------------------------------------------- %
% % Outputs: E and nu                                                      %
% % ---------------------------------------------------------------------- %
% %   E     = the eccentric anomaly at time t                              %
% %   nu    = the true anomaly at time t                                   %                            
% % ---------------------------------------------------------------------- %
% % IMPORTANT: the units of a, mu, t0, and t must be consistent            %
% % ---------------------------------------------------------------------- %

%Calculate the initial eccentric anomaly using inital true anomaly
E0 = 2*atan2(sqrt(1-e)*sin(nu0/2),sqrt(1+e)*cos(nu0/2));
%Set the number of Newton iterations to perform
Iter = 60;
%Define the orbital period
tau = 2*pi*sqrt((a^3)/mu);
%Find the current time since last periapsis crossing
tp = t0 - (sqrt(a^3/mu)*(E0-e*sin(E0)));
%Calculate number of periapsis crossings so far
k = floor((t - tp)/tau);
%Calculate the C constant at the current time
C = -sqrt(mu/(a^3))*(t-t0) + 2*pi*k - (E0 - e*sin(E0));
%Calculate mean anomaly
E = E0 - e*sin(E0);
%Find the root of Kepler's equation using fixed point iteration
for ii = 1:Iter
    %Calculate eccentic anomaly using mean anomaly
    E = e*sin(E) - C;
end
%Calculate the terminal true anomaly using terminal eccentric anomaly 
nu = 2*atan2(sqrt(1+e)*sin(E/2),sqrt(1-e)*cos(E/2));
end

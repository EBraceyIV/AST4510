function nu = rootFinder(f,t0,nu0,t,p,e,mu,N)
% --------------------------------------------------------------%
% ------------------------ ROOTFINDER.M ------------------------%
% --------------------------------------------------------------%
% Find the root of a function of the form                       %
%               G(nu) = F(nu) - (t-t0) = 0                      %
% where F(nu) is given as                                       %
%                                                               %
%                           /nu                                 %
%                          /                                    %
%                          |                                    %
%             F(nu)   =    | f(q)dq                             %
%                          |                                    %
%                          /                                    %
%                         /nu0                                  %
%                                                               %
% The inputs and outputs of this function are as follows:       %
% Inputs:                                                       %
%     f: integrand of function                                  %
%   nu0: true anomaly at current time                           %
%    t0: current time                                           %
%     t: terminal time                                          %
%     p: semi-latus rectum                                      %
%     e: eccentricity                                           %
%    mu: gravitational parameter                                %
%     N: number of Legendre-Gauss points                        %
% Output:                                                       %
%    nu: true anomaly at time t                                 %
% --------------------------------------------------------------%
N_It = 20; %Define the number of Newton iterations to run
nu = nu0 + 5*pi/180;  %Define guess for next true anomaly (+5 degrees)
%Run Newton's method for N_It iterations
for ii = 1:N_It
    %Calculate the time change between the two true anomalies
    [x, w] = GaussPointsWeights(nu0,nu,N);
    F = w'*f(x);
    %Calculate the values needed to use Newton's method
    G = F-(t-t0);
    g = f(nu);
    %Calculate the new nu
    nu = nu - G/g;
end
end
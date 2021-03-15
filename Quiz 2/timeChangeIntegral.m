function deltat = timeChangeIntegral(f,nu1,nu2,c,N)
% --------------------------------------------------------------% 
% This function employs Legendre-Gauss quadrature to compute an % 
% approximation to                                              %
%                              /nu2                             %
%                             /                                 % 
%                             |                                 % 
%          deltat = t2 - t1 = | f(nu,c)dnu                      % 
%                             |                                 % 
%                             /                                 % 
%                            /nu1                               % 
% where f(nu) is a function that used to define the change in   % 
% The inputs and outputs of this function are as follows:       %
% Inputs:                                                       %
%   f    	 = a handle to the function to be integrated        %
%   nu1      = lower integration limit                          %
%            = initial true anomaly (rad)                       %
%   nu2      = upper integration limit                          %
%            = terminal true anomaly (rad)                      %
%   c        = a column vector of constants                     %
%   N        = number of Gauss points & weights                 %
%              used to approximate the integral                 %
% Output:                                                       %
%   deltat   = Gauss quadrature approximation of                %
%               deltat, where deltat = t2 - t1 is the           %
%               time change from nu1 to nu2                     %
% --------------------------------------------------------------%
%Obtain the points and weights for Gauss quadrature
[x,w] = GaussPointsWeights(nu1,nu2,N); 
%Calculate the F vector running the f function at the Gauss points
Fv = f(x,c);
%Calculate the deltat by weighting each value of the F vector
deltat = w.'*Fv;
end
function [x,w] = GaussPointsWeights(a,b,N)
% --------------------------------------------------------------%
% This function computes the N Gauss quadrature points and the  %
% corresponding N Gauss quadrature weights on [a,b].            %
% Inputs:                                                       %
%   a = lower limit of integral to be evaluated                 %
%   b = upper limit of integral to be evaluated                 %
%   N = Number of Gauss Points and Weights to Be Computed.      %
% Output:                                                       %
%   x = Column Vector with the N Gauss Points on [a,b].         %
%   w = Column Vector with the N Gauss Quadrature Weights.      %
% --------------------------------------------------------------%
%         This code uses the Golub-Welsch algorithm.            %
% --------------------------------------------------------------%

beta = .5./sqrt(1-(2*(1:N-1)).^(-2)); % 3-term recurrence coeffs
T = diag(beta,1) + diag(beta,-1);     % Jacobi matrix
[V,D] = eig(T);                       % Eigenvalue decomposition
x= diag(D); [x,i] = sort(x);          % Gauss points on [-1,+1]
w= (2*V(1,i).^2).';                   % Gauss weights on [-1,+1]
x= (b-a)/2*x + (b+a)/2;               % Gauss points on [a,b]
w= (b-a)/2*w;                         % Gauss weights on [a,b]
end
 function K = KEppan(a,b,h)
%  Computes the Epanechnikov kernel for each
%  value in X with a bandwidth H. 
%  The Epanechnikov kernel is defined as:
%  K(x-x0) = 0.75 * (1 - (x-x0)^2) for |x| <= 1
%  K(x) = 0 otherwise

% Calculate x / h
x= a-b;
x = x / h;

% Compute the Epanechnikov kernel
K = zeros(size(x));
indices = abs(x) <= 1;
K(indices) = 0.75 * (1 - x(indices).^2);

end
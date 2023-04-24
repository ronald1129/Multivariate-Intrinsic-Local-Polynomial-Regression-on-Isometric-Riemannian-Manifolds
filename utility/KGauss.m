      function [K] = KGauss(a,b,h)
        [dd,~] = size(a);
         K = exp(-0.5*FNorm(a-b)/h^2)/(sqrt(2*pi)*h);
    end
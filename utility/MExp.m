function [E] = MExp(A)
       E = diag(exp(diag(A)));
    end
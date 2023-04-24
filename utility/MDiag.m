 function [D] = MDiag(A)
        [d,~] = size(A);
        D = eye(d).*A;
    end
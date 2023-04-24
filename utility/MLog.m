    function [L] = MLog(A)
       L = diag(log(diag(A)));
    end
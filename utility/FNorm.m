    function [FN] = FNorm(A)
    % Compute Frobenis norm of A
        FN=trace(A*A');
    end

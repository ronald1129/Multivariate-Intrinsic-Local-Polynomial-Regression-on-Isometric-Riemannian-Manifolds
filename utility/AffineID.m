 % Affine Invariant distance
        function [D] = AffineID(A,B)
            D= sqrt(FNorm(logm(A\B)));
        end
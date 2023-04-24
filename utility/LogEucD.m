% Log Euclidean distance
        function [D] = LogEucD(A,B)
            D = sqrt(FNorm(logm(A)-logm(B)));
        end
        
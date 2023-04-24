function [D] = LogCholD(A,B)
            LA = chol(A)';
            LB = chol(B)';
            L_delta1 = MLowp(LA-LB);
            L_delta2 = MLog(MDiag(LA)) - MLog(MDiag(LB));
            D = sqrt(FNorm(L_delta1) + FNorm(L_delta2));
        end
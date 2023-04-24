function [Inv] = FInv(A)
        Inv=diag(diag(A).^(-1));
    end
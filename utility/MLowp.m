function [L] = MLowp(A)
        [d,~] = size(A);
        L = tril(ones(d)-eye(d)).*A;
    end
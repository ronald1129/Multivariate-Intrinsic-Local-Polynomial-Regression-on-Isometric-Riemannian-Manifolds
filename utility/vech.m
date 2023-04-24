function [Vec]=vech(A)
        [d,~]=size(A);
        Vec=diag(A);
        for i=1:d
            for j=1:d
                if i>j
                    Vec=[Vec;A(i,j)];
                end
            end
        end
    end
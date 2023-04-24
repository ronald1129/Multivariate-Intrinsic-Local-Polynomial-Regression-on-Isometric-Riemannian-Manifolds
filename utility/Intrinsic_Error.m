function [Error,MError] = Intrinsic_Error(A,B,Distance)

[~,~,N] = size(A);
Error = zeros(1,N);

for j = 1:N

   if strcmp(Distance,'LogChol')
        error = LogCholD(A(:,:,j),B(:,:,j));
    elseif strcmp(Distance,'LogEuc')
        error = LogEucD(A(:,:,j), B(:,:,j));
    elseif strcmp(Distance,'Forb')
        error = FrobD(A(:,:,j), B(:,:,j));    
    else 
        error = AffineID(A(:,:,j),B(:,:,j));
   end

   Error(j) = error;
end
MError = mean(Error);
end


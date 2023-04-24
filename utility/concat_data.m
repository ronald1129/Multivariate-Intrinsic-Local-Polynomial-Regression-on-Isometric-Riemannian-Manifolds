function [C,labels]=concat_data(A,B)
[d,~,Na]=size(A);
[d,~,Nb]=size(B);
C=zeros(d,d,Na+Nb);
for j=1:Na+Nb
    if j<Na+1
        C(:,:,j)=A(:,:,j);
    else 
        C(:,:,j)=B(:,:,j-Na);
    end
end
end


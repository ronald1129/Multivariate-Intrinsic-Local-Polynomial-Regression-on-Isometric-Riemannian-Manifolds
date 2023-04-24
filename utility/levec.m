function [L,V]=levec(A,M)
[d,~,N]=size(A);
L=zeros(d,d,N);
nd=d*(d+1)/2;
V=zeros(nd,N);
for j=1:N
    L(:,:,j)=logm(A(:,:,j));
    V(:,j)=M.Vec(L(:,:,j));
end
end

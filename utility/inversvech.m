function [M,L]=inversvech(v)
%this functi,on its the inverse operator of half vec
v=v';
N=max(size(v));
n=floor((-1+sqrt(1+8*N))/2);
M=diag(v(1:n))/2;
l=1;
k=2;
for j=1:N-n
   if k<n+1
        for i=1:k-1
          M(k,i)=v(1,n+l);
          l=l+1;
        end
   end
   k=k+1;
end
L=M;
M=M+M';
end
%
% v=rand(20);
% v=(v+v')/2;
% r=inversvech(HalfVec(v)')
% v
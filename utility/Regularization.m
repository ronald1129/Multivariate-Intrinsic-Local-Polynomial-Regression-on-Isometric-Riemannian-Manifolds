function [Xreg,X]=Regularization(X)
X=(X+X')/2;

emin=min(eig(X));
p=0;
delta=0.01;
[d,~]=size(X);
L=trace(X)/d;

Xreg=X;
%------------------------------------
while emin<0&& p<1
    Xreg=(1-p)*X+p*L*eye(size(X));
    emin=min(eig(Xreg));
    p=p+delta;
end
Xreg =Xreg+0.000001*eye(size(Xreg));
%------------------------------------
% if lplot==1
%     figure(1)
%       surf(abs(X));
%     figure(2)
%       surf(abs(Xreg));
% end

end
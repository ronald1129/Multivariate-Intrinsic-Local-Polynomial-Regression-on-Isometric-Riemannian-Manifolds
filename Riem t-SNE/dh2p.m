function [P]=dh2p(D,H,perp,tol)
% Input 
        % Geodesic distance matrix: D nxn
        % Geodesic distance matrix: H nxn
        % Data space dimensionality: dim
        % Desired perplexity:perp
% Output
        % nxn matrix of conditional probabilities
%-------------------------------------------------
[N,~]=size(D);
P=zeros(N,N);
bin_search_steps=100;
if ~exist("tol")|| isempty(tol)
    tol=10^(-4);
end
perp_tol=tol;
for i=1:N
    t_min=-10^(10);
    t_max=10^(10);
    t=1;
    for l=1:bin_search_steps
        row_sum=0;
        for j=1:N
            P(i,j)=H(i,j)*exp(-D(i,j)^2/(2*t));
            row_sum=row_sum+P(i,j);
        end
        entropy=0;
        for j=1:N
            P(i,j)=P(i,j)/row_sum;
            if P(i,j)~=0.0
                entropy=entropy+P(i,j)*log2(P(i,j));
            end
        end
        entropy_diff=-entropy-log2(perp);
        if abs(entropy_diff)<=perp_tol
            break 
        end 
        if entropy_diff<0
            t_min=t;
            if t_max==10^(16)
                t=2*t;
            else
                t=(t+t_max)/2;
            end
        else 
            t_max=t;
            if t_min==10^(16)
                t=t/2;
            else
                t=(t+t_min)/2;
            end
        end
    end
end
end





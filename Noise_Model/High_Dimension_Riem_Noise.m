function [DataNoise]=High_Dimension_Riem_Noise(Data,noise_level)

% Data its on the HPD Manifold

[d,~,N] = size(Data);

% Data with noise

DataNoise = zeros(d,d,N);

% Parameters
%s=rand(d*(d+1)/2);
%sigma = Regularization(s+s')/2;
sigma = noise_level*eye(d*(d+1)/2);
mu = zeros(1,d*(d+1)/2);
epsilon = mvnrnd(mu,sigma,N)';

for j = 1:N

    G=chol(Data(:,:,j));
    DataNoise(:,:,j)=Regularization(G*expm(inversvech(epsilon(:,j)))*G');
    
end

end



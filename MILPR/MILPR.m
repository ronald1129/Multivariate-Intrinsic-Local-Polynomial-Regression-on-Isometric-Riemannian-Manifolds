function [S_x] = MILPR(Sr,xv,x0,k,h,f,invf,K_Type,fdata)
%% --------
% Kernel

% Check the kernel type
if K_Type == 'KEppan'
    Kern = @(a1,a2,hs) KEppan(a1,a2,hs);
elseif K_Type == 'KGauss'
    Kern = @(a1,a2,hs) KGauss(a1,a2,hs);
end

%% ------
% Dimensions

% Covariate dimension
p = size(x0,1);
% Manifold and sample dimensions
[n,~,N] = size(Sr);
% Regularization parameter
lamb = 10^(-3);

% Calculate the polynomial dimension
if p == 1
    pdim = k + 1;
else
    pdim = (p^(k+1)-1)/(p-1);
end

%% ----
% Regression matrices

% Initialize the response matrix
Y = [];

% If function data is not provided, evaluate the function at the sample points
if isempty(fdata)
    for i = 1:N
        Y = [Y, f(Sr(:,:,i))];
    end
% Otherwise, use the provided function data
elseif ~isempty(fdata)
    for i = 1:N
        Y = [Y, fdata(:,:,i)];
    end
end

% Initialize the design matrix and weight matrix
X = [];
W = [];                

% Initialize the local regression function handle
fS_x = @(x) 0;

% Loop over the sample points to calculate the weight and design matrices
for i = 1:N 
    % Calculate the weight for the current sample point
    W = [W, Kern(xv(:,i),x0,h)];
    
    % Calculate the design matrix for the current sample point
    Xcol_i = [];
    for j = 1: k+1
        Xcol_i = [Xcol_i; Kron_power(xv(:,i)-x0,j-1)];
    end
    X = [X, Xcol_i];
end

% Construct the diagonal weight matrix
W = diag(W);

%% Regression

% Calculate the coefficient vector
theta = W * X' * inv(X * W * X' + lamb^2 * eye(pdim));

% Construct the coefficient tensor
Theta = Box_Product(eye(n),theta,ones(1,N),p.^(0:1:k));

% Calculate the response function coefficients
beta = Y * Theta;

%%
% Local Polynomial

% Initialize the local polynomial function handle
for j = 1:k+1
     if p ==1 
         j_start = 1 + n*(j-1);
         j_end   = n*j;
     else
         j_start = 1 + n*(p^(j-1)-1)/(p-1);
         j_end   = n*(p^j-1)/(p-1);
     end
     fS_x = @(x) fS_x(x) + beta(:,j_start:j_end) * kron(eye(n),Kron_power(x-x0, j-1));
end

% Construct the final function handle
S_x = @(x) invf(fS_x(x));

% Define the Kronecker power function
function [x_power] = Kron_power(xx,pow)
if pow == 0
x_power = 1;
else 
    x_power = kron(Kron_power(xx,pow-1),xx);
end

end

end












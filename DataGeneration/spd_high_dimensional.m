function [Data,x] = spd_high_dimensional(cdim,mdim,sdim)

% initialize data matrix and vector of parameters
Data = zeros(mdim,mdim,sdim);
x =zeros(cdim,sdim);

% generate random parameter vectors
for j = 1:cdim
    x(j,:) = linspace(-rand(1),rand(1),sdim);
end

% create symbolic variable y
y = sym('y');

% generate matrices S_ij for each pair (i,j) of indices
for i = 1:mdim
    for j = 1:mdim
        if i < j
            % generate random beta vector for off-diagonal entries
            beta_ij = rand(1,cdim);
            % define S_ij as a function of y
            S{i,j} = @(y) beta_ij*y/sqrt(FNorm(beta_ij));
        elseif i == j 
            % use constant beta vector for diagonal entries
            beta_ij = 0.5*ones(1,cdim);
            S{i,j} = @(y) beta_ij*y/sqrt(FNorm(beta_ij));
        end
    end
end

% evaluate data matrix entries for each value of m
for m = 1:sdim
    for i = 1:mdim
        for j = 1:mdim
           if i<=j
                % evaluate S_ij for the given parameter vector x(:,m)
                Data(i,j,m) = (S{i,j}(x(:,m)));
           end
        end
    end
    % make the data matrix symmetric by averaging with its transpose
    Data(:,:,m) = (Data(:,:,m)+Data(:,:,m)')/2;
    % apply the matrix exponential to each slice of the data matrix
    Data(:,:,m) = expm(Data(:,:,m));
end

end

% helper function to compute the Frobenius norm of a vector
function norm = FNorm(x)
    norm = sqrt(sum(x.^2));
end

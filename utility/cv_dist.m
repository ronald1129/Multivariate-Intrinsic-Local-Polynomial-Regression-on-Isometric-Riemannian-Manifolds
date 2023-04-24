function [d] = cv_dist(Sr,xv,x0,k,h,f,invf,kernt,Distance,fdata,eval_option)

[dim, dim, N] = size(Sr);
[dimc, ~] = size(xv);

% Preallocate memory for test data
Sr_test = zeros(dim, dim, N-1);
xv_test = zeros(dimc, N-1);

% Define a counter for summing up the errors
error_sum = 0;
Sr2 = Sr;

if ~isempty(fdata)
    Sr = fdata;
    end


% Loop over each fold in the cross-validation
for j = 1:N
    % Split the data into training and validation set
    training_indices = setdiff(1:N, j);
    Sr_train = Sr(:, :, training_indices);
    xv_train = xv(:, training_indices);
    
    % Compute mean of the training data
    x0_eval = xv(:,j);
    
    % Run the regression on the training data
    Reg = eval_MILPR(Sr_train, xv_train, x0_eval, k, h, f, invf, kernt,Sr_train,eval_option);
    
    % Compute the error for the validation data
    if strcmp(Distance, 'LogChol')
        error = LogCholD(Reg, Sr2(:, :, j))^2;
    elseif strcmp(Distance,'LogEuc')
        error = LogEucD(Reg, Sr2(:, :, j))^2;
    elseif strcmp(Distance, 'AffineI')
        error = AffineID(Reg, Sr2(:, :, j))^2;
    end
    
    % Add the error to the sum
    error_sum = error_sum + error;
end

% Compute the average error over all the folds
d = (error_sum / N);



end

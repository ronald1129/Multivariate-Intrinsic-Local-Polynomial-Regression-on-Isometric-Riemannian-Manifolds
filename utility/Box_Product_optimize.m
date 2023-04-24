function [C] = Box_Product_optimize(A,B,part1,part2)
%Bib

%[1] - MathWorks. (n.d.). Kronecker product of matrices - MATLAB kron. 
%      Retrieved from https://www.mathworks.com/help/matlab/ref/kron.html
%[2] - 

% Get the number of rows and columns of partition vectors
N = length(part1);
M = length(part2);

% Initialize the resultant matrix
C = [];

% Loop over the rows of the partition vector part1
for i = 1:N
    row_i = [];
    % Loop over the columns of the partitin vector part2
    for j = 1:M
        B_ij = block_element(i,j,B,part1,part2);
        % Calculate the Kronecker product of A and B_ij
        row_i = [row_i, kron(A,B_ij)];
    end
    % Add the row_i to the resultant matrix C
    C = [C;row_i];
end

end

% Helper function to get the block element of matrix B
function [B_xy] = block_element(x,y,B,part1,part2)
    % Calculate the indices for the block element
    yy_0 = sum(part2(1:y)) + 1 - part2(y);
    yy_1 = sum(part2(1:y));

    xx_0 = sum(part1(1:x)) + 1 - part1(x);
    xx_1 = sum(part1(1:x));

    % Get the block element from the matrix B
    B_xy = B(xx_0:xx_1,yy_0:yy_1);
end
function [C] =Box_Product(A,B,part1,part2)
N = length(part1);
M = length(part2);
C=[];
for i = 1:N
    row_i = [];
    for j = 1:M
        B_ij = block_element(i,j,B,part1,part2);
        row_i= [row_i, kron(A,B_ij)];
    end
    C = [C;row_i];
end

    function [B_xy] = block_element(x,y,B,part1,part2)

        yy_0 = sum(part2(1:y)) + 1 - part2(y);
        yy_1 = sum(part2(1:y));

        xx_0 = sum(part1(1:x)) + 1 - part1(x);
        xx_1 = sum(part1(1:x));

        B_xy = B(xx_0:xx_1,yy_0:yy_1);

    end
end



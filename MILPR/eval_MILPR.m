function [S_eval] = eval_MILPR(Sr,xv,x_star,k,h,f,invf,K_Type,fdata,option)

N =  size(x_star,2);
[d,~,~] = size(Sr);
S_eval = zeros(d,d,N);
x0 = mean(x_star,2);


    if strcmp(option, 'mean')

        Reg =  MILPR(Sr,xv,x0,k,h,f,invf,K_Type,fdata);
        for j = 1:N
             S_eval(:,:,j) = Reg(xv(:,j));
        end
            
    elseif strcmp(option, 'exhaustive')

        for j = 1:N
            Reg =  MILPR(Sr,xv,x_star(:,j),k,h,f,invf,K_Type,fdata);
            S_eval(:,:,j) = Reg(xv(:,j));
        end

    end

end

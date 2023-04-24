 function [S,x] = gen_functional_spd(Cdim,Sdim,Rdim)
        x =zeros(Cdim,Sdim);
        for j = 1:Cdim
             x(j,:) = linspace(-rand(1),rand(1),Sdim);
        end
        Sfun = @(x) expm(kron(([-0.1*(x+0.1),0.2*(x+0.1),sin(0.75*x);...
    0.2*(x+0.1),0.6*(x+0.1),-0.4*(x+0.1);...
    sin(0.75*x),-0.4*(x+0.1),0.5*(x+0.1)]),eye(Rdim)));
        S = zeros(3*Rdim,3*Rdim,Sdim);
        s = sum(x,1);
        for i = 1:Sdim
            S(:,:,i) = Sfun(s(i));
        end
    end
function [Manifold]=CholeskyManifold(dim)
% This function creates the Cholesky manifold with dimension dim
%------------------
% Symbolic Variables 
%Information on the Cholesky Manifold
x = sym('x');   
y = sym('y');   
t = sym('t');   
X = sym('X');   
Y = sym('Y');   
T = sym('T');   
NSample = sym('NSample'); 
%---------------------------------------
%----------------------------------------------------------
%% Computational convenient expression of the basic operations in the Manifold
%
Manifold.Diag = @(x) MDiag(x);    % Diag operator
Manifold.Lowp = @(x) MLowp(x);    % Strictly low part 
Manifold.Exp  = @(x) MExp(x);     % Exp of Diag matrices
Manifold.Log  = @(x) MLog(x);     % Log of Diag Matrices
Manifold.FNorm = @(x) FNorm(x);   % F norm
Manifold.KGauss = @(x,y,t) KGauss(x,y,t);
Manifold.KEppan = @(x,y,t) epan_kernel(x,y,t);
Manifold.Vec=@(x) HalfVec(x);
%----------------------------------------------------------
%----------------------------------------------------------
%% Riemannian operator exp, log maps, parallel transport, 
% Observation: x is the reference point on the Cholesky manifold
%Exp Map at x
Manifold.RiemOp.Exp = @(x,t) RiemExp(x,t);
%Log Map at x
Manifold.RiemOp.Log = @(x,y) RiemLog(x,y);
%Parallel Transport for x to y
Manifold.RiemOp.Parallel = @(x,y,t) RiemParll(x,y,t);
%Geodesic distance x to y 
Manifold.RiemOp.Dist = @(x,y) RiemDist(x,y);
%Tangent space product
Manifold.RiemOp.TanProd = @(x,y,t) tanspace_prod(x,y,t);
%Generate Base in tangent space on the Cholesky manifold 
Manifold.RiemOp.TanBase = @(x) TanBase(x);
% Generate metric tensor at point x with tangent base y
Manifold.RiemOp.MetricTensor2 = @(x) metrictensor2(x);
%-----------------------------------------------------------
%-----------------------------------------------------------
%% Basic Statistic
%Generate Random Sample 
Manifold.BasicStat.RandomSample = @(x,NSample) GenRandomSample(x,NSample);
%Kacher Mean  for the data X
Manifold.BasicStat.KM_Cholesky = @(X) Kacher_Mean_Cholesky(X);
Manifold.BasicStat.KM_HPD = @(X) Kacher_Mean_Cholesky(X)*Kacher_Mean_Cholesky(X)';
%ERROR X and Y
Manifold.BasicStat.Error.RMS = @(X,Y) Error_RMS(X,Y);
Manifold.BasicStat.Error.RMSLOG = @(X,Y) Error_RMSLOG(X,Y);
Manifold.BasicStat.Error.AGM = @(X,Y) Error_AGM(X,Y);
Manifold.BasicStat.Error.AGM_LOG=@(X,Y) Error_AGM_LE(X,Y);
Manifold.BasicStat.Error.AGM_E=@(X,Y) Error_AGM_E(X,Y);
Manifold.BasicStat.Error.AGM_AI=@(X,Y) Error_AGM_AI(X,Y);
%Generate Regression problem
%This function creates a regression problem  
Manifold.BasicStat.RegressionProblem.Data = @(x,y,t) gen_functional_spd(x,y,t);
% Manifold.BasicStat.RegressionProblem.Kernel = ()
%------------------------------------------------------------
%------------------------------------------------------------
%% Data manipulation
%Transform X data on the Cholesky manifold to the HPD
Manifold.Manipulate.Chol2HPD = @ (X) Chol2HPD(X);
%Transform X data on the HPD manifold to the CholeskyM
Manifold.Manipulate.HPD2Chol = @ (X) HPD2Chol(X);
%Transform X data on the HPD manifold to the Euclidean log 
Manifold.Manipulate.HPD2Log = @ (X) HPD2Log(X);
%Delete information at point Xi
Manifold.Manipulate.DeleteInfo = @(X,x,t) eraseinfo(X,x,t);
%Riemannian Vectorization on the HPD at the Kacher Mean
Manifold.Manipulate.RiemVec = @ (X) RiemannianVec(X);
%Riemannian projection into the tangent space at the Kacher mean
Manifold.Manipulate.RiemProj_KM = @ (X) Project_Karcher(X);;
%Riemannian Vectorization on the HPD on the Identity
Manifold.Manipulate.RiemVecIdentity = @ (X) RiemannianVecIdentity(X);
%Riemannian projection into the tangent space at the identity
Manifold.Manipulate.RiemProj_Identity = @ (X,Y) Parallel2Identity(X,Y);
% ------------------------------------------------------------
%------------------------------------------------------------
%% Auxx functions
    %-------
    function [D] = MDiag(A)
        [d,~] = size(A);
        D = eye(d).*A;
    end
    %-------
    function [L] = MLowp(A)
        [d,~] = size(A);
        L = tril(ones(d)-eye(d)).*A;
    end
    %-------
    function [E] = MExp(A)
       E = diag(exp(diag(A)));
    end
    %-------
    function [L] = MLog(A)
       L = diag(log(diag(A)));
    end
    %-------
    function [FN] = FNorm(A)
        FN=trace(A*A');
    end
    %-------
    function [Fp] = Fprod(A,B)
        Fp=trace(A*B');
    end
    %-------
    function [Inv] = FInv(A)
        Inv=diag(diag(A).^(-1));
    end
    %-------
    function [Vec]=HalfVec(A)
        [d,~]=size(A);
        Vec=diag(A);
        for i=1:d
            for j=1:d
                if i>j
                    Vec=[Vec;A(i,j)];
                end
            end
        end
    end
    %-------
%% Riemannian Function 
   
    %------- 
    % inner product in the target space 
    
    function [p]=tanspace_prod(a,b,t)
        p=Fprod(MLowp(a),MLowp(b))+Fprod(FInv(MDiag(t))*MDiag(a),FInv(MDiag(t))*MDiag(b));
    end
    %-------
    
    function [E] = RiemExp(a,b)
        E = MLowp(a + b) + MDiag(a)*MExp(MDiag(b)*FInv(MDiag(a)));
    end
    
    %-------
    
    function [E] = RiemLog(a,b)
       E = MLowp(b - a) + MDiag(a)*MLog(MDiag(b)*FInv(MDiag(a)));
    end
    
    %-------
    function [P] = RiemParll(a,b,tt)
       P = MLowp(tt) + MDiag(b)*MDiag(tt)*FInv(MDiag(a));
    end
    
    %-------
    function [D] = RiemDist(a,b)
        D = sqrt(FNorm(MLowp(a - b)) + FNorm(MLog(MDiag(a)) - MLog(MDiag(b))));
    end
    %-------
    function [E] = TanBase(mdim)
        dim2 = mdim*(mdim+1)/2;
        position = cell(1,dim2);
        E = zeros(mdim,mdim,dim2);
        k = 1;
        for i = 1:mdim
            for j = 1:i
                position{k} = [i,j];
                k = k+1;
            end
        end
        for i = 1:dim2
             pos = position{i};
             ii = pos(1);
             jj = pos(2);
             E(ii,jj,i) = 1;
        end
    end
   %-------
     
   % Generate Metrix tensor 
 
    function [G,d] = metrictensor2(a)
        [mdim,~] = size(a);
        nd = mdim*(mdim+1)/2;
        G = zeros(nd,nd);
        d = 1;
        tic
        for i = 1:nd
            if  i <= mdim
                G(i,i) = a(i,i)^(-2);
                d = d*G(i,i);
            else 
                G(i,i) = 1;
            end 
        end
        toc
    end
%% Auxx Statistical func
    %-------
    function [kt] = Kacher_Mean_Cholesky(A)
       [d,~,N] = size(A);
       ktl = zeros(d);
       ktd = zeros(d);
       tic
       for j = 1:N
           ktl = ktl + MLowp(A(:,:,j));
           ktd = ktd +  MLog(MDiag(A(:,:,j)));
       end
       ktl = ktl/N;
       ktd = ktd/N;
       kt = ktl + MExp(ktd);
       toc
    end
    %-------
    function [SData] = GenRandomSample(d,N)
             SData = zeros(d,d,N);
             for j = 1:N
                 s = rand(d);
                 SData(:,:,j) = MLowp(s) + abs(MDiag(s));
             end
    end
    %-------
    %-------
    function [E,R] = Error_AGM(A,B)
        [~,~,N] = size(A);
        E = 0;
        R=zeros(1,N);
        for j = 1:N
            R(j) = RiemDist(A(:,:,j),B(:,:,j));
            E = E+R(j)^2;
        end
        E = E/N;
    end
    
  %---------
    function  D = LogEucDist(A,B)
               D = sqrt(FNorm(logm(A)-logm(B)));
    end
    
    function [E,R] = Error_AGM_LE(A,B)
        [~,~,N] = size(A);
        E = 0;
        R=zeros(1,N);
        for j = 1:N
            R(j) = LogEucDist(A(:,:,j),B(:,:,j));
            E = E+R(j)^2;
        end
        E = E/N;
    end
    
     function [E,R] = Error_AGM_E(A,B)
        [~,~,N] = size(A);
        E = 0;
        R=zeros(1,N);
        for j = 1:N
            R(j) = sqrt(FNorm(A(:,:,j)-B(:,:,j)));
            E = E+R(j)^2;
        end
        E = E/N;
    end
%-------
   
    
        
    %% Data Manipulation Auxx
    %-------
    function [PA] = Chol2HPD(A)
        [d,~,N] = size(A);
        PA = zeros(d,d,N);
        for j = 1:N
            PA(:,:,j) = A(:,:,j)*A(:,:,j)';
        end
    end
    %-------
    function [A] = HPD2Chol(PA)
        [d,~,N] = size(PA);
        A = zeros(d,d,N);
        for j = 1:N
            PA(:,:,j)=Regularization(PA(:,:,j));
            A(:,:,j) = chol(PA(:,:,j))';
        end   
       
    end
    %-------
    function [A] = HPD2Log(PA)
        [d,~,N] = size(PA);
        A = zeros(d,d,N);
        for j = 1:N
            A(:,:,j) = logm(PA(:,:,j));
        end   
       
    end
    %Riemannian vectorization at the Karcher mean
    function [Tan,K] = RiemannianVec(A)
        % A is on Cholesky Manifold
        % for A in HPD manifold use first HPD2CHol;
        tic
        K = Kacher_Mean_Cholesky(A);
        [d,~,N] = size(A);
        Ndim = floor(d*(d+1)/2);
        Tan = zeros(Ndim,N); % Tangent data at Karcker Mean on HPD
        for j = 1:N
            Tj = RiemLog(K,A(:,:,j));
            Tan(:,j) = HalfVec(K*Tj'+Tj*K');
        end
        toc
    end
   % Projection of the data on the Cholesky Manifold on the tangent space
   % at the Karcher mean  on the Choelsky manifold
    function [Tan,K] = Project_Karcher(A)
        tic
        K = Kacher_Mean_Cholesky(A);
        [d,~,N] = size(A);
        Ndim = d;
        Tan = zeros(Ndim,Ndim,N); % Tangent data at Karcker Mean on HPD
        for j = 1:N
            Tan(:,:,j) = RiemLog(K,A(:,:,j));
        end
        toc
    end
   %Transport the Tangent space from the Karcher mean to the tangent space
   %at the identity
    function [TanI] = Parallel2Identity(Tan,K)
        % input data Tan its is in the tangent space of the Karcher Mean
        %Mean on the Cholesky manifold
        [d,~] = size(K);
        I = eye(d);
        [~,~,N] = size(Tan);
        Ndim = floor(d*(d+1)/2);
        TanI = zeros(Ndim,N);
       
        for j = 1:N
            Tj = RiemParll(K,I,Tan(:,:,j));
            TanI(:,j) = HalfVec(Tj'+Tj);
        end
       
    end
   
    % Riemannian Vec at identity
    function [TanI] = RiemannianVecIdentity(A)
        % input A is on the Cholesky Manifold 
        [Tan,K] = Project_Karcher(A);
        TanI = Parallel2Identity(Tan,K);
    end
   
    %--------
  
    function [S,x] = gen_functional_spd(Cdim,Sdim,Rdim)
        x =zeros(Cdim,Sdim);
        for j = 1:Cdim
             x(j,:) = linspace(-rand(1),and(1),Sdim);
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
    %----------
    function [K] = KGauss(a,b,h)
        [dd,~] = size(a);
         K = exp(-0.5*FNorm(a-b)/h^2)/(sqrt(2*pi)*h);
    end
    function u = KEppan(a,b,h)
    % u=max(eps,3/4*(1-sum(u.^2,3)));
    u = max(0,3/4*(1-FNorm(a-b)/h^2));
    end
    

end

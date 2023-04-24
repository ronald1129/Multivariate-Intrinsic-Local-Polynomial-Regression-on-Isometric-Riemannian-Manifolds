function [s,v,t,E,T,Estimation] = simualtion_milpr(Cdim,Mdim,Rdim,Simulation_dim,Sample_Dim,eval_option,noise_level)


M=CholeskyManifold(3);
y=sym('y');
N = Sample_Dim;

EL = zeros(1,Simulation_dim);
EC = zeros(1,Simulation_dim);
EA = zeros(1,Simulation_dim);

tL = zeros(1,Simulation_dim);
tC = zeros(1,Simulation_dim);
tA = zeros(1,Simulation_dim);

t = zeros(1,3);
s = zeros(1,3);
v = zeros(1,3);

options = optimset('Display','iter');

for i=1:Simulation_dim

    [S,xv]= spd_high_dimensional(Cdim,Mdim,N); %M.BasicStat.RegressionProblem.Data(Cdim,N,Mdim/3);
    
    % Noise model
 
    Sr=High_Dimension_Riem_Noise(S,noise_level); 
    %--------------------------------,

    x0 = mean(xv,2);
    sigma = 300;%max(std(xv'));
    
  
    t=tic;

    LC = @(t) MLowp(chol(t)') + MLog(MDiag(chol(t)'));
    iLC = @(t) MLowp(t) +  expm(MDiag(t));
    invLC =@(t) (iLC(t)*iLC(t)');

    LCdata = zeros(Mdim,Mdim,N);


    for j = 1:N
        LCdata(:,:,j) = LC(Sr(:,:,j));
    end

    
    %Bandwidth selection using LOOCV
      cv_LC = @(h) cv_dist(Sr,xv,x0,Rdim,h,LC,invLC,'KGauss','AffineI',LCdata,eval_option);
      hfC = fminbnd(cv_LC,0,3*sigma,options);
      %hfC=1000;
    
    % Regression

    Sp = eval_MILPR(Sr,xv,xv,Rdim,hfC,LC,invLC,'KGauss',LCdata,eval_option);
    
    tC(i)=toc(t);

       
    %%
    clear t
    t=tic;
    MA = positive_definite_karcher_mean(Sr);
   LMA = logm(MA);
   MA1 = expm(-0.5*LMA);
    MA2 = expm(0.5*LMA);

    AA = @(t) MA2*logm(MA1*t*MA1)*MA2;
    invAA =@(t) MA2*expm(MA1*t*MA1)*MA2;

    AAdata = zeros(Mdim,Mdim,N);
    for j = 1:N
        AAdata(:,:,j) =  AA(Sr(:,:,j));
    end

    
    % Bandwidth selection using LOOCV
      
      cv_AA = @(h) cv_dist(Sr,xv,x0,Rdim,h,AA,invAA,'KGauss','AffineI',AAdata,eval_option);
      hfA = fminbnd(cv_AA,0,3*sigma,options);
      %hfA=1000;

   
    % Regression

    SpA = eval_MILPR(Sr,xv,xv,Rdim,hfA,AA,invAA,'KGauss',AAdata,eval_option);

    tA(i)=toc(t);
    

    %%

    clear t
    t=tic;

    %ML = LogEuc(Sr);
    LL =@(t) logm(t);
    invLL =@(t) expm(t);

    LLdata = zeros(Mdim,Mdim,N);
    for j = 1:N
        LLdata(:,:,j) = LL(Sr(:,:,j));
    end


    % Bandwidth selection using LOOCV
      
      cv_LL = @(h) cv_dist(Sr,xv,x0,Rdim,h,LL,invLL,'KGauss','AffineI',LLdata,eval_option);
      hfL = fminbnd(cv_LL,0,3*sigma,options);
      %hfL=1000;

    % Regression
    SpL = eval_MILPR(Sr,xv,xv,Rdim,hfL,LL,invLL,'KGauss',LLdata,eval_option);

    tL(i)=toc(t);
    
    

    %%
    
    [~,EL(i)] = Intrinsic_Error(S,SpL,'AffineI');
    [~,EC(i)] = Intrinsic_Error(S,Sp,'AffineI');
    [~,EA(i)] = Intrinsic_Error(S,SpA,'AffineI');

end

E=[EC;EL;EA];
T=[tC;tL;tA];

% Error
 s = mean(E,2)';

%Variance

%  v(1)=std(EC);
%  v(2)=std(EL);
%  v(3)=std(EA);

%Time
 t = mean(T,2)';

% Estimated Data

Estimation.LC = Sp;
Estimation.LL = SpL;
Estimation.LA = SpA;
Estimation.True = S;
Estimation.Noise = Sr;


end




%% Figure 4 A

M =CholeskyManifold;

Cdim = 1;   % Num of covariates
Rdim = 1;   % Regression degree 0 or 1
pMdim = 1;  % Manifold dimension 3xpMdim
N = 50; %100;  % Sample size by test 
Ntest = 1; % Number of test 

boxplots = 1;
%y = sym('y');
l1 = sym('l1');
l2 = sym('l2');
Mdim = 3*pMdim;
eval_option = 'exhaustive';
noise_level = 0.5; % var


%Simulated Regression

[error,~,~,~,~,Estimation] = simualtion_milpr(Cdim,Mdim,Rdim,Ntest,N,eval_option,noise_level);
S = Estimation.True;
Sr = Estimation.Noise;
LLCf = Estimation.LC;
LLLf = Estimation.LL;
LLAf = Estimation.LA;

%% Ellipsoidal Representation 
if pMdim ==1 && N < 101
    delta = 1;
        
    subplot(5,1,1)
    plot_spd_eig(S,delta);title('Ture SPD');
    axis off
    
    subplot(5,1,2)
    plot_spd_eig(Sr,delta);title('SPD + Noise')
    axis off
    subplot(5,1,3)
    
    plot_spd_eig(LLCf,delta);title('Log - Chol')
    axis off
    subplot(5,1,4)
    
    plot_spd_eig(LLLf,delta);title('Log - Eucl')
    axis off
    subplot(5,1,5)
    
    plot_spd_eig(LLAf,delta);title('AI')
    axis off
end
    
error


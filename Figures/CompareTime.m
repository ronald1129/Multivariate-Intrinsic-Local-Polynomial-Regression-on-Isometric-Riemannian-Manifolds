
% parameters
pMdim = 1;
Cdim = 1;   % Num of covariates
Mdim = 3*pMdim;
Rdim = 1;   % Regression degree 0 or 1

Simulation_dim = 100;
Sample_Dim = 100; 

eval_option = 'mean';
methods = {'LC','LE','AI'};

% run methods
[s,v,t,E,T] = simualtion_milpr(Cdim,Mdim,Rdim,Simulation_dim,Sample_Dim,eval_option,0.5);
E = E';
T = T';

%% plot results
tbl = table(E(:,1),E(:,2),E(:,3),'VariableNames',methods);

set(gcf,'color',[0.8 0.9 0.8]);
al_goodplot(tbl.LC,1,[],[],[],[],[],1);
al_goodplot(tbl.LE,2,[],[],[],[],[],1);
al_goodplot(tbl.AI,3,[],[],[],[],[],1);

xticks([1 2 3])
xticklabels(methods)
ylabel('RSME');
title('RMSE for all methods')

figure(2)
tbl = table(T(:,1),T(:,2),T(:,3),'VariableNames',methods);

set(gcf,'color',[0.8 0.9 0.8]);
al_goodplot(tbl.LC,1,[],[],[],[],[],1);
al_goodplot(tbl.LE,2,[],[],[],[],[],1);
al_goodplot(tbl.AI,3,[],[],[],[],[],1);

xticks([1 2 3])
xticklabels(methods)
ylabel('CPU-Time');
title('bla')



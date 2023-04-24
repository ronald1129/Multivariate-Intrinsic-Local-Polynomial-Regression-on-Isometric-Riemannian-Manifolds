%% Figure 5 B

% Clearing the workspace and closing any open figures
clear all
close all

% Define the dimensions of Mdim and Cdim
Mdim = 3:3:18;
Cdim = 1:1:6;

% Initializing matrices to hold the results of simulations
EC = zeros(6,6);
TC = zeros(6,6);
EL = zeros(6,6);
TL = zeros(6,6);
EA = zeros(6,6);
TA = zeros(6,6);

% Initializing simulation_stage and option variables
simulation_stage = 1;
option = 'mean';

% Loop through Cdim and Mdim
for i = 1:6
for j = 1:6
    
    % Run the simulation_milpr function with given parameters and store the results
    [e,v,t,~,~,~] = simualtion_milpr(Cdim(i),Mdim(j),1,2,10,option,0.5);
    
    % Assign the results to the appropriate matrix
    EC(i,j) = e(1);
    TC(i,j) = t(1);
    EL(i,j) = e(2);
    TL(i,j) = t(2);
    EA(i,j) = e(3);
    TA(i,j) = t(3);
    
    % Update the simulation stage counter
    simulation_stage = simulation_stage + 1;
end
end

% Define the limits for the color scale of the surfaces
zlimE0 = min([min(min(log10(EC))),min(min(log10(EL))),min(min(log10(EA)))]);
zlimE1 = max([max(max(log10(EC))),max(max(log10(EL))),max(max(log10(EA)))]);
zlimT0 = min([min(min(log10(TC))),min(min(log10(TL))),min(min(log10(TA)))]);
zlimT1 = max([max(max(log10(TC))),max(max(log10(TL))),max(max(log10(TA)))]);
zlimTO0 = min([min(min(log10(EC.*TC)/(zlimE1 + zlimT1))),min(min(log10(EL.*TL)/(zlimE1 + zlimT1))),min(min(log10(EA.*TA)/(zlimE1 + zlimT1)))]);
zlimTO1 = max([max(max(log10(EC.*TC)/(zlimE1 + zlimT1))),max(max(log10(EA.*TA)/(zlimE1 + zlimT1)))]);

% Plotting the results 
%LogC
subplot(3,3,1)
surf(Cdim,Mdim,log10(EC),'FaceColor','texturemap','EdgeColor','none');
xlabel('CDim')
ylabel('MDim')
caxis([zlimE0 zlimE1])
zlim([zlimE0 zlimE1])
ylim([3 18 ])
xlim([1 6])
colorbar()
 colormap(turbo)


subplot(3,3,2)
surf(Cdim,Mdim,log10(TC),'FaceColor','texturemap','EdgeColor','none');
xlabel('CDim')
ylabel('MDim')
caxis([zlimT0 zlimT1])
zlim([zlimT0 zlimT1])
ylim([3 18 ])
xlim([1 6])
colorbar()
 colormap(turbo)



subplot(3,3,3)
surf(Cdim,Mdim,log10(EC.*TC)/(zlimE1+zlimT1),'FaceColor','texturemap','EdgeColor','none');
xlabel('CDim')
ylabel('MDim')
zlim([zlimTO0 1])
ylim([3 18 ])
xlim([1 6])
colorbar()
 colormap(turbo)
 caxis([zlimTO0  zlimTO1]);


%LogEuc

subplot(3,3,4)
surf(Cdim,Mdim,log10(EL),'FaceColor','texturemap','EdgeColor','none');
xlabel('CDim')
ylabel('MDim')
caxis([zlimE0 zlimE1])
zlim([zlimE0 zlimE1])
ylim([3 18 ])
xlim([1 6])
colorbar()


colormap(turbo)

subplot(3,3,5)
surf(Cdim,Mdim,log10(TL),'FaceColor','texturemap','EdgeColor','none');
xlabel('CDim')
ylabel('MDim')
caxis([zlimT0 zlimT1])
zlim([zlimT0 zlimT1])
ylim([3 18 ])
xlim([1 6])
colorbar()
colormap(turbo)


subplot(3,3,6)
surf(Cdim,Mdim,log10(EL.*TL)/(zlimE1+zlimT1),'FaceColor','texturemap','EdgeColor','none');
xlabel('CDim')
ylabel('MDim')
caxis([zlimTO0 zlimTO1])
zlim([zlimTO0 1])
ylim([3 18 ])
xlim([1 6])
colorbar() 
colormap(turbo)
caxis([zlimTO0  zlimTO1]);

%AffineI

subplot(3,3,7)
surf(Cdim,Mdim,log10(EA),'FaceColor','texturemap','EdgeColor','none');
xlabel('CDim')
ylabel('MDim')
caxis([zlimE0 zlimE1])
zlim([zlimE0 zlimE1])
ylim([3 18 ])
xlim([1 6])
colorbar()
colormap(turbo)

subplot(3,3,8)
surf(Cdim,Mdim,log10(TA),'FaceColor','texturemap','EdgeColor','none');
xlabel('CDim')
ylabel('MDim')
caxis([zlimT0 zlimT1]);
zlim([zlimT0 zlimT1]);
ylim([3 18 ])
xlim([1 6])
colorbar()
 colormap(turbo)


subplot(3,3,9)
surf(Cdim,Mdim,log10(EA.*TA)/(zlimE1+zlimT1),'FaceColor','texturemap','EdgeColor','none');
xlabel('CDim')
ylabel('MDim')
zlim([zlimTO0 1])
ylim([3 18 ])
xlim([1 6])
colorbar()
 colormap(turbo)
caxis([zlimTO0  zlimTO1]);

%% Figure 5 A


M = 6;
C = 6;


Error = [reshape(EC(1:C,1:M),[M*C,1]);reshape(EL(1:C,1:M),[M*C,1]);reshape(EA(1:C,1:M),[M*C,1])];
Time = [reshape(TC(1:C,1:M),[M*C,1]);reshape(TL(1:C,1:M),[M*C,1]);reshape(TA(1:C,1:M),[M*C,1])];
Error = Error./max(Error);
Time = Time./max(Time);
legend2=[1*ones(C*M,1);2*ones(C*M,1);3*ones(C*M,1)];
hold on
gscatter(Error,Time,legend2)
xlim([0 1])
ylim([0 1])



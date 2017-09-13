close('all')
%% User Input
% Basis b
% 1 = Boyd
% 2 = Thew
b=1;

closefig=0; % Close figures after completion?

simIndex=1; % Choose index for detailed report of this simulation

% Choose Active MaximumLikelihood Solvers
ga=1;
fmin=1;
lstsqr=1;

solverdetails=[61,54,34]; % Choose simulations for which detailed solver analysis is done


outputFolder = 'OutputSolverTest100_HighNoise_WNM08/';

if (exist(outputFolder) == false)
    mkdir(outputFolder);
end

%% Complete tomographic set: Measurement basis

% Boyd or Thew
switch b
    case 1 % Boyd
        chi={[1;0],[0;1],1/sqrt(2)*([1;0]+[0;1]),1/sqrt(2)*([1;0]+[0;1]*exp(1i*pi/2.0))};
        basis='Boyd';
    case 2 % Thew
        phase=exp(2i*pi/3.0);
        chi={[1;0],1/sqrt(2)*([1;0]+[0;1]),1/sqrt(2)*([1;0]+[0;1]*phase),1/sqrt(2)*([1;0]+[0;1]*phase^2)};
        %chi={[1;0],[0;1],1/sqrt(2)*([1;0]+[0;1]),1/sqrt(2)*([1;0]+[0;1]*exp(2i*pi/3.0))};
        basis='Thew';
end

% Kronecker product of measurment basis
Chi=cell(16,1);
for ii=1:4
    for jj=1:4
       Chi{4*(ii-1)+jj}=kron(chi{jj},chi{ii});
    end
end

%% Extended Tomographic Set: Basis for NS projection
chi{5}=1/sqrt(2)*([1;0]-[0;1]);
chi{6}=1/sqrt(2)*([-1;1i]);

ChiProj=cell(36,1);
for ii=1:6
   for j=1:6
      ChiProj{6*(ii-1)+j}=kron(chi{ii},chi{j}); 
   end
end

% MM(x,a) resp. MM(y,b)
MM=cell(3,2);
MM{1,1}=chi{1};
MM{1,2}=chi{2};
MM{2,1}=chi{3};
MM{2,2}=chi{5};
MM{3,1}=chi{4};
MM{3,2}=chi{6};

%% Initial State
plotRho(rho_Ideal,[outputFolder,'rho_targetState.pdf'],'\textbf{Target State $\rho_{ideal}$}');

%% Check Signaling of Noisy and Projected Behaviours
nonsignalingList(:,1) = num2cell(cellfun(@checkNonSignaling,{P_Noise{:,1}})); % Noisy Behaviour
nonsignalingList(:,2) = num2cell(cellfun(@checkNonSignaling,{P_Noise{:,2}})); % Projected Behaviour

%% Calculate Traces
for i=1:8
    traceRhoList(:,i) = cellfun(@(A) trace(A), {rho_Noise{:,i}});
    traceRho2List(:,i) = cellfun(@(A) trace(A^2), {rho_Noise{:,i}}); 
end

%% Calculate Delta P
delta_P(:,1)=cellfun(@(A) A-P_Ideal,{P_Noise{:,1}}, 'UniformOutput', false);
delta_P(:,2)=cellfun(@(A) A-P_Ideal,{P_Noise{:,2}}, 'UniformOutput', false);
delta_P_Vector=cell2mat(cellfun(@(A) reshape(A-P_Ideal,36,1), {P_Noise{:,1}}, 'UniformOutput', false));
delta_P_Vector_Proj=cell2mat(cellfun(@(A) reshape(A-P_Ideal,36,1), {P_Noise{:,2}}, 'UniformOutput', false));

%% Calculate Fidelity
for i = 1:8
    fidelityList(:,i)=(cellfun(@(A) Fidelity(A,rho_Ideal),{rho_Noise{:,i}}));
end
    

%% Calculate Pideal + (1-x)dP
deltaP = P_Noise{simIndex,1}-P_Noise{simIndex,2}; % P_Noise - P_Proj_LI
x=linspace(0,1,11);


for i = 1:11
    P_List{i} =  P_Noise{simIndex,2} + (1-x(i))*deltaP;
    rhoList{i}= qst_linearinversion(Chi,P_QST_Selection(P_List{i}));
    fid(i) = Fidelity(rhoList{i},rho_Ideal);
    
end

%%

%% Plot Reconstruced States
plotRho(rho_Noise{simIndex,1},[outputFolder,'rho_noiseState_li.pdf'],{'\textbf{Reconstructed State $\rho_{Noise}$}','Linear Inversion'});
plotRho(rho_Noise{simIndex,2},[outputFolder,'rho_projState_li.pdf'],{'\textbf{Reconstructed State after Projection $\rho_{Proj}$}','Linear Inversion'});

plotRho(rho_Noise{simIndex,3},[outputFolder,'rho_noiseState_ml_fmin.pdf'],{'\textbf{Reconstructed State $\rho_{Noise}$}','Maximum Likelihood'});
plotRho(rho_Noise{simIndex,4},[outputFolder,'rho_projState_ml_fmin.pdf'],{'\textbf{Reconstructed State after Projection $\rho_{Proj}$}','Maximum Likelihood'});

plotRho(rho_Noise{simIndex,5},[outputFolder,'rho_noiseState_ml_ga.pdf'],{'\textbf{Reconstructed State $\rho_{Noise}$}','Maximum Likelihood'});
plotRho(rho_Noise{simIndex,6},[outputFolder,'rho_projState_ml_ga.pdf'],{'\textbf{Reconstructed State after Projection $\rho_{Proj}$}','Maximum Likelihood'});

plotRho(rho_Noise{simIndex,7},[outputFolder,'rho_noiseState_ml_lstsqr.pdf'],{'\textbf{Reconstructed State $\rho_{Noise}$}','Maximum Likelihood'});
plotRho(rho_Noise{simIndex,8},[outputFolder,'rho_projState_ml_lstsqr.pdf'],{'\textbf{Reconstructed State after Projection $\rho_{Proj}$}','Maximum Likelihood'});


plotRho(rho_Noise{simIndex,1}-rho_Ideal,[outputFolder,'rho_deltaNoiseState_li.pdf'],{'\textbf{Difference to Target State before Projection $\rho_{Noise}-\rho_{Ideal}$}','Linear Inversion'});
plotRho(rho_Noise{simIndex,2}-rho_Ideal,[outputFolder,'rho_deltaProjState_li.pdf'],{'\textbf{Difference to Target State after Projection $\rho_{Proj}-\rho_{Ideal}$}','Maximum Likelihood'});
plotRho(rho_Noise{simIndex,3}-rho_Ideal,[outputFolder,'rho_deltaNoiseState_ml.pdf'],{'\textbf{Difference to Target State before Projection $\rho_{Noise}-\rho_{Ideal}$}','Maximum Likelihood'});
plotRho(rho_Noise{simIndex,4}-rho_Ideal,[outputFolder,'rho_deltaProjState_ml.pdf'],{'\textbf{Difference to Target State after Projection $\rho_{Proj}-\rho_{Ideal}$}','Maximum Likelihood'});

%% Plot Behaviour
plotBehaviour(P_Ideal,P_Noise, simIndex,[outputFolder,'behaviour.pdf'],['\textbf{Behaviour $P(ab| xy)$ of Simulation ',num2str(simIndex),'}']);

P = cell2mat(cellfun(@(A) {reshape(A,36,1)}, {P_Noise{:,1}})) % Behaviour with Noise
P_Noise_mean{1,1}=mean(P,2)
P = cell2mat(cellfun(@(A) {reshape(A,36,1)}, {P_Noise{:,2}})) % Behaviour with Noise after Projection
P_Noise_mean{1,2}=mean(P,2)

plotBehaviour(P_Ideal, P_Noise_mean, 1, [outputFolder,'behaviour_mean.pdf'],'\textbf{Mean Behaviour $P(ab| xy)$}')

%% Check Poisson Noise
figure;
pIndex=10;
countFactor = 100;
maxCount = 1;
lambda_Target = P_Ideal(pIndex)*countFactor/ maxCount;
%hist(abs(delta_P_Matrix(pIndex,:))')
nBins = max(abs(delta_P_Vector(pIndex,:))')-min(abs(delta_P_Vector(pIndex,:))')+1;
histfit(abs(delta_P_Vector(pIndex,:))'*countFactor,10,'poisson')
pd = fitdist(abs(delta_P_Vector(pIndex,:))','poisson');
h=legend(['$\lambda_{fit}=',num2str(pd.lambda), '$ $\lambda_{target}=', num2str(lambda_Target),'$'])
set(h,'Interpreter','latex')

for jj=1:36 
    pd = fitdist(abs(delta_P_Vector(jj,:))','poisson');
    pd_proj = fitdist(abs(delta_P_Vector_Proj(jj,:))','poisson');
    var_fit(jj)=pd.lambda*countFactor;
    var_fit_proj(jj)=pd_proj.lambda;

    var_P(jj)=var(delta_P_Vector(jj,:)*countFactor);
    var_P_proj(jj)=var(delta_P_Vector_Proj(jj,:)*countFactor);
end

fig=figure;
plot(var_P)
hold on
plot(var_P_proj)
plot(var_fit,'--')
plot(var_fit_proj,'--')
plot(reshape(P_Ideal,36,1),'LineWidth',2)
legend('Poisson Fit', 'Poisson Fit after NS Projection','Variance Func','Variance Func after NS Projection','P Ideal')
%saveas(fig,[outputFolder,'noise.pdf']);
print('-fillpage',[outputFolder, 'noise'],'-dpdf')

close(fig)

mean_delta_P{1}=zeros(size(delta_P{1,1}));
for ii=1:nMeasurement
    mean_delta_P{1} = mean_delta_P{1} + delta_P{ii,1};
end
mean_delta_P{1}=mean_delta_P{1} / nMeasurement;

%% Plot Fidelity and P element of NS
fig=figure('units','normalized','outerposition',[0 0 0.8 0.8]);
subplot(2,1,1)
plot(real((fidelityList(:,1))),'b')
hold on
plot(real((fidelityList(:,3))),'b--')
plot(real((fidelityList(:,5))),'b:')
plot(real((fidelityList(:,7))),'b-.')
xlabel('# simulation','interpreter','latex')
ylabel('Fidelity $\mathcal{F}$','interpreter','latex')

yyaxis right
ylabel('$P(ab|xy) \in \mathcal{NS}$','interpreter','latex')
set(gca,'ColorOrder',[0 0 1;1 0 1])
set(gca,'FontSize',12)
plot(cell2mat(nonsignalingList),'d')
ylim([-0.25 1.25])
yticks([0 1])
plt = gca;
plt.YAxis(2).Color='k';


h=legend('$\mathcal{F}(\rho_{Noise},\rho_{Ideal})$ Linear Inversion',...
    '$\mathcal{F}(\rho_{Noise},\rho_{Ideal})$ Maximum Likelihood fminsearch',...
    '$\mathcal{F}(\rho_{Noise},\rho_{Ideal})$ Maximum Likelihood ga',...
    '$\mathcal{F}(\rho_{Noise},\rho_{Ideal})$ Maximum Likelihood lstsqr',...
    '$P_{Noise} \in \mathcal{NS}$', '$P_{Proj} \in \mathcal{NS}$');
legend('Location','bestoutside')
set(h,'Interpreter','latex')

subplot(2,1,2)
plot(real((fidelityList(:,2))),'r')
hold on
plot(real((fidelityList(:,4))),'r--')
plot(real((fidelityList(:,6))),'r:')
plot(real((fidelityList(:,8))),'r-.')
xlabel('# simulation','interpreter','latex')
ylabel('Fidelity $\mathcal{F}$','interpreter','latex')


xlim([1,size(fidelityList,1)])


yyaxis right
ylabel('$P(ab|xy) \in \mathcal{NS}$','interpreter','latex')
set(gca,'ColorOrder',[0 0 1;1 0 1])
set(gca,'FontSize',12)
plot(cell2mat(nonsignalingList),'d')
ylim([-0.25 1.25])
yticks([0 1])
plt = gca;
plt.YAxis(2).Color='k';


h=legend('$\mathcal{F}(\rho_{Proj},\rho_{Ideal})$ Linear Inversion',...
    '$\mathcal{F}(\rho_{Proj},\rho_{Ideal})$ Maximum Likelihood fminsearch',...
    '$\mathcal{F}(\rho_{Proj},\rho_{Ideal})$ Maximum Likelihood ga',...
    '$\mathcal{F}(\rho_{Proj},\rho_{Ideal})$ Maximum Likelihood lstsqr',...
    '$P_{Noise} \in \mathcal{NS}$', '$P_{Proj} \in \mathcal{NS}$');
legend('Location','bestoutside')
set(h,'Interpreter','latex')



fig.PaperOrientation='landscape';
%saveas(fig,[outputFolder,'fidelity.pdf']);
print('-fillpage',[outputFolder, 'fidelity'],'-dpdf')

%close(fig)

%% Plot Fidelity Histogram

fig=figure()
edges=linspace(0.9,1.05,20);
subplot(2,1,1)
h1=histogram(real(fidelityList(:,1)),edges) % Noisy with LI
set(gca,'YLim',[0,length(fidelityList)*1.1]) 
hold on
histogram(real(fidelityList(:,5)),edges) % Noisy with ML
ax = gca;
ax.ColorOrderIndex = 1;
plot([mean(fidelityList(:,1)), mean(fidelityList(:,1))],get(gca,'YLim'),'LineStyle','--')
plot([mean(fidelityList(:,5)), mean(fidelityList(:,5))],get(gca,'YLim'),'LineStyle','--')
h=legend('$\mathcal{F}(\rho_{Noise},\rho_{Ideal})$ Linear Inversion','$\mathcal{F}(\rho_{Noise},\rho_{Ideal})$ Maximum Likelihood ga',...
    'average','average');
legend('Location','bestoutside')
set(h,'Interpreter','latex')

subplot(2,1,2)
histogram(real(fidelityList(:,2)),edges) % Proj with LI
set(gca,'YLim',[0,length(fidelityList)*1.1]) 
hold on
histogram(real(fidelityList(:,6)),edges) % Proj with ML
ax = gca;
ax.ColorOrderIndex = 1;
plot([mean(fidelityList(:,2)), mean(fidelityList(:,2))],get(gca,'YLim'),'LineStyle','--')
plot([mean(fidelityList(:,6)), mean(fidelityList(:,6))],get(gca,'YLim'),'LineStyle','--')
h=legend('$\mathcal{F}(\rho_{Proj},\rho_{Ideal})$ Linear Inversion','$\mathcal{F}(\rho_{Proj},\rho_{Ideal})$ Maximum Likelihood ga',...
    'average','average');
legend('Location','bestoutside')
set(h,'Interpreter','latex')

fig.PaperOrientation='landscape';
print('-fillpage',[outputFolder, 'fidelityHistogramm'],'-dpdf')

%% Plot Fidelity of Simulation 1
fig=figure('units','normalized','outerposition',[0 0 0.8 0.8]);

plot(x,real((fid)),'o--')
grid on
title(['Fidelity $\mathcal{F}(\rho_x,\rho_{ideal} )$ with $P_x = x P_{proj} + (1-x) P_{noise}$'],'interpreter','latex')
xlabel('x','interpreter','latex')
ylabel('Fidelity $\mathcal{F}$','interpreter','latex')
%xlim([1,size(fidelityList,1)])

set(gca,'FontSize',12)

fig.PaperOrientation='landscape';
print('-fillpage',[outputFolder, 'fidelity_sim', num2str(simIndex)],'-dpdf')
%saveas(fig,[outputFolder,'fidelity_sim', num2str(simIndex),'.pdf']);
close(fig)

%% Plot Trace
fig=figure();%'units','normalized','outerposition',[0 0 0.8 0.8]);
%plot(cell2mat(traceRhoList))
hold on
h=plot(traceRho2List)
set(h,{'LineStyle'},{'-';'-';'--';'--';'-.';'-.';':';':'})
set(h,{'Color'},{'r';'b';'r';'b';'r';'b';'r';'b'})
xlabel('# simulation','interpreter','latex')
ylabel('Trace','interpreter','latex')
xlim([1,size(fidelityList,1)])
grid on
% h=legend('Tr($\rho_{Noise}$) Linear Inversion','Tr($\rho_{Proj}$) Linear Inversion',...
%     'Tr($\rho_{Noise}$) Maximum Likelihood','Tr($\rho_{Proj}$) Maximum Likelihood',...
%     'Tr($\rho_{Noise}^2$)','Tr($\rho_{Proj}^2$)');
h=legend('Tr($\rho_{Noise}^2$) Linear Inversion','Tr($\rho_{Proj}^2$) Linear Inversion',...
    'Tr($\rho_{Noise}^2$) Maximum Likelihood fmin','Tr($\rho_{Proj}^2$) Maximum Likelihood fmin',...
    'Tr($\rho_{Noise}^2$) Maximum Likelihood ga','Tr($\rho_{Proj}^2$) Maximum Likelihood ga',...
    'Tr($\rho_{Noise}^2$) Maximum Likelihood lstsqr','Tr($\rho_{Proj}^2$) Maximum Likelihood lstsqr');
%legend('Location','bestoutside')
set(h,'Interpreter','latex')
fig.PaperOrientation='landscape';
fig.PaperPositionMode = 'auto';
%saveas(fig,[outputFolder,'trace.pdf']);
print('-fillpage',[outputFolder, 'trace'],'-dpdf')
%close(fig)

%% Plot Time
fig=figure();
avgrTime = mean(tElapsed(:,3:8),1);
plot(tElapsed(:,3:8))
title('Evaluation time for different Maximum Likelihood Solvers','interpreter','latex')
xlabel('# simulation','interpreter','latex')
ylabel('Time [sec]','interpreter','latex')
h=legend(['fminsearch on $\rho_{Noise}$, Average: ',num2str(avgrTime(1)),' sec'],...
    ['fminsearch on $\rho_{Proj}$, Average: ',num2str(avgrTime(2)),' sec'],...
    ['genetic algorithm on $\rho_{Noise}$, Average: ',num2str(avgrTime(3)),' sec'],...
    ['genetic algorithm on $\rho_{Proj}$, Average: ',num2str(avgrTime(4)),' sec'],...
    ['lstsqr on $\rho_{Noise}$, Average: ',num2str(avgrTime(5)),' sec'],...
    ['lstsqr on $\rho_{Proj}$, Average: ',num2str(avgrTime(6)),' sec']);
set(h,'Interpreter','latex')
fig.PaperOrientation='landscape';
print('-fillpage',[outputFolder, 'timesolverml'],'-dpdf')

%% Generate Convergence Plots
if (fmin==1 || ga==1 || lstsqr==1)
    formatSpec = '%d\t %f\t %s\n';

    for i=1:length(solverdetails)
        fig=figure();
        if fmin==1
            subplot(2,1,1)
            filename = 'solver_fmin_noise_';
            fileID = fopen([outputFolder, filename, num2str(i),'.txt'],'r');
            data = textscan(fileID,formatSpec);
            fclose(fileID);
            semilogy(data{:,2},'DisplayName','fmin Noise')
            hold on
            subplot(2,1,2)
            filename = 'solver_fmin_proj_';
            fileID = fopen([outputFolder, filename, num2str(i),'.txt'],'r');
            data = textscan(fileID,formatSpec);
            fclose(fileID);
            semilogy(data{:,2},'DisplayName','fmin Proj')
            hold on
        end
        if ga==1
            subplot(2,1,1)
            filename = 'solver_ga_noise_';
            fileID = fopen([outputFolder, filename, num2str(i),'.txt'],'r');
            data = textscan(fileID,formatSpec,'HeaderLines',1);
            fclose(fileID);
            semilogy(data{:,2},'DisplayName','GA Noise')
            hold on
            subplot(2,1,2)
            filename = 'solver_ga_proj_';
            fileID = fopen([outputFolder, filename, num2str(i),'.txt'],'r');
            data = textscan(fileID,formatSpec,'HeaderLines',1);
            fclose(fileID);
            semilogy(data{:,2},'DisplayName','GA Proj')
            hold on
        end
        if lstsqr==1
            subplot(2,1,1)
            filename = 'solver_lstsqr_noise_';
            fileID = fopen([outputFolder, filename, num2str(i),'.txt'],'r');
            data = textscan(fileID,formatSpec);
            fclose(fileID);
            semilogy(data{:,2},'DisplayName','LSTSQR Noise')
            hold on
            subplot(2,1,2)
            filename = 'solver_lstsqr_proj_';
            fileID = fopen([outputFolder, filename, num2str(i),'.txt'],'r');
            data = textscan(fileID,formatSpec);
            fclose(fileID);
            semilogy(data{:,2},'DisplayName','LSTSQR Proj')
        end
        
        subplot(2,1,1)
        title(['Convergence in Simulation #',num2str(solverdetails(i))])

        legend('show')
        ylim([10e-5 10e-3])
        xlim([0 2000])
        xlabel('# Iterations','interpreter','latex')
        ylabel('min($\mathcal{L}$)','interpreter','latex')
        
        subplot(2,1,2)
        legend('show')
        ylim([0 0.005])
        xlim([0 2000])
        xlabel('# Iterations','interpreter','latex')
        ylabel('min($\mathcal{L}$)','interpreter','latex')
        
        fig.PaperOrientation='landscape';
        print('-fillpage',[outputFolder, 'convergence_sim_',num2str(solverdetails(i))],'-dpdf')
    end
end

%%
if(closefig==1)
    close('all')
end
%autoArrangeFigures()

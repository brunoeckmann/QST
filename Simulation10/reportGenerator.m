close('all')
%% User Input
% Basis b
% 1 = Boyd
% 2 = Thew
b=1;

closefig=0;
simIndex=1; % Choose index for detailed report of this simulation


outputFolder='OutputSolverTest5/';

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
traceRhoList(:,1) = cellfun(@(A) trace(A), {rho_Noise{:,1}}); % rho Noise LinearInversion
traceRhoList(:,2) = cellfun(@trace, {rho_Noise{:,2}}); % rho Proj LinearInversion
traceRhoList(:,3) = cellfun(@trace, {rho_Noise{:,3}}); % rho Noise MaximumLikelihood
traceRhoList(:,4) = cellfun(@trace, {rho_Noise{:,4}}); % rho Proj MaximumLikelihood
traceRho2List(:,1) = cellfun(@(A) trace(A^2), {rho_Noise{:,1}}); % rho Noise LinearInversion
traceRho2List(:,2) = cellfun(@(A) trace(A^2), {rho_Noise{:,2}}); % rho Proj LinearInversion
traceRho2List(:,3) = cellfun(@(A) trace(A^2), {rho_Noise{:,3}}); % rho Noise MaximumLikelihood
traceRho2List(:,4) = cellfun(@(A) trace(A^2), {rho_Noise{:,4}}); % rho Proj MaximumLikelihood


%% Calculate Delta P
delta_P(:,1)=cellfun(@(A) A-P_Ideal,{P_Noise{:,1}}, 'UniformOutput', false);
delta_P(:,2)=cellfun(@(A) A-P_Ideal,{P_Noise{:,2}}, 'UniformOutput', false);
delta_P_Vector=cell2mat(cellfun(@(A) reshape(A-P_Ideal,36,1), {P_Noise{:,1}}, 'UniformOutput', false));
delta_P_Vector_Proj=cell2mat(cellfun(@(A) reshape(A-P_Ideal,36,1), {P_Noise{:,2}}, 'UniformOutput', false));

%% Calculate Fidelity
fidelityList(:,1)=(cellfun(@(A) Fidelity(A,rho_Ideal),{rho_Noise{:,1}}));
fidelityList(:,2)=(cellfun(@(A) Fidelity(A,rho_Ideal),{rho_Noise{:,2}}));
fidelityList(:,3)=(cellfun(@(A) Fidelity(A,rho_Ideal),{rho_Noise{:,3}}));
fidelityList(:,4)=(cellfun(@(A) Fidelity(A,rho_Ideal),{rho_Noise{:,4}}));

    

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

plotRho(rho_Noise{simIndex,3},[outputFolder,'rho_noiseState_ml.pdf'],{'\textbf{Reconstructed State $\rho_{Noise}$}','Maximum Likelihood'});
plotRho(rho_Noise{simIndex,4},[outputFolder,'rho_projState_ml.pdf'],{'\textbf{Reconstructed State after Projection $\rho_{Proj}$}','Maximum Likelihood'});

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
saveas(fig,[outputFolder,'noise.pdf']);

mean_delta_P{1}=zeros(size(delta_P{1,1}));
for ii=1:nMeasurement
    mean_delta_P{1} = mean_delta_P{1} + delta_P{ii,1};
end
mean_delta_P{1}=mean_delta_P{1} / nMeasurement;

%% Plot Fidelity and P element of NS
fig=figure('units','normalized','outerposition',[0 0 0.8 0.8]);
plot(real((fidelityList(:,1))),'b')
hold on
plot(real((fidelityList(:,3))),'b--')
plot(real((fidelityList(:,2))),'r')
plot(real((fidelityList(:,4))),'r--')
xlabel('# simulation','interpreter','latex')
ylabel('Fidelity $\mathcal{F}$','interpreter','latex')
xlim([1,size(fidelityList,1)])

yyaxis right
ylabel('$P(ab|xy) \in \mathcal{NS}$','interpreter','latex')
set(gca,'ColorOrder',[0 0 1;1 0 1])
set(gca,'FontSize',12)

plot(cell2mat(nonsignalingList),'d')
ylim([-0.25 1.25])
h=legend('$\mathcal{F}(\rho_{Noise},\rho_{Ideal})$ Linear Inversion','$\mathcal{F}(\rho_{Noise},\rho_{Ideal})$ Maximum Likelihood',...
    '$\mathcal{F}(\rho_{Proj},\rho_{Ideal})$ Linear Inversion','$\mathcal{F}(\rho_{Proj},\rho_{Ideal})$ Maximum Likelihood',...
    '$P_{Noise} \in \mathcal{NS}$', '$P_{Proj} \in \mathcal{NS}$');
legend('Location','bestoutside')
set(h,'Interpreter','latex')
plt = gca;
plt.YAxis(2).Color='k';
fig.PaperOrientation='landscape';
saveas(fig,[outputFolder,'fidelity.pdf']);

%% Plot Fidelity Histogram

edges=linspace(0.9,1.05,20);
h1=histogram(real(fidelityList(:,1)),edges)
set(gca,'YLim',[0,length(fidelityList)*1.1])
hold on
histogram(real(fidelityList(:,3)),edges)
histogram(real(fidelityList(:,2)),edges)

histogram(real(fidelityList(:,4)),edges)

h=legend('$\mathcal{F}(\rho_{Noise},\rho_{Ideal})$ Linear Inversion','$\mathcal{F}(\rho_{Noise},\rho_{Ideal})$ Maximum Likelihood',...
    '$\mathcal{F}(\rho_{Proj},\rho_{Ideal})$ Linear Inversion','$\mathcal{F}(\rho_{Proj},\rho_{Ideal})$ Maximum Likelihood',...
    '$P_{Noise} \in \mathcal{NS}$', '$P_{Proj} \in \mathcal{NS}$');
legend('Location','bestoutside')
set(h,'Interpreter','latex')
ax = gca;
ax.ColorOrderIndex = 1;
plot([mean(fidelityList(:,1)), mean(fidelityList(:,1))],get(gca,'YLim'),'LineStyle','--')
plot([mean(fidelityList(:,2)), mean(fidelityList(:,2))],get(gca,'YLim'),'LineStyle','--')
plot([mean(fidelityList(:,3)), mean(fidelityList(:,3))],get(gca,'YLim'),'LineStyle','--')
plot([mean(fidelityList(:,4)), mean(fidelityList(:,4))],get(gca,'YLim'),'LineStyle','--')


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
saveas(fig,[outputFolder,'fidelity_sim', num2str(simIndex),'.pdf']);

%% Plot Trace
fig=figure('units','normalized','outerposition',[0 0 0.8 0.8]);
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
hh=legend('Tr($\rho_{Noise}^2$) Linear Inversion','Tr($\rho_{Proj}^2$) Linear Inversion',...
    'Tr($\rho_{Noise}^2$) Maximum Likelihood','Tr($\rho_{Proj}^2$) Maximum Likelihood');
%legend('Location','bestoutside')
set(hh,'Interpreter','latex')
fig.PaperOrientation='landscape';
saveas(fig,[outputFolder,'trace.pdf']);

%%
%autoArrangeFigures()
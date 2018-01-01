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
lsqnonlin=1;

solverdetails=[1,100]; % Choose simulations for which detailed solver analysis is done


outputFolder = 'OutputSolverTest100_Noise100_WNM08/';

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

%% Extended Tomographic Set: Basis for NS regularization
chi{5}=1/sqrt(2)*([1;0]-[0;1]);
chi{6}=1/sqrt(2)*([-1;1i]);

ChiReg=cell(36,1);
for ii=1:6
   for j=1:6
      ChiReg{6*(ii-1)+j}=kron(chi{ii},chi{j}); 
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
plotRho(rho_WNM,[outputFolder,'rho_targetStateWNM.pdf'],['\textbf{Target State $\rho_{WNM}$ White Noise $\lambda=$',num2str(lambda),'}']);


%% Check Signaling of Noisy and Regularized Behaviours
nonsignalingList(:,1) = num2cell(cellfun(@checkNonSignaling,{P_Noise{:,1}})); % Noisy Behaviour
nonsignalingList(:,2) = num2cell(cellfun(@checkNonSignaling,{P_Noise{:,2}})); % regected Behaviour

%% Calculate Traces
for i=1:8
    traceRhoList(:,i) = cellfun(@(A) trace(A), {rho_Noise{:,i}});
    traceRho2List(:,i) = cellfun(@(A) trace(A^2), {rho_Noise{:,i}}); 
end

%% Calculate Delta P
delta_P(:,1)=cellfun(@(A) A-P_Ideal,{P_Noise{:,1}}, 'UniformOutput', false);
delta_P(:,2)=cellfun(@(A) A-P_Ideal,{P_Noise{:,2}}, 'UniformOutput', false);
delta_P_Vector=cell2mat(cellfun(@(A) reshape(A-P_Ideal,36,1), {P_Noise{:,1}}, 'UniformOutput', false));
delta_P_Vector_reg=cell2mat(cellfun(@(A) reshape(A-P_Ideal,36,1), {P_Noise{:,2}}, 'UniformOutput', false));

%% Calculate Fidelity
ll = [1,2];
if fmin==1
    ll = [ll,3,4];
end
if ga==1
    ll = [ll,5,6];
end
if lsqnonlin==1
    ll=[ll,7,8];
end

for i = ll
    fidelityList(:,i)=(cellfun(@(A) Fidelity(A,rho_Ideal),{rho_Noise{:,i}})); % Fidelity with respect to 
    %superfidelityList(:,i)=(cellfun(@(A) SuperFidelity(A,rho_Ideal),{rho_Noise{:,i}})); % Fidelity with respect to 

    %fidelityList(:,i)=(cellfun(@(A)
    %Fidelity(A,rho_WNM),{rho_Noise{:,i}})); %Fidelity with respect to
    %rho_WNM
end
    

%% Calculate Pideal + (1-x)dP
deltaP = P_Noise{simIndex,1}-P_Noise{simIndex,2}; % P_Noise - P_reg_LI
x=linspace(0,1,11);


for i = 1:11
    P_List{i} =  P_Noise{simIndex,2} + (1-x(i))*deltaP;
    rhoList{i}= qst_linearinversion(Chi,P_QST_Selection(P_List{i}));
    fid(i) = Fidelity(rhoList{i},rho_Ideal);
    
end

%%

%% Plot Reconstruced States
plotRho(rho_Noise{simIndex,1},[outputFolder,'rho_noiseState_li.pdf'],{'\textbf{Reconstructed State $\rho_{Noise}$}','Linear Inversion'});
plotRho(rho_Noise{simIndex,2},[outputFolder,'rho_regState_li.pdf'],{'\textbf{Reconstructed State after Regularization $\rho_{reg}$}','Linear Inversion'});

plotRho(rho_Noise{simIndex,3},[outputFolder,'rho_noiseState_ml_fmin.pdf'],{'\textbf{Reconstructed State $\rho_{Noise}$}','Maximum Likelihood'});
plotRho(rho_Noise{simIndex,4},[outputFolder,'rho_regState_ml_fmin.pdf'],{'\textbf{Reconstructed State after Regularization $\rho_{reg}$}','Maximum Likelihood'});

plotRho(rho_Noise{simIndex,5},[outputFolder,'rho_noiseState_ml_ga.pdf'],{'\textbf{Reconstructed State $\rho_{Noise}$}','Maximum Likelihood'});
plotRho(rho_Noise{simIndex,6},[outputFolder,'rho_regState_ml_ga.pdf'],{'\textbf{Reconstructed State after Regularization $\rho_{reg}$}','Maximum Likelihood'});

plotRho(rho_Noise{simIndex,7},[outputFolder,'rho_noiseState_ml_lsqnonlin.pdf'],{'\textbf{Reconstructed State $\rho_{Noise}$}','Maximum Likelihood'});
plotRho(rho_Noise{simIndex,8},[outputFolder,'rho_regState_ml_lsqnonlin.pdf'],{'\textbf{Reconstructed State after Regularization $\rho_{reg}$}','Maximum Likelihood'});


plotRho(rho_Noise{simIndex,1}-rho_Ideal,[outputFolder,'rho_deltaNoiseState_li.pdf'],{'\textbf{Difference to Target State before Regularization $\rho_{Noise}-\rho_{Ideal}$}','Linear Inversion'});
plotRho(rho_Noise{simIndex,2}-rho_Ideal,[outputFolder,'rho_deltaregState_li.pdf'],{'\textbf{Difference to Target State after Regularization $\rho_{reg}-\rho_{Ideal}$}','Maximum Likelihood'});
plotRho(rho_Noise{simIndex,3}-rho_Ideal,[outputFolder,'rho_deltaNoiseState_ml.pdf'],{'\textbf{Difference to Target State before Regularization $\rho_{Noise}-\rho_{Ideal}$}','Maximum Likelihood'});
plotRho(rho_Noise{simIndex,4}-rho_Ideal,[outputFolder,'rho_deltaregState_ml.pdf'],{'\textbf{Difference to Target State after Regularization $\rho_{reg}-\rho_{Ideal}$}','Maximum Likelihood'});

%% Plot Behaviour
plotBehaviour(P_Ideal,P_Noise, simIndex,[outputFolder,'behaviour.pdf'],['\textbf{Behaviour $P(ab| xy)$ of Simulation ',num2str(simIndex),'}']);

P = cell2mat(cellfun(@(A) {reshape(A,36,1)}, {P_Noise{:,1}})) % Behaviour with Noise
P_Noise_mean{1,1}=mean(P,2)
P = cell2mat(cellfun(@(A) {reshape(A,36,1)}, {P_Noise{:,2}})) % Behaviour with Noise after Regularization
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
    pd_reg = fitdist(abs(delta_P_Vector_reg(jj,:))','poisson');
    var_fit(jj)=pd.lambda*countFactor;
    var_fit_reg(jj)=pd_reg.lambda;

    var_P(jj)=var(delta_P_Vector(jj,:)*countFactor);
    var_P_reg(jj)=var(delta_P_Vector_reg(jj,:)*countFactor);
end

fig=figure;
plot(var_P)
hold on
plot(var_P_reg)
plot(var_fit,'--')
plot(var_fit_reg,'--')
plot(reshape(P_Ideal,36,1),'LineWidth',2)
legend('Poisson Fit', 'Poisson Fit after NS Regularization','Variance Func','Variance Func after NS Regularization','P Ideal')
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

plot(real((fidelityList(:,1))),'b','DisplayName','$\mathcal{F}(\rho_{Noise},\rho_{Ideal})$ Linear Inversion')
hold on


if fmin==1
    plot(real((fidelityList(:,3))),'b--','DisplayName','$\mathcal{F}(\rho_{Noise},\rho_{Ideal})$ Maximum Likelihood fminsearch')
end
if ga==1
    plot(real((fidelityList(:,5))),'b:','DisplayName','$\mathcal{F}(\rho_{Noise},\rho_{Ideal})$ Maximum Likelihood ga')
end
if lsqnonlin==1
    plot(real((fidelityList(:,7))),'b-.','DisplayName','$\mathcal{F}(\rho_{Noise},\rho_{Ideal})$ Maximum Likelihood lsqnonlin')
end
xlabel('# simulation','interpreter','latex')
ylabel('Fidelity $\mathcal{F}$','interpreter','latex')

yyaxis right
ylabel('$P(ab|xy) \in \mathcal{NS}$','interpreter','latex')
set(gca,'ColorOrder',[0 0 1;1 0 1])
set(gca,'FontSize',12)
plot(cell2mat(nonsignalingList(:,1)),'d','DisplayName','$P_{Noise} \in \mathcal{NS}$');
ylim([-0.25 1.25])
yticks([0 1])
plt = gca;
plt.YAxis(2).Color='k';

h=legend('show')
legend('Location','bestoutside')
set(h,'Interpreter','latex')

subplot(2,1,2)
plot(real((fidelityList(:,2))),'r','DisplayName','$\mathcal{F}(\rho_{reg},\rho_{Ideal})$ Linear Inversion')
hold on
if fmin==1
    plot(real((fidelityList(:,4))),'r--','DisplayName','$\mathcal{F}(\rho_{reg},\rho_{Ideal})$ Maximum Likelihood fminsearch')
end
if ga==1
    plot(real((fidelityList(:,6))),'r:','DisplayName','$\mathcal{F}(\rho_{reg},\rho_{Ideal})$ Maximum Likelihood ga')
end
if lsqnonlin==1
    plot(real((fidelityList(:,8))),'r-.','DisplayName','$\mathcal{F}(\rho_{reg},\rho_{Ideal})$ Maximum Likelihood lsqnonlin')
end
xlabel('# simulation','interpreter','latex')
ylabel('Fidelity $\mathcal{F}$','interpreter','latex')
xlim([1,size(fidelityList,1)])
yyaxis right
ylabel('$P(ab|xy) \in \mathcal{NS}$','interpreter','latex')
set(gca,'ColorOrder',[0 0 1;1 0 1])
set(gca,'FontSize',12)
plot(cell2mat(nonsignalingList(:,2)),'d','DisplayName','$P_{reg} \in \mathcal{NS}$');ylim([-0.25 1.25])
yticks([0 1])
plt = gca;
plt.YAxis(2).Color='k';
h=legend('show')
legend('Location','bestoutside')
set(h,'Interpreter','latex')



fig.PaperOrientation='landscape';
%saveas(fig,[outputFolder,'fidelity.pdf']);
print('-fillpage',[outputFolder, 'fidelity'],'-dpdf')

%close(fig)

%% Plot Fidelity Histogram
str_list={'fminsearch','genetic algorithm','lsqnonlin'};

for i=1:3
    fig=figure()
    edges=linspace(0,2,300);
    subplot(2,1,1)
    histogram(real(fidelityList(:,1)),edges) % Noisy with LI
    %set(gca,'YLim',[0,length(fidelityList)*1.1]) 
    hold on
    histogram(real(fidelityList(:,2)),edges) % reg with LI
    ax = gca;
    ax.ColorOrderIndex = 1;
    ylabel('# States')
     h=legend('$\mathcal{F}(\rho_{Noise},\rho_{Ideal})$ Linear Inversion','$\mathcal{F}(\rho_{reg},\rho_{Ideal})$ Linear Inversion');
    legend('Location','northwest')
    set(h,'Interpreter','latex')
    plot([mean(fidelityList(:,1)), mean(fidelityList(:,1))],get(gca,'YLim'),'LineStyle','--','LineWidth',2)
    plot([mean(fidelityList(:,2)), mean(fidelityList(:,2))],get(gca,'YLim'),'LineStyle','--','LineWidth',2)
    plot([0.8 0.8],get(gca,'YLim'),'k--','LineWidth',2)

    xlabel('Fidelity w.r.t. $\rho_{target}$','interpreter','latex')
    xlim([0.6 1])
    subplot(2,1,2)
    %set(gca,'YLim',[0,length(fidelityList)*1.1]) 
    if fmin==1 && i==1
        h1=histogram(real(fidelityList(:,3)),edges); % Noisy with ML
        hold on
        h2=histogram(real(fidelityList(:,4)),edges); % reg with ML
        ax = gca;
        ax.ColorOrderIndex = 1;
        plot([mean(fidelityList(:,3)), mean(fidelityList(:,3))],get(gca,'YLim'),'LineStyle','--','LineWidth',2)
        plot([mean(fidelityList(:,4)), mean(fidelityList(:,4))],get(gca,'YLim'),'LineStyle','--','LineWidth',2)
        h=legend([h1 h2],{['$\mathcal{F}(\rho_{Noise},\rho_{Ideal})$ ML ',str_list{i}],['$\mathcal{F}(\rho_{reg},\rho_{Ideal})$ ML ',str_list{i}]});

    end
    
    if ga==1 && i==2
        h1=histogram(real(fidelityList(:,5)),edges); % Noisy with ga
        hold on
        h2=histogram(real(fidelityList(:,6)),edges); % reg with ga
            ax = gca;
        ax.ColorOrderIndex = 1;
        plot([mean(fidelityList(:,5)), mean(fidelityList(:,5))],get(gca,'YLim'),'LineStyle','--','LineWidth',2)
        plot([mean(fidelityList(:,6)), mean(fidelityList(:,6))],get(gca,'YLim'),'LineStyle','--','LineWidth',2)
        h=legend([h1 h2],{['$\mathcal{F}(\rho_{Noise},\rho_{Ideal})$ ML ',str_list{i}],['$\mathcal{F}(\rho_{reg},\rho_{Ideal})$ ML ',str_list{i}]});

    end
    
    if lsqnonlin==1 && i==3
        h1=histogram(real(fidelityList(:,7)),edges); % Noisy with nonlinlsq
        hold on
        h2=histogram(real(fidelityList(:,8)),edges); % reg with nonlinlsq
        ax = gca;
        ax.ColorOrderIndex = 1;    
        plot([mean(fidelityList(:,7)), mean(fidelityList(:,7))],get(gca,'YLim'),'LineStyle','--','LineWidth',2)
        plot([mean(fidelityList(:,8)), mean(fidelityList(:,8))],get(gca,'YLim'),'LineStyle','--','LineWidth',2)
        h=legend([h1 h2],{['$\mathcal{F}(\rho_{Noise},\rho_{Ideal})$ ML ',str_list{i}],['$\mathcal{F}(\rho_{reg},\rho_{Ideal})$ ML ',str_list{i}]});

    end
   
    plot([0.8 0.8],get(gca,'YLim'),'k--','LineWidth',2)
    xlim([0.6 1])

    
    set(h,'Interpreter','latex')
    set(h,'Location','northwest')
    fig.PaperOrientation='landscape';
    ylabel('# States')
    xlabel('Fidelity w.r.t. $\rho_{target}$','interpreter','latex')
    print('-fillpage',[outputFolder, 'fidelityHistogramm_',str_list{i}],'-dpdf')
end
    
%% Plot Fidelity of Simulation 1
fig=figure('units','normalized','outerposition',[0 0 0.8 0.8]);

plot(x,real((fid)),'o--')
grid on
title(['Fidelity $\mathcal{F}(\rho_x,\rho_{ideal} )$ with $P_x = x P_{reg} + (1-x) P_{noise}$'],'interpreter','latex')
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
% h=legend('Tr($\rho_{Noise}$) Linear Inversion','Tr($\rho_{reg}$) Linear Inversion',...
%     'Tr($\rho_{Noise}$) Maximum Likelihood','Tr($\rho_{reg}$) Maximum Likelihood',...
%     'Tr($\rho_{Noise}^2$)','Tr($\rho_{reg}^2$)');
h=legend('Tr($\rho_{Noise}^2$) Linear Inversion','Tr($\rho_{reg}^2$) Linear Inversion',...
    'Tr($\rho_{Noise}^2$) Maximum Likelihood fmin','Tr($\rho_{reg}^2$) Maximum Likelihood fmin',...
    'Tr($\rho_{Noise}^2$) Maximum Likelihood ga','Tr($\rho_{reg}^2$) Maximum Likelihood ga',...
    'Tr($\rho_{Noise}^2$) Maximum Likelihood lsqnonlin','Tr($\rho_{reg}^2$) Maximum Likelihood lsqnonlin');
%legend('Location','bestoutside')
set(h,'Interpreter','latex')
fig.PaperOrientation='landscape';
fig.PaperPositionMode = 'auto';
%saveas(fig,[outputFolder,'trace.pdf']);
print('-fillpage',[outputFolder, 'trace'],'-dpdf')
%close(fig)

%% Histrogram Trace
str_list={'fminsearch','genetic algorithm','lsqnonlin'};
x0=0.2;
x1=1.05;
for i=1:3
    fig=figure()
    edges=linspace(0,2,300);
    subplot(2,1,1)
    histogram(real(traceRho2List(:,1)),edges) % Noisy with LI
    %set(gca,'YLim',[0,length(fidelityList)*1.1]) 
    hold on
    histogram(real(traceRho2List(:,2)),edges) % reg with LI
    ax = gca;
    ax.ColorOrderIndex = 1;
    ylabel('# States')
    h=legend('$\rho^{noise}$ Linear Inversion','$\rho^{reg}$ Linear Inversion');
    legend('Location','northwest')
    set(h,'Interpreter','latex')
    plot([mean(traceRho2List(:,1)), mean(traceRho2List(:,1))],get(gca,'YLim'),'LineStyle','--','LineWidth',2)
    plot([mean(traceRho2List(:,2)), mean(traceRho2List(:,2))],get(gca,'YLim'),'LineStyle','--','LineWidth',2)
    plot([1 1],get(gca,'YLim'),'k--','LineWidth',2)

    xlabel('Tr[$\rho^2$]','interpreter','latex')
    xlim([x0 x1])
    subplot(2,1,2)
    %set(gca,'YLim',[0,length(fidelityList)*1.1]) 
    if fmin==1 && i==1
        h1=histogram(real(traceRho2List(:,3)),edges); % Noisy with ML
        hold on
        h2=histogram(real(traceRho2List(:,4)),edges); % reg with ML
        ax = gca;
        ax.ColorOrderIndex = 1;
        plot([mean(traceRho2List(:,3)), mean(traceRho2List(:,3))],get(gca,'YLim'),'LineStyle','--','LineWidth',2)
        plot([mean(traceRho2List(:,4)), mean(traceRho2List(:,4))],get(gca,'YLim'),'LineStyle','--','LineWidth',2)
        h=legend([h1 h2],{['$\rho^{noise}$ ML ',str_list{i}],['$\rho^{reg}$ ML ',str_list{i}]});

    end
    
    if ga==1 && i==2
        h1=histogram(real(traceRho2List(:,5)),edges); % Noisy with ga
        hold on
        h2=histogram(real(traceRho2List(:,6)),edges); % reg with ga
            ax = gca;
        ax.ColorOrderIndex = 1;
        plot([mean(traceRho2List(:,5)), mean(traceRho2List(:,5))],get(gca,'YLim'),'LineStyle','--','LineWidth',2)
        plot([mean(traceRho2List(:,6)), mean(traceRho2List(:,6))],get(gca,'YLim'),'LineStyle','--','LineWidth',2)
        h=legend([h1 h2],{['$\rho^{noise}$ ML ',str_list{i}],['$\rho^{reg}$ ML ',str_list{i}]});

    end
    
    if lsqnonlin==1 && i==3
        h1=histogram(real(traceRho2List(:,7)),edges); % Noisy with nonlinlsq
        hold on
        h2=histogram(real(traceRho2List(:,8)),edges); % reg with nonlinlsq
        ax = gca;
        ax.ColorOrderIndex = 1;    
        plot([mean(traceRho2List(:,7)), mean(traceRho2List(:,7))],get(gca,'YLim'),'LineStyle','--','LineWidth',2)
        plot([mean(traceRho2List(:,8)), mean(traceRho2List(:,8))],get(gca,'YLim'),'LineStyle','--','LineWidth',2)
        h=legend([h1 h2],{['$\rho^{noise}$ ML ',str_list{i}],['$\rho^{reg}$ ML ',str_list{i}]});

    end
   
    plot([1 1],get(gca,'YLim'),'k--','LineWidth',2)
    xlim([x0 x1])

    
    set(h,'Interpreter','latex')
    set(h,'Location','northwest')
    fig.PaperOrientation='landscape';
    ylabel('# States')
    xlabel('Tr[$\rho^2$]','interpreter','latex')
    print('-fillpage',[outputFolder, 'traceHistogramm_',str_list{i}],'-dpdf')
end

%% Histrogram Eigenvalues
str_list={'fminsearch','genetic algorithm','lsqnonlin'};

for i=1:3
    eigs_Noise=[];
    eigs_NS=[];
    for m=1:nMeasurement
        eigs_Noise = [eigs_Noise; eig(rho_Noise{m,1})];
        eigs_NS = [eigs_NS; eig(rho_Noise{m,2})];
    end
    
    x0=min(eigs_Noise);  
    x1=max(eigs_Noise)+0.1;
    
    fig=figure()
    edges=linspace(-1,1,200);
    subplot(2,1,1)
    histogram(eigs_Noise,edges) % Noisy with LI
    %set(gca,'YLim',[0,length(fidelityList)*1.1]) 
    hold on
    histogram(eigs_NS,edges) % reg with LI
    ax = gca;
    ax.ColorOrderIndex = 1;
    ylabel('# Eigenvalues')
    h=legend('$\rho^{noise}$ Linear Inversion','$\rho^{reg}$ Linear Inversion');
    legend('Location','northwest')
    set(h,'Interpreter','latex')
    plot([0 0],get(gca,'YLim'),'k--','LineWidth',2)
    plot([1 1],get(gca,'YLim'),'k--','LineWidth',2)
    yll=get(gca,'YLim');

    xlabel('Eigenvalues','interpreter','latex')
    xlim([x0 x1])
    subplot(2,1,2)
    %set(gca,'YLim',[0,length(fidelityList)*1.1]) 
    if fmin==1 && i==1
        eigs_Noise=[];
        eigs_NS=[];
        for m=1:nMeasurement
            eigs_Noise = [eigs_Noise; eig(rho_Noise{m,3})];
            eigs_NS = [eigs_NS; eig(rho_Noise{m,4})];
        end
        h1=histogram(eigs_Noise,edges); % Noisy with ML
        hold on
        h2=histogram(eigs_NS,edges); % reg with ML
        ax = gca;
        ax.ColorOrderIndex = 1;
        h=legend([h1 h2],{['$\rho^{noise}$ ML ',str_list{i}],['$\rho^{reg}$ ML ',str_list{i}]});

    end
    
    if ga==1 && i==2
        eigs_Noise=[];
        eigs_NS=[];
        for m=1:nMeasurement
            eigs_Noise = [eigs_Noise; eig(rho_Noise{m,5})];
            eigs_NS = [eigs_NS; eig(rho_Noise{m,6})];
        end
        h1=histogram(eigs_Noise,edges); % Noisy with ML
        hold on
        h2=histogram(eigs_NS,edges); % reg with ML
        ax = gca;
        ax.ColorOrderIndex = 1;
        h=legend([h1 h2],{['$\rho^{noise}$ ML ',str_list{i}],['$\rho^{reg}$ ML ',str_list{i}]});

    end
    
    if lsqnonlin==1 && i==3        
        eigs_Noise=[];
        eigs_NS=[];
        for m=1:nMeasurement
            eigs_Noise = [eigs_Noise; eig(rho_Noise{m,7})];
            eigs_NS = [eigs_NS; eig(rho_Noise{m,8})];
        end
        h1=histogram(eigs_Noise,edges); % Noisy with ML
        hold on
        h2=histogram(eigs_NS,edges); % reg with ML
        ax = gca;
        ax.ColorOrderIndex = 1;    
        h=legend([h1 h2],{['$\rho^{noise}$ ML ',str_list{i}],['$\rho^{reg}$ ML ',str_list{i}]});

    end
    
    plot([0 0],get(gca,'YLim'),'k--','LineWidth',2)
    plot([1 1],get(gca,'YLim'),'k--','LineWidth',2)
    xlim([x0 x1])
    ylim(yll)

    
    set(h,'Interpreter','latex')
    set(h,'Location','northwest')
    fig.PaperOrientation='landscape';
    ylabel('# Eigenvalues')
    xlabel('Eigenvalues','interpreter','latex')
    print('-fillpage',[outputFolder, 'eigenvaluesHistogramm_',str_list{i}],'-dpdf')
end

%% Plot Solver Time
fig=figure();
avgrTime = mean(tElapsed(:,3:8),1);
plot(tElapsed(:,3:8))
title('Evaluation time for different Maximum Likelihood Solvers','interpreter','latex')
xlabel('# simulation','interpreter','latex')
ylabel('Time [sec]','interpreter','latex')
h=legend(['fminsearch on $\rho_{Noise}$, Average: ',num2str(avgrTime(1)),' sec'],...
    ['fminsearch on $\rho_{reg}$, Average: ',num2str(avgrTime(2)),' sec'],...
    ['genetic algorithm on $\rho_{Noise}$, Average: ',num2str(avgrTime(3)),' sec'],...
    ['genetic algorithm on $\rho_{reg}$, Average: ',num2str(avgrTime(4)),' sec'],...
    ['lsqnonlin on $\rho_{Noise}$, Average: ',num2str(avgrTime(5)),' sec'],...
    ['lsqnonlin on $\rho_{reg}$, Average: ',num2str(avgrTime(6)),' sec']);
set(h,'Interpreter','latex')
fig.PaperOrientation='landscape';
print('-fillpage',[outputFolder, 'timesolverml'],'-dpdf')

%% Histrogram Plot Solver Time
fig=figure();
subplot(2,1,1)
if fmin==1
h1=histogram(tElapsed(:,3),'DisplayName','fminsearch $\rho_{noise}$')
h1.BinWidth=20;
hold on
end
if ga==1
h2=histogram(tElapsed(:,5),'DisplayName','genetic algorithm $\rho_{noise}$')
h2.BinWidth=20;
hold on
end
if lsqnonlin==1
h3=histogram(tElapsed(:,7),'DisplayName','lsqnonlin $\rho_{noise}$')
h3.BinWidth=20;
end
ylabel('run','interpreter','latex')
xlabel('evaluation time [sec]','interpreter','latex')
l=legend('show')
set(l,'interpreter','latex')
xl=get(gca,'xlim')
yl=get(gca,'ylim')
subplot(2,1,2)
if fmin==1
h1=histogram(tElapsed(:,4),'DisplayName','fminsearch $\rho_{reg}$')
h1.BinWidth=10;
hold on
end
if ga==1
h2=histogram(tElapsed(:,6),'DisplayName','genetic algorithm $\rho_{reg}$')
h2.BinWidth=10;
hold on
end
if lsqnonlin==1
h3=histogram(tElapsed(:,8),'DisplayName','lsqnonlin $\rho_{reg}$')
h3.BinWidth=10;
end
l=legend('show')
set(l,'interpreter','latex')
%set(gca,'xlim',xl)
%set(gca,'ylim',yl)
ylabel('run','interpreter','latex')
xlabel('evaluation time [sec]','interpreter','latex')
fig.PaperOrientation='landscape';
print('-fillpage',[outputFolder, 'time_histogram'],'-dpdf')

%% Generate Convergence Plots
if (fmin==1 || ga==1 || lsqnonlin==1)
    formatSpec = '%d\t %f\t %d\t %s';

    for i=1:length(solverdetails)
        fig=figure();
        if fmin==1
            subplot(2,1,1)
            filename = 'solver_fmin_noise_';
            fileID = fopen([outputFolder, filename, num2str(solverdetails(i)),'.txt'],'r');
            data = textscan(fileID,formatSpec);
            fclose(fileID);
            semilogy(data{3},data{:,2},'DisplayName','fminsearch')
            hold on
            subplot(2,1,2)
            filename = 'solver_fmin_proj_';
            fileID = fopen([outputFolder, filename, num2str(solverdetails(i)),'.txt'],'r');
            data = textscan(fileID,formatSpec);
            fclose(fileID);
            semilogy(data{3},data{:,2},'DisplayName','fminsearch')
            hold on
        end
        if ga==1
            subplot(2,1,1)
            filename = 'solver_ga_noise_';
            fileID = fopen([outputFolder, filename, num2str(solverdetails(i)),'.txt'],'r');
            data = textscan(fileID,formatSpec,'HeaderLines',1);
            fclose(fileID);
            semilogy(data{3},data{:,2},'DisplayName','genetic algorithm')
            hold on
            subplot(2,1,2)
            filename = 'solver_ga_proj_';
            fileID = fopen([outputFolder, filename, num2str(solverdetails(i)),'.txt'],'r');
            data = textscan(fileID,formatSpec,'HeaderLines',1);
            fclose(fileID);
            semilogy(data{3},data{2},'DisplayName','genetic algorithm')
            hold on
        end
        if lsqnonlin==1
            subplot(2,1,1)
            filename = 'solver_lsqnonlin_noise_';
            fileID = fopen([outputFolder, filename, num2str(solverdetails(i)),'.txt'],'r');
            data = textscan(fileID,formatSpec);
            fclose(fileID);
            semilogy(data{3},data{:,2},'DisplayName','lsqnonlin')
            hold on
            subplot(2,1,2)
            filename = 'solver_lsqnonlin_proj_';
            fileID = fopen([outputFolder, filename, num2str(solverdetails(i)),'.txt'],'r');
            data = textscan(fileID,formatSpec);
            fclose(fileID);
            semilogy(data{3},data{:,2},'DisplayName','lsqnonlin')
        end
        
        subplot(2,1,1)
        title(['Convergence in Simulation #',num2str(solverdetails(i)),' QST Maximum Likelihood \rho_{Noise}'])

        legend('show')
        %ylim([10e-5 10e-3])
        %xlim([0 2000])
        xlabel('# function evaluations','interpreter','latex')
        ylabel('log[min($\mathcal{L}$)]','interpreter','latex')
        
        subplot(2,1,2)
        title(['Convergence in Simulation #',num2str(solverdetails(i)),' QST Maximum Likelihood \rho_{reg}'])
        legend('show')
        %lim([0 0.005])
        %xlim([0 2000])
        xlabel('# function evaluations','interpreter','latex')
        ylabel('log[min($\mathcal{L}$)]','interpreter','latex')
        
        fig.PaperOrientation='landscape';
        print('-fillpage',[outputFolder, 'convergence_sim_',num2str(solverdetails(i))],'-dpdf')
    end
end

%% Histogram Convergence Residuals
% Read residuals
formatSpec = '%d\t %f\t %d\t %s';
res=[];
for i=1:nMeasurement
    if fmin==1
        filename = 'solver_fmin_noise_';
        fileID = fopen([outputFolder, filename, num2str(i),'.txt'],'r');
        data = textscan(fileID,formatSpec);
        fclose(fileID);
        res(i,1)=data{2}(end);
        filename = 'solver_fmin_proj_';
        fileID = fopen([outputFolder, filename, num2str(i),'.txt'],'r');
        data = textscan(fileID,formatSpec);
        fclose(fileID);
        res(i,2)=(data{2}(end));
    end
    
    if ga==1
        filename = 'solver_ga_noise_';
        fileID = fopen([outputFolder, filename, num2str(i),'.txt'],'r');
        data = textscan(fileID,formatSpec,'HeaderLines',1);
        fclose(fileID);
        res(i,3)=(data{2}(end));
        filename = 'solver_ga_proj_';
        fileID = fopen([outputFolder, filename, num2str(i),'.txt'],'r');
        data = textscan(fileID,formatSpec,'HeaderLines',1);
        fclose(fileID);
        res(i,4)=(data{2}(end));
    end
    
    if lsqnonlin==1
        filename = 'solver_lsqnonlin_noise_';
        fileID = fopen([outputFolder, filename, num2str(i),'.txt'],'r');
        data = textscan(fileID,formatSpec);
        fclose(fileID);
        res(i,5)=(data{2}(end));
        filename = 'solver_lsqnonlin_proj_';
        fileID = fopen([outputFolder, filename, num2str(i),'.txt'],'r');
        data = textscan(fileID,formatSpec);
        fclose(fileID);
        res(i,6)=(data{2}(end));
    end
end

fig=figure();
fig.Position = [10   10   700   500];
subplot(2,1,1)
if fmin==1
h1=histogram(log10(res(:,1)),'DisplayName','fminsearch')
h1.BinWidth=0.4;
hold on
end
if ga==1
h2=histogram(log10(res(:,3)),'DisplayName','genetic algorithm')
h2.BinWidth=0.4;
hold on
end
if lsqnonlin==1
h3=histogram(log10(res(:,5)),'DisplayName','lsqnonlin')
h3.BinWidth=0.4;
end
title('before regularization')
ylabel('\# samples','interpreter','latex')
l=legend('show')
set(l,'interpreter','latex','box','off','location','north')
xl=get(gca,'xlim')
yl=get(gca,'ylim')
box
grid

subplot(2,1,2)
if fmin==1
h1=histogram(log10(res(:,2)),'DisplayName','fminsearch')
h1.BinWidth=2;
hold on
end
if ga==1
h2=histogram(log10(res(:,4)),'DisplayName','genetic algorithm')
h2.BinWidth=2;
hold on
end
if lsqnonlin==1
h3=histogram(log10(res(:,6)),'DisplayName','lsqnonlin')
h3.BinWidth=2;
end
l=legend('show');

set(l,'interpreter','latex','box','off','location','north')
%set(gca,'xlim',xl)
%set(gca,'ylim',yl)
title('after regularization')
ylabel('\# samples','interpreter','latex')
xlabel('log[min($L$)]','interpreter','latex')
fig.PaperOrientation='landscape';
print('-fillpage',[outputFolder, 'convergence_residual_histogram'],'-dpdf')
set(gcf, 'Renderer', 'opengl')
set(gca,'FontName','Helvetica');
box
grid
tightfig;
%print([outputFolder, 'convergence_residual_histogram.eps'],'-depsc2')
saveas(gcf,[outputFolder, 'convergence_residual_histogram'],'epsc')


%% 

%%
if(closefig==1)
    close('all')
end
%autoArrangeFigures()

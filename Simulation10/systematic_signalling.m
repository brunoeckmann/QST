close('all')
clear
fclose('all');
rng default; %reset random number generator
poolobj = gcp('nocreate');
delete(poolobj);
%% How does a systematic signalling source in the setup affect the distance to the NS polytope

% Choose Basis b
% 1 = Boyd
% 2 = Thew
b=1;

% Choose NS Regularization Method
ns_method = 'liang';
%ns_method = 'alberto';

% Choose Number of Simulations
nMeasurement=1000;
nAlpha=20;

% Choose White Noise Model lambda factor
lambda = 0.52;

% Choose Noise Factor
countFactor1=100; % >0, if countFactor is low, Noise is higher. Good values are around 10-200
countFactor2=200;
maxCount=1;

%% Complete tomographic set: Measurement basis

% Boyd or Thew
switch b
    case 1 % Boyd
        chi={[1;0],[0;1],1/sqrt(2)*([1;0]+[0;1]),1/sqrt(2)*([1;0]+[0;1]*exp(1i*pi/2.0))};
        basis='Boyd';
    case 2 % Thew
        phase=exp(2i*pi/3.0);
        chi={[1;0],1/sqrt(2)*([1;0]+[0;1]),1/sqrt(2)*([1;0]+[0;1]*phase),1/sqrt(2)*([1;0]+[0;1]*phase^2)};
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
    
% Wave function
alpha = 1;
beta = 1;
psi = alpha*kron([1;0],[1;0]) + beta*kron([0;1],[0;1]);
psi = psi*1/sqrt(alpha^2+beta^2); % Normalization
rho_Ideal = psi*psi'; % Density matrix of Ideal State

% Add white noise
lambda=0.8;
rho_Ideal=rho_Ideal * lambda + (1-lambda)/4*eye(4);

% BEHAVIOUR of Ideal State
P_Ideal=behaviour(MM, rho_Ideal);

% Behaviour of state with systematic signalling
sig_coeff=linspace(0,0.99,nAlpha);


%% Generate Systematic Signaling Behaviour
P_sys_sig=P_Ideal;

P_sys_sig(:,:,1,1) = [1 0;0 0];

%% Start Simulation
P_Noise=cell(length(sig_coeff),nMeasurement,4);



for m=1:length(sig_coeff)
    rng default;
    for ii=1:nMeasurement
        fprintf('\n#### \t Start Sim.\t %i von\t %i \t #### \n', ii, nMeasurement);

        P = P_Noise(m,ii,:);

        [P{1}, Sig] = P_Noise_Poisson(P_Ideal, maxCount,countFactor1);
        [P{2}, Sig] = P_Noise_Poisson(P_Ideal, maxCount,countFactor2);


        % Apply Noise on each Measurement
        %[P{1}, Sigma{ii,1}] = P_Noise_Poisson(P_Ideal, maxCount,countFactor);
                
        %P_SYS_SIG = P{1}*(1-sig_coeff(m)) + sig_coeff(m)*P_sys_sig;
        
        %[P{2}, Sigma{ii,1}] = P_Noise_Poisson(P_SYS_SIG, maxCount,countFactor);
        P{3} = P_Noise_Poisson(P_Ideal*(1-sig_coeff(m)) + sig_coeff(m)*P_sys_sig, maxCount,countFactor1);
        P{4} = P_Noise_Poisson(P_Ideal*(1-sig_coeff(m)) + sig_coeff(m)*P_sys_sig, maxCount,countFactor2);
        P_Noise(m,ii,:) = P;

    end
    check(m) = checkNormalization(P{2});

end
%%
for m=1:length(sig_coeff)
    % Calculate average Behaviour
    for j=1:4
        P = cell2mat(cellfun(@(A) {reshape(A,36,1)}, {P_Noise{m,:,j}})); % Behaviour with Noise after Regularization
        P_mean{m,j}=reshape(mean(P,2),2,2,3,3);
    end
end

%% Calculate distances from P_mean to NS Polytope
for m=1:length(sig_coeff)
    for j=1:4
        P_Reg{m,j}=FindNA_KL_PENLAB(P_mean{m,j},[]);
    end
end

%% Calculate Kullback-Leibler Divergence
for i=1:nMeasurement
   D_KL_stat1(i) = kullback_leibler_divergence(FindNA_KL_PENLAB(P_Noise{1,i,1},[],0),P_Noise{1,i,1});
   D_KL_stat2(i) = kullback_leibler_divergence(FindNA_KL_PENLAB(P_Noise{1,i,2},[],0),P_Noise{1,i,2});
end

for m=1:length(sig_coeff)
    for i=1:nMeasurement
       D_KL_sys1(m,i) = kullback_leibler_divergence(FindNA_KL_PENLAB(P_Noise{m,i,3},[],0),P_Noise{m,i,3});
       D_KL_sys2(m,i) = kullback_leibler_divergence(FindNA_KL_PENLAB(P_Noise{m,i,4},[],0),P_Noise{m,i,4});
    end
end

%%
% Std-Deviation for kullback-leibler divergence    
stdDev_stat1 = std(D_KL_stat1(:));
stdDev_stat2 = std(D_KL_stat2(:));

D_KL_mean_stat1=mean(D_KL_stat1(:));
D_KL_mean_stat2=mean(D_KL_stat2(:));
fprintf('<D(reg,stat)> / sigma = %f\t Noise = %i\n', D_KL_mean_stat1/stdDev_stat1,countFactor1)
fprintf('<D(reg,stat)> / sigma = %f\t Noise = %i\n', D_KL_mean_stat2/stdDev_stat2,countFactor2)

for m=1:length(sig_coeff)
    stdDev_sys1(m) = std(D_KL_sys1(m,:)); % Systematic signalling
    stdDev_sys2(m) = std(D_KL_sys2(m,:)); % Systematic signalling
    D_KL_mean_sys1(m)=mean(D_KL_sys1(m,:));
    D_KL_mean_sys2(m)=mean(D_KL_sys2(m,:));
    fprintf('<D(reg,sys)> / sigma = %f\t Noise = %i\n', D_KL_mean_sys1(m)/stdDev_sys1(m),countFactor1)
    fprintf('<D(reg,sys)> / sigma = %f\t Noise = %i\n', D_KL_mean_sys1(m)/stdDev_sys1(m),countFactor2)

end





%% Plot
f=figure();
set(gcf, 'Position', get(0, 'Screensize'));
plot(sig_coeff(1:end-1),D_KL_mean_sys1(1:end-1)./stdDev_sys1(1:end-1),'DisplayName','$P_\alpha^{stat+sys}=\alpha P^{sys} + (1-\alpha)P^{stat}$, Low Noise','Marker','s','LineWidth',1)
hold on
plot(sig_coeff(1:end-1),D_KL_mean_sys2(1:end-1)./stdDev_sys1(1:end-1),'DisplayName','$P_\alpha^{stat+sys}=\alpha P^{sys} + (1-\alpha)P^{stat}$, High Noise','Marker','s','LineWidth',1)
xlabel('$\alpha$','interpreter','latex')
yl=ylabel('$\frac{D_{KL}(P_{reg}^\alpha||P_\alpha^{stat+sys})}{\sigma_\alpha}$','interpreter','latex');
set(gca,'FontSize',12)
h=legend('show');
set(h,'interpreter','latex','Location','NorthWest')
%set(get(gca,'YLabel'),'rotation',0)
set(get(gca,'YLabel'),'FontSize',16)
set(get(gca,'XLabel'),'FontSize',16)
h1 = get(yl, 'position');
%set(yl, 'position', [-0.1 h1(2) h1(3)])
set(gca,'FontSize',20)
saveas(gcf,'OutputSystematicSignaling/dkl_alpha_plot_100_200.fig')
saveas(gcf,'OutputSystematicSignaling/dkl_alpha_plot_100_200','eps')

%% Plot Behaviour
figure()
subplot(2,1,1)
plot(reshape(P_mean{1},36,1),'b','DisplayName','$\langle P_{stat}\rangle$')
hold on
plot(reshape(P_Ideal,36,1),'k','DisplayName','$P^{ideal}$')
plot(reshape(P_Reg{1},36,1),'r--','DisplayName','$P^{reg}$')
h=legend('show');
set(h,'interpreter','latex')
title('Behaviour with statistic signalling')

subplot(2,1,2)
plot(reshape(P_Noise{10,1,3},36,1),'b','DisplayName','$P^{stat+sys}_{\alpha=0.47}$ low noise')
hold on
plot(reshape(P_Noise{10,1,4},36,1),'b','DisplayName','$P^{stat+sys}_{\alpha=0.47}$ high noise')
plot(reshape(P_Ideal,36,1),'k','DisplayName','$P_{ideal}$')
plot(reshape(P_Reg{2},36,1),'r--','DisplayName','$P_{reg}$')
h=legend('show');
set(h,'interpreter','latex')
title('Example Behaviour with systematic signalling')

%% Histogramm Plot
figure
set(gcf, 'Position', get(0, 'Screensize'));
edges = [0:0.01:1];
subplot(2,1,1)
histogram(D_KL_sys1(1,:),edges,'DisplayName','$D_{KL}(P_{reg}^\alpha,P_{\alpha=0}^{stat+sys})$','Normalization','probability')
hold on
histogram(D_KL_sys1(11,:),edges,'DisplayName','$D_{KL}(P_{reg}^\alpha,P_{\alpha=0.47}^{stat+sys})$','Normalization','probability')
histogram(D_KL_sys1(end-1,:),edges,'DisplayName','$D_{KL}(P_{reg}^\alpha,P_{\alpha=0.95}^{stat+sys})$','Normalization','probability')
h=legend('show');
set(h,'interpreter','latex')   
xlim([0 0.6])
ylim([0 0.4])
xlabel('Kullback-Leibler Divergence','interpreter','latex')
ylabel('Probability','interpreter','latex')
title(['\textbf{count factor } $\mathcal{N}=',num2str(countFactor1),'$'],'interpreter','latex')
set(gca,'FontSize',20)


subplot(2,1,2)
histogram(D_KL_sys2(1,:),edges,'DisplayName','$D_{KL}(P_{reg}^\alpha,P_{\alpha=0}^{stat+sys})$','Normalization','probability')
hold on
histogram(D_KL_sys2(11,:),edges,'DisplayName','$D_{KL}(P_{reg}^\alpha,P_{\alpha=0.47}^{stat+sys})$','Normalization','probability')
histogram(D_KL_sys2(end-1,:),edges,'DisplayName','$D_{KL}(P_{reg}^\alpha,P_{\alpha=0.95}^{stat+sys})$','Normalization','probability')
h=legend('show');
set(h,'interpreter','latex')   
xlim([0 0.6])
ylim([0 0.4])
xlabel('Kullback-Leibler Divergence','interpreter','latex')
ylabel('Probability','interpreter','latex')
title(['\textbf{count factor } $\mathcal{N}=',num2str(countFactor2),'$'],'interpreter','latex')
set(gca,'FontSize',20)

saveas(gcf,'OutputSystematicSignaling/histogram_Dkl_100_200.fig')
saveas(gcf,'OutputSystematicSignaling/histogram_Dkl_100_200','eps')
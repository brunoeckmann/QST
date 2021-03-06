close('all')
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

% Choose White Noise Model lambda factor
lambda = 0.8;

% Choose Noise Factor
countFactor=100; % >0, if countFactor is low, Noise is higher. Good values are around 10-200
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


%% Generate Artificial Signaling Behaviour
p_sig=zeros(2,2,3,3);

p_sig(:,:,1,1) = [0.5 0;0 0.5];
p_sig(:,:,1,2) = [0.5 0;0 0.5];
p_sig(:,:,1,3) = [1 0;0 0];
p_sig(:,:,2,1) = [0.5 0;0 0.5];
p_sig(:,:,2,2) = [0.5 0;0 0.5];
p_sig(:,:,2,3) = [0.5 0;0 0.5];
p_sig(:,:,3,1) = [0.5 0;0 0.5];
p_sig(:,:,3,2) = [0.5 0;0 0.5];
p_sig(:,:,3,3) = [0.5 0;0 0.5];


%% Initial State
    
% Wave function
alpha = 1;
beta = 1;
psi = alpha*kron([1;0],[1;0]) + beta*kron([0;1],[0;1]);
psi = psi*1/sqrt(alpha^2+beta^2); % Normalization
rho_Ideal = psi*psi'; % Density matrix of Ideal State

% BEHAVIOUR of Ideal State
P_Ideal=behaviour(MM, rho_Ideal);

% Behaviour of state with systematic signalling
sig_coeff=linspace(0,1,3);



%% Start Simulation
P_Noise=cell(length(sig_coeff),nMeasurement,2);

for m=1:length(sig_coeff)
    P_SIG = P_Ideal*(1-sig_coeff(m)) + sig_coeff(m)*p_sig;
    check = checkNormalization(P_SIG)


    for ii=1:nMeasurement
        fprintf('\n#### \t Start Sim.\t %i von\t %i \t #### \n', ii, nMeasurement);

        P = P_Noise(m,ii,:);



        % Apply Noise on each Measurement
        [P{1}, Sigma{ii,1}] = P_Noise_Poisson(P_Ideal, maxCount,countFactor);
        [P{2}, Sigma{ii,1}] = P_Noise_Poisson(P_SIG, maxCount,countFactor);

        P_Noise(m,ii,:) = P;

    end
end
%%
for m=1:length(sig_coeff)
    % Calculate average Behaviour
    for j=1:2
        P = cell2mat(cellfun(@(A) {reshape(A,36,1)}, {P_Noise{m,:,j}})); % Behaviour with Noise after Regularization
        P_mean{m,j}=reshape(mean(P,2),2,2,3,3);
    end
    
end

%% Calculate distances from P_mean to NS Polytope
for m=1:length(sig_coeff)
    P_Reg{m,1}=FindNA_KL_PENLAB(P_mean{m,1},[]);
    P_Reg{m,2}=FindNA_KL_PENLAB(P_mean{m,2},[]);
end

%% Calculate Kullback-Leibler Divergence
for m=1:length(sig_coeff)
    for i=1:nMeasurement
       D_KL(m,i,1) = kullback_leibler_divergence(P_Reg{m,1},P_Noise{m,i,1});
       D_KL(m,i,2) = kullback_leibler_divergence(P_Reg{m,2},P_Noise{m,i,2});
    end

stdDev(m,1) = std(D_KL(m,:,1));
stdDev(m,2) = std(D_KL(m,:,2));

D_KL_mean(m,1)=kullback_leibler_divergence(P_Reg{m,1},P_mean{m,1});
D_KL_mean(m,2)=kullback_leibler_divergence(P_Reg{m,2},P_mean{m,2});

fprintf('<D(reg,stat)> / sigma = %f\n', D_KL_mean(m,1)/stdDev(m,1))
fprintf('<D(reg,syst)> / sigma = %f\n', D_KL_mean(m,2)/stdDev(m,2))

end

%% Plot
figure()
plot(D_KL_mean(:,1)./stdDev(:,1))

%% Plot Behaviour
figure()
subplot(2,1,1)
plot(reshape(P_mean{1},36,1),'b','DisplayName','$\langle P_{stat}\rangle$')
hold on
plot(reshape(P_Ideal,36,1),'k','DisplayName','$P_{ideal}$')
plot(reshape(P_Reg{1},36,1),'r--','DisplayName','$P_{reg}$')
h=legend('show');
set(h,'interpreter','latex')
title('Behaviour with statistic signalling')

subplot(2,1,2)
plot(reshape(P_mean{2},36,1),'b','DisplayName','$\langle P_{sys}\rangle$')
hold on
plot(reshape(P_Ideal,36,1),'k','DisplayName','$P_{ideal}$')
plot(reshape(P_Reg{2},36,1),'r--','DisplayName','$P_{reg}$')
h=legend('show');
set(h,'interpreter','latex')
title('Behaviour with systematic signalling')

%% Histogramm Plot
figure
histogram(D_KL(:,1),'DisplayName','$D_{KL}(p_{reg},p_{stat})$')
hold on
histogram(D_KL(:,2),'DisplayName','$D_{KL}(p_{reg},p_{sys})$')
h=legend('show');
set(h,'interpreter','latex')

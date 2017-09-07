close('all')
%% User Input

% Choose Basis b
% 1 = Boyd
% 2 = Thew
b=1;

% Choose Number of Simulations
nMeasurement=2;

% Create Output file of command window
diary('log_sequential.out')
diary on

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

%plotRho(rho_Ideal,[outputFolder,'rho_targetState.pdf'],'\textbf{Target State $\rho_{ideal}$}');

n = countSignal(Chi,rho_Ideal);
P=behaviour(MM, rho_Ideal); 

fprintf('\n#### \t Ideal State QST \t #### \n');
rho_Ideal_qst=qst_maximumlikelihood(rho_Ideal,Chi,n);
fprintf('\n#### \t done \t #### \n');


% BEHAVIOUR of Ideal State
P_Ideal=behaviour(MM, rho_Ideal);
lambda = 0.5;

%rho_WNM = rho_Ideal * lambda + (1-lambda)/4*eye(4);
%P_WNM = behaviour(MM, rho_WNM);
%[rho_WNM_QST, M_WNM] = qst_linearinversion(Chi,P_QST_Selection(P_WNM));


%% Start Simulation
P_Noise=cell(nMeasurement,4);
M_Noise=cell(nMeasurement,2);
rho_Noise=cell(nMeasurement,4);
ii=1;

fprintf('\n Starting sequential Simulations... \n');

for ii=1:nMeasurement
    fprintf('\n#### \t Start Sim.\t %i von\t %i \t #### \n', ii, nMeasurement);

    P = P_Noise(ii,:);
    M = M_Noise(ii,:);
    r = rho_Noise(ii,:);
     
    % Apply Noise on each Measurement
    countFactor=100;
    maxCount=1;
    [P{1}, Sigma{ii,1}] = P_Noise_Poisson(P_Ideal, maxCount,countFactor);

    % Project on NS Polytop
    P{2} = nonSignaling(P{1});


    % QST Linear Inversion: Reconstruct Density Matrix
    tic;
    [r{1}, M{1}] = qst_linearinversion(Chi,P_QST_Selection(P{1}));
    t1 = toc;
    tic;
    [r{2}, M{2}] = qst_linearinversion(Chi,P_QST_Selection(P{2}));
    t2 = toc;
    
    % QST Maximum Likelihood: Reconstruct Density Matrix
    tic;
    r{3} = qst_maximumlikelihood(r{2},Chi,P_QST_Selection(P{1}));
    t3 = toc;
    tic;
    r{4} = qst_maximumlikelihood(r{2},Chi,P_QST_Selection(P{2}));
    t4 = toc;
    
    % Save data
    P_Noise(ii,:) = P;
    M_Noise(ii,:) = M;
    rho_Noise(ii,:) = r;
    
    % Output
    fprintf('\n#### \t Summary Sim.\t %i of %i \t #### \n', ii, nMeasurement);
    fprintf('\t QST Linear Inversion for Noisy Behaviour in %f Seconds\n',t1);
    fprintf('\t QST Linear Inversion for Projected Behaviour in %f Seconds\n',t2);
    fprintf('\t QST Maximum Likelihood for Noisy Behaviour in %f Seconds\n',t3);
    fprintf('\t QST Maximum Likelihood for Projected Behaviour in %f Seconds\n',t4);
    fprintf('#### \t \t \t #### \t \t \t \t ####\n');
    
end
diary off




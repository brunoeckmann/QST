close('all')
%% User Input

% Choose Basis b
% 1 = Boyd
% 2 = Thew
b=1;

% Choose Number of Simulations
nMeasurement=10;

% Create Output file of command window
diary
diary('log.txt')
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

rho_Ideal_qst=qst_maximumlikelihood(rho_Ideal,Chi,n);

% BEHAVIOUR of Ideal State
P_Ideal=behaviour(MM, rho_Ideal);
lambda = 0.5;

%rho_WNM = rho_Ideal * lambda + (1-lambda)/4*eye(4);
%P_WNM = behaviour(MM, rho_WNM);
%[rho_WNM_QST, M_WNM] = qst_linearinversion(Chi,P_QST_Selection(P_WNM));


%% Start Simulation
Behaviour_P=cell(nMeasurement,4);
M_Noise=cell(nMeasurement,2);
rho_Noise=cell(nMeasurement,4);

parfor ii=1:nMeasurement
    P = Behaviour_P(ii,:);
    M = M_Noise(ii,:);
    r = rho_Noise(ii,:);
     
    fprintf('\n#### \t Sim.\t %i von\t %i \t #### \n', ii, nMeasurement);

    % Apply Noise on each Measurement
    countFactor=100;
    maxCount=1;
    fprintf('\t Apply Noise to Ideal State\n');
    [P{1}, Sigma{ii,1}] = P_Noise_Poisson(P_Ideal, maxCount,countFactor);

    % Project on NS Polytop
    fprintf('\t Project Behaviour to Local Polytope\n');
    P{2} = nonSignaling(P{1});


    % QST Linear Inversion: Reconstruct Density Matrix
    fprintf('\t Start QST Linear Inversion for Noisy Behaviour\n');
    tic;
    [r{1}, M{1}] = qst_linearinversion(Chi,P_QST_Selection(P{1}));
    t = toc;
    fprintf('\t Elapsed time: %f \n',t);
    fprintf('\t Start QST Linear Inversion for Projected Behaviour\n');
    tic;
    [r{2}, M{2}] = qst_linearinversion(Chi,P_QST_Selection(P{2}));
    t = toc;
    fprintf('\t Elapsed time: %f \n',t);
    
    % QST Maximum Likelihood: Reconstruct Density Matrix
    fprintf('\t Start QST Maximum Likelihood for Noisy Behaviour\n');
    tic;
    r{3} = qst_maximumlikelihood(r{2},Chi,P_QST_Selection(P{1}));
    t = toc;
    fprintf('\t Elapsed time: %f \n',t);
    fprintf('\t Start QST Maximum Likelihood for Projected Behaviour\n');
    tic;
    r{4} = qst_maximumlikelihood(r{2},Chi,P_QST_Selection(P{2}));
    t = toc;
    fprintf('\t Elapsed time: %f \n',t);
    
    Behaviour_P(ii,:) = P;
    M_Noise(ii,:) = M;
    rho_Noise(ii,:) = r;
    
    
    fprintf('#### \t #### \t ####\n');
    
    

end






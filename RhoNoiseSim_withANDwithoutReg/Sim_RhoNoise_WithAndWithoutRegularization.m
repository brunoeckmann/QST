% Simulation RhoNoise with and without Regularization
% Units:    [fsec, mum]         
% Date:     07.12.2017
% Author:   Bruno Eckmann, Sacha Schwarz
%
% Abstact:  As a function of the total sample number N, the density matrix
%           reconstruction of noisy relative frequencies as well as of 
%           regularized noisy relative frequencies is performed by means of
%           LIN or ML as reconstruction scheme.
%
% Comments: -> Code is written for mixed two-qubit states as target state
%           -> Boyd Basis is investigated
%           -> SOLVER is only fmin investigated
% 
%**************************************************************************
clear all
close all
clc
scrsz = get(0,'ScreenSize');

set(0,'defaulttextinterpreter','latex')

%% USER INPUT
%--------------------------------------------------------------------------
% total number of simulation samples
nMeasurement = 100;

% mixing parameter for Werner state
lambda = 0.52;

% noise factor
N_tot=5120; % >0, if countFactor is low, Noise is higher. Good values are around 10-200
maxCount=1;

%% DEFINE MIXED TWO_QUBIT STATE
%--------------------------------------------------------------------------
% state dimension yields corresponding Bell scenario
d = 2;
numbInputs = 2*d-1;
numbOfOutputs = d*(numbInputs-d);

% pure quantum state
ket0 = [1 0];
ket1 = [0 1];
psi = kron(ket0,ket0) + kron(ket1,ket1);
psiNorm = psi*psi';
psi = sqrt(1/psiNorm)*psi;
RhoPure = psi'*psi;

% Werner state
RhoTarget = lambda*RhoPure+(1-lambda)/d^2*eye(d^2);

%% DEFINING OUTPUT ARRAY DIMENSIONS
%--------------------------------------------------------------------------
P_Noise=cell(nMeasurement,3);
M_Noise=cell(nMeasurement,2);
rho_Noise=cell(nMeasurement,9);
tElapsed=zeros(nMeasurement,8);

%% DEFINE TOMOGRAPHIC COMPLETE SET AND LOCAL MEASUREMENT BASIS
%--------------------------------------------------------------------------
% define Boyd basis
chi={[1;0], ...
     [0;1], ...
     1/sqrt(2)*([1;0]+[0;1]), ...
     1/sqrt(2)*([1;0]+[0;1]*exp(1i*pi/2.0))};
 
% Kronecker product of measurment basis
Chi=cell(16,1);
for ii=1:4
    for jj=1:4
       Chi{4*(ii-1)+jj}=kron(chi{jj},chi{ii});
    end
end

% Extended Tomographic Set: Basis for NS projection
chi{5}=1/sqrt(2)*([1;0]-[0;1]);
chi{6}=1/sqrt(2)*([-1;1i]);

ChiProj=cell(36,1);
for ii=1:6
   for j=1:6
      ChiProj{6*(ii-1)+j}=kron(chi{ii},chi{j}); 
   end
end

% define matrix for simulation
MM=cell(3,2);
MM{1,1}=chi{1};
MM{1,2}=chi{2};
MM{2,1}=chi{3};
MM{2,2}=chi{5};
MM{3,1}=chi{4};
MM{3,2}=chi{6};

%% START SIMULATION
%--------------------------------------------------------------------------
% initial behavior
P_WNM = behaviour(MM, RhoTarget);

%fprintf('\n Starting parallel Simulations... \n');
%poolobj = parpool;
%ppm = ParforProgMon('Simulation Progress ', nMeasurement);

for ii=1:nMeasurement
    fprintf('\n#### \t Start Sim.\t %i von\t %i \t #### \n', ii, nMeasurement);

    P = P_Noise(ii,:);
    M = M_Noise(ii,:);
    r = rho_Noise(ii,:);
    t = tElapsed(ii,:);

    % Apply Noise on each Measurement
    [P{1}, Sigma{ii,1}] = P_Noise_Poisson(P_WNM, maxCount,N_tot);

    % Regularization on NS Polytop
    [P{2} flag prob] = FindNA_KL_PENLAB(P{1},[]);

    % Regularization to Q
    [P{3} flag prob] = FindNA_KL_PENLAB(P{1},qvx.di.HierarchyLevel.local(1));
    r{9} = qst_linearinversion(Chi,P_QST_Selection(P{3}));
    
    % QST Linear Inversion: Reconstruct Density Matrix
    tic;
    [r{1}, M{1}] = qst_linearinversion(Chi,P_QST_Selection(P{1}));
    t(1) = toc;
    tic;
    [r{2}, M{2}] = qst_linearinversion(Chi,P_QST_Selection(P{2}));
    t(2) = toc; 

    formatSpec = '%i\t %e\t %i\t %s\n';

    % QST Maximum Likelihood fminsearch: Reconstruct Density Matrix
    fprintf('\n#### \t ML\t %i von\t %i \t #### \n', ii, nMeasurement);
    tic;
    [r{3},Output] = qst_maximumlikelihood_fmin(r{2},Chi,P_QST_Selection(P{1}));
    t(3) = toc;
    tic;
    [r{4},Output] = qst_maximumlikelihood_fmin(r{2},Chi,P_QST_Selection(P{2}));
    t(4) = toc;

    % Save data
    P_Noise(ii,:) = P;
    M_Noise(ii,:) = M;
    rho_Noise(ii,:) = r;
    tElapsed(ii,:) = t;

    % Output
    fprintf('\n#### \t Summary Sim.\t %i of %i \t #### \n', ii, nMeasurement);
    fprintf('\t QST Linear Inversion for Noisy Behaviour in %f Seconds\n',t(1));
    fprintf('\t QST Linear Inversion for Projected Behaviour in %f Seconds\n',t(2));
    fprintf('\t QST Maximum Likelihood fminsearch for Noisy Behaviour in %f Seconds\n',t(3));
    fprintf('\t QST Maximum Likelihood fminsearch for Projected Behaviour in %f Seconds\n',t(4));
    fprintf('\t QST Maximum Likelihood ga for Noisy Behaviour in %f Seconds\n',t(5));
    fprintf('\t QST Maximum Likelihood ga for Projected Behaviour in %f Seconds\n',t(6));
    fprintf('\t QST Maximum Likelihood lstsqr for Noisy Behaviour in %f Seconds\n',t(7));
    fprintf('\t QST Maximum Likelihood lstsqr for Projected Behaviour in %f Seconds\n',t(8));
    fprintf('#### \t \t \t #### \t \t \t \t ####\n');
    
    %ppm.increment(); 

end
%delete(poolobj)
diary off

%% SAVE RHO_NOISE TO MAT
%--------------------------------------------------------------------------
% if(exist('RESULTS','dir') ~= 7)
%     mkdir RESULTS
% end
%matName = sprintf('RESULTS/rho_Noise_%d.mat',nMeasurement);
save('Rho_Pure')
save('Rho_Target')
matName = sprintf('rho_Noise_R%d_N%d.mat',nMeasurement,N_tot);
if N_tot < 1000
    matName = sprintf('rho_Noise_R%d_N0%d.mat',nMeasurement,N_tot);
end
if N_tot < 100
    matName = sprintf('rho_Noise_R%d_N00%d.mat',nMeasurement,N_tot);
end
if N_tot < 10
    matName = sprintf('rho_Noise_R%d_N000%d.mat',nMeasurement,N_tot);
end        

save(matName,'rho_Noise')







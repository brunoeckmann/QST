close('all')
fclose('all');
poolobj = gcp('nocreate');
delete(poolobj);
%% User Input

% Choose Basis b
% 1 = Boyd
% 2 = Thew
b=1;

% Choose Number of Simulations
nMeasurement=100;

% Choose White Noise Model lambda factor
lambda = 0.8;

% Choose Noise Factor
countFactor=100; % >0, if countFactor is low, Noise is higher. Good values are around 10-200
maxCount=1;

% Choose Active MaximumLikelihood Solvers
ga=1;
fmin=1;
lstsqr=1;

% Choose if Solver Convergence Output is created
solverOutput=1;
OutputFolder = 'OutputSolverTest100_Noise100_WNM08/';

% Choose Logfile
logfile = [OutputFolder,'log.txt'];


%% Activate Logfile
if(exist(OutputFolder)==0)
    mkdir(OutputFolder)
end

if (exist(logfile))
    delete(logfile);
end
diary(logfile)
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
%P=behaviour(MM, rho_Ideal); 

fprintf('\n#### \t Ideal State QST \t #### \n');
tic
[rho_Ideal_qst, Output_Ideal]=qst_maximumlikelihood_lstsqr(rho_Ideal,Chi,n);
toc
fprintf('\n#### \t done \t #### \n');


% BEHAVIOUR of Ideal State
P_Ideal=behaviour(MM, rho_Ideal);

% White Noise Model
rho_WNM = rho_Ideal * lambda + (1-lambda)/4*eye(4);
P_WNM = behaviour(MM, rho_WNM);
fprintf('\n#### \t Ideal State with White Noise QST \t #### \n');
tic
[rho_WNM_qst, Output_WNM] = qst_maximumlikelihood_lstsqr(rho_WNM,Chi,P_QST_Selection(P_WNM));
toc
fprintf('\n#### \t done \t #### \n');


%% Start Simulation
P_Noise=cell(nMeasurement,2);
M_Noise=cell(nMeasurement,2);
rho_Noise=cell(nMeasurement,8);
tElapsed=zeros(nMeasurement,8);

fprintf('\n Starting parallel Simulations... \n');
%poolobj = parpool;
ppm = ParforProgMon('Simulation Progress', nMeasurement);

parfor ii=1:nMeasurement
    fprintf('\n#### \t Start Sim.\t %i von\t %i \t #### \n', ii, nMeasurement);

    P = P_Noise(ii,:);
    M = M_Noise(ii,:);
    r = rho_Noise(ii,:);
    t = tElapsed(ii,:);
    
    % Apply Noise on each Measurement
    [P{1}, Sigma{ii,1}] = P_Noise_Poisson(P_WNM, maxCount,countFactor);

    % Project on NS Polytop
    P{2} = nonSignaling(P{1});

    
    % QST Linear Inversion: Reconstruct Density Matrix
    tic;
    [r{1}, M{1}] = qst_linearinversion(Chi,P_QST_Selection(P{1}));
    t(1) = toc;
    tic;
    [r{2}, M{2}] = qst_linearinversion(Chi,P_QST_Selection(P{2}));
    t(2) = toc;
    
    formatSpec = '%i\t %e\t %i\t %s\n';
    
    if fmin==1
        
        % QST Maximum Likelihood fminsearch: Reconstruct Density Matrix
        tic;
        [r{3},Output] = qst_maximumlikelihood_fmin(r{2},Chi,P_QST_Selection(P{1}));
        fileID = fopen([OutputFolder, 'solver_fmin_noise_',num2str(ii),'.txt'],'w');
        O=Output.';
        fprintf(fileID,formatSpec, O{:,:});
        fclose(fileID);
        t(3) = toc;
        tic;
        [r{4},Output] = qst_maximumlikelihood_fmin(r{2},Chi,P_QST_Selection(P{2}));
        fileID = fopen([OutputFolder, 'solver_fmin_proj_',num2str(ii),'.txt'],'w');
        O=Output.';
        fprintf(fileID,formatSpec, O{:,:});
        fclose(fileID);
        t(4) = toc;
    end
    
    if ga==1
        % QST Maximum Likelihood ga: Reconstruct Density Matrix
        tic;
        [r{5},Output] = qst_maximumlikelihood_ga(r{2},Chi,P_QST_Selection(P{1}));
        fileID = fopen([OutputFolder, 'solver_ga_noise_',num2str(ii),'.txt'],'w');
        O=Output.';
        fprintf(fileID,formatSpec, O{:,:});
        fclose(fileID);
        t(5) = toc;
        tic;
        [r{6},Output] = qst_maximumlikelihood_ga(r{2},Chi,P_QST_Selection(P{2}));
        fileID = fopen([OutputFolder, 'solver_ga_proj_',num2str(ii),'.txt'],'w');
        O=Output.';
        fprintf(fileID,formatSpec, O{:,:});
        fclose(fileID);
        t(6) = toc;
    end
    
    if lstsqr==1
        % QST Maximum Likelihood lstsqr: Reconstruct Density Matrix
        tic;
        [r{7},Output] = qst_maximumlikelihood_lstsqr(r{2},Chi,P_QST_Selection(P{1}));
        fileID = fopen([OutputFolder, 'solver_lstsqr_noise_',num2str(ii),'.txt'],'w');
        O=Output.';
        fprintf(fileID,formatSpec, O{:,:});        
        fclose(fileID);
        t(7) = toc;
        tic;
        [r{8},Output] = qst_maximumlikelihood_lstsqr(r{2},Chi,P_QST_Selection(P{2}));
        fileID = fopen([OutputFolder, 'solver_lstsqr_proj_',num2str(ii),'.txt'],'w');
        O=Output.';
        fprintf(fileID,formatSpec, O{:,:});
        fclose(fileID);
        t(8) = toc;
    end
    
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
    
    ppm.increment(); 

end
delete(poolobj)
diary off

% %% Test
% P = P_Noise(ii,:);
% M = M_Noise(ii,:);
% r = rho_Noise(ii,:);
% t = tElapsed(ii,:);
% 
% tic
% [A,B] = qst_maximumlikelihood_gs(r{2},Chi,P_QST_Selection(P{2}));
% fileID = fopen([OutputFolder, 'solver_gs_proj_',num2str(ii),'.txt'],'w');
% O=Output.';
% fprintf(fileID,formatSpec, O{:,:});
% fclose(fileID);
% toc



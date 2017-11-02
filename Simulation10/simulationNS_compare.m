%%
% Comparison betweeen NS Subtraction of Alberto and Liang

close('all')
fclose('all');
rng default; %reset random number generator
poolobj = gcp('nocreate');
delete(poolobj);
%% User Input

% Choose Basis b
% 1 = Boyd
% 2 = Thew
b=1;

% Choose NS Projection Method
ns_method = 'liang';
%ns_method = 'alberto';

% Choose Number of Simulations
nMeasurement=5;

% Choose White Noise Model lambda factor
lambda = 1;

% Choose Noise Factor
countFactor=100; % >0, if countFactor is low, Noise is higher. Good values are around 10-200
maxCount=1;

% Choose Active MaximumLikelihood Solvers
ga=1;
fmin=1;
lstsqr=1;

% Choose if Solver Convergence Output is created
solverOutput=1;
OutputFolder = 'OutputLiangTest_1/';

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

%% Simulation
for i=1:nMeasurement
    [P{1,i}, Sigma] = P_Noise_Poisson(P_WNM, maxCount,countFactor);
    P{2,i} = nonSignaling(P{1,i});
    P{3,i} = FindNA_KL_PENLAB(P{1,i},[]);

    
    fig=figure(i)
    set(fig,'units','normalized')
    set(fig,'outerposition',[0 0 0.8 0.8]);
    
    plot(reshape(P_Ideal,36,1),'r-o','DisplayName','Ideal')    
    hold on
    plot(reshape(P{1,i},36,1),'k-o','DisplayName','Noise')    

    plot(reshape(P{2,i},36,1),'c--d','DisplayName','Alberto')
    plot(reshape(P{3,i},36,1),'b--d','DisplayName','Liang')
    legend('show')
end

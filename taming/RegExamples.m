% Generate random signaling observations
nParties = 2;
nInputs = 3;
nOutputs = 2;
scenario = qvx.di.Scenario.homogenous(nParties, nInputs, nOutputs);

PObs = qvx.di.Correlations.randomSignaling(scenario, 'Pabxy');

% PObs has size = [nOutputs,nOutputs...,nInputs,nInputs...], indices (a,b...,x,y...)

%% Nearest nonsignaling / min Kullback-Leibler distance to NS polytope

% Nearest nonsignaling approximation, KL divergence, PENLAB solver

level = []; % set to [] to use the nonsignaling set
PNNS_PL = FindNA_KL_PENLAB(PObs, level);

%% Nearest quantum approximation

% with PENLAB
level = qvx.di.HierarchyLevel.local(1); % almost quantum set
PNQA = FindNA_KL_PENLAB(PObs, level);



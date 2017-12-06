% Simulation RhoNoise with and without Regularization
% Units:    [fsec, mum]         
% Date:     06.12.2017
% Author:   Bruno Eckmann, Sacha Schwarz
%
% Abstact:  As a function of the total sample number N, the density matrix
%           reconstruction of noisy relative frequencies as well as of 
%           regularized noisy relative frequencies is performed by means of
%           LIN or ML as reconstruction scheme.
%
% Comments: -> Code is written for mixed two-qubit states as target state
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
Ntot = 1000;

% mixing parameter for Werner state
lambda = 0.52;

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
rhoPure = psi'*psi;

% Werner state
RhoTarget = lambda*rhoPure+(1-lambda)/d^2*eye(d^2);

%% DEFINING OUTPUT ARRAY DIMENSIONS
%--------------------------------------------------------------------------
% behaviors
Pvec = zeros(numbInputs^2*numbOfOutputs^2);
Pvec_Noise = zeros(numbInputs^2*numbOfOutputs^2);
Pvec_Noise_Reg = zeros(numbInputs^2*numbOfOutputs^2);

% reconstructed density matrices
Rho_Noise_LIN = zeros(d^4);
Rho_Noise_ML = zeros(d^4);
Rho_Noise_Reg_LIN = zeros(d^4);
Rho_Noise_Reg_ML = zeros(d^4);

%% SIMULATION
%--------------------------------------------------------------------------



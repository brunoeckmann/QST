% Evaluate Fidelity in Large Statistics Limit
% Units:    [fsec, mum]         
% Date:     07.12.2017
% Author:   Sacha Schwarz
%
% Abstact:  Since in the large statistics limit the fidelity must be 
%           approximated by every reconstruction scheme. We show that with 
%           out method this is done more accurately.
%
% Comments: -> 
% 
%**************************************************************************
clear all
close all
clc
scrsz = get(0,'ScreenSize');

set(0,'defaulttextinterpreter','latex')

%% USER INPUT
%--------------------------------------------------------------------------


%% INPUT DATA
%--------------------------------------------------------------------------
load('Rho_Pure.mat')
load('Rho_Target.mat')
[mFile,Path] = uigetfile('*.mat','Pick rho_Noise Mat-Files in folder RESULTS','MultiSelect','on')

%% CALCULATE FIDELTIY
%--------------------------------------------------------------------------
F_eff = real(Fidelity(RhoTarget, RhoPure))
F = [];
numbOfFpoints = length(mFile);
if iscell(mFile) == 0
    numbOfFpoints = 1;
end
    
for i = 1:numbOfFpoints
    if numbOfFpoints == 1
        load(mFile)
    else
        load(mFile{i})
    end
    nMeasurements = size(rho_Noise,1);
    for k = 1:4
        for kk = 1:nMeasurements
            F(kk,k) = real(Fidelity(rho_Noise{kk,k}, RhoTarget));
        end
        F_mean(k) = mean(F(:,k));
        F_std(k) = std(F(:,k))/sqrt(nMeasurements);
    end
end
F_mean
F_std








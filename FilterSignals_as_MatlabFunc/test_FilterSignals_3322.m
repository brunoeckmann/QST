% test_FilterSignals_3322
% =======================
% Units:   [fsec, mum]         
% Date:    25.06.2017
% Author:  Sacha Schwarz
%
% Comments: -> 
% 
%**************************************************************************
clear all
close all
clc
scrsz = get(0,'ScreenSize');

%%
% USER INPUT
%--------------------------------------------------------------------------
M = 3;  % number of measurements
W = 2;  % number of outcomes

%%
% INPUT DATA FILE
%--------------------------------------------------------------------------
P = dlmread("TestBehaviour_3322/behaviour_noise.txt");
P_in = reshape(P, [M M W W]);

% === OLD WAY OF IMPORTING 1D ARRAY TO 4D ARRAY ===
% for a = 1:M
%     for b = 1:M
%         for r = 1:W
%             for s = 1:W
%                 P_in(r,s,a,b) = P(s+(r-1)*W+(b-1)*W*W+(a-1)*W*W*M);
%             end
%         end
%     end
% end
% ===

%%
% RUN FILTER_SIGNALS AND PLOT RESULT
%--------------------------------------------------------------------------
P_out = filter_signals_mex(P_in, W, M);
squeeze(sum(sum(P_out,2),1));

%%
% FROM 4D ARRAY TO 1D FOR TESTING
%--------------------------------------------------------------------------
P_out_vec = reshape(P_out, [M*M*W*W 1])








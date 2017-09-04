%% Nonsignaling Projection
% function calls c++ programm written by ????
% Author: Bruno Eckmann, breckman@phys.ethz.ch

% Note: Only works for N=1 !!!

% Input
% Behaviour P(ab|xy) is a (W x W x M x M) 4dim Array
% N: Number of Behaviours (only works for N=1)
% M: Number of Inputs
% W: Number of Measurement Outcomes

% Output
% Projected Behaviour P(ab|xy) on NS polytop

function P_proj=nonSignaling(P)
    filename = 'behaviours.txt';
    N=1;
    W = size(P,1);
    M = size(P,3);
    
    % Ordering P to be used in filter_signals.exe
    % x,y = 1...M
    % a,b = 1...W
    % filter_signals will read the array as:
    % for x...
    %     for y...
    %         for a...
    %             for b...
    
    Psave=[];
    
    for n=1:N
       for x = 1:M
           for y = 1:M
               for a = 1:W
                   for b = 1:W
                       Psave((M*W*W)*(x-1)+(W*W)*(y-1)+W*(a-1)+b,n)=P(a,b,x,y);                  
                   end
               end
           end
       end
    end
    
    % Save P(ab|xy) to file
    dlmwrite(filename,Psave)

    % Call NS Projection
    if isunix==false
        command = 'filter_signals.exe behaviours.txt';
    else
        command = './filterRun behaviours.txt';
    end
    

    [status,cmdout] = system(command);
    
    % Read output
    out=textscan(cmdout,'%f');
    out=out{1}(3:M*M*W*W+2);

    % Store output vector in 4dim array P(ab|xy)
    for n=1:N
       for x = 1:M
           for y = 1:M
               for a = 1:W
                   for b = 1:W
                       P_proj(a,b,x,y)=out((M*W*W)*(x-1)+(W*W)*(y-1)+W*(a-1)+b);
                   end
               end
           end
       end
    end

end
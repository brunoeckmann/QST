function [Pnoisy, sigma]=P_Noise_Poisson(P_Ideal, maxCountRate, countFactor)
% P_NOISE_POISSON
% Adds Poisson distributed noise to Behaviour P_Ideal.
% Input: P-Matrix (Behaviour), maxCountRate, countFactor
% Output: P-Matrix with Noise, the Behaviour is normalized.

% from dark count measurement: lambda=1.8966


    %deltaP=poissrnd(lambda,size(Behaviour));
    %Pnoisy = Behaviour+deltaP;
    
    %P = uint16(Behaviour*scaleFactor); %uint16 needed for propper use of imnoise
    
    % Apply Poisson Noise
    %Pnoisy = double(imnoise(P,'poisson'))/scaleFactor;
    
    % P = P + dP where dP is Poisson distributed
    
    sigma = sqrt(P_Ideal*countFactor)./ sqrt(maxCountRate);
    lambda = sigma.^2;
    
    deltaP = poissrnd(lambda, size(P_Ideal))/countFactor*maxCountRate;
    Pnoisy = P_Ideal+deltaP;
    
    % Normalize: Sum(P(:,:,x,y)) != 1
    [A,B,X,Y]=size(P_Ideal);
    for x=1:X
        for y=1:Y
            N=sum(reshape(Pnoisy(:,:,x,y),A*B,1)); % Normalization factor
            Pnoisy(:,:,x,y) = Pnoisy(:,:,x,y) / N;
        end
    end
end

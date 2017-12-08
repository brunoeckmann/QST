function [Pnoisy, sigma]=P_Noise_Poisson(P_Ideal, maxCountRate, countFactor)
% P_NOISE_POISSON
% Apply Poisson distributed noise to Behaviour P_Ideal.
% Input: P-Matrix (Behaviour), maxCountRate, countFactor
% Output: P-Matrix with Noise, the Behaviour is normalized.

% from dark count measurement: lambda=1.8966

    sigma = sqrt(P_Ideal*countFactor)./ sqrt(maxCountRate);
    lambda = sigma.^2;
    
    Nnoisy = poissrnd(lambda, size(P_Ideal));
    % Normalize: Sum(P(:,:,x,y)) != 1
    [A,B,X,Y]=size(P_Ideal);
    for x=1:X
        for y=1:Y
            N=sum(reshape(Nnoisy(:,:,x,y),A*B,1)); % Normalization factor
            if N < 10e-8
                Pnoisy(:,:,x,y) = Nnoisy(:,:,x,y);
            else
                Pnoisy(:,:,x,y) = Nnoisy(:,:,x,y) / N;
            end
        end
    end
end

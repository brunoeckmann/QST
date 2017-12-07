%% Gauss Noise
% symmetry = 0: only positiv scattering
% symmetry = 1: pos and neg. scattering
function Pnoisy=P_Noise_Gauss(Behaviour, sigma, symmetry)
    if symmetry==0
        Pnoisy = Behaviour+abs(normrnd(0,sigma,size(Behaviour)));
    else
        Pnoisy = normrnd(Behaviour,sigma,size(Behaviour));    
    end
end
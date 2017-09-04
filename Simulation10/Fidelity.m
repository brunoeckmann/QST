%% Fidelity between two density matrices

function F=Fidelity(rho, rhoTarget)
    F = trace(sqrtm(sqrtm(rho)*rhoTarget*sqrtm(rho)));  
end
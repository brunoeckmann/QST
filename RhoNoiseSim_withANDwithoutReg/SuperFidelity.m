%% SuperFidelity between two density matrices

function F=SuperFidelity(rho, rhoTarget)
    F = trace(rho*rhoTarget) + sqrt(1-trace(rho^2))* sqrt(1-trace(rhoTarget^2));
end
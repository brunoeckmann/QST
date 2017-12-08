function [Iopt C] = FindNonsignalingMaximumC(scenario, space, I)
% Finds the nonsignaling maximum of the given inequality
%
% Requires CVX
    cvx_begin
    variable C(length(I)) % the correlation vector is a variable
    C == qvx.di.cvx.NonsignalingSet(scenario, space)
    maximize I * C
    cvx_end
    Iopt = I * C;
end

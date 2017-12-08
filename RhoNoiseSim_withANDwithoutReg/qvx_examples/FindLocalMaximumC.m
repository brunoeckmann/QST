function [Iopt C] = FindLocalMaximumC(scenario, space, I)
% Finds the local maximum of the given inequality
%
% Requires CVX
    cvx_begin
    variable C(length(I)) % variable: the correlation vector
    C == qvx.di.cvx.LocalSet(scenario, space) % C is in the local set
    maximize I * C % maximize
    cvx_end
    Iopt = I * C;
end

function P = LocalSet(scenario, space)
% LocalSet Polytope representing the set of local correlations
%
% LocalSet(scenario, space) represents the local set for the given
% scenario, in the given correlation space.
%
% space: default value 'P'
    if nargin < 2
        space = 'P';
    end
    [M D] = qvx.di.LocalModels.localStrategies(scenario, space);
    cvx_begin set
    variable P(size(M, 1))
    variable q(size(M, 2))
    q >= 0
    D * P == M * q
    sum(q) == 1
    cvx_end
end

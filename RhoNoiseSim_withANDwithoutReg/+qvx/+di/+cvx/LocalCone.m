function P = LocalCone(scenario, space)
    if nargin < 2
        space = 'P';
    end
    [M D] = qvx.di.LocalModels.localStrategies(scenario, space);
    cvx_begin set
    variable P(size(M, 1))
    variable q(size(M, 2))
    q >= 0
    D * P == M * q
    cvx_end
end

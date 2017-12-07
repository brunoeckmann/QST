function C = NonsignalingCone(scenario, space)
% NonsignalingCone Describes the conic extension of the nonsignaling polytope
    if nargin < 2
        space = 'P';
    end
    d = qvx.di.Correlations.dimension(scenario, space);
    cvx_begin set
    variable C(d)
    if isequal(space, 'P')
        [Aeq beq] = qvx.di.Correlations.nonsignalingLinearConstraintsP(scenario);
        Aeq * C == beq
        C >= 0
    elseif isequal(space, 'SG') || isequal(space, 'SC')
        C(~scenario.nonsignalingMask) == 0
        [M D] = scenario.conversionMatrix('P', space);
        M*C >= 0
    elseif isequal(space, 'NG') || isequal(space, 'NC')
        [M D] = scenario.conversionMatrix('P', space);
        M*C >= 0
    end
    cvx_end
end

function C = NonsignalingSet(scenario, space)
% NonsignalingSet Polytope representing the set of nonsignaling correlations
%
% Nonsignaling(scenario, space) represents the nonsignaling set for the given
% scenario, in the given correlation space.
%
% space: default value 'P'
    if nargin < 2
        space = 'P';
    end
    d = qvx.di.Correlations.dimension(scenario, space);
    cvx_begin set
    variable C(d)
    if isequal(space, 'P')
        [Aeq beq] = qvx.di.Correlations.nonsignalingConstraintsP(scenario);
        Aeq * C == beq
        C >= 0
    elseif isequal(space, 'SG') || isequal(space, 'SC')
        C(~scenario.nonsignalingMask) == 0
        [M D] = scenario.conversionMatrix('P', space);
        M*C >= 0
        C(1) == 1
    elseif isequal(space, 'NG') || isequal(space, 'NC')
        [M D] = scenario.conversionMatrix('P', space);
        M*C >= 0
        C(1) == 1
    end
    cvx_end
end

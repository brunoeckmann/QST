function [Iopt C] = FindQuantumMaximumC(scenario, space, I, hierarchyLevel)
% FindQuantumMaximumC Computes the Tsirelson bound of an inequality
%
% FindQuantumMaximumC(scenario, space, I, hierarchyLevel) returns an
% upper bound on the quantum violation
%
% scenario - Scenario considered
% space    - Correlation space (internally, 'NG' is used)
% I        - Inequality coefficients
% hierarchyLevel - Level and type of hierarchy used
%
% This function works internally in the Collins-Gisin representation,
% conversion is done as required.
%
% Requires CVX
    IG = qvx.di.Expressions.convert(scenario, 'NG', space, I);
    cvx_begin
    variable PG(length(IG)) % Correlation vector in Collins-Gisin form
    PG == qvx.di.cvx.QuantumSupersetCG(scenario, hierarchyLevel)
    PG(1) == 1 % normalization, as the moment matrix is not normalized by default
    maximize IG * PG
    cvx_end
    Iopt = IG * PG;
    if nargout > 1
        C = qvx.di.Correlations.convert(scenario, space, 'NG', PG);
    end
end

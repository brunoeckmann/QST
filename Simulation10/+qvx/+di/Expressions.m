classdef Expressions
% Expressions Collection of methods that work on Bell expressions/inequalities
%
% Expressions are row vectors that represents linear maps acting on
% correlations, such as the l.h.s. of Bell inequalities, in a
% correlation space (see also qvx.di.Spaces)
    methods (Static)
        function res = convert(scenario, toSpace, fromSpace, expression)
            expression = expression(:)';
            [M D] = scenario.conversionMatrix(fromSpace, toSpace); % reversed order compared to
                                                                   % correlations
            res = expression * M / D;
        end
        function res = one(scenario, space)
            [M D] = scenario.conversionMatrix('NC', space);
            expression = zeros(1, size(M, 1));
            expression(1) = 1;
            res = expression * M / D;
        end
    end
end

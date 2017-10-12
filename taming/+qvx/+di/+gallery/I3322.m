classdef I3322
    methods (Static)
        function S = scenario
            S = qvx.di.Scenario.homogenous(2,3,2);
        end
        function I = expression(space)
            if nargin < 1
                space = 'P';
            end
            % expression taken from arXiv:1006.3032v1
        %        id A1 A2 A3
            I = [ 0  0 -1  0  % id
                 -1  1  1 -1  % B1
                 -2  1  1  1  % B2
                  0 -1  1  0];% B3
            I = I(:)';
            I = qvx.di.Expressions.convert(qvx.di.gallery.I3322.scenario, space, 'NG', I);
        end
    end
end

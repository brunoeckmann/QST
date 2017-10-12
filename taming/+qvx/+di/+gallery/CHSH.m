classdef CHSH
    methods (Static)
        function S = scenario
            S = qvx.di.Scenario.homogenous(2,2,2);
        end
        function I = expression(space)
            if nargin < 1
                space = 'P';
            end 
            I = [0 0 0
                 0 1 1
                 0 1 -1];
            I = I(:)';
            I = qvx.di.Expressions.convert(qvx.di.gallery.CHSH.scenario, space, 'NC', I);
        end
        function P = PRBox(space)
            if nargin < 1
                space = 'P';
            end
            P = [1 0 0
                 0 1 1
                 0 1 -1];
            P = P(:);
            P = qvx.di.Correlations.convert(qvx.di.gallery.CHSH.scenario, space, 'NC', P);
        end
    end
end

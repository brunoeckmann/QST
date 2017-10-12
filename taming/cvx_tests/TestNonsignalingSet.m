classdef TestNonsignalingSet < matlab.unittest.TestCase
    properties (TestParameter)
        scenario = {qvx.di.gallery.CHSH.scenario};
        %                    qvx.di.gallery.CGLMP.scenario(3) ...
        %                    qvx.di.gallery.CGLMP.scenario(4) ...
        %                    qvx.di.gallery.CGLMP.scenario(5) ...
        %                    qvx.di.gallery.I3322.scenario};
        expression = {qvx.di.gallery.CHSH.expression('NG')};
        %                      qvx.di.gallery.CGLMP.expression(3, 'NG') ...
        %                      qvx.di.gallery.CGLMP.expression(4, 'NG') ...
        %                      qvx.di.gallery.CGLMP.expression(5, 'NG') ...
        %                      qvx.di.gallery.I3322.expression('NG')};
        maximum = {4};
                  
    end
    methods (Test, ParameterCombination='sequential')
        function testSize(testCase, scenario, expression, maximum)
            cvx_begin
            variable PG(length(expression)) % Correlation vector in Collins-Gisin form
            PG == qvx.di.cvx.NonsignalingSet(scenario, 'NG');
            maximize expression * PG
            cvx_end
            Iopt = expression * PG;
            tol = 1e-4;
            assert(abs(Iopt - maximum) < tol);
        end
    end
end

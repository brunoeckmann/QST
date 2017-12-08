classdef TestNPAHierarchy < matlab.unittest.TestCase
    properties (TestParameter)
        scenario = {qvx.di.gallery.CHSH.scenario ...
                    qvx.di.gallery.CGLMP.scenario(3) qvx.di.gallery.CGLMP.scenario(3) ...
                    qvx.di.gallery.CGLMP.scenario(4) qvx.di.gallery.CGLMP.scenario(4) ...
                    qvx.di.gallery.CGLMP.scenario(5) qvx.di.gallery.CGLMP.scenario(5) ...
                    qvx.di.gallery.I3322.scenario qvx.di.gallery.I3322.scenario ...
                    qvx.di.gallery.I3322.scenario qvx.di.gallery.I3322.scenario};
        expression = {qvx.di.gallery.CHSH.expression('NG') ...
                      qvx.di.gallery.CGLMP.expression(3, 'NG') qvx.di.gallery.CGLMP.expression(3, 'NG') ...
                      qvx.di.gallery.CGLMP.expression(4, 'NG') qvx.di.gallery.CGLMP.expression(4, 'NG') ...
                      qvx.di.gallery.CGLMP.expression(5, 'NG') qvx.di.gallery.CGLMP.expression(5, 'NG') ...
                      qvx.di.gallery.I3322.expression('NG') qvx.di.gallery.I3322.expression('NG') ...
                      qvx.di.gallery.I3322.expression('NG') qvx.di.gallery.I3322.expression('NG')};
        level = {qvx.di.HierarchyLevel.NPA(1) ...
                 qvx.di.HierarchyLevel.NPAplus(1) qvx.di.HierarchyLevel.almostQuantum ... 
                 qvx.di.HierarchyLevel.NPAplus(1) qvx.di.HierarchyLevel.almostQuantum ... 
                 qvx.di.HierarchyLevel.NPAplus(1) qvx.di.HierarchyLevel.almostQuantum ... 
                 qvx.di.HierarchyLevel.NPAplus(1) qvx.di.HierarchyLevel.almostQuantum ...
                 qvx.di.HierarchyLevel.NPA(2) qvx.di.HierarchyLevel.NPA(3)};
        maximum = {2*sqrt(2) ...
                   3.1547 2.9149 ...
                   3.2126 2.9727 ...
                   3.2997 3.0157 ...
                   0.36603 0.25147 ...
                   0.25094 0.25088};
    end
    methods (Test, ParameterCombination='sequential')
        function testSize(testCase, scenario, expression, level, ...
                          maximum)
            cvx_begin
            variable PG(length(expression)) % Correlation vector in Collins-Gisin form
            PG == qvx.di.cvx.QuantumSupersetCG(scenario, level);
            PG(1) == 1 % normalization, as the moment matrix is not normalized by default
            maximize expression * PG
            cvx_end
            Iopt = expression * PG;
            tol = 1e-4;
            assert(abs(Iopt - maximum) < tol);
        end
    end
end

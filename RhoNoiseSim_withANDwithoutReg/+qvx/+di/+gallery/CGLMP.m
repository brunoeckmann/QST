classdef CGLMP
    methods (Static)
        function S = scenario(d)
            S = qvx.di.Scenario.homogenous(2,2,d);
        end
        function I = expression(d, space)
            if nargin < 2
                space = 'P';
            end
            I = zeros(d, 2, d, 2);
            for k = 0:floor(d/2)-1
                coeff = 1 - 2*k/(d-1);
                for a = 1:d
                    for b = 1:d
                        cond = @(n) mod(4*d + n, d) == 0;
                        if cond(a-b-k)
                            I(a,1,b,1) = I(a,1,b,1) + coeff;
                        end
                        if cond(b-a-k-1)
                            I(a,2,b,1) = I(a,2,b,1) + coeff;
                        end
                        if cond(a-b-k)
                            I(a,2,b,2) = I(a,2,b,2) + coeff;
                        end
                        if cond(b-a-k)
                            I(a,1,b,2) = I(a,1,b,2) + coeff;
                        end
                        if cond(a-b+k+1)
                            I(a,1,b,1) = I(a,1,b,1) - coeff;
                        end
                        if cond(b-a+k)
                            I(a,2,b,1) = I(a,2,b,1) - coeff;
                        end
                        if cond(a-b+k+1)
                            I(a,2,b,2) = I(a,2,b,2) - coeff;
                        end
                        if cond(b-a+k+1)
                            I(a,1,b,2) = I(a,1,b,2) - coeff;
                        end
                    end
                end
            end
            I = permute(I, [3 4 1 2]);
            I = I(:);
            I = qvx.di.Expressions.convert(qvx.di.gallery.CGLMP.scenario(d), space, 'P', I);
        end
    end
end

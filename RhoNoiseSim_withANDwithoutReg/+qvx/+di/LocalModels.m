classdef LocalModels
% LocalModels - methods for handling local models
    methods (Static)
        function [M D] = localStrategiesForParty(party, space)
        % Returns a matrix with the deterministic strategies
        %
        % space - optional argument specifying the correlation
        %         space, by default 'P'
        %
        % Outputs a matrix with one deterministic strategy vector
        % per column
        %
        % (Note: the output obeys the rational convention) 
            if nargin < 2
                space = 'P';
            end
            switch space
              case 'P'
                n = length(party.inputs);
                d = sum(party.inputs);
                p = prod(party.inputs);
                D = 1;
                M = zeros(d, p);
                for c = 1:p
                    r = 1;
                    t = c - 1;
                    for j = 1:n
                        o = party.inputs(j);
                        alpha = mod(t, o);
                        t = floor(t / o);
                        M(r + alpha, c) = 1;
                        r = r + o;
                    end
                end
              otherwise
                [M1 D1] = qvx.di.LocalModels.localStrategiesForParty(party, 'P');
                [M2 D2] = party.conversionMatrix(space, 'P');
                M = M2 * M1;
                D = D2 * D1;
            end
            if nargout < 2
                M = M/D;
            end
        end
        function [M D] = localStrategies(scenario, space)
            import qvx.di.LocalModels
            if nargin < 2
                space = 'P';
            end
            M = 1;
            D = 1;
            for i = 1:length(scenario.parties)
                [Mp Dp] = LocalModels.localStrategiesForParty(scenario.parties{i}, space);
                M = kron(Mp, M);
                D = Dp * D;
            end
            if nargout < 2
                M = M/D;
            end
        end        
    end
end

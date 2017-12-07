classdef Correlations
% Correlations Collection of methods that work on correlation vectors
% 
% Correlations are column vectors that represent probabilities in a correlation
% space (see also qvx.di.Spaces)
%
% Use the 'describeOrder' method to get an explicit enumeration of the coefficients.
    methods (Static)
        
        % Generate points
        function prob = uniformlyRandom(scenario, space)
            import qvx.di.Correlations
            if nargin < 2
                space = 'P';
            end
            prob = 1;
            for i = 1:length(scenario.parties)
                prob = kron(Correlations.uniformlyRandomForParty(scenario.parties{i}, 'P'), prob);
            end
            prob = Correlations.convert(scenario, space, 'P', prob);
        end

        function prob = randomLocal(scenario, space)
        % Correlations.randomLocal Generates random local correlations
            M = qvx.di.LocalModels.localStrategies(scenario, space);
            q = rand(size(M, 2), 1);
            q = q / sum(q);
            prob = M*q;
        end
        
        function prob = randomNonsignaling(scenario, space)
        % Correlations.randomNonsignaling Generates random nonsignaling correlations
            import qvx.di.Correlations;
            dNG = Correlations.dimension(scenario, 'NG');
            NG = rand(dNG, 1);
            NG(1) = 1;
            NG0 = Correlations.uniformlyRandom(scenario, 'NG');
            P = Correlations.convert(scenario, 'P', 'NG', NG);
            % v P + (1-v) P0 = 0
            while min(P) < 0
                NG = (NG + NG0)/2;
                P = Correlations.convert(scenario, 'P', 'NG', NG);
            end
            prob = Correlations.convert(scenario, space, 'NG', NG);
        end
        
        function prob = randomSignaling(scenario, space)
        % Correlations.randomSignaling Generates random signaling correlations
            import qvx.di.Correlations;
            assert(~isequal(space, 'NG') && ~isequal(space, 'NC'));
            d = Correlations.dimension(scenario, 'SC');
            C = rand(d, 1);
            C(scenario.balancedNormalizationMask) = 0; % balance normalization
            C(1) = 1; % fix overall normalization to 1
            C0 = Correlations.uniformlyRandom(scenario, 'SC');
            P = Correlations.convert(scenario, 'P', 'SC', C);
            % v P + (1-v) P0 = 0
            while min(P) < 0
                C = (C + C0)/2;
                P = Correlations.convert(scenario, 'P', 'SC', C);
            end
            prob = Correlations.convert(scenario, space, 'SC', C);            
        end
        
        % Methods for correlation vectors
        function d = dimensionForParty(party, space)
            switch space
              case {'P', 'SG', 'SC'}
                d = sum(party.inputs);
              case {'NG', 'NC'}

                d = sum(party.inputs) - length(party.inputs) + 1;
            end
        end
        function d = dimension(scenario, space)
            d = prod(cellfun(@(x) qvx.di.Correlations.dimensionForParty(x, space), scenario.parties));
        end
        function prob = uniformlyRandomForParty(party, space)
            if nargin < 2
                space = 'P';
            end
            prob = [];
            for i = 1:length(party.inputs)
                prob = [prob
                        ones(party.inputs(i),1)/party.inputs(i)];
            end
            prob = party.conversionMatrix(space, 'P') * prob;
        end
        function P = fromPabxy(Pabxy)
        % Converts a (2*n)-dim. tensor P(a,b,..,x,y,..) to a vector of correlations in P space
            s = size(Pabxy);
            n = ndims(Pabxy)/2;
            order = [fliplr(1:n)
                     fliplr(n+1:2*n)];
            reord = permute(Pabxy, order(:)');
            P = reord(:);
        end
        function Pabxy = toPabxy(scenario, P)
        % Converts a vector of correlations in P space to a (2*n)-dim. tensor P(a,b,..,x,y,..)
            n = scenario.nParties;
            assert(scenario.hasHomogenousParties);
            outs = cellfun(@(p) p.inputs(1), scenario.parties);
            outs = outs(:)';
            ins = cellfun(@(p) p.nInputs, scenario.parties);
            ins = ins(:)';
            s = [fliplr(outs)
                 fliplr(ins)];
            s = s(:)';
            P1 = reshape(P, s);
            Pabxy = permute(P1, [fliplr(1:2:2*n-1) fliplr(2:2:2*n)]);
        end
        function res = convert(scenario, toSpace, fromSpace, correlation)
        % Converts the 'correlation' vector in the space 'fromSpace' to a vector
        % in the space 'toSpace' for the given 'scenario'
            if isequal(fromSpace, 'Pabxy')
                fromSpace = 'P';
                correlation = qvx.di.Correlations.fromPabxy(correlation);
            end
            if isequal(toSpace, 'Pabxy')
                toSpace = 'P';
                toPabxy = true;
            else
                toPabxy = false;
            end
            correlation = correlation(:);
            [M D] = scenario.conversionMatrix(toSpace, fromSpace);
            res = M * correlation / D;
            if toPabxy
                res = qvx.di.Correlations.toPabxy(scenario, res);
            end
        end
        function [Aeq beq] = nonsignalingLinearConstraintsP(scenario)
        % Returns the linear constraints obeyed by (even subnormalized) nonsignaling correlations
            [M D] = scenario.conversionMatrix('SG', 'P');
            mask = scenario.nonsignalingMask;
            Aeq = qvx.util.matrowgcd(M(~mask, :));
            if nargout > 1
                beq = zeros(size(Aeq, 1), 1);
            end
        end
        function [Aeq beq] = overallNormalizationConstraintP(scenario)
        % Returns the general normalization constraint
            [M D] = scenario.conversionMatrix('SG', 'P');
            Aeq = M(1,:);
            beq = D;
            [Aeq beq] = qvx.util.matvecrowgcd(Aeq, beq);
        end
        function [Aeq beq] = nonsignalingConstraintsP(scenario)
        % Returns all the equality constraints obeyed by normalized, nonsignaling correlations
            [A1 b1] = scenario.overallNormalizationConstraintP;
            [A2 b2] = scenario.nonsignalingLinearConstraintsP;
            Aeq = [A1
                   A2];
            beq = [b1
                   b2];
        end
        function [M D] = nonsignalingProjectionMatrix(scenario, space)
        % Returns the matrix that projects into the nonsignaling subspace
        %
        % Optionally, a space ('SG', 'SC') can be provided instead of the default 'P'
            if nargin < 2
                space = 'P';
            end
            switch space
              case 'P'
                [M1 D1] = scenario.conversionMatrix('P', 'NG');
                [M2 D2] = scenario.conversionMatrix('NG', 'P');
              case 'SG'
                [M1 D1] = scenario.conversionMatrix('SG', 'NG');
                [M2 D2] = scenario.conversionMatrix('NG', 'SG');
              case 'SC'
                [M1 D1] = scenario.conversionMatrix('SC', 'NC');
                [M2 D2] = scenario.conversionMatrix('NC', 'SC');
            end
            M = M1 * M2;
            D = D1 * D2;
            if nargout < 2
                M = M / D;
            end
        end
        function proj = nonsignalingProjection(scenario, space, coeffs)
        % nonsignalingProjection Projects a correlation vector in the nonsignaling subspace
        % 
        % Correlations.nonsignalingProjection(scenario, space, x) where:
        % - scenario is the scenario considered
        % - space is the correlation space type (P, SG or SC)
        % - coeffs is a correlation vector
            import qvx.di.Correlations;
            coeffs = coeffs(:);
            switch space
              case {'P', 'Pabxy'}
                if isequal(space, 'Pabxy')
                    coeffs = Correlations.fromPabxy(coeffs);
                end
                ng = Correlations.convert(scenario, 'NG', 'P', coeffs);
                proj = Correlations.convert(scenario, 'P', 'NG', ng);
                if isequal(space, 'Pabxy')
                    proj = Correlations.toPabxy(scenario, proj);
                end
              case 'SG'
                ng = Correlations.convert(scenario, 'NG', 'SG', coeffs);
                proj = Correlations.convert(scenario, 'SG', 'NG', ng);
              case 'SC'
                ng = Correlations.convert(scenario, 'NC', 'SC', coeffs);
                proj = Correlations.convert(scenario, 'SC', 'NC', ng);
            end
        end
        function strings = describeOrder(scenario, space)
        % describe Returns a list of strings describing the ordering of coefficients
            [O I] = scenario.outputInputIndices(space);
            strings = cell(size(O,1), 1);
            switch space
              case 'SC'
                error('Not implemented')
              case 'SG'
                error('Not implemented')
              case 'NC'
                for i = 1:size(O,1)
                    str = '<';
                    for j = 1:size(O,2)
                        switch O(i,j)
                          case 0
                          case 1
                            str = [str char('A'+j-1) num2str(I(i,j))];
                          case 2
                            error('Not implemented for nonbinary correlators')
                        end
                    end
                    str = [str '>'];
                    strings{i} = str;
                end
              otherwise
                for i = 1:size(O,1)
                    str = 'P';
                    if any(O(i,:) == 0)
                        for j = 1:size(O,2)
                            if O(i,j) ~= 0
                                str = [str char('A'+j-1)];
                            end
                        end
                    end
                    str = [str '('];
                    sep = '';
                    for j = 1:size(O,2)
                        if O(i,j) ~= 0
                            str = [str sep num2str(O(i,j))];
                            sep = ',';
                        end
                    end
                    str = [str '|'];
                    sep = '';
                    for j = 1:size(I,2)
                        if I(i,j) ~= 0
                            str = [str sep num2str(I(i,j))];
                            sep = ',';
                        end
                    end
                    str = [str ')'];
                    strings{i} = str;
                end
            end
        end
    end
end

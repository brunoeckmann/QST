classdef (Sealed) Scenario
% Scenario Describes a Bell scenario
    properties
        parties = {};
    end
    methods
        % Constructors
        function S = Scenario(parties)
            S.parties = parties;
        end
        % General methods
        function b = hasHomogenousParties(S)
            b = all(cellfun(@(p) p.isHomogenous, S.parties));
        end
        function b = isHomogenous(S)
            if ~all(cellfun(@(p) p.isHomogenous, S.parties))
                b = false;
            else
                b = true;
                for i = 2:length(S.parties)
                    if ~isequal(S.parties{1}.inputs, S.parties{i}.inputs)
                        b = false;
                    end
                end
            end
        end
        function n = nParties(S)
        % scenario.nParties Returns the number of parties.
            n = length(S.parties);
        end
        function str = toString(S)
            str = '{';
            sep = '';
            for i = 1:length(S.parties)
                str = [str sep S.parties{i}.toString];
                sep = ' ';
            end
            str = [str '}'];
        end
        function str = disp(S)
            disp(S.toString);
        end
        % Functions returning matrices
        function mask = productMask(S, normC, corrC, sigC)
            mask = 1;
            for i = 1:length(S.parties)
                mask = kron(S.parties{i}.mask(normC, corrC, sigC), mask);
            end
        end
        function mask = nonsignalingMask(S)
            mask = S.productMask(1, 1, 0) == 1;
        end
        function mask = balancedNormalizationMask(S)
            mask = S.productMask(1, 0, 2) > 1;
        end
        function [O I] = outputInputIndices(S, space)
            [O I] = S.parties{1}.outputInputIndices(space);
            for i = 2:length(S.parties)
                [Op Ip] = S.parties{i}.outputInputIndices(space);
                newO = [];
                newI = [];
                nO = size(O,1);
                for k = 1:length(Op)
                    newO = [newO
                            Op(k)*ones(nO,1) O];
                    newI = [newI
                            Ip(k)*ones(nO,1) I];
                end
                O = newO;
                I = newI;
            end
        end
        function [M D] = conversionMatrix(S, to, from)
            M = 1;
            D = 1;
            for i = 1:length(S.parties)
                [Mp Dp] = S.parties{i}.conversionMatrix(to, from);
                M = kron(Mp, M);
                D = D * Dp;
            end
            if nargout < 2
                M = M / D;
            end
        end
    end
    methods (Static)
        function S = homogenous(nParties, nInputs, nOutputs)
            party = qvx.di.Party.homogenous(nInputs, nOutputs);
            parties = cell(1, nParties);
            parties(1,:) = {party};
            S = qvx.di.Scenario(parties);
        end
        function S = fromPabxy(Pabxy)
        % Scenario.fromPabxy(Pabxy) Constructs a scenario from a tensor shape
        %
        % For example, Scenario.fromPabxy(ones(2,2,3,3)/4) returns the scenario of I3322
            s = size(Pabxy);
            n = ndims(Pabxy)/2;
            outs = s(1:n);
            ins = s(n+1:2*n);
            parties = cell(1, n);
            for i = 1:n
                parties{i} = qvx.di.Party.homogenous(ins(i), outs(i));
            end
            S = qvx.di.Scenario(parties);
        end
    end
end

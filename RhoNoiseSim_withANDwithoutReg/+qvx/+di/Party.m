classdef (Sealed) Party
% Party Describes a party in a Bell scenario
    properties
        inputs = []; % 1xn vector representing the number of outputs per input
    end
    methods
        function P = Party(inputs)
        % Constructs a party from a 1xn vector of positive integers
            P.inputs = inputs;
        end
        function n = nInputs(P)
        % party.nInputs Returns the number of inputs
            n = length(P.inputs);
        end
        function b = isHomogenous(P)
        % isHomogenous Returns whether the party is homogenous
            b = all(P.inputs == P.inputs(1));
        end
        function str = toString(P)
        % toString Returns a string representation of the party
            str = '[';
            sep = '';
            for i = 1:length(P.inputs)
                str = [str sep num2str(P.inputs(i))];
                sep = ' ';
            end
            str = [str ']'];
        end
        function disp(P)
            disp(P.toString);
        end
        function [O I] = outputInputIndices(P, space)
            switch space
              case 'SC'
                error('Not implemented for signaling correlators')
              case 'SG'
                error('Not implemented for signaling Collins-Gisin')
              case 'P'
                O = [];
                I = [];
                for i = 1:length(P.inputs)
                    O = [O
                         (1:P.inputs(i))'];
                    I = [I
                         i*ones(P.inputs(i), 1)];
                end
              case 'NG'
                O = [0];
                I = [0];
                for i = 1:length(P.inputs)
                    O = [O
                         (1:P.inputs(i)-1)'];
                    I = [I
                         i*ones(P.inputs(i)-1, 1)];
                end
              case 'NC'
                [O I] = P.outputInputIndices('NG');
            end
        end
        function mask = mask(P, normC, corrC, sigC)
            d = sum(P.inputs);
            n = length(P.inputs);
            c = d - n;
            mask = [normC corrC*ones(1, c) sigC*ones(1, n-1)];
        end
        function [M D] = conversionMatrix(P, to, from)
        % matrix Returns the conversion matrix between correlation spaces
        % 
        % M = party.matrix(to, from) returns the conversion such
        %     that y = M * x has y in the 'to' space and x in the
        %     'from' space
        %
        % Correlation spaces are given as strings (see qvx.di.Spaces for help)
        % 
        % Example: to project x, written as a probability version,
        % in the nonsignaling subspace, compute: 
        %
        % p.conversionMatrix('P', 'NG') * p.conversionMatrix('NG','P')
        %
        % If only a single return argument is asked, the function
        % returns the requested double matrix, however in
        % floating-point precision.
        % 
        % If two return arguments are asked, the function returns
        % the matrix as a numerator/denominator pair, where the
        % numerator is an integer matrix and the denominator an
        % integer, such that num/den is the requested matrix.
            select = [to '_' from];
            switch select
              case 'NC_P'
                [M1 D1] = P.conversionMatrix('NC', 'SC');
                [M2 D2] = P.conversionMatrix('SC', 'P');
                M = M1 * M2;
                D = D1 * D2;
              case 'P_NC'
                [M1 D1] = P.conversionMatrix('P', 'SC');
                [M2 D2] = P.conversionMatrix('SC', 'NC');
                M = M1 * M2;
                D = D1 * D2;
              case 'SC_P'
                d = sum(P.inputs);
                n = length(P.inputs);
                D = n;
                M = zeros(d, d);
                M(1, :) = 1;
                r = 2;
                c = 1;
                for i = 1:n
                    o = P.inputs(i);
                    [q qD] = qvx.di.QMatrix(o).matrix;
                    assert(qD == 1);
                    M(r:r+o-2,c:c+o-1) = q(2:end,:) * n;
                    r = r + o - 1;
                    c = c + o;
                end
                [q qD] = qvx.di.QMatrix(n).matrix;
                assert(qD == 1);
                for i = 2:n
                    c = 1;
                    for j = 1:n
                        o = P.inputs(j);
                        M(r, c:c+o-1) = q(i,j) * n;
                        c = c + o;
                    end
                    r = r + 1;
                end
              case 'P_SC'
                d = sum(P.inputs);
                n = length(P.inputs);
                [qi qiD] = qvx.di.QMatrix(n).matrixInverse;
                l = 1;
                for i = 1:n
                    o = P.inputs(i);
                    l = lcm(l, o*qiD);
                    [qo qoD] = qvx.di.QMatrix(o).matrixInverse;
                    l = lcm(l, qoD);
                end
                M = zeros(d, d);
                r = 1;
                c = 2;
                for i = 1:n
                    o = P.inputs(i);
                    [qo qoD] = qvx.di.QMatrix(o).matrixInverse;
                    M(r:r+o-1, 1) = l/o;
                    M(r:r+o-1, c:c+o-2) = qo(:,2:end) * (l/qoD);
                    r = r + o;
                    c = c + o - 1;
                end
                for j = 2:n
                    r = 1;
                    for i = 1:n
                        o = P.inputs(i);
                        M(r:r+o-1, c) = qi(i, j) * (l/o/qiD);
                        r = r + o;
                    end
                    c = c + 1;
                end
                D = l;
              case 'NG_P'
                [M1 D1] = P.conversionMatrix('NG', 'SG');
                [M2 D2] = P.conversionMatrix('SG', 'P');
                M = M1 * M2;
                D = D1 * D2;
              case 'P_NG'
                [M1 D1] = P.conversionMatrix('P', 'SG');
                [M2 D2] = P.conversionMatrix('SG', 'NG');
                M = M1 * M2;
                D = D1 * D2;
              case 'SC_NC'
                [M D] = P.conversionMatrix('SG', 'NG');
              case 'NC_SC'
                [M D] = P.conversionMatrix('NG', 'SG');
              case 'SG_NG'
                [M D] = P.conversionMatrix('NG', 'SG');
                M = M';
              case 'NG_SG'
                d = sum(P.inputs);
                n = length(P.inputs);
                c = d - n + 1;
                M = zeros(c, d);
                M(1:c, 1:c) = eye(c);
                D = 1;
              case 'NG_NG'
                d = sum(P.inputs);
                n = length(P.inputs);
                c = d - n + 1;
                M = eye(c, c);
                D = 1;
              case 'NC_NC'
                [M D] = P.conversionMatrix('NG', 'NG');
              case 'P_P'
                M = eye(sum(P.inputs));
                D = 1;
              case 'SC_SC'
                [M D] = P.conversionMatrix('P', 'P');
              case 'SG_SG'
                [M D] = P.conversionMatrix('P', 'P');
              case 'P_SG'
                n = length(P.inputs);
                [Mqinv Dqinv] = qvx.di.QMatrix(n).matrixInverse;
                l = P.inputs(1);
                for i = 2:n
                    l = lcm(l, P.inputs(i));
                end
                D = l * Dqinv;
                d = sum(P.inputs);
                M = zeros(d, d);
                N = zeros(d, d);
                c = 2;
                r = 1;
                for i = 1:n
                    o = P.inputs(i);
                    M(r + o - 1, 1) = D;
                    M(r + o - 1, c:c+o-2) = -D;
                    M(r:r+o-2, c:c+o-2) = eye(o - 1, o - 1) * D;
                    r = r + o;
                    c = c + o - 1;
                end
                for j = 2:n
                    r = 1;
                    for i = 1:n
                        o = P.inputs(i);
                        M(r:r+o-1, c) = Mqinv(i, j)*(l/o);
                        r = r + o;
                    end
                    c = c + 1;
                end
              case 'SG_P'
                n = length(P.inputs);
                d = sum(P.inputs);
                l = P.inputs(1);
                for i = 2:n
                    l = lcm(l, P.inputs(i));
                end
                M = zeros(d, d);
                M(1, :) = l;
                D = l*n;
                c = 1;
                r = 2;
                for i = 1:n
                    o = P.inputs(i);
                    a = -(n - 1) * (l / o); % -(n-1)/(o * n)
                    b = l / o; % 1/(o * n)
                    M(r:r+o-2,:) = b;
                    for dr = 0:o-2
                        for dc = 0:o-1
                            if dr == dc
                                M(r+dr, c+dc) = a + D;
                            else
                                M(r+dr, c+dc) = a;
                            end
                        end
                    end
                    r = r + o - 1;
                    c = c + o;
                end
                [Mq Dq] = qvx.di.QMatrix(n).matrix;
                assert(Dq == 1);
                for i = 2:n
                    c = 1;
                    for j = 1:n
                        o = P.inputs(j);
                        M(r, c:c+o-1) = Mq(i,j) * D;
                        c = c + o;
                    end
                    r = r + 1;
                end
              otherwise
                assert(~isequal(to, 'P'));
                assert(~isequal(from, 'P'));
                [M1 D1] = P.conversionMatrix(to, 'P');
                [M2 D2] = P.conversionMatrix('P', from);
                M = M1 * M2;
                D = D1 * D2;
            end
            if nargout < 2
                M = M/D;
            end
        end
    end
    methods (Static)
        function P = homogenous(nInputs, nOutputs)
            P = qvx.di.Party(ones(1, nInputs) * nOutputs);
        end
    end
end

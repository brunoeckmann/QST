classdef CGMM
% CGMM Implementation of the NPA algebra for Collins-Gisin projectors
%
% Data structure:
% Products of Collins-Gisin operators are stored in a cell array
% one cell per party
% In each cell, a n x 2 integer array, where the first column contains
% the input indices and the second column contains the output indices
    methods (Static)
        function res = reduceParty(x)
        % takes a n x 2 integer array as input and reduces the product of operators
        % according to the algebraic rules
        % if the expression contains neighboring projectors with different outputs
        % for the same input, returns [], otherwise returns a n' x 2 integer array
            res = [];
            n = size(x, 1);
            if n == 0
                res = x;
                return
            end
            lastI = x(1, 1);
            lastO = x(1, 2);
            for k = 2:size(x, 1)
                i = x(k, 1);
                o = x(k, 2);
                if i == lastI
                    if o ~= lastO
                        res = [];
                        return
                    else
                        % do nothing, we compact both operators
                    end
                else
                    if lastI ~= 0
                        res = [res
                               lastI lastO];
                    end
                    lastI = i;
                    lastO = o;
                end
            end
            res = [res
                   lastI lastO];
        end
        function res = transpose(x)
        % transposes a product of CG operators
        % input: a 1 x n cell array where n is the number of partiesx
            n = length(x);
            res = cell(1, n);
            for i = 1:n
                res{i} = flipud(x{i});
            end
        end
        function res = product(x, y)
        % computes the product of x and y, where x and y are 1 x n cell
        % arrays representing products of CG operators
            n = length(x);
            assert(n == length(y));
            res = cell(1, n);
            for i = 1:n
                res{i} = qvx.di.internal.CGMM.reduceParty([x{i}; y{i}]);
                if isequal(res{i}, [])
                    res = [];
                    return
                end
            end
        end
        function ops = opsForParty(party, level)
            ops = {zeros(0,2)};
            for l = 1:level
                newOps = {zeros(0,2)};
                for k = 1:length(ops)
                    op = ops{k};
                    if size(op, 1) ~= 0
                        lastI = op(end, 1);
                    else
                        lastI = 0;
                    end
                    for i = 1:length(party.inputs)
                        if i ~= lastI
                            for o = 1:party.inputs(i)-1
                                newOps{1, end+1} = [op
                                                    i o];
                            end
                        end
                    end
                end
                ops = newOps;
            end
        end
        function res = localLevelAccumulate(list, parties, level)
            import qvx.di.internal.CGMM;
            if length(parties) == 0
                res = list;
                return
            end
            ops = CGMM.opsForParty(parties{1}, level);
            newList = {};
            for i = 1:length(list)
                prev = list{i};
                for j = 1:length(ops)
                    newList{1, end+1} = {prev{:} ops{j}};
                end
            end
            res = CGMM.localLevelAccumulate(newList, parties(2:end), level);
        end
        function str = text(x)
        % returns the text representation of the provided product of CG operators
        % input: a 1 x n cell array
            n = length(x);
            str = '';
            for i = 1:n
                if size(x{i}, 1) ~= 0
                    io = x{i};
                    for r = 1:size(io, 1)
                        str = [str sprintf('%c%d_%d', 'A'+i-1, io(r,1), io(r,2))];
                    end
                end
            end
        end
        function l = opLength(x)
            l = sum(cellfun(@(pty) size(pty, 1), x));
        end
        function MM = buildFromVector(scenario, vector, addNonNegativity)
            import qvx.di.MomentMatrix;
            import qvx.di.internal.CGMM;
            d = length(vector);
            MB = qvx.di.internal.MomentMatrixBuilder(d);
            n = length(scenario.parties);
            [O I] = scenario.outputInputIndices('NC');
            outputDim = size(O, 1);
            % register the moments for the output space
            for i = 1:outputDim
                moment = cell(1, n);
                for j = 1:n
                    if O(i,j) == 0
                        moment{j} = zeros(0, 2);
                    else
                        moment{j} = [I(i,j) O(i,j)];
                    end
                end
                assert(MB.getMoment(CGMM.text(moment)) == i);
            end
            % construct the index matrix
            for r = 1:d
                for c = r:d
                    res = CGMM.product(CGMM.transpose(vector{r}), vector{c});
                    if ~isequal(res, [])
                        resT = CGMM.transpose(res);
                        res = CGMM.text(res);
                        resT = CGMM.text(resT);
                        MB.setIndexMatrix(r, c, res, resT);
                    end
                end
            end
            if addNonNegativity
                [M D] = scenario.conversionMatrix('P', 'NG');
                nFreeMoments = MB.lastMoment - size(M, 2);
                iA = -[M zeros(size(M, 1), nFreeMoments)];
                ib = zeros(size(M, 1), 1);
                for r = 1:size(iA, 1)
                    g = iA(r, 1);
                    for c = 2:size(iA, 2)
                        g = gcd(g, iA(r,c));
                    end
                    if g ~= 1
                        iA(r, :) = iA(r, :)/g;
                    end
                end
            else
                iA = [];
                ib = [];
            end
            MM = MomentMatrix(MB.indexMatrix, MB.moments, MB.momentMap, outputDim, iA, ib);
        end
        function MM = buildLocalLevel(scenario, level)
            import qvx.di.internal.CGMM;
            vector = CGMM.localLevelAccumulate({{}}, scenario.parties, level);
            MM = CGMM.buildFromVector(scenario, vector, false);
        end
        function MM = buildNPALevel(scenario, level, addNonNegativity)
        % Constructs the NPA variant
        % addNonNegativity specifies if nonnegativity constraints are added
        % on top, regardless on the provided level
            import qvx.di.internal.CGMM;
            vector = CGMM.localLevelAccumulate({{}}, scenario.parties, level);
            lens = cellfun(@(x) CGMM.opLength(x), vector);
            vector = vector(lens <= level);
            MM = CGMM.buildFromVector(scenario, vector, addNonNegativity);
        end
        function MM = build(scenario, hierarchyLevel)
            switch hierarchyLevel.hierarchyType
              case 'local'
                MM = qvx.di.internal.CGMM.buildLocalLevel(scenario, hierarchyLevel.level);
              case 'NPA'
                MM = qvx.di.internal.CGMM.buildNPALevel(scenario, hierarchyLevel.level, hierarchyLevel.addNonNegativity);
            end
        end
    end
end

classdef MomentMatrix
% MomentMatrix Describes a convex set represented using moment matrix constraints
    properties
        outputDim = 0; % dimension of the represented convex set 
        
        indexMatrix = []; % indices of the moments appearing in the moment matrix
                          % indexMatrix is a m x m integer matrix
        moments = {}; % 1 x m cell array containing the moment string representations
        momentMap = []; % map from string to integer, mapping the string representation of a
                        % moment, or its transpose, to the relevant index
        iA = []; % optional linear inequality constraints
        ib = []; % such that iA * momentsValues(:) <= ib
    end
    methods
        function MM = MomentMatrix(indexMatrix, moments, momentMap, outputDim, iA, ib)
        % should not be called by the user
            MM.indexMatrix = indexMatrix;
            MM.moments = moments;
            MM.momentMap = momentMap;
            MM.outputDim = outputDim;
            MM.iA = iA;
            MM.ib = ib;
        end
        function mat = matrixForMoment(M, momentIndex)
            mat = (M.indexMatrix == momentIndex);
        end
        function n = nMoments(M)
            n = length(M.moments);
        end
        function n = nFreeMoments(M)
            n = M.nMoments - M.outputDim;
        end
    end
    methods (Static)
        function MM = quantumCollinsGisin(scenario, hierarchyLevel)
            MM = qvx.di.internal.CGMM.build(scenario, hierarchyLevel);
        end
    end
end

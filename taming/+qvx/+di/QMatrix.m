classdef (Sealed) QMatrix
    properties
        n = [];
    end
    methods
        function Q = QMatrix(n);
            Q.n = n;
        end
        function [M D] = matrix(Q)
            n = Q.n;
            M = zeros(n, n);
            M(1, :) = 1;
            for r = 2:n
                M(r, r-1:r) = [1 -1];
            end
            D = 1;
        end
        function [M D] = matrixInverse(Q)
            n = Q.n;
            M = zeros(n, n);
            D = n;
            M(:,1) = 1;
            for c = 2:n
                M(1:c-1,c) = (n-c+1);
                M(c:n,c) = (-c+1);
            end
            if nargout == 1
                M = M / D;
            end
        end
    end
end

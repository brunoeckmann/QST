classdef MomentMatrixBuilder < handle
    properties
        indexMatrix;
        moments;
        momentMap;
        lastMoment;
    end
    methods
        function M = MomentMatrixBuilder(d)
            M.indexMatrix = zeros(d, d);
            M.moments = {};
            M.momentMap = containers.Map('KeyType', 'char', 'ValueType', 'double');
            M.lastMoment = 0;
        end
        function setIndexMatrix(M, r, c, moment, varargin)
            moment = M.getMoment(moment, varargin{:});
            M.indexMatrix(r, c) = moment;
            M.indexMatrix(c, r) = moment;
        end
        function index = getMoment(M, moment, varargin)
            if M.momentMap.isKey(moment)
                index = M.momentMap(moment);
            else
                M.lastMoment = M.lastMoment + 1;
                M.momentMap(moment) = M.lastMoment;
                for i = 1:length(varargin)
                    M.momentMap(varargin{i}) = M.lastMoment;
                end
                index = M.lastMoment;
                M.moments{1, end+1} = moment;
            end
        end
    end
end

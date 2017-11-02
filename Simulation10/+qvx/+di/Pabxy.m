% Pabxy - Probability distributions represented by (2*n)-dimensional tensors
%
% When a Pabxy object is accepted or returned by a function, it corresponds to
% a (2*n)-dimension tensor, where n is the number of parties. The first n indices
% correspond to the outputs of the joint distribution, and the next n indices correspond
% to the inputs.
%
% For example, a bipartite distribution P(ab|xy) is indexed as P(a,b,x,y).
%
% This representation can be used only when the parties are homogenous (same number
% of outputs for each input).

% Spaces - Correlation spaces
%
% The following strings can be used to specify a space.
%
% Below, we describe, for each space, the enumeration of
% coefficients that applies to a single party. For scenarios
% involving multiple parties, the enumeration of coefficients
% is such that a product distribution C_AB = kron(C_A, C_B).
% 
% 
% 'P'  - full probabilities P(a|x), enumerated by increasing
%      first 'a' and then 'x', such that all coefficients
%      corresponding to the same input are grouped together
%
%      Some functions also accept the 'Pabxy' variant, where
%      the vector is provided as a (2*n)-dimension tensor
%      where n is the number of parties (see Pabxy.m)
%
% 'NG' - nonsignaling Collins-Gisin, where the first
%      coefficient is the normalization 1, and then the
%      next coefficients correspond to the coefficients of
%      P(a|x), where the last output (with 'a' maximal) is
%      omitted for each input. Coefficients are ordered
%      by first increasing 'a' and then 'x'
%
% 'SG' - signaling Collins-Gisin, where the first
%      coefficients are exactly the same as in the 'NG'
%      form, and the additional coefficients are always
%      zero when the distributions are properly normalized
%      and nonsignaling
%
% 'NC' - nonsignaling correlators, where the first coefficient
%      is the normalization 1, and the subsequent coefficients 
%      correspond (for binary outputs) to the usual binary 
%      correlators, one coefficient per input.
%      When the number of outputs is not binary, we consider
%      instead the generalization described in
%      Rosset 2014, J. Phys. A., 'Classifying 50 years of
%      Bell inequalities'.
%
% 'SC' - signaling correlators, where the first coefficients
%      correspond to the nonsignaling correlators (see NC), and
%      then coefficients corresponding to the signaling subspace.

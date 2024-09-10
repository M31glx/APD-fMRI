function [T_pos,T_neg] = transitivity_wu_sign(W)
%TRANSITIVITY_WU_SIGN    Transitivity
%
%   [T_pos,T_neg] = transitivity_wu_sign(W);
%
%   Transitivity is the ratio of 'triangles to triplets' in the network.
%   (A classical version of the clustering coefficient).
%
%   Input:      W       weighted undirected connection matrix
%
%   Output:     T       transitivity scalar
%
%   Note:      All weights must be between 0 and 1.
%              This may be achieved using the weight_conversion.m function,
%              W_nrm = weight_conversion(W, 'normalize');
%
%   Reference: Rubinov M, Sporns O (2010) NeuroImage 52:1059-69
%              based on Onnela et al. (2005) Phys Rev E 71:065103
%
%
%   Mika Rubinov, UNSW/U Cambridge, 2010-2015

%   Modification history:
%   2010: Original
%   2015: Expanded documentation
%   2015: Now computed for positive and negative weights separately and 
%         automatically enforces that weights must be between 0 and 1 (JMS)

n = length(W);                                    %number of nodes

if max(abs(W(:)))>1
    W = W/max(abs(W(:)));
end

W_pos            = threshold_absolute(W,0);
W_pos(1:n+1:end) = 0;
W_neg            = threshold_absolute(W*-1,0);
W_neg(1:n+1:end) = 0;

K_pos    = sum(W_pos~=0,2);
K_neg    = sum(W_neg~=0,2);

cyc3_pos = diag((W_pos.^(1/3))^3);
cyc3_neg = diag((W_neg.^(1/3))^3);

T_pos    = sum(cyc3_pos)./sum((K_pos.*(K_pos-1)));       %transitivity
T_neg    = sum(cyc3_neg)./sum((K_neg.*(K_neg-1)));       %transitivity

function [v_pos,v_neg] = eigenvector_centrality_und_sign(CIJ)
%EIGENVECTOR_CENTRALITY_UND_SIGN   Spectral measure of centrality
%
%   [v_pos,v_neg] = eigenvector_centrality_und_sign(CIJ)
%
%   Eigenector centrality is a self-referential measure of centrality:
%   nodes have high eigenvector centrality if they connect to other nodes
%   that have high eigenvector centrality. The eigenvector centrality of
%   node i is equivalent to the ith element in the eigenvector 
%   corresponding to the largest eigenvalue of the adjacency matrix.
%
%   Inputs:           CIJ,  binary/weighted undirected adjacency matrix.
%
%   Outputs: v_pos, v_neg,  eigenvectors associated with the largest
%                           eigenvalues of the adjacency matrix CIJ for 
%                           positive and negative weights separately.
%
%   Reference: Newman, MEJ (2002). The mathematics of networks.
%
%   Contributors:
%   Xi-Nian Zuo, Chinese Academy of Sciences, 2010
%   Rick Betzel, Indiana University, 2012
%   Mika Rubinov, University of Cambridge, 2015

%   MODIFICATION HISTORY
%   2010/2012: original (XNZ, RB)
%   2015: ensure the use of leading eigenvector (MR)
%   2016: computed for positive and negative weights separately (JMS)


n = length(CIJ);
CIJ_pos = threshold_absolute(CIJ,0);
CIJ_pos(1:n+1:end) = 0;
CIJ_neg = threshold_absolute(CIJ*-1,0);
CIJ_neg(1:n+1:end) = 0;

if n < 1000
    [V_pos,D_pos] = eig(CIJ_pos);
    [V_neg,D_neg] = eig(CIJ_neg);
else
    [V_pos,D_pos] = eigs(sparse(CIJ_pos));
    [V_neg,D_neg] = eigs(sparse(CIJ_neg));
end

[~,idx_pos] = max(diag(D_pos));
[~,idx_neg] = max(diag(D_neg));

ec_pos = abs(V_pos(:,idx_pos));
ec_neg = abs(V_neg(:,idx_neg));

v_pos = reshape(ec_pos,length(ec_pos),1);
v_neg = reshape(ec_neg,length(ec_neg),1);
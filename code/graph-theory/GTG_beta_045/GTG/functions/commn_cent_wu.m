function [CC_pos,CC_neg] = commn_cent_wu(W,Ci)
%COMMN_CENT_WU_SIGN     Commn Centrality
%
%   [CC_pos,CC_neg] = commn_cent_wu(W,Ci);
%
%   Commn centrality is a weighted combination of intra-module-strength and
%   extra-module-strength. This reflects the extent to which a node is an 
%   important bridge between modules. Adapted from Gupta et al. to allow 
%   weighted/signed networks (if a binary matrix is input, the output will 
%   reflect the original formula). 
%
%   Inputs:     W,        undirected binary or weighted connection matrix
%                         (if weighted, must have positive weights, 
%                         negative weights optional)
%
%   Output:     CC_pos, commn centrality from positive weights
%               CC_neg, commn centrality from negative weights
%
%   Reference: Gupta N, Singh A, Cherifi H. arXiv:1601.07108v1
%
%
%   Jeff Spielberg, Boston University

%   Modification History:
%   February 2016: Original

n              = size(W,1);
W(1:(n+1):end) = 0;
W_pos          = W.*(W>0);
W_neg          = -W.*(W<0);

% Ensure that module labels are continuous
comms = unique(Ci);
ncs = length(comms);
if ~all(comms(:)==(1:ncs)')
    for cu = 1:ncs
        Ci(comms(cu)) = cu;
    end
end

% Preallocate space
Kc_pos            = zeros(n,ncs);
Kc_neg            = zeros(n,ncs);
intra_cln_pos     = zeros(n,1);
intra_cln_neg     = zeros(n,1);
extra_cln_pos     = zeros(n,1);
extra_cln_neg     = zeros(n,1);
intra_cl_pos      = zeros(ncs,1);
intra_cl_neg      = zeros(ncs,1);
extra_cl_pos      = zeros(ncs,1);
extra_cl_neg      = zeros(ncs,1);
max_intra_cln_pos = zeros(ncs,1);
max_intra_cln_neg = zeros(ncs,1);
max_extra_cln_pos = zeros(ncs,1);
max_extra_cln_neg = zeros(ncs,1);
CC_pos            = zeros(n,1);
CC_neg            = zeros(n,1);

% Map community affiliation
Gc_pos = (W_pos~=0)*diag(Ci);
Gc_neg = (W_neg~=0)*diag(Ci);

% Calculate node strength for connections to each module
for cu = 1:ncs
    Kc_pos(:,cu) = sum(W_pos.*(Gc_pos==cu),2);
    Kc_neg(:,cu) = sum(W_neg.*(Gc_neg==cu),2);
end

for cu = 1:ncs
    % Calculate intra-strength for each node
    intra_cln_pos(Ci==cu)  = Kc_pos(Ci==cu,cu);
    intra_cln_neg(Ci==cu)  = Kc_neg(Ci==cu,cu);
    
    % Calculate extra-strength for each node
    extra_cln_pos(Ci==cu) = sum(Kc_pos(Ci==cu,:),2)-Kc_pos(Ci==cu,cu);
    extra_cln_neg(Ci==cu) = sum(Kc_neg(Ci==cu,:),2)-Kc_neg(Ci==cu,cu);
    
    % Calculate intra-strength for each community
    intra_cl_pos(cu)  = sum(intra_cln_pos(Ci==cu))/2;
    intra_cl_neg(cu)  = sum(intra_cln_neg(Ci==cu))/2;
    
    % Calculate extra-strength for each community
    extra_cl_pos(cu) = sum(extra_cln_pos(Ci==cu));
    extra_cl_neg(cu) = sum(extra_cln_neg(Ci==cu));
    
    % Find maximum in-strength for each community
    max_intra_cln_pos(cu)  = max(intra_cln_pos(Ci==cu));
    max_intra_cln_neg(cu)  = max(intra_cln_neg(Ci==cu));
    
    % Find maximum extra-strength for each community
    max_extra_cln_pos(cu) = max(extra_cln_pos(Ci==cu));
    max_extra_cln_neg(cu) = max(extra_cln_neg(Ci==cu));
end

% Compute ratio of intra-strength to total strength for each community
mu_c_pos = intra_cl_pos./(intra_cl_pos+extra_cl_pos);
mu_c_neg = intra_cl_neg./(intra_cl_neg+extra_cl_neg);

% Calculate commn centrality
for cn = 1:n
    CC_pos(cn) = ((1+mu_c_pos(Ci(cn)))*intra_cln_pos(cn)) + ((1-mu_c_pos(Ci(cn)))*(((extra_cln_pos(cn)/max_extra_cln_pos(Ci(cn)))*max_intra_cln_pos(Ci(cn)))^2));
    CC_neg(cn) = ((1+mu_c_neg(Ci(cn)))*intra_cln_neg(cn)) + ((1-mu_c_neg(Ci(cn)))*(((extra_cln_neg(cn)/max_extra_cln_neg(Ci(cn)))*max_intra_cln_neg(Ci(cn)))^2));
end

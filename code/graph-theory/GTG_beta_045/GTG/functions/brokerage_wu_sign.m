function [bkg_pos,bkg_neg,bkg_tot_pos,bkg_tot_neg] = brokerage_wu_sign(W)
%BROKERAGE_WU_SIGN     Brokerage
%
%   [bkg_pos,bkg_neg,bkg_tot_pos,bkg_tot_neg] = brokerage_wu_sign(W);
%
%   Brokerage is a node-strength-weighted version of betweeness centrality.
%   This measure reflects the extent to which a node controls shortcut
%   'bridges.' Adapted from Everett et al. to allow weighted/signed 
%   networks. 
%
%   Inputs:     W,        undirected connection matrix with positive and
%                         possibly negative weights
%
%   Output:     bkg_pos,     brokerage from positive weights
%               bkg_neg,     brokerage from negative weights
%               bkg_pos_tot, total brokerage from positive weights
%               bkg_neg_tot, total brokerage from negative weights
%
%   Reference: Everett M, Valente TW. (2016). Social Networks 44, 202-208.
%
%
%   Jeff Spielberg, Boston University

%   Modification History:
%   February 2016: Original

n              = size(W,1);
W(1:(n+1):end) = 0;
W_pos          = W.*(W>0);
W_neg          = -W.*(W<0);

% Compute node betweeness
bkg_pos(:,1) = betweenness_wei(weight_conversion(W_pos,'lengths')).*2;
bkg_neg(:,1) = betweenness_wei(weight_conversion(W_neg,'lengths')).*2;

% Compute node strengths
S_pos(:,1) = strengths_und(W_pos);
S_neg(:,1) = strengths_und(W_neg);

% Compute brokerage
bkg_pos(bkg_pos>0) = (bkg_pos(bkg_pos>0)+n-1)./S_pos(bkg_pos>0);
bkg_neg(bkg_neg>0) = (bkg_neg(bkg_neg>0)+n-1)./S_neg(bkg_neg>0);

% Compute total brokerage for the network
bkg_tot_pos = sum(bkg_pos);
bkg_tot_neg = sum(bkg_neg);

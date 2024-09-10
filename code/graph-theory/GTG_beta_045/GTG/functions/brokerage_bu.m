function [bkg,bkg_tot] = brokerage_bu(W)
%BROKERAGE_WU     Brokerage
%
%   [bkg_pos,bkg_tot_pos] = brokerage_wu(W);
%
%   Brokerage is a degree-weighted version of betweeness centrality.
%   This measure reflects the extent to which a node controls shortcut
%   'bridges.' 
%
%   Inputs:     W,       undirected binary connection matrix
%
%   Output:     bkg,     brokerage
%               bkg_tot, total brokerage
%
%   Reference: Everett M, Valente TW. (2016). Social Networks 44, 202-208.
%
%
%   Jeff Spielberg, Boston University

%   Modification History:
%   February 2016: Original

n              = size(W,1);
W(1:(n+1):end) = 0;

% Compute node betweeness
bkg(:,1) = betweenness_bin(weight_conversion(W,'lengths')).*2;

% Compute degrees
S(:,1) = degrees_und(W);

% Compute brokerage
bkg(bkg>0) = (bkg(bkg>0)+n-1)./S(bkg>0);

% Compute total brokerage for the network
bkg_tot = sum(bkg);

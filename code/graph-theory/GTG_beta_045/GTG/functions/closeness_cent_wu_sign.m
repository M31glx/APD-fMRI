function [CC_pos,CC_neg] = closeness_cent_wu_sign(W,flag)

% Author: Jeffrey M. Spielberg (jspielb2@gmail.com)
% Version: 05.20.15
% 
% Calculates the closeness centrality of a node. 
% 
% USAGE:
% 
% W      = weighted (signed) connectivity matrix
% 
% CC_pos = closeness centrality values for each node based on positive
%                  weights
% 
% CC_neg = closeness centrality values for each node based on negative
%                  weights
% 
% WARNING: This is a beta version. There no known bugs, but only limited 
% testing has been perfomed. This software comes with no warranty (even the
% implied warranty of merchantability or fitness for a particular purpose).
% Therefore, USE AT YOUR OWN RISK!!!
%
% Copyleft 2014-2015. Software can be modified and redistributed, but 
% modifed, redistributed versions must have the same rights.

if ~exist('flag','var')
    flag = 1;
end

n              = size(W,1);
W(1:(n+1):end) = 0;

L     = weight_conversion(W,'lengths');
L_pos = L.*(L>0);
L_neg = -L.*(L<0);

D_pos = distance_wei(L_pos);
D_neg = distance_wei(L_neg);

if any(isnan(L_pos))
    CC_pos = NaN(n,1);
else
    switch flag
        case 1 % Use arithmetic mean
            CC_pos = (1./sum(D_pos,2))./(n-1);
        case 2 % Use harmonic mean (may work better in disconnected matrices)
            invD_pos              = 1./D_pos;
            invD_pos(1:(n+1):end) = 0;   % Replace infs with 0
            CC_pos                = (sum(invD_pos,2))./(n-1);
        case 3 % Dangalchev (2006)'s formula
            invD_pos              = 1./(2.^D_pos);
            invD_pos(1:(n+1):end) = 0;   % Replace infs with 0
            CC_pos                = (sum(invD_pos,2))./(n-1);
    end
end

if any(isnan(L_neg))
    CC_neg = NaN(n,1);
else
    switch flag
        case 1 % Use arithmetic mean
            CC_neg = (1./sum(D_neg,2))./(n-1);
        case 2 % Use harmonic mean (may work better in disconnected matrices)
            invD_neg              = 1./D_neg;
            invD_neg(1:(n+1):end) = 0;   % Replace infs with 0
            CC_neg                = (sum(invD_neg,2))./(n-1);
        case 3 % Dangalchev (2006)'s formula
            invD_neg              = 1./(2.^D_neg);
            invD_neg(1:(n+1):end) = 0;   % Replace infs with 0
            CC_neg                = (sum(invD_neg,2))./(n-1);
    end
end


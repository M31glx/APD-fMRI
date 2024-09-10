function [outmats] = reignin_graph_out(inmats,maxval)

% Winsorize (reign-in) outlier values in a set of connectivity matrices.
% Specify the maximum value in terms of standard deviations. However, 
% returned matrices will be in original units (instead of in standard 
% deviation units). 

% Author: Jeffrey M. Spielberg (jspielb2@gmail.com)
% Version: 03.23.15
% 
% WARNING: This is a beta version. There no known bugs, but only limited 
% testing has been perfomed. This software comes with no warranty (even the
% implied warranty of merchantability or fitness for a particular purpose).
% Therefore, USE AT YOUR OWN RISK!!!
%
% Copyleft 2014-2016. Software can be modified and redistributed, but 
% modifed, redistributed versions must have the same rights

maxval = abs(maxval);
minval = -maxval;

zmats = zscore(inmats,0,3);
high_outliers = logical(zmats > maxval);
low_outliers  = logical(zmats < minval);

newhighvals = repmat(mean(inmats,3),[1,1,size(inmats,3)])+(repmat(std(inmats,0,3),[1,1,size(inmats,3)]).*maxval);
newlowvals = repmat(mean(inmats,3),[1,1,size(inmats,3)])+(repmat(std(inmats,0,3),[1,1,size(inmats,3)]).*minval);

outmats = inmats;
outmats(high_outliers) = newhighvals(high_outliers);
outmats(low_outliers)  = newlowvals(low_outliers);

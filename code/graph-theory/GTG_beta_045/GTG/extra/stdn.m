function out = stdn(in)

% Calculate standard deviation across all dimensions
%
% Author: Jeffrey M. Spielberg (jspielb2@gmail.com)
% Version: 12.16.15
% 
% WARNING: This is a beta version. There no known bugs, but only limited 
% testing has been perfomed. This software comes with no warranty (even the
% implied warranty of merchantability or fitness for a particular purpose).
% Therefore, USE AT YOUR OWN RISK!!!
%
% Copyleft 2014-2016. Software can be modified and redistributed, but 
% modifed, redistributed versions must have the same rights

out = std(in(:));

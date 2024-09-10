function out_ts = rescale_ts(in_ts,maxval,minval)

% Rescale vector to the specified minimum and maximum value
% 
% Author: Jeffrey M. Spielberg (jspielb2@gmail.com)
% Version: 01.07.16
% 
% WARNING: This is a beta version. There no known bugs, but only limited 
% testing has been perfomed. This software comes with no warranty (even the
% implied warranty of merchantability or fitness for a particular purpose).
% Therefore, USE AT YOUR OWN RISK!!!
%
% Copyleft 2014-2016. Software can be modified and redistributed, but 
% modifed, redistributed versions must have the same rights

if nargin<3
    maxval = 1;
    minval = 0;
end

% Set new range
out_ts = in_ts.*((maxval-minval)./range(in_ts));

% Set new 0 point
out_ts = out_ts-min(out_ts)+minval;

function [r] = corrn(X1,X2,mask)

% Pearson correlation coefficient for matrices. 
% 
% Author: Jeffrey M. Spielberg (jspielb2@gmail.com)
% Version: 01.30.16
% 
% WARNING: This is a beta version. There no known bugs, but only limited 
% testing has been perfomed. This software comes with no warranty (even the
% implied warranty of merchantability or fitness for a particular purpose).
% Therefore, USE AT YOUR OWN RISK!!!
%
% Copyleft 2014-2016. Software can be modified and redistributed, but 
% modifed, redistributed versions must have the same rights

if any(size(X1) ~= size(X2))
    error('X1 and X2 must be the same size');
end

if ~isa(X1,'double')
    X1 = double(X1);
end
if ~isa(X2,'double')
    X2 = double(X2);
end

if nargin==2
    mask = true(size(X1));
elseif numel(size(mask))==numel(size(X1)) && size(mask)==size(X1)
    mask = logical(mask);
elseif ndims(mask)<ndims(X1)
    dims = size(X1);
    mask = logical(repmat(mask,[ones(1,ndims(mask)),dims(end)]));
end

X1 = X1(mask)-mean(X1(mask));
X2 = X2(mask)-mean(X2(mask));
X11prod  = X1.*X1;
X22prod  = X2.*X2;
X12prod  = X1.*X2;
r = sum(X12prod(:))/sqrt(sum(X11prod(:))*sum(X22prod(:)));

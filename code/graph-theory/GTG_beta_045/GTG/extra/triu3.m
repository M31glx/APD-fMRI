function out = triu3(mat,k)

% Return upper triangular matrix for 3-dimensional matrix (upper matrix is
% from the first two dimensions). 

% Author: Jeffrey M. Spielberg (jspielb2@gmail.com)
% Version: 01.24.16
%
% WARNING: This is a beta version. There no known bugs, but only limited
% testing has been perfomed. This software comes with no warranty (even the
% implied warranty of merchantability or fitness for a particular purpose).
% Therefore, USE AT YOUR OWN RISK!!!
%
% Copyleft 2014-2016. Software can be modified and redistributed, but 
% modifed, redistributed versions must have the same rights. 

if ~exist('k','var')
    k = 0;
end

out = zeros(size(mat));
for d = 1:size(mat,3)
    out(:,:,d) = triu(mat(:,:,d),k);
end


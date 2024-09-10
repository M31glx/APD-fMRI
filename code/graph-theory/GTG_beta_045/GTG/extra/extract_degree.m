function [deg,outlabels] = extract_degree(W,inlabels)

% Author: Jeffrey M. Spielberg (jspielb2@gmail.com)
% Version: 03.15.16
% 
% WARNING: This is a beta version. There no known bugs, but only limited 
% testing has been perfomed. This software comes with no warranty (even the
% implied warranty of merchantability or fitness for a particular purpose).
% Therefore, USE AT YOUR OWN RISK!!!
%
% Copyleft 2014-2016. Software can be modified and redistributed, but 
% modifed, redistributed versions must have the same rights


W                    = full(W)~=0;
W(1:size(W,1)+1:end) = 0;
if maxn(abs(W-W'))>1e-1
    W = W+W';
end
[deg,i]   = sort(sum(W,2),'descend');
outlabels = inlabels(i);
if nargout==1
    deg = [outlabels,num2cell(deg)];
end

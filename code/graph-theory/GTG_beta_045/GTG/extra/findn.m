function varargout = findn(X)

% Find within an n-dimensional matrix. Input should be the desired logical
% relationship (the same usage as find.m). Up to 20 dimensions are
% accepted. 
%
% Example:  [x,y,z] = findx(W==45)
% 
% 
% Author: Jeffrey M. Spielberg (jspielb2@gmail.com)
% Version: 10.27.15
% 
% WARNING: This is a beta version. There no known bugs, but only limited 
% testing has been perfomed. This software comes with no warranty (even the
% implied warranty of merchantability or fitness for a particular purpose).
% Therefore, USE AT YOUR OWN RISK!!!
%
% Copyleft 2014-2016. Software can be modified and redistributed, but 
% modifed, redistributed versions must have the same rights

if nargout==1
    varargout{1} = find(X==1);
elseif nargout==2
    [varargout{1},varargout{2}] = find(X==1);
elseif nargout==3
    [varargout{1},varargout{2},varargout{3}] = ind2sub(size(X),find(X==1));
elseif nargout==4
    [varargout{1},varargout{2},varargout{3},varargout{4}] = ind2sub(size(X),find(X==1));
elseif nargout==5
    [varargout{1},varargout{2},varargout{3},varargout{4},varargout{5}] = ind2sub(size(X),find(X==1));
elseif nargout==6
    [varargout{1},varargout{2},varargout{3},varargout{4},varargout{5}, ...
        varargout{6}] = ind2sub(size(X),find(X==1));
elseif nargout==7
    [varargout{1},varargout{2},varargout{3},varargout{4},varargout{5}, ...
        varargout{6},varargout{7}] = ind2sub(size(X),find(X==1));
elseif nargout==8
    [varargout{1},varargout{2},varargout{3},varargout{4},varargout{5}, ...
        varargout{6},varargout{7},varargout{8}] = ind2sub(size(X),find(X==1));
elseif nargout==9
    [varargout{1},varargout{2},varargout{3},varargout{4},varargout{5}, ...
        varargout{6},varargout{7},varargout{8},varargout{9}] = ind2sub(size(X),find(X==1));
elseif nargout==10
    [varargout{1},varargout{2},varargout{3},varargout{4},varargout{5}, ...
        varargout{6},varargout{7},varargout{8},varargout{9}, ...
        varargout{10}] = ind2sub(size(X),find(X==1));
elseif nargout==11
    [varargout{1},varargout{2},varargout{3},varargout{4},varargout{5}, ...
        varargout{6},varargout{7},varargout{8},varargout{9}, ...
        varargout{10},varargout{11}] = ind2sub(size(X),find(X==1));
elseif nargout==12
    [varargout{1},varargout{2},varargout{3},varargout{4},varargout{5}, ...
        varargout{6},varargout{7},varargout{8},varargout{9}, ...
        varargout{10},varargout{11},varargout{12}] = ind2sub(size(X),find(X==1));
elseif nargout==13
    [varargout{1},varargout{2},varargout{3},varargout{4},varargout{5}, ...
        varargout{6},varargout{7},varargout{8},varargout{9}, ...
        varargout{10},varargout{11},varargout{12},varargout{13}] = ind2sub(size(X),find(X==1));
elseif nargout==14
    [varargout{1},varargout{2},varargout{3},varargout{4},varargout{5}, ...
        varargout{6},varargout{7},varargout{8},varargout{9}, ...
        varargout{10},varargout{11},varargout{12},varargout{13}, ...
        varargout{14}] = ind2sub(size(X),find(X==1));
elseif nargout==15
    [varargout{1},varargout{2},varargout{3},varargout{4},varargout{5}, ...
        varargout{6},varargout{7},varargout{8},varargout{9}, ...
        varargout{10},varargout{11},varargout{12},varargout{13}, ...
        varargout{14},varargout{15}] = ind2sub(size(X),find(X==1));
elseif nargout==16
    [varargout{1},varargout{2},varargout{3},varargout{4},varargout{5}, ...
        varargout{6},varargout{7},varargout{8},varargout{9}, ...
        varargout{10},varargout{11},varargout{12},varargout{13}, ...
        varargout{14},varargout{15},varargout{16}] = ind2sub(size(X),find(X==1));
elseif nargout==17
    [varargout{1},varargout{2},varargout{3},varargout{4},varargout{5}, ...
        varargout{6},varargout{7},varargout{8},varargout{9}, ...
        varargout{10},varargout{11},varargout{12},varargout{13}, ...
        varargout{14},varargout{15},varargout{16},varargout{17}] = ind2sub(size(X),find(X==1));
elseif nargout==18
    [varargout{1},varargout{2},varargout{3},varargout{4},varargout{5}, ...
        varargout{6},varargout{7},varargout{8},varargout{9}, ...
        varargout{10},varargout{11},varargout{12},varargout{13}, ...
        varargout{14},varargout{15},varargout{16},varargout{17}, ...
        varargout{18}] = ind2sub(size(X),find(X==1));
elseif nargout==19
    [varargout{1},varargout{2},varargout{3},varargout{4},varargout{5}, ...
        varargout{6},varargout{7},varargout{8},varargout{9}, ...
        varargout{10},varargout{11},varargout{12},varargout{13}, ...
        varargout{14},varargout{15},varargout{16},varargout{17}, ...
        varargout{18},varargout{19}] = ind2sub(size(X),find(X==1));
elseif nargout==20
    [varargout{1},varargout{2},varargout{3},varargout{4},varargout{5}, ...
        varargout{6},varargout{7},varargout{8},varargout{9}, ...
        varargout{10},varargout{11},varargout{12},varargout{13}, ...
        varargout{14},varargout{15},varargout{16},varargout{17}, ...
        varargout{18},varargout{19},varargout{20}] = ind2sub(size(X),find(X==1));
end

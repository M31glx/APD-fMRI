function [mapping,consistency] = match_grp(grp1,grp2)

% Find mapping between two group assignments with the largest overlap. 
%
% Usage: [mapping,consistency] = match_grp(grp1,grp2)
% 
% Inputs: grp1         = First group assignment for each entry
%         grp2         = Second group assignment for each entry
% 
% Outputs: mapping     = Mapping between groups
%          consistency = How similar are the two groupings, based on the
%                        computed mapping (i.e., # of entries that
%                        overlap/total # of entries). Range = 0:1.
% 
% 
% Author: Jeffrey M. Spielberg (jspielb2@gmail.com)
% Version: 03.09.16
% 
% WARNING: This is a beta version. There no known bugs, but only limited 
% testing has been perfomed. This software comes with no warranty (even the
% implied warranty of merchantability or fitness for a particular purpose).
% Therefore, USE AT YOUR OWN RISK!!!
%
% Copyleft 2014-2016. Software can be modified and redistributed, but 
% modifed, redistributed versions must have the same rights

vals1 = unique(grp1);
vals2 = unique(grp2);
count = 1;
olp   = zeros(3,length(vals1)*length(vals2));

for v1 = 1:length(vals1)
    for v2 = 1:length(vals2)
        olp(1,count) = vals1(v1);
        olp(2,count) = vals2(v2);
        olp(3,count) = sum((grp1==vals1(v1)).*(grp2==vals2(v2)));
        count = count+1;
    end
end

if v1<v2
    mapping = zeros(3,length(vals1));
    for v = 1:v1
        h = find(olp(3,:)==max(olp(3,:)),1);
        mapping(:,v) = olp(:,h);
        olp(:,olp(2,:)==mapping(2,v)) = [];
        olp(:,olp(3,:)==mapping(3,v)) = [];
    end
else
    mapping = zeros(3,length(vals2));
    for v = 1:v2
        h = find(olp(3,:)==max(olp(3,:)),1);
        mapping(:,v) = olp(:,h);
        olp(:,olp(2,:)==mapping(2,v)) = [];
        olp(:,olp(3,:)==mapping(3,v)) = [];
    end
end

consistency = sum(mapping(3,:))/numel(grp1);

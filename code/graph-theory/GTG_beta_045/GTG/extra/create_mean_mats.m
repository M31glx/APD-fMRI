function [full_mean_mat,grp_mean_mat] = create_mean_mats(Z,varargin)

% Calculate mean matrices from a set of matrices using fisher r-to-z
% transformations. Mean matrices for subgroups within the larger set can
% also be calculated by either specifying a grouping variable or specifying
% the number of groups/participants per group if no specific grouping is
% desired. If the grp_mean_mat output is specified, but no grouping
% information provided, the script will default to groups of approximately
% 50. 
% 
% Author: Jeffrey M. Spielberg (jspielb2@gmail.com)
% Version: 07.29.14
% 
% WARNING: This is a beta version. There are no known bugs, but only limited 
% testing has been perfomed. This software comes with no warranty (even the
% implied warranty of merchantability or fitness for a particular purpose).
% Therefore, USE AT YOUR OWN RISK!!!
%
% Copyleft 2014-2016. Software can be modified and redistributed, but 
% modifed, redistributed versions must have the same rights.

if nargin > 1
    grouping_var = varargin{1};
    contin       = varargin{2};
    if nargin > 3
        switch varargin{3};
            case '-num_grps'
                num_grps     = varargin{4};
            case '-sub_per_grp'
                subs_per_grp = varargin{4};
        end
    else
        subs_per_grp = 50;
    end
end

fish_z                                       = atanh(Z);
full_mean_mat                                = atan(mean(fish_z,3));
full_mean_mat(1:size(full_mean_mat,1)+1:end) = 0;

if exist('grouping_var','var')
    if isempty(contin) || var(contin) == 0 
        group_mem = grouping_var;
    else
        if exist('num_grps','var')
            subs_per_grp = round(size(Z,3)/num_grps);
        else
            num_grps     = round(size(Z,3)/subs_per_grp);
            subs_per_grp = round(size(Z,3)/num_grps);
        end
        
        [~,ind] = sort(grouping_var);
        fish_z  = fish_z(:,:,ind);
        
        for grp = 1:num_grps
            if grp ~= num_grps
                group_mem((((grp-1)*subs_per_grp)+1):(grp*subs_per_grp)) = grp;
            else
                group_mem((((grp-1)*subs_per_grp)+1):size(Z,3))          = grp;
            end
        end
    end
    grps         = unique(group_mem);
    grp_mean_mat = zeros([size(full_mean_mat),length(grps)]);
    for grp = 1:length(grps)
        temp_conmat                              = atan(mean(fish_z(:,:,logical(group_mem == grps(grp))),3));
        temp_conmat(1:size(temp_conmat,1)+1:end) = 0;
        grp_mean_mat(:,:,grp)                    = temp_conmat;
    end
end
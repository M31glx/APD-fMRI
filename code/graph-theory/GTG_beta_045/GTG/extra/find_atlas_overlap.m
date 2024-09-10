function overlap_info = find_atlas_overlap(at1,at2,at2_names,sim_measure,thr1,two_labels,thr2)

% Find overlap between sets of ROI atlases, and label at1 ROIs with 1 or 2 
% labels from the at2 ROI(s) with the greatest overlap. 
%
% Inputs:  at1:         First atlas. 
%          at2:         Second atlas (for which you already have labels). 
%          at2_names:   Labels for each ROI in at2. 
%          sim_measure: Measure used to determine 'overlap'. Default is an
%                       asymmetric dice coefficient, but any measure can be
%                       used that is present in simbin.m.
%          two_labels:  Specify whether to use the top 1 or 2 overlapping
%                       labels. 
% 
% 
% 
% Author: Jeffrey M. Spielberg (jspielb2@gmail.com)
% Version: 03.07.16
% 
% WARNING: This is a beta version. There no known bugs, but only limited 
% testing has been perfomed. This software comes with no warranty (even the
% implied warranty of merchantability or fitness for a particular purpose).
% Therefore, USE AT YOUR OWN RISK!!!
%
% Copyleft 2014-2016. Software can be modified and redistributed, but 
% modifed, redistributed versions must have the same rights

if ~exist('thr1','var')
    thr1 = 0;
end
if ~exist('two_labels','var')
    two_labels = 0;
elseif ~exist('thr2','var')
    thr2 = 0;
end
if ~exist('sim_measure','var')
    sim_measure = 'diceasym1';
end


at1                       = double(at1);
at2                       = double(at2);
at1_labels                = unique(at1);
at1_labels(at1_labels==0) = [];
overlap_info              = cell(length(at1_labels),2);

if mod(size(at1,1),2)
    mid_ROIs = unique(at1(ceil(size(at1,1)/2),:,:));
else
    mid_ROIs = unique(at1(floor(size(at1,1)/2):ceil(size(at1,1)/2),:,:));
end

for curr_lab = 1:length(at1_labels)
    curr_at1_mask                 = at1==at1_labels(curr_lab);
    overlap_info{curr_lab,2}      = at1_labels(curr_lab);
    curr_overlap                  = curr_at1_mask.*at2;
    overlap_labs                  = unique(curr_overlap);
    overlap_labs(overlap_labs==0) = [];
    hval                          = 0;
    if two_labels
        hval2 = 0;
        overlap_info2 = cell(size(overlap_info));
    end
    
    for curr_overlap_lab = 1:length(overlap_labs)
        curr_at2_mask = at2==overlap_labs(curr_overlap_lab);
        sim           = simbin(curr_at1_mask,curr_at2_mask,sim_measure);
        if sim>hval && sim>thr1
            overlap_info{curr_lab,1} = at2_names{overlap_labs(curr_overlap_lab)};
            hval                     = sim;
        elseif two_labels && sim>hval2 && sim>thr2
            overlap_info2{curr_lab,1} = at2_names{overlap_labs(curr_overlap_lab)};
            hval2                     = sim;
        end
    end
    if two_labels
        if ~isempty(overlap_info2{curr_lab,1})
            overlap_info{curr_lab,1} = [overlap_info{curr_lab,1},'__',overlap_info2{curr_lab,1}];
        end
    end
    
    bal = sumn(curr_at1_mask(1:floor(size(curr_at1_mask,1)/2),:,:))/sumn(curr_at1_mask);
    if isempty(overlap_info{curr_lab,1})
        if bal>0.5
            overlap_info{curr_lab,1} = 'R_NONE';
        else
            overlap_info{curr_lab,1} = 'L_NONE';
        end
    end
    
    if any(ismember(mid_ROIs,at1_labels(curr_lab))) && bal>0.3 && bal<0.7
        overlap_info{curr_lab,1}(1:2) = 'M_';
    end
end

ulabs = unique(overlap_info(:,1));
for ulab = 1:length(ulabs)
    if sum(ismember(overlap_info(:,1),ulabs{ulab})) > 1
        ulab_locs = find(ismember(overlap_info(:,1),ulabs{ulab}));
        for lab = 1:length(ulab_locs)
            overlap_info{ulab_locs(lab),1} = [overlap_info{ulab_locs(lab),1},'_',num2str(lab)];
        end
    end
end
